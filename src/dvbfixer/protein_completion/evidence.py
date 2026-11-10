"""Private manifest-last cohort evidence, deliberately disjoint from public bundles.

Machine-owned canonical JSON and complete deterministic replay are required.
Hashes establish integrity, not authorship or scientific qualification.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.bundle_commit import commit_manifest_last
from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes
from dvbfixer.protein_completion.cohort import (
    COHORT_REVISION,
    MAX_TARGETS,
    CohortPlan,
    CohortResult,
    MaterializedCohort,
    audit_materialized_cohort,
    complete_cohort,
    materialize_cohort,
    plan_cohort,
)
from dvbfixer.protein_completion.source import SourceCompletionError

KIND = "experimental-protein-completion-cohort-evidence"
MAX_ARTIFACT_BYTES = 64*1024*1024
MAX_TOTAL_BYTES = 128*1024*1024
_BASE = {"source.pdb", "targets.json", "plan.json", "attempts.json", "audit.json"}
_ACCEPTED = _BASE | {"candidate.pdb", "candidate.dat"}
_SAFE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}\Z")


@dataclass(frozen=True)
class CohortEvidence:
    status: Literal["accepted", "rejected"]
    stage: Literal["raw", "rounded"]
    plan: CohortPlan
    result: CohortResult
    materialized: MaterializedCohort | None


def _json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                       allow_nan=False)+"\n").encode("utf-8")


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _prepare(source: bytes, plan: CohortPlan, result: CohortResult) -> tuple[CohortEvidence, dict[str, bytes]]:
    # Never trust result status, attempt audits, selected indices or omissions.
    if not isinstance(plan, CohortPlan) or not isinstance(result, CohortResult):
        raise SourceCompletionError("original cohort plan/result required; other policies need separate evidence contracts")
    authoritative = plan_cohort(source, tuple(r.residue for r in plan.requests))
    if _json(asdict(plan)) != _json(asdict(authoritative)):
        raise SourceCompletionError("cohort evidence plan differs from canonical source replay")
    expected = complete_cohort(source, authoritative)
    if _json(asdict(result)) != _json(asdict(expected)):
        raise SourceCompletionError("retained cohort result differs from deterministic replay")
    status: Literal["accepted", "rejected"] = "rejected"
    stage: Literal["raw", "rounded"] = "raw"
    refusal: str | None = result.reason
    materialized = None
    if result.selected is not None:
        stage = "rounded"
        try:
            materialized = materialize_cohort(source, plan, result.selected)
        except SourceCompletionError as exc:
            # Only a known scientific rounded-audit rejection is evidence of a
            # scientific refusal. I/O, programming errors and other failures
            # must propagate, not be relabelled as failed geometry.
            if not str(exc).startswith("rounded joint PDB rejected:"):
                raise
            refusal = str(exc)
        else:
            status, refusal = "accepted", None
    audit = {"status": status, "stage": stage, "refusal": refusal,
             "source_sha256": _sha(source), "plan_digest": plan.digest,
             "generated": [asdict(a) for a in materialized.generated] if materialized else [],
             "public_profile_enabled": False, "scientific_qualification": False}
    payloads = {"source.pdb": source, "targets.json": _json([asdict(r.residue) for r in plan.requests]),
                "plan.json": _json(asdict(plan)), "attempts.json": _json(asdict(result)), "audit.json": _json(audit)}
    if materialized is not None:
        payloads["candidate.pdb"] = materialized.pdb
        payloads["candidate.dat"] = _json(materialized.dat_record().to_dict())
    if (any(not 0 < len(data) <= MAX_ARTIFACT_BYTES for data in payloads.values())
            or sum(map(len, payloads.values())) > MAX_TOTAL_BYTES):
        raise SourceCompletionError("cohort evidence exceeds bounded artifact budget")
    return CohortEvidence(status, stage, plan, result, materialized), payloads


def _index(evidence: CohortEvidence, payloads: dict[str, bytes]) -> bytes:
    return _json({"schema_version": 1, "kind": KIND, "policy": COHORT_REVISION,
                  "status": evidence.status, "source_sha256": _sha(payloads["source.pdb"]),
                  "plan_digest": evidence.plan.digest,
                  "public_profile_enabled": False, "scientific_qualification": False,
                  "artifacts": {name: {"sha256": _sha(data), "bytes": len(data)} for name, data in payloads.items()}})


def _targets(data: bytes) -> tuple[ComponentInstanceRef, ...]:
    raw = json.loads(data)
    if not isinstance(raw, list) or not 1 <= len(raw) <= MAX_TARGETS:
        raise SourceCompletionError("bounded target array required")
    result = []
    keys = {"model", "chain_id", "sequence_number", "insertion_code", "alternate_location", "occurrence"}
    for row in raw:
        if (not isinstance(row, dict) or set(row) != keys or type(row["model"]) is not int
                or type(row["occurrence"]) is not int
                or any(not isinstance(row[k], str) for k in keys-{"model", "occurrence"})):
            raise SourceCompletionError("strict exact target fields required")
        result.append(ComponentInstanceRef(row["model"], row["chain_id"], row["sequence_number"],
                                           row["insertion_code"], row["alternate_location"], row["occurrence"]))
    return tuple(result)


def read_cohort_evidence(directory: Path, *, expected_source_sha256: str | None = None) -> CohortEvidence:
    """Require a committed index and replay source planning, all attempts and output.

    This reader accepts only its own experimental kind. Public diffusion and
    atom-completion bundle readers remain unchanged and reject these indices.
    """
    if directory.is_symlink() or not directory.is_dir():
        raise SourceCompletionError("cohort evidence requires a regular directory")
    root = directory.resolve(strict=True)
    raw_index = read_private_file_bytes(root, "bundle.json", max_bytes=64_000)
    index = json.loads(raw_index)
    if (not isinstance(index, dict) or index.get("kind") != KIND or type(index.get("schema_version")) is not int
            or index.get("schema_version") != 1 or index.get("policy") != COHORT_REVISION
            or index.get("status") not in ("accepted", "rejected")
            or index.get("public_profile_enabled") is not False or index.get("scientific_qualification") is not False):
        raise SourceCompletionError("unsupported/private-only cohort evidence index")
    names = _ACCEPTED if index["status"] == "accepted" else _BASE
    artifacts = index.get("artifacts")
    if (not isinstance(artifacts, dict) or set(artifacts) != names
            or {p.name for p in root.iterdir()} != names | {"bundle.json"}):
        raise SourceCompletionError("cohort evidence must contain exactly its flat declared artifacts")
    payloads = {}
    total = 0
    for name in sorted(names):
        row = artifacts[name]
        if (not isinstance(row, dict) or set(row) != {"sha256", "bytes"}
                or type(row["bytes"]) is not int or not 0 < row["bytes"] <= MAX_ARTIFACT_BYTES
                or not isinstance(row["sha256"], str) or re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is None):
            raise SourceCompletionError("invalid cohort artifact reference")
        total += row["bytes"]
        if total > MAX_TOTAL_BYTES:
            raise SourceCompletionError("cohort evidence exceeds total artifact budget")
        data = read_private_file_bytes(root, name, max_bytes=row["bytes"])
        if len(data) != row["bytes"] or _sha(data) != row["sha256"]:
            raise SourceCompletionError("cohort artifact digest/size mismatch")
        payloads[name] = data
    source = payloads["source.pdb"]
    if expected_source_sha256 is not None and _sha(source) != expected_source_sha256:
        raise SourceCompletionError("cohort source differs from externally expected SHA")
    plan = plan_cohort(source, _targets(payloads["targets.json"]))
    replay = complete_cohort(source, plan)
    evidence, expected_payloads = _prepare(source, plan, replay)
    if payloads != expected_payloads or raw_index != _index(evidence, expected_payloads):
        raise SourceCompletionError("cohort evidence differs from independent canonical replay")
    if evidence.materialized is not None:
        # Canonical .dat decoding remains owned by DatRecord; the full bytes above
        # also check derived totals/extra fields that a tolerant decoder can drop.
        dat = DatRecord.from_json(payloads["candidate.dat"].decode("utf-8"))
        if dat.to_dict() != evidence.materialized.dat_record().to_dict():
            raise SourceCompletionError("cohort sidecar generated ownership mismatch")
        chosen = evidence.result.selected
        if chosen is None:
            raise SourceCompletionError("accepted evidence has no selected cohort")
        audit_materialized_cohort(source, plan, chosen, payloads["candidate.pdb"])
    return evidence


def _durable_write(path: Path, data: bytes) -> None:
    with path.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())


def write_cohort_evidence(source: bytes, plan: CohortPlan, result: CohortResult, destination: Path) -> Path:
    """Preflight/replay before staging; exclusive manifest-last no-overwrite commit.

    Existing/crash-left paths are never adopted or deleted. Failed staging is
    left for diagnosis; generic commit rollback may remove its own transferred
    entries, never caller inputs. Original run/source files are not mutated.
    """
    if (not _SAFE_NAME.fullmatch(destination.name) or destination.parent.is_symlink()
            or destination.name in {".", ".."}):
        raise SourceCompletionError("explicit safe destination and regular parent required")
    destination = destination.parent.resolve(strict=True)/destination.name
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    evidence, payloads = _prepare(source, plan, result)
    staging = Path(tempfile.mkdtemp(prefix=".protein-cohort-staging-", dir=destination.parent))
    for name, data in payloads.items():
        _durable_write(staging/name, data)
    _durable_write(staging/"bundle.json", _index(evidence, payloads))
    read_cohort_evidence(staging, expected_source_sha256=_sha(source))
    commit_manifest_last(staging, destination)
    descriptor = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return destination
