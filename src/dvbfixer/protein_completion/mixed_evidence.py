"""Distinct private mixed-policy evidence with admission/raw/rounded replay.

The complete source derives modes, partners, attempts and global stereo anew.
Hashes are integrity checks, not authentication or scientific qualification.
"""
from __future__ import annotations

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
from dvbfixer.protein_completion.cohort import CohortResult
from dvbfixer.protein_completion.evidence import _durable_write, _json, _sha, _targets
from dvbfixer.protein_completion.mixed_cohort import (
    MIXED_POLICY,
    MaterializedMixedCohort,
    MixedCohortPlan,
    audit_materialized_mixed_cohort,
    complete_mixed_cohort,
    materialize_mixed_cohort,
    plan_mixed_cohort,
)
from dvbfixer.protein_completion.source import SourceCompletionError

KIND = "experimental-protein-completion-mixed-evidence"
MAX_ARTIFACT_BYTES = 64*1024*1024
MAX_TOTAL_BYTES = 128*1024*1024
_BASE = {"source.pdb", "targets.json", "audit.json"}
_PLANNED = _BASE | {"plan.json", "attempts.json"}
_ACCEPTED = _PLANNED | {"candidate.pdb", "candidate.dat"}


@dataclass(frozen=True)
class MixedEvidence:
    status: Literal["accepted", "rejected"]
    stage: Literal["admission", "raw", "rounded"]
    plan: MixedCohortPlan | None
    result: CohortResult | None
    materialized: MaterializedMixedCohort | None
    refusal: str | None


def _prepare(source: bytes, targets: tuple[ComponentInstanceRef, ...], *,
             plan: MixedCohortPlan | None = None,
             result: CohortResult | None = None) -> tuple[MixedEvidence, dict[str, bytes]]:
    if not isinstance(source, bytes) or not 0 < len(source) <= 16*1024*1024:
        raise SourceCompletionError("bounded nonempty original source bytes required")
    if (not isinstance(targets, tuple) or not 1 <= len(targets) <= 4
            or any(not isinstance(t, ComponentInstanceRef) for t in targets)):
        raise SourceCompletionError("immutable exact target tuple required")
    target_bytes = _json([asdict(t) for t in targets])
    # Strict bounded transport before scientific admission; duplicate/absent
    # targets and unsupported chemistry remain replayable admission refusals.
    _targets(target_bytes)
    if (plan is None) != (result is None):
        raise SourceCompletionError("supply both retained plan and result or neither")
    if plan is not None and (not isinstance(plan, MixedCohortPlan) or not isinstance(result, CohortResult)):
        raise SourceCompletionError("mixed plan/result required; old policy cannot opt in")
    status: Literal["accepted", "rejected"] = "rejected"
    stage: Literal["admission", "raw", "rounded"] = "admission"
    authoritative = None
    replay = None
    materialized = None
    refusal = None
    try:
        authoritative = plan_mixed_cohort(source, targets)
    except ValueError as exc:
        # Planning is a pure admission function. No I/O or generation occurs
        # here; runtime, programming and system errors are not swallowed.
        if plan is not None:
            raise SourceCompletionError("retained mixed plan no longer admits source") from exc
        refusal = str(exc)
    else:
        stage = "raw"
        replay = complete_mixed_cohort(source, authoritative)
        if plan is not None:
            if result is None:
                raise SourceCompletionError("retained mixed result required")
            if (_json(asdict(plan)) != _json(asdict(authoritative))
                    or _json(asdict(result)) != _json(asdict(replay))):
                raise SourceCompletionError("retained mixed plan/result differs from canonical replay")
        refusal = replay.reason
        if replay.selected is not None:
            stage = "rounded"
            try:
                materialized = materialize_mixed_cohort(source, authoritative, replay.selected)
            except SourceCompletionError as exc:
                if not str(exc).startswith("rounded mixed PDB rejected:"):
                    raise
                refusal = str(exc)
            else:
                status, refusal = "accepted", None
    evidence = MixedEvidence(status, stage, authoritative, replay, materialized, refusal)
    audit = {"status": status, "stage": stage, "refusal": refusal,
             "source_sha256": _sha(source), "plan_digest": authoritative.digest if authoritative else None,
             "rounded_audit": asdict(materialized.audit) if materialized else None,
             "public_profile_enabled": False, "scientific_qualification": False}
    payloads = {"source.pdb": source, "targets.json": target_bytes, "audit.json": _json(audit)}
    if authoritative is not None and replay is not None:
        payloads.update({"plan.json": _json(asdict(authoritative)), "attempts.json": _json(asdict(replay))})
    if materialized is not None:
        payloads.update({"candidate.pdb": materialized.pdb, "candidate.dat": _json(materialized.dat_record().to_dict())})
    if (any(not 0 < len(data) <= MAX_ARTIFACT_BYTES for data in payloads.values())
            or sum(map(len, payloads.values())) > MAX_TOTAL_BYTES):
        raise SourceCompletionError("mixed evidence exceeds bounded artifact budget")
    return evidence, payloads


def _index(evidence: MixedEvidence, payloads: dict[str, bytes]) -> bytes:
    return _json({"schema_version": 1, "kind": KIND, "policy": MIXED_POLICY,
                  "status": evidence.status, "stage": evidence.stage,
                  "source_sha256": _sha(payloads["source.pdb"]),
                  "plan_digest": evidence.plan.digest if evidence.plan else None,
                  "public_profile_enabled": False, "scientific_qualification": False,
                  "artifacts": {name: {"sha256": _sha(data), "bytes": len(data)} for name, data in payloads.items()}})


def read_mixed_evidence(directory: Path, *, expected_source_sha256: str | None = None) -> MixedEvidence:
    """Verify exact private inventory, then repeat admission/attempts/output audits."""
    if directory.is_symlink() or not directory.is_dir():
        raise SourceCompletionError("mixed evidence requires a regular directory")
    root = directory.resolve(strict=True)
    raw_index = read_private_file_bytes(root, "bundle.json", max_bytes=64_000)
    index = json.loads(raw_index)
    if (not isinstance(index, dict) or index.get("kind") != KIND or type(index.get("schema_version")) is not int
            or index.get("schema_version") != 1 or index.get("policy") != MIXED_POLICY
            or index.get("status") not in ("accepted", "rejected")
            or index.get("stage") not in ("admission", "raw", "rounded")
            or index.get("public_profile_enabled") is not False or index.get("scientific_qualification") is not False
            or (index["status"] == "accepted" and index["stage"] != "rounded")):
        raise SourceCompletionError("unsupported/private-only mixed evidence index")
    names = _ACCEPTED if index["status"] == "accepted" else _BASE if index["stage"] == "admission" else _PLANNED
    artifacts = index.get("artifacts")
    if (not isinstance(artifacts, dict) or set(artifacts) != names
            or {p.name for p in root.iterdir()} != names | {"bundle.json"}):
        raise SourceCompletionError("mixed evidence requires exactly its flat declared artifacts")
    payloads = {}
    total = 0
    for name in sorted(names):
        row = artifacts[name]
        if (not isinstance(row, dict) or set(row) != {"sha256", "bytes"}
                or type(row["bytes"]) is not int or not 0 < row["bytes"] <= MAX_ARTIFACT_BYTES
                or not isinstance(row["sha256"], str) or re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is None):
            raise SourceCompletionError("invalid mixed artifact reference")
        total += row["bytes"]
        if total > MAX_TOTAL_BYTES:
            raise SourceCompletionError("mixed evidence exceeds total artifact budget")
        data = read_private_file_bytes(root, name, max_bytes=row["bytes"])
        if len(data) != row["bytes"] or _sha(data) != row["sha256"]:
            raise SourceCompletionError("mixed artifact digest/size mismatch")
        payloads[name] = data
    source = payloads["source.pdb"]
    if expected_source_sha256 is not None and _sha(source) != expected_source_sha256:
        raise SourceCompletionError("mixed source differs from externally expected SHA")
    evidence, expected = _prepare(source, _targets(payloads["targets.json"]))
    if payloads != expected or raw_index != _index(evidence, expected):
        raise SourceCompletionError("mixed evidence differs from independent canonical replay")
    if evidence.materialized is not None:
        dat = DatRecord.from_json(payloads["candidate.dat"].decode("utf-8"))
        if dat.to_dict() != evidence.materialized.dat_record().to_dict():
            raise SourceCompletionError("mixed sidecar generated ownership mismatch")
        if evidence.plan is None or evidence.result is None or evidence.result.selected is None:
            raise SourceCompletionError("accepted mixed evidence missing plan or selected result")
        audit_materialized_mixed_cohort(source, evidence.plan, evidence.result.selected, payloads["candidate.pdb"])
    return evidence


def write_mixed_evidence(source: bytes, targets: tuple[ComponentInstanceRef, ...], destination: Path, *,
                         plan: MixedCohortPlan | None = None, result: CohortResult | None = None) -> Path:
    """Exclusive manifest-last commit; admission refusals preserve original bytes.

    Supplied plan/result are optional retained evidence, not scientific authority.
    Preflight precedes staging; existing paths are never adopted or overwritten.
    """
    if (re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}", destination.name) is None
            or destination.parent.is_symlink() or destination.name in {".", ".."}):
        raise SourceCompletionError("explicit safe mixed destination and regular parent required")
    destination = destination.parent.resolve(strict=True)/destination.name
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    evidence, payloads = _prepare(source, targets, plan=plan, result=result)
    staging = Path(tempfile.mkdtemp(prefix=".mixed-cohort-staging-", dir=destination.parent))
    for name, data in payloads.items():
        _durable_write(staging/name, data)
    _durable_write(staging/"bundle.json", _index(evidence, payloads))
    read_mixed_evidence(staging, expected_source_sha256=_sha(source))
    commit_manifest_last(staging, destination)
    descriptor = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return destination
