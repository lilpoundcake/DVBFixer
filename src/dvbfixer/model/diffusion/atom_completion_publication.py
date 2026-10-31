"""Manifest-last experimental completion bundles, disjoint from public bundles."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import ClassVar

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
from dvbfixer.model.diffusion.atom_completion_evidence import (
    CompletionEvidenceAudit,
    audit_completion_evidence,
)
from dvbfixer.model.diffusion.atom_completion_trace import CompletionSamplerTrace
from dvbfixer.model.diffusion.atom_completion_validation import audit_completion_candidate
from dvbfixer.model.diffusion.bundle_commit import commit_manifest_last
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    DiffusionContractError,
    _JsonContract,
)
from dvbfixer.model.diffusion.private_artifact import (
    read_private_artifact_bytes,
    read_private_file_bytes,
)
from dvbfixer.model.diffusion.scope import _parse_coordinate_records

_PAYLOAD_NAMES = {
    "request": "request.json", "source": "source.pdb", "motif": "motif.pdb", "reader": "reader.npz",
    "native": "native.npz", "raw": "raw.pdb", "refined": "candidate.pdb", "engine_view": "engine-view.json",
    "preflight": "preflight.json", "refinement": "refinement.json", "progress": "progress.json",
    "adapter_source": "adapter-source.py",
}


@dataclass(frozen=True, slots=True)
class CompletionBundleIndex(_JsonContract):
    trace: ArtifactReference
    dat: ArtifactReference
    raw_validation: ArtifactReference
    refined_validation: ArtifactReference
    kind: str = "experimental-atom-completion-bundle"
    scientific_qualification: bool = False
    public_profile_enabled: bool = False

    _schema_version: ClassVar[int] = 1

    def __post_init__(self) -> None:
        super(CompletionBundleIndex, self).__post_init__()
        if (type(self.schema_version) is not int or self.kind != "experimental-atom-completion-bundle"
                or self.scientific_qualification is not False or self.public_profile_enabled is not False):
            raise DiffusionContractError("completion bundle cannot claim public qualification")
        if len({a.path for a in (self.trace, self.dat, self.raw_validation, self.refined_validation)}) != 4:
            raise DiffusionContractError("completion bundle index artifact paths must be distinct")


def completion_dat_record(request: AtomCompletionRequest, source: bytes, candidate: bytes) -> DatRecord:
    """Only requested missing atoms enter downstream restraint ownership."""
    audit_completion_candidate(source, candidate, request)
    records, _, _ = _parse_coordinate_records(candidate.decode("utf-8").splitlines())
    by_atom = {record.identity: record for record in records}
    result = DatRecord(description="Experimental cc94 atom completion; generated-only refinement; not qualified")
    for atom in request.generated_atoms:
        record = by_atom[atom]
        result.added_atoms.append({"chain": atom.chain, "resid": atom.residue_number, "icode": atom.insertion_code,
                                   "resname": record.residue_name, "atom": atom.atom_name, "element": record.element})
        key = f"{atom.chain}/{record.residue_name}{atom.residue_number}{atom.insertion_code}"
        summary = result.residue_summary.setdefault(key, {"heavy": 0, "hydrogen": 0})
        summary["heavy"] += 1
    return result


def publish_completion_bundle(trace: CompletionSamplerTrace, workspace: Path, destination: Path) -> Path:
    """Audit before reservation, snapshot bounded evidence, commit readiness last.

    Existing/crash-left destinations are never adopted or removed. Failed staging
    is retained for diagnosis; original run artifacts are never mutated/deleted.
    """
    if destination.name in {"", ".", ".."} or destination.parent.is_symlink():
        raise ValueError("completion destination must have a regular parent and explicit name")
    destination = destination.parent.resolve(strict=True) / destination.name
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    audit = audit_completion_evidence(trace, workspace)
    snapshot = CompletionSamplerTrace.from_dict({**trace.to_dict(),
        **{role: asdict(ArtifactReference(_PAYLOAD_NAMES[role], artifact.sha256))
           for role, artifact in trace.artifacts().items()}})
    payloads = {_PAYLOAD_NAMES[role]: data for role, data in audit.payloads.items()}
    payloads["sampler-trace.json"] = snapshot.to_json().encode("utf-8")
    dat = completion_dat_record(audit.request, audit.payloads["source"], audit.payloads["refined"])
    payloads["candidate.dat"] = _json_bytes(dat.to_dict())
    payloads["raw-validation.json"] = _json_bytes(asdict(audit.raw_validation))
    payloads["refined-validation.json"] = _json_bytes(asdict(audit.refined_validation))
    index = CompletionBundleIndex(1, _reference("sampler-trace.json", payloads), _reference("candidate.dat", payloads),
                                  _reference("raw-validation.json", payloads), _reference("refined-validation.json", payloads))
    staging = Path(tempfile.mkdtemp(prefix=".completion-staging-", dir=destination.parent))
    for name, data in payloads.items():
        _durable_write(staging / name, data)
    _durable_write(staging / "bundle.json", index.to_json().encode("utf-8"))
    # Verify exactly the bytes about to be committed, not the live runner tree.
    read_completion_bundle(staging)
    commit_manifest_last(staging, destination)
    descriptor = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return destination


def read_completion_bundle(directory: Path) -> CompletionEvidenceAudit:
    """Readiness requires the committed index plus every independently valid payload."""
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("completion bundle must be a regular directory")
    root = directory.resolve(strict=True)
    index = CompletionBundleIndex.from_json(read_private_file_bytes(root, "bundle.json", max_bytes=64_000).decode("utf-8"))
    trace = CompletionSamplerTrace.from_json(read_private_artifact_bytes(root, index.trace, max_bytes=64_000).decode("utf-8"))
    paths = [a.path for a in trace.artifacts().values()] + [a.path for a in (
        index.trace, index.dat, index.raw_validation, index.refined_validation)] + ["bundle.json"]
    if (len(set(paths)) != len(paths) or any(Path(path).name != path for path in paths)
            or {p.name for p in root.iterdir()} != set(paths)):
        raise DiffusionContractError("completion bundle must contain exactly its distinct flat declared artifacts")
    audit = audit_completion_evidence(trace, root)
    dat = DatRecord.from_json(read_private_artifact_bytes(root, index.dat, max_bytes=2_000_000).decode("utf-8"))
    expected = completion_dat_record(audit.request, audit.payloads["source"], audit.payloads["refined"])
    if dat.to_dict() != expected.to_dict():
        raise DiffusionContractError("completion .dat does not match exact generated-only ownership")
    for artifact, validation in ((index.raw_validation, audit.raw_validation), (index.refined_validation, audit.refined_validation)):
        data = read_private_artifact_bytes(root, artifact, max_bytes=2_000_000)
        if json.loads(data) != json.loads(_json_bytes(asdict(validation))):
            raise DiffusionContractError("completion bundled validation differs from independent recomputation")
    return audit


def _json_bytes(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False).encode("utf-8")


def _reference(name: str, payloads: dict[str, bytes]) -> ArtifactReference:
    return ArtifactReference(name, hashlib.sha256(payloads[name]).hexdigest())


def _durable_write(path: Path, data: bytes) -> None:
    with path.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
