"""Explicit private VAL/THR source path with global THR-CB and rounded readback.

Default source/cohort/evidence entry points retain their original capability.
No RDKit, observed-coordinate repair or filesystem publication is introduced.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, replace
from typing import Literal

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.pdb_component_reconstruction import PdbAtom, _component_ref
from dvbfixer.protein_completion.branched import named_beta_volume
from dvbfixer.protein_completion.contracts import Addition, Candidate, CompletionRequest
from dvbfixer.protein_completion.source import (
    MaterializedCompletion,
    SourceCompletionError,
    _global_ca,
    _insert_candidates,
    _layout,
    _plan_source,
    _readback_additions,
)
from dvbfixer.protein_completion.validation import audit_candidate

BRANCHED_SOURCE_POLICY = "paired-val-thr-source-global-thr-beta-v1"
_TARGETS = frozenset({"VAL", "THR"})


@dataclass(frozen=True)
class ThrBetaCheck:
    residue: ComponentInstanceRef
    state: Literal["observed_R", "target_pending"]
    volume_angstrom3: float | None


@dataclass(frozen=True)
class BranchedSourcePlan:
    request: CompletionRequest
    source_thr_beta: tuple[ThrBetaCheck, ...]
    policy: str = BRANCHED_SOURCE_POLICY

    def __post_init__(self) -> None:
        if (not isinstance(self.request, CompletionRequest) or self.request.residue_name not in _TARGETS
                or self.policy != BRANCHED_SOURCE_POLICY or not isinstance(self.source_thr_beta, tuple)
                or any(not isinstance(c, ThrBetaCheck) for c in self.source_thr_beta)):
            raise SourceCompletionError("explicit immutable branched-source plan required")

    @property
    def digest(self) -> str:
        payload = {"policy": self.policy, "request_digest": self.request.digest,
                   "source_thr_beta": [asdict(c) for c in self.source_thr_beta]}
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                         allow_nan=False).encode()).hexdigest()


@dataclass(frozen=True)
class BranchedSourceAudit:
    plan_digest: str
    generated: tuple[Addition, ...]
    thr_beta: tuple[ThrBetaCheck, ...]


@dataclass(frozen=True)
class MaterializedBranchedSource:
    completion: MaterializedCompletion
    audit: BranchedSourceAudit

    @property
    def pdb(self) -> bytes:
        return self.completion.pdb

    def dat_record(self) -> DatRecord:
        return self.completion.dat_record()


def _thr_beta_checks(atoms: list[PdbAtom], *, pending_target: ComponentInstanceRef | None,
                      pending_targets: frozenset[ComponentInstanceRef] = frozenset()) -> tuple[ThrBetaCheck, ...]:
    residues: dict[ComponentInstanceRef, list[PdbAtom]] = {}
    for atom in atoms:
        if atom.residue_name == "THR":
            residues.setdefault(_component_ref(atom), []).append(atom)
    checks = []
    for residue, group in residues.items():
        xyz = {a.name: a.position for a in group}
        if not {"CA", "CB", "OG1", "CG2"} <= xyz.keys():
            if residue == pending_target or residue in pending_targets:
                checks.append(ThrBetaCheck(residue, "target_pending", None))
                continue  # The independent target/mask admission already ran.
            raise SourceCompletionError(f"undefined global THR CB: {residue}")
        volume = named_beta_volume(xyz["CA"], xyz["CB"], xyz["OG1"], xyz["CG2"])
        if volume <= .5:
            raise SourceCompletionError(f"non-R/nearplanar global THR CB: {residue}")
        checks.append(ThrBetaCheck(residue, "observed_R", volume))
    return tuple(checks)


def plan_branched_source(source: bytes, residue: ComponentInstanceRef) -> BranchedSourcePlan:
    """All source rules plus global observed THR-CB; only target may be incomplete."""
    request = _plan_source(source, residue, target_capability=_TARGETS)
    _, atoms, _ = _layout(source)
    checks = _thr_beta_checks(atoms, pending_target=residue)
    return BranchedSourcePlan(request, checks)


def _rebind(source: bytes, plan: BranchedSourcePlan) -> CompletionRequest:
    if not isinstance(plan, BranchedSourcePlan):
        raise SourceCompletionError("branched-source plan required")
    expected = plan_branched_source(source, plan.request.residue)
    if plan != expected or plan.digest != expected.digest:
        raise SourceCompletionError("branched source/context/stereo/policy binding mismatch")
    return expected.request


def audit_branched_materialized_source(source: bytes, plan: BranchedSourcePlan,
                                       candidate: Candidate, output: bytes) -> BranchedSourceAudit:
    request = _rebind(source, plan)
    audit = audit_candidate(request, candidate)
    if not audit.accepted:
        raise SourceCompletionError(f"raw branched candidate rejected: {audit.findings}")
    added = _readback_additions(source, output, candidate.additions)
    final = audit_candidate(request, replace(candidate, additions=added), pdb_roundtrip=True)
    if not final.accepted:
        raise SourceCompletionError(f"rounded branched PDB rejected: {final.findings}")
    _, atoms, _ = _layout(output)
    _global_ca(atoms, allow_missing_cb=None)
    checks = _thr_beta_checks(atoms, pending_target=None)
    return BranchedSourceAudit(plan.digest, added, checks)


def materialize_branched_source(source: bytes, plan: BranchedSourcePlan,
                                candidate: Candidate) -> MaterializedBranchedSource:
    """One in-memory insert; all original bytes unchanged and no incomplete beta centre."""
    request = _rebind(source, plan)
    audit = audit_candidate(request, candidate)
    if not audit.accepted:
        raise SourceCompletionError(f"raw branched candidate rejected: {audit.findings}")
    output = _insert_candidates(source, (request,), (candidate,))
    readback = audit_branched_materialized_source(source, plan, candidate, output)
    completion = MaterializedCompletion(output, request.source_sha256, hashlib.sha256(output).hexdigest(),
                                         request.digest, readback.generated)
    return MaterializedBranchedSource(completion, readback)
