"""Explicit private singleton VAL/THR source policy; paired policies are unchanged.

The observed distal partner, full source and global THR-CB checks are digest-bound.
No file publishing, empirical library or observed-coordinate repair is introduced.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, replace

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.protein_completion.anchored_branch import audit_anchored_branch
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.branched_source import (
    BranchedSourceAudit,
    ThrBetaCheck,
    _thr_beta_checks,
)
from dvbfixer.protein_completion.contracts import Candidate, CompletionRequest, Observation
from dvbfixer.protein_completion.source import (
    MaterializedCompletion,
    SourceCompletionError,
    _global_ca,
    _insert_candidates,
    _layout,
    _plan_source,
    _readback_additions,
)
from dvbfixer.protein_completion.validation import _admission_reason

ANCHORED_SOURCE_POLICY = "singleton-val-thr-source-fixed-partner-global-thr-beta-v1"
_TARGETS = frozenset({"VAL", "THR"})


def _anchored_admission(request: CompletionRequest) -> str | None:
    return _admission_reason(request, anchored_branch=True)


@dataclass(frozen=True)
class AnchoredSourcePlan:
    request: CompletionRequest
    fixed_partner: Observation
    source_thr_beta: tuple[ThrBetaCheck, ...]
    policy: str = ANCHORED_SOURCE_POLICY

    def __post_init__(self) -> None:
        if (not isinstance(self.request, CompletionRequest) or self.request.residue_name not in _TARGETS
                or _anchored_admission(self.request) is not None or not isinstance(self.fixed_partner, Observation)
                or self.policy != ANCHORED_SOURCE_POLICY or not isinstance(self.source_thr_beta, tuple)
                or any(not isinstance(c, ThrBetaCheck) for c in self.source_thr_beta)):
            raise SourceCompletionError("explicit immutable anchored-source plan required")
        branch = BRANCHES[self.request.residue_name]
        partner = branch.second if self.request.missing_atoms[0] == branch.first else branch.first
        if (self.fixed_partner.atom.component != self.request.residue or self.fixed_partner.atom.name != partner
                or self.fixed_partner not in self.request.observations):
            raise SourceCompletionError("exact observed distal partner required")

    @property
    def digest(self) -> str:
        payload = {"policy": self.policy, "request_digest": self.request.digest,
                   "fixed_partner": asdict(self.fixed_partner),
                   "source_thr_beta": [asdict(c) for c in self.source_thr_beta]}
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                         allow_nan=False).encode()).hexdigest()


@dataclass(frozen=True)
class MaterializedAnchoredSource:
    completion: MaterializedCompletion
    audit: BranchedSourceAudit

    @property
    def pdb(self) -> bytes:
        return self.completion.pdb

    def dat_record(self) -> DatRecord:
        return self.completion.dat_record()


def plan_anchored_source(source: bytes, residue: ComponentInstanceRef) -> AnchoredSourcePlan:
    request = _plan_source(source, residue, target_capability=_TARGETS, validate_admission=_anchored_admission)
    _, atoms, _ = _layout(source)
    checks = _thr_beta_checks(atoms, pending_target=residue)
    branch = BRANCHES[request.residue_name]
    partner_name = branch.second if request.missing_atoms[0] == branch.first else branch.first
    partner = next(o for o in request.observations if o.atom.component == residue and o.atom.name == partner_name)
    return AnchoredSourcePlan(request, partner, checks)


def _rebind(source: bytes, plan: AnchoredSourcePlan) -> CompletionRequest:
    if not isinstance(plan, AnchoredSourcePlan):
        raise SourceCompletionError("anchored-source plan required; paired plans cannot opt in")
    expected = plan_anchored_source(source, plan.request.residue)
    if plan != expected or plan.digest != expected.digest:
        raise SourceCompletionError("anchored source/partner/context/stereo/policy binding mismatch")
    return expected.request


def audit_anchored_materialized_source(source: bytes, plan: AnchoredSourcePlan,
                                      candidate: Candidate, output: bytes) -> BranchedSourceAudit:
    request = _rebind(source, plan)
    raw = audit_anchored_branch(request, candidate)
    if not raw.accepted:
        raise SourceCompletionError(f"raw anchored candidate rejected: {raw.findings}")
    added = _readback_additions(source, output, candidate.additions)
    rounded = audit_anchored_branch(request, replace(candidate, additions=added), pdb_roundtrip=True)
    if not rounded.accepted:
        raise SourceCompletionError(f"rounded anchored PDB rejected: {rounded.findings}")
    _, atoms, _ = _layout(output)
    _global_ca(atoms, allow_missing_cb=None)
    checks = _thr_beta_checks(atoms, pending_target=None)
    return BranchedSourceAudit(plan.digest, added, checks)


def materialize_anchored_source(source: bytes, plan: AnchoredSourcePlan,
                               candidate: Candidate) -> MaterializedAnchoredSource:
    request = _rebind(source, plan)
    audit = audit_anchored_branch(request, candidate)
    if not audit.accepted:
        raise SourceCompletionError(f"raw anchored candidate rejected: {audit.findings}")
    output = _insert_candidates(source, (request,), (candidate,))
    checked = audit_anchored_materialized_source(source, plan, candidate, output)
    completion = MaterializedCompletion(output, request.source_sha256, hashlib.sha256(output).hexdigest(),
                                         request.digest, checked.generated)
    return MaterializedAnchoredSource(completion, checked)
