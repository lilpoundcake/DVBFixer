"""New private mixed-mask cohort; old source/cohort/evidence contracts are unchanged.

One fixed source, explicit per-target mode, <=81 joint combinations and complete
global THR-CB readback. No empirical data, publishing or production routing.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, replace
from typing import Literal

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.pdb_component_reconstruction import _component_ref
from dvbfixer.protein_completion.anchored_branch import (
    audit_anchored_branch,
    complete_anchored_branch,
)
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.branched_source import ThrBetaCheck, _thr_beta_checks
from dvbfixer.protein_completion.cohort import (
    MAX_COMBINATIONS,
    MAX_TARGETS,
    CohortResult,
    MaterializedCohort,
    _enumerate_joint,
)
from dvbfixer.protein_completion.contracts import (
    MAX_CONTEXT_ATOMS,
    Addition,
    Audit,
    Candidate,
    CompletionRequest,
    Observation,
)
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import (
    SourceCompletionError,
    _global_ca,
    _insert_candidates,
    _layout,
    _plan_source,
    _readback_additions,
)
from dvbfixer.protein_completion.validation import _admission_reason, audit_candidate

MIXED_POLICY = "mixed-canonical-paired-anchored-source-global-thr-beta-v1"
Mode = Literal["standard", "paired", "anchored"]
_STANDARD = frozenset({"ALA", "SER", "CYS"})
_CAPABILITY = _STANDARD | frozenset({"VAL", "THR"})


def _mode(request: CompletionRequest) -> Mode:
    if request.residue_name in _STANDARD:
        return "standard"
    if request.residue_name not in _CAPABILITY:
        raise SourceCompletionError("unsupported mixed target")
    branch = BRANCHES.get(request.residue_name)
    if branch is None:
        raise SourceCompletionError("unsupported mixed target")
    if len(request.missing_atoms) == 1 and request.missing_atoms[0] in {branch.first, branch.second}:
        return "anchored"
    return "paired"


def _admission(request: CompletionRequest) -> str | None:
    if not request.missing_atoms:
        return "mixed_targets_require_missing_atoms"
    return _admission_reason(request, anchored_branch=_mode(request) == "anchored")


@dataclass(frozen=True)
class MixedTarget:
    request: CompletionRequest
    mode: Mode
    fixed_partner: Observation | None = None

    def __post_init__(self) -> None:
        if (not isinstance(self.request, CompletionRequest) or self.mode != _mode(self.request)
                or _admission(self.request) is not None):
            raise SourceCompletionError("source-derived mixed mode and valid mask required")
        if self.mode == "anchored":
            branch = BRANCHES[self.request.residue_name]
            partner = branch.second if self.request.missing_atoms[0] == branch.first else branch.first
            if (not isinstance(self.fixed_partner, Observation)
                    or self.fixed_partner.atom.component != self.request.residue
                    or self.fixed_partner.atom.name != partner
                    or self.fixed_partner not in self.request.observations):
                raise SourceCompletionError("mixed anchored mode requires its exact observed partner")
        elif self.fixed_partner is not None:
            raise SourceCompletionError("non-anchored mixed mode cannot claim a partner")


@dataclass(frozen=True)
class MixedCohortPlan:
    targets: tuple[MixedTarget, ...]
    source_thr_beta: tuple[ThrBetaCheck, ...]
    policy: str = MIXED_POLICY

    def __post_init__(self) -> None:
        if (not isinstance(self.targets, tuple) or not 1 <= len(self.targets) <= MAX_TARGETS
                or any(not isinstance(t, MixedTarget) for t in self.targets) or self.policy != MIXED_POLICY
                or not isinstance(self.source_thr_beta, tuple)
                or any(not isinstance(c, ThrBetaCheck) for c in self.source_thr_beta)):
            raise SourceCompletionError("bounded immutable mixed plan and exact policy required")
        if (len({t.request.residue for t in self.targets}) != len(self.targets)
                or len({t.request.source_sha256 for t in self.targets}) != 1):
            raise SourceCompletionError("mixed targets must be unique and share one source SHA")

    @property
    def digest(self) -> str:
        payload = {"policy": self.policy, "targets": [
            {"request_digest": t.request.digest, "mode": t.mode,
             "fixed_partner": asdict(t.fixed_partner) if t.fixed_partner is not None else None} for t in self.targets],
            "source_thr_beta": [asdict(c) for c in self.source_thr_beta]}
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


@dataclass(frozen=True)
class MixedCohortAudit:
    plan_digest: str
    generated: tuple[Addition, ...]
    targets: tuple[Audit, ...]
    thr_beta: tuple[ThrBetaCheck, ...]


@dataclass(frozen=True)
class MaterializedMixedCohort:
    completion: MaterializedCohort
    audit: MixedCohortAudit

    @property
    def pdb(self) -> bytes:
        return self.completion.pdb

    def dat_record(self) -> DatRecord:
        record = self.completion.dat_record()
        record.description = "Private mixed sidechain completion; generated atoms only"
        return record


def plan_mixed_cohort(source: bytes, targets: tuple[ComponentInstanceRef, ...]) -> MixedCohortPlan:
    if (not isinstance(targets, tuple) or not 1 <= len(targets) <= MAX_TARGETS
            or any(not isinstance(t, ComponentInstanceRef) for t in targets) or len(set(targets)) != len(targets)):
        raise SourceCompletionError("one to four distinct immutable mixed targets required")
    _, atoms, _ = _layout(source)
    order = list(dict.fromkeys(_component_ref(a) for a in atoms))
    if any(t not in order for t in targets):
        raise SourceCompletionError("exact mixed target absent")
    chosen = frozenset(targets)
    planned = []
    for target in order:
        if target not in chosen:
            continue
        request = _plan_source(source, target, other_missing_cb=chosen-{target},
                               target_capability=_CAPABILITY, validate_admission=_admission)
        mode = _mode(request)
        partner = None
        if mode == "anchored":
            branch = BRANCHES[request.residue_name]
            name = branch.second if request.missing_atoms[0] == branch.first else branch.first
            partner = next(o for o in request.observations if o.atom.component == target and o.atom.name == name)
        planned.append(MixedTarget(request, mode, partner))
    total = sum(len(t.request.missing_atoms) for t in planned)
    if any(len(t.request.observations)+total-len(t.request.missing_atoms) > MAX_CONTEXT_ATOMS for t in planned):
        raise SourceCompletionError("mixed joint context exceeds bounded atom budget")
    checks = _thr_beta_checks(atoms, pending_target=None, pending_targets=chosen)
    return MixedCohortPlan(tuple(planned), checks)


def _rebind(source: bytes, plan: MixedCohortPlan) -> None:
    if not isinstance(plan, MixedCohortPlan):
        raise SourceCompletionError("mixed plan required; older plans cannot opt in")
    expected = plan_mixed_cohort(source, tuple(t.request.residue for t in plan.targets))
    if plan != expected or plan.digest != expected.digest:
        raise SourceCompletionError("mixed source/mode/partner/context/stereo binding mismatch")


def _audit(target: MixedTarget, candidate: Candidate, *, rounded: bool = False) -> Audit:
    if target.mode == "anchored":
        return audit_anchored_branch(target.request, candidate, pdb_roundtrip=rounded)
    return audit_candidate(target.request, candidate, pdb_roundtrip=rounded)


def _joint(plan: MixedCohortPlan, candidates: tuple[Candidate, ...], *, rounded: bool = False) -> tuple[Audit, ...]:
    if not isinstance(candidates, tuple) or len(candidates) != len(plan.targets):
        raise SourceCompletionError("exact immutable mixed candidate tuple required")
    audits = []
    for i, (target, candidate) in enumerate(zip(plan.targets, candidates, strict=True)):
        others = tuple(Observation(a.atom, a.element, a.xyz) for j, c in enumerate(candidates) if j != i for a in c.additions)
        request = replace(target.request, observations=target.request.observations+others)
        contextual = replace(target, request=request)
        audits.append(_audit(contextual, replace(candidate, request_digest=request.digest), rounded=rounded))
    return tuple(audits)


def _raw(plan: MixedCohortPlan, candidates: tuple[Candidate, ...]) -> None:
    if not isinstance(candidates, tuple) or len(candidates) != len(plan.targets):
        raise SourceCompletionError("exact immutable mixed candidate tuple required")
    if not all(_audit(t, c).accepted for t, c in zip(plan.targets, candidates, strict=True)):
        raise SourceCompletionError("raw mixed candidate rejected")
    if not all(a.accepted for a in _joint(plan, candidates)):
        raise SourceCompletionError("joint mixed generated geometry rejected")


def complete_mixed_cohort(source: bytes, plan: MixedCohortPlan) -> CohortResult:
    _rebind(source, plan)
    results = tuple(complete_anchored_branch(t.request) if t.mode == "anchored" else complete_sidechain(t.request)
                    for t in plan.targets)
    return _enumerate_joint(results, lambda candidates: _joint(plan, candidates), plan.digest,
                            max_combinations=MAX_COMBINATIONS, completed_reason="mixed_joint_geometry_control")


def audit_materialized_mixed_cohort(source: bytes, plan: MixedCohortPlan,
                                   candidates: tuple[Candidate, ...], output: bytes) -> MixedCohortAudit:
    _rebind(source, plan)
    _raw(plan, candidates)
    added = _readback_additions(source, output, tuple(a for c in candidates for a in c.additions))
    rounded = tuple(replace(c, additions=tuple(a for a in added if a.atom.component == t.request.residue))
                    for t, c in zip(plan.targets, candidates, strict=True))
    audits = _joint(plan, rounded, rounded=True)
    if not all(a.accepted for a in audits):
        raise SourceCompletionError(f"rounded mixed PDB rejected: {audits}")
    _, atoms, _ = _layout(output)
    _global_ca(atoms, allow_missing_cb=None)
    checks = _thr_beta_checks(atoms, pending_target=None)
    return MixedCohortAudit(plan.digest, added, audits, checks)


def materialize_mixed_cohort(source: bytes, plan: MixedCohortPlan,
                             candidates: tuple[Candidate, ...]) -> MaterializedMixedCohort:
    _rebind(source, plan)
    _raw(plan, candidates)
    requests = tuple(t.request for t in plan.targets)
    output = _insert_candidates(source, requests, candidates)
    checked = audit_materialized_mixed_cohort(source, plan, candidates, output)
    completion = MaterializedCohort(output, hashlib.sha256(source).hexdigest(), hashlib.sha256(output).hexdigest(),
                                    plan.digest, checked.generated, tuple((r.residue, r.residue_name) for r in requests))
    return MaterializedMixedCohort(completion, checked)
