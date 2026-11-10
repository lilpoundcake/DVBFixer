"""Bounded private multi-target staggered controls with joint clash/readback gates.

Not coupled energetic packing, empirical rotamer qualification or a public API.
All targets share unchanged source bytes; generated atoms enter every other
target's independent audit before any accepted output is returned.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections.abc import Callable
from dataclasses import dataclass, replace
from typing import Literal

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import AddedAtom, DatRecord
from dvbfixer.pdb_component_reconstruction import _component_ref
from dvbfixer.protein_completion.contracts import (
    MAX_CONTEXT_ATOMS,
    Addition,
    Audit,
    Candidate,
    CompletionRequest,
    CompletionResult,
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
from dvbfixer.protein_completion.validation import audit_candidate

MAX_TARGETS = 4
MAX_COMBINATIONS = 81
COHORT_REVISION = "canonical-sidechain-cohort-staggered-v1"


@dataclass(frozen=True)
class CohortPlan:
    requests: tuple[CompletionRequest, ...]  # canonical source encounter order

    def __post_init__(self) -> None:
        if (not isinstance(self.requests, tuple) or not 1 <= len(self.requests) <= MAX_TARGETS
                or any(not isinstance(r, CompletionRequest) for r in self.requests)):
            raise SourceCompletionError("bounded immutable completion requests required")
        if (len({r.residue for r in self.requests}) != len(self.requests)
                or len({r.source_sha256 for r in self.requests}) != 1):
            raise SourceCompletionError("cohort requires unique targets and one source SHA")

    @property
    def digest(self) -> str:
        payload = {"revision": COHORT_REVISION, "requests": [r.digest for r in self.requests]}
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@dataclass(frozen=True)
class JointAttempt:
    indices: tuple[int, ...]  # indices into each retained single-target attempt list
    audits: tuple[Audit, ...]

    @property
    def accepted(self) -> bool:
        return bool(self.audits) and all(a.accepted for a in self.audits)


@dataclass(frozen=True)
class CohortResult:
    status: Literal["completed", "rejected"]
    plan_digest: str
    reason: str
    targets: tuple[CompletionResult, ...]
    attempts: tuple[JointAttempt, ...] = ()
    selected_index: int | None = None

    @property
    def selected(self) -> tuple[Candidate, ...] | None:
        if self.selected_index is None:
            return None
        indices = self.attempts[self.selected_index].indices
        return tuple(r.attempts[i].candidate for r, i in zip(self.targets, indices, strict=True))


@dataclass(frozen=True)
class MaterializedCohort:
    pdb: bytes
    source_sha256: str
    output_sha256: str
    plan_digest: str
    generated: tuple[Addition, ...]
    residue_names: tuple[tuple[ComponentInstanceRef, str], ...]

    def dat_record(self) -> DatRecord:
        names = dict(self.residue_names)
        added: list[AddedAtom] = [
            {"chain": a.atom.component.chain_id, "resid": a.atom.component.sequence_number,
             "icode": a.atom.component.insertion_code, "resname": names[a.atom.component],
             "atom": a.atom.name, "element": a.element} for a in self.generated
        ]
        return DatRecord(description="Private joint staggered completion; generated atoms only", added_atoms=added)


def plan_cohort(source: bytes, targets: tuple[ComponentInstanceRef, ...]) -> CohortPlan:
    if (not isinstance(targets, tuple) or not 1 <= len(targets) <= MAX_TARGETS
            or any(not isinstance(t, ComponentInstanceRef) for t in targets) or len(set(targets)) != len(targets)):
        raise SourceCompletionError("one to four distinct immutable exact targets required")
    _, atoms, _ = _layout(source)
    order = list(dict.fromkeys(_component_ref(a) for a in atoms))
    if any(t not in order for t in targets):
        raise SourceCompletionError("exact cohort target absent")
    chosen = frozenset(targets)
    requests = tuple(_plan_source(source, target, other_missing_cb=chosen-{target})
                     for target in order if target in chosen)
    if any(not r.missing_atoms for r in requests):
        raise SourceCompletionError("cohort targets must each have missing atoms")
    total_new = sum(len(r.missing_atoms) for r in requests)
    if any(len(r.observations)+total_new-len(r.missing_atoms) > MAX_CONTEXT_ATOMS for r in requests):
        raise SourceCompletionError("joint context exceeds bounded atom budget")
    return CohortPlan(requests)


def _rebind(source: bytes, plan: CohortPlan) -> None:
    if not isinstance(plan, CohortPlan) or not isinstance(plan.requests, tuple):
        raise SourceCompletionError("immutable source-bound cohort plan required")
    if plan_cohort(source, tuple(r.residue for r in plan.requests)) != plan:
        raise SourceCompletionError("cohort source/context/order binding mismatch")


def _joint_audits(plan: CohortPlan, candidates: tuple[Candidate, ...], *,
                  rounded: bool = False) -> tuple[Audit, ...]:
    if not isinstance(candidates, tuple) or len(candidates) != len(plan.requests):
        raise SourceCompletionError("exact immutable candidate cohort required")
    audits = []
    for i, (request, candidate) in enumerate(zip(plan.requests, candidates, strict=True)):
        # Source-only candidate binding is checked before this contextual audit.
        # The ephemeral digest explicitly includes the OTHER new atoms; it is
        # never substituted for the original source-bound candidate in evidence.
        others = tuple(Observation(a.atom, a.element, a.xyz)
                       for j, c in enumerate(candidates) if j != i for a in c.additions)
        contextual = replace(request, observations=request.observations+others)
        bound = replace(candidate, request_digest=contextual.digest)
        audits.append(audit_candidate(contextual, bound, pdb_roundtrip=rounded))
    return tuple(audits)


def _audit_raw(plan: CohortPlan, candidates: tuple[Candidate, ...]) -> None:
    if not isinstance(candidates, tuple) or len(candidates) != len(plan.requests):
        raise SourceCompletionError("exact immutable candidate cohort required")
    for request, candidate in zip(plan.requests, candidates, strict=True):
        audit = audit_candidate(request, candidate)
        if not audit.accepted:
            raise SourceCompletionError(f"raw cohort candidate rejected: {audit.findings}")
    if not all(a.accepted for a in _joint_audits(plan, candidates)):
        raise SourceCompletionError("joint generated geometry rejected")


def complete_cohort(source: bytes, plan: CohortPlan) -> CohortResult:
    """Enumerate <=81 source-only-feasible combinations; retain all joint audits.

    No empirical library parameter exists in this control-only cohort slice.
    Selection maximizes capped worst clearance; stable enumeration breaks ties.
    """
    _rebind(source, plan)
    results = tuple(complete_sidechain(r) for r in plan.requests)
    return _enumerate_joint(results, lambda candidates: _joint_audits(plan, candidates), plan.digest,
                            max_combinations=MAX_COMBINATIONS)


def _enumerate_joint(results: tuple[CompletionResult, ...],
                      audit: Callable[[tuple[Candidate, ...]], tuple[Audit, ...]],
                      plan_digest: str, *, max_combinations: int,
                      completed_reason: str = "joint_geometry_control") -> CohortResult:
    """Shared bounded control search; callers own admission, modes and source audit."""
    pools = [tuple(i for i, a in enumerate(r.attempts) if a.audit.accepted) for r in results]
    if any(not pool for pool in pools):
        return CohortResult("rejected", plan_digest, "individual_target_failed", results)
    count = 1
    for pool in pools:
        count *= len(pool)
    if count > max_combinations:
        return CohortResult("rejected", plan_digest, "combination_budget_exceeded", results)
    attempts = []
    for indices in itertools.product(*pools):
        candidates = tuple(r.attempts[i].candidate for r, i in zip(results, indices, strict=True))
        attempts.append(JointAttempt(indices, audit(candidates)))
    passing = [i for i, a in enumerate(attempts) if a.accepted]

    def score(index: int) -> float:
        return min(min(0., a.minimum_nonbonded_clearance_angstrom or 0.) for a in attempts[index].audits)

    selected = max(passing, key=score) if passing else None
    return CohortResult("completed" if passing else "rejected", plan_digest,
                        completed_reason if passing else "all_joint_candidates_failed",
                        results, tuple(attempts), selected)


def audit_materialized_cohort(source: bytes, plan: CohortPlan, candidates: tuple[Candidate, ...],
                              output: bytes) -> tuple[Addition, ...]:
    _rebind(source, plan)
    _audit_raw(plan, candidates)
    added = _readback_additions(source, output, tuple(a for c in candidates for a in c.additions))
    rounded = tuple(replace(c, additions=tuple(a for a in added if a.atom.component == r.residue))
                    for r, c in zip(plan.requests, candidates, strict=True))
    audits = _joint_audits(plan, rounded, rounded=True)
    if not all(a.accepted for a in audits):
        raise SourceCompletionError(f"rounded joint PDB rejected: {audits}")
    _, atoms, _ = _layout(output)
    _global_ca(atoms, allow_missing_cb=None)
    return added


def materialize_cohort(source: bytes, plan: CohortPlan,
                       candidates: tuple[Candidate, ...]) -> MaterializedCohort:
    """All-or-nothing in-memory insertion/readback; no partial accepted output."""
    _rebind(source, plan)
    _audit_raw(plan, candidates)
    output = _insert_candidates(source, plan.requests, candidates)
    added = audit_materialized_cohort(source, plan, candidates, output)
    return MaterializedCohort(output, hashlib.sha256(source).hexdigest(), hashlib.sha256(output).hexdigest(),
                              plan.digest, added, tuple((r.residue, r.residue_name) for r in plan.requests))
