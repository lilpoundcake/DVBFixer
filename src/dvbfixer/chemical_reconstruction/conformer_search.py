"""Bounded nonempirical fragment search; retained-context ranking without references."""
from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product
from typing import Literal

from dvbfixer.domain.structure_identity import ExactAtomRef

from .candidate_audit import audit_candidate
from .construction import propose
from .contracts import ReconstructionRequest, _number, _sha, digest_record
from .peptide_state import LocalPeptideAuthority
from .proposals import Candidate, CandidateAudit, GeometryPlan, ProposalResult, Strategy
from .steric_ranking import (
    StericContext,
    StericScore,
    compile_steric_context,
    radii_digest,
    score_sterics,
)
from .torsion_geometry import (
    ConformerError,
    TorsionTopology,
    TorsionTraceAudit,
    audit_torsion_offsets,
    compile_fragment_torsions,
    rotate_fragments,
)


@dataclass(frozen=True, kw_only=True)
class ConformerSearchPlan(GeometryPlan):
    base_plan: GeometryPlan
    base_candidate: Candidate
    authority_digest: str
    topology: TorsionTopology
    radii_digest: str
    grid_offsets: tuple[float, ...] = (0., -120., 120.)
    max_combinations: int = 256

    def __post_init__(self) -> None:
        super().__post_init__()
        _sha(self.authority_digest)
        _sha(self.radii_digest)
        if (self.strategy != Strategy.TORSION or type(self.base_plan) is not GeometryPlan
                or self.base_plan.strategy != Strategy.INTERNAL or self.base_plan.request_digest != self.request_digest
                or self.base_candidate.request_digest != self.request_digest or self.base_candidate.plan_digest != self.base_plan.digest
                or self.base_candidate.strategy != Strategy.INTERNAL or self.placements != self.base_plan.placements
                or self.policy != self.base_plan.policy or not isinstance(self.topology, TorsionTopology)
                or not isinstance(self.grid_offsets, tuple) or any(not _number(v) for v in self.grid_offsets)
                or self.grid_offsets != (0., -120., 120.) or type(self.max_combinations) is not int or self.max_combinations != 256):
            raise ConformerError("conformer-plan-binding-or-policy")
        if len(self.grid_offsets)**len(self.topology.rotors) > self.max_combinations:
            raise ConformerError("conformer-grid-budget-exceeded")

    @property
    def digest(self) -> str:
        return digest_record(("chemical-conformer-search-plan-v1", self))


@dataclass(frozen=True)
class ConformerAttempt:
    index: int
    offsets: tuple[float, ...]
    reason: str
    raw: Candidate | None = None
    raw_audit: CandidateAudit | None = None
    raw_trace: TorsionTraceAudit | None = None
    raw_sterics: StericScore | None = None
    rounded: Candidate | None = None
    rounded_audit: CandidateAudit | None = None
    rounded_trace: TorsionTraceAudit | None = None
    rounded_sterics: StericScore | None = None

    @property
    def selectable(self) -> bool:
        return bool(self.raw is not None and self.rounded is not None
                    and self.raw_audit and self.raw_audit.passed and self.rounded_audit and self.rounded_audit.passed
                    and self.raw_trace and self.raw_trace.passed and self.rounded_trace and self.rounded_trace.passed
                    and self.raw_sterics is not None and self.rounded_sterics is not None)


@dataclass(frozen=True)
class ConformerSearchResult:
    status: str
    reason: str
    baseline: ProposalResult
    plan: ConformerSearchPlan | None = None
    attempts: tuple[ConformerAttempt, ...] = ()
    selected_index: int | None = None

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def native_conformation_qualified(self) -> Literal[False]:
        return False


def _rounded(candidate: Candidate) -> Candidate:
    additions = []
    for o in candidate.additions:
        fields = tuple(f"{x:8.3f}" for x in o.position)
        if any(len(field) != 8 for field in fields):
            raise ConformerError("conformer-PDB-coordinate-range")
        additions.append(replace(o, position=(float(fields[0]), float(fields[1]), float(fields[2]))))
    return replace(candidate, additions=tuple(additions))


def _attempt_measurements(request: ReconstructionRequest, plan: ConformerSearchPlan, candidate: Candidate,
                           offsets: tuple[float, ...], context: StericContext, *, rounded: bool
                           ) -> tuple[CandidateAudit, TorsionTraceAudit, StericScore]:
    return (audit_candidate(request, plan, candidate),
            audit_torsion_offsets(request, plan.base_candidate, candidate, plan.topology, offsets, rounded=rounded),
            score_sterics(request, plan, candidate, context))


def selected_attempt_index(attempts: tuple[ConformerAttempt, ...]) -> int | None:
    selectable = [a for a in attempts if a.selectable]
    if not selectable:
        return None

    def key(a: ConformerAttempt) -> tuple[float, float, int]:
        assert a.rounded_sterics is not None
        return a.rounded_sterics.sum_squared_overlap_angstrom2, a.rounded_sterics.worst_overlap_angstrom, a.index

    return min(selectable, key=key).index


def search_conformers(request: ReconstructionRequest, base_plan: GeometryPlan, authority: LocalPeptideAuthority,
                       roots: tuple[ExactAtomRef, ...], *, expected_baseline_signature: str | None = None) -> ConformerSearchResult:
    """Enumerate the frozen grid; native-reference input is deliberately absent."""
    baseline = propose(request, base_plan)
    if baseline.status != "declared-geometry-passed" or baseline.raw is None or baseline.rounded is None:
        return ConformerSearchResult("base-refused", "baseline-raw-rounded-geometry-not-passed", baseline)
    if expected_baseline_signature is not None:
        _sha(expected_baseline_signature)
        if digest_record(baseline.rounded.additions) != expected_baseline_signature:
            return ConformerSearchResult("baseline-refused", "frozen-baseline-coordinate-signature-mismatch", baseline)
    try:
        topology = compile_fragment_torsions(request, authority, roots)
        context = compile_steric_context(request)
        plan = ConformerSearchPlan(request.digest, Strategy.TORSION, base_plan.policy, base_plan.placements,
                                  base_plan=base_plan, base_candidate=baseline.raw, authority_digest=authority.digest,
                                  topology=topology, radii_digest=context.radii_digest)
    except ValueError as exc:
        return ConformerSearchResult("unsupported", str(exc), baseline)
    attempts = []
    for index, offsets in enumerate(product(plan.grid_offsets, repeat=len(plan.topology.rotors))):
        raw = rounded_candidate = None
        raw_audit = raw_trace = raw_sterics = rounded_audit = rounded_trace = rounded_sterics = None
        try:
            additions = rotate_fragments(request, authority, plan.base_candidate, plan.topology, offsets)
            raw = Candidate(request.digest, plan.digest, Strategy.TORSION, additions)
            raw_audit, raw_trace, raw_sterics = _attempt_measurements(request, plan, raw, offsets, context, rounded=False)
            rounded_candidate = _rounded(raw)
            rounded_audit, rounded_trace, rounded_sterics = _attempt_measurements(request, plan, rounded_candidate, offsets, context, rounded=True)
            reason = "measured-grid-attempt"
        except ValueError as exc:
            reason = f"construction-or-measurement-failed:{exc}"
        attempts.append(ConformerAttempt(index, offsets, reason, raw, raw_audit, raw_trace, raw_sterics,
                                         rounded_candidate, rounded_audit, rounded_trace, rounded_sterics))
    frozen = tuple(attempts)
    selected = selected_attempt_index(frozen)
    return ConformerSearchResult("selected" if selected is not None else "rejected", "steric-grid-not-native-conformation",
                                 baseline, plan, frozen, selected)


def audit_search_result(request: ReconstructionRequest, authority: LocalPeptideAuthority,
                         result: ConformerSearchResult) -> tuple[str, ...]:
    """Readback/selection audit without proposing or rotating ANY coordinates."""
    plan = result.plan
    if plan is None:
        return ("conformer-no-search-plan-to-audit",)
    errors = set()
    if (plan.request_digest != request.digest or plan.authority_digest != authority.digest
            or plan.radii_digest != radii_digest() or result.baseline.raw != plan.base_candidate
            or result.baseline.status != "declared-geometry-passed"
            or result.baseline.rounded != _rounded(plan.base_candidate)):
        errors.add("conformer-search-binding")
    if (not audit_candidate(request, plan.base_plan, plan.base_candidate).passed
            or result.baseline.rounded is None or not audit_candidate(request, plan.base_plan, result.baseline.rounded).passed):
        errors.add("conformer-invalid-baseline")
    if compile_fragment_torsions(request, authority, plan.topology.roots) != plan.topology:
        errors.add("conformer-topology-binding")
    context = compile_steric_context(request)
    expected = tuple(product(plan.grid_offsets, repeat=len(plan.topology.rotors)))
    if len(result.attempts) != len(expected):
        errors.add("conformer-incomplete-grid-records")
    for index, (attempt, offsets) in enumerate(zip(result.attempts, expected)):
        if attempt.index != index or attempt.offsets != offsets:
            errors.add("conformer-grid-index-or-offset-binding")
        if attempt.raw is None:
            if attempt.reason == "measured-grid-attempt" or attempt.selectable:
                errors.add("conformer-unaccounted-construction-failure")
            continue
        raw = _attempt_measurements(request, plan, attempt.raw, offsets, context, rounded=False)
        if raw != (attempt.raw_audit, attempt.raw_trace, attempt.raw_sterics):
            errors.add("conformer-raw-audit-or-score-mismatch")
        if attempt.rounded is None:
            if attempt.reason == "measured-grid-attempt":
                errors.add("conformer-missing-rounded-record")
            continue
        if attempt.rounded != _rounded(attempt.raw):
            errors.add("conformer-rounded-transport")
        rounded_values = _attempt_measurements(request, plan, attempt.rounded, offsets, context, rounded=True)
        if rounded_values != (attempt.rounded_audit, attempt.rounded_trace, attempt.rounded_sterics):
            errors.add("conformer-rounded-audit-or-score-mismatch")
    selected = selected_attempt_index(result.attempts)
    if result.selected_index != selected or result.status != ("selected" if selected is not None else "rejected"):
        errors.add("conformer-selection-mismatch")
    return tuple(sorted(errors))
