"""Source-authenticated saved numerical readback and independent branch selection."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .contracts import Unit, _sha, _tuple, digest_record
from .cycle_audit import CycleGeometryAudit, CyclePairAudit, measure_cycle_pair
from .cycle_residuals import CycleNumericalError, residual_values, squared_residual, vector_points
from .cycle_solve_plan import CycleSolvePlan, compile_cycle_solve_plan
from .cycle_solver import CycleBranchResult, CycleConstructionResult, cycle_environment


@dataclass(frozen=True)
class CycleBranchAudit:
    index: int
    initial: CycleGeometryAudit
    pair: CyclePairAudit | None
    selectable: bool
    reason: str


@dataclass(frozen=True)
class CycleConstructionAudit:
    plan_digest: str
    result_digest: str
    branches: tuple[CycleBranchAudit, ...]
    selected_index: int | None
    status: Literal["selected-declared-geometry-only", "rejected", "unavailable", "aborted"]

    def __post_init__(self) -> None:
        _sha(self.plan_digest)
        _sha(self.result_digest)
        _tuple(self.branches, CycleBranchAudit, "cycle branch audits", maximum=8)
        if ((self.selected_index is not None) != (self.status == "selected-declared-geometry-only")
                or self.selected_index is not None and (type(self.selected_index) is not int or not 0 <= self.selected_index < 8)
                or self.status not in {"selected-declared-geometry-only", "rejected", "unavailable", "aborted"}):
            raise ValueError("cycle-selection-audit-shape")

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def H_parent_approved(self) -> Literal[False]:
        return False

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False


def _same(actual: object, expected: object, code: str) -> None:
    if digest_record(actual) != digest_record(expected):
        raise ValueError(code)


def _audit_trace(plan: CycleSolvePlan, branch: CycleBranchResult) -> None:
    """Verify transport/objective/counter consistency; NEVER rerun the optimizer.

    Counts are observed execution metadata, not a proof that an optimizer followed
    a particular trajectory. Trusted artifact hashes/code pins bind that evidence;
    independent physical acceptance below does not rely on counts or objective.
    """
    expected = plan.branches[branch.index]
    if branch.branch_digest != expected.digest:
        raise ValueError("cycle-trace-plan-binding")
    if not expected.ready:
        if (branch.status != "refused" or branch.reason != expected.reason or branch.stats is not None
                or branch.returned_x is not None or branch.raw or branch.rounded or branch.initial_residuals is not None
                or branch.final_residuals is not None or branch.optimizer_residual_calls or branch.external_residual_calls
                or branch.solver_calls or branch.last_finite_x is not None):
            raise ValueError("cycle-trace-refusal-binding")
        return
    if branch.status == "refused" or branch.external_residual_calls < 1:
        raise ValueError("cycle-trace-admitted-branch-accounting")
    zero = tuple(0. for _ in range(3*len(plan.initializer.problem.coverage.variables)))
    try:
        initial = residual_values(plan, branch.index, zero)
        squared_residual(initial)
    except CycleNumericalError as exc:
        if (branch.status != "numerical-failed" or branch.reason != str(exc) or branch.initial_residuals is not None
                or branch.optimizer_residual_calls or branch.solver_calls or branch.returned_x is not None or branch.stats is not None
                or branch.raw or branch.rounded or branch.final_residuals is not None or branch.last_finite_x is not None
                or branch.external_residual_calls != 1):
            raise ValueError("cycle-trace-initial-failure-binding") from exc
        return
    _same(branch.initial_residuals, initial, "cycle-trace-initial-objective")
    if branch.last_finite_x is None or branch.last_finite_residuals is None:
        raise ValueError("cycle-trace-missing-partial-evidence")
    partial = residual_values(plan, branch.index, branch.last_finite_x)
    squared_residual(partial)
    _same(branch.last_finite_residuals, partial, "cycle-trace-partial-objective")
    if branch.status == "unavailable":
        if (not branch.reason.startswith("cycle-solver-dependency-unavailable:") or branch.solver_calls
                or branch.optimizer_residual_calls or branch.stats is not None or branch.returned_x is not None
                or branch.raw or branch.rounded or branch.final_residuals is not None or branch.external_residual_calls != 1):
            raise ValueError("cycle-trace-unavailable-binding")
        _same(branch.last_finite_x, zero, "cycle-trace-unavailable-partial")
        return
    if branch.solver_calls != 1:
        raise ValueError("cycle-trace-one-solver-call-required")
    if branch.stats is not None and branch.optimizer_residual_calls < branch.stats.nfev:
        raise ValueError("cycle-trace-solver-evaluation-count")
    if branch.optimizer_residual_calls > plan.settings.max_optimizer_calls:
        if branch.status != "numerical-failed" or branch.reason != "cycle-optimizer-residual-call-budget":
            raise ValueError("cycle-trace-budget-breach-not-retained")
    if branch.returned_x is None:
        if (branch.status not in {"numerical-failed", "aborted"} or branch.raw or branch.rounded
                or branch.final_residuals is not None or branch.external_residual_calls != 1):
            raise ValueError("cycle-trace-incomplete-output-status")
        return
    try:
        points = vector_points(plan, branch.index, branch.returned_x)
    except CycleNumericalError as exc:
        if branch.status != "numerical-failed" or branch.reason != str(exc) or branch.raw or branch.rounded:
            raise ValueError("cycle-trace-invalid-return-binding") from exc
        return
    start = plan.initializer.starts[branch.index]
    if tuple(o.identity for o in branch.raw) != tuple(o.identity for o in start.additions):
        raise ValueError("cycle-trace-returned-mask-or-order")
    for actual, original in zip(branch.raw, start.additions, strict=True):
        _same((actual.position, actual.element, actual.unit, actual.occupancy),
              (points[actual.identity], original.element, Unit.ANGSTROM, 1.), "cycle-trace-displacement-transport")
    if branch.stats is None or branch.external_residual_calls != 2:
        raise ValueError("cycle-trace-final-accounting")
    try:
        final = residual_values(plan, branch.index, branch.returned_x)
        squared_residual(final)
    except CycleNumericalError as exc:
        if branch.status != "numerical-failed" or branch.reason != str(exc) or branch.final_residuals is not None:
            raise ValueError("cycle-trace-final-failure-binding") from exc
        return
    if (branch.status != "coordinates-returned" or branch.reason != "finite-numerical-return-not-geometry-approved"
            or branch.optimizer_residual_calls < 1):
        raise ValueError("cycle-trace-complete-output-status")
    _same(branch.final_residuals, final, "cycle-trace-final-objective")
    _same(branch.last_finite_x, branch.returned_x, "cycle-trace-final-partial-vector")
    _same(branch.last_finite_residuals, final, "cycle-trace-final-partial-objective")


def measure_cycle_construction(plan: CycleSolvePlan, prepared: bytes,
                               result: CycleConstructionResult) -> CycleConstructionAudit:
    """Full source authentication/physical acceptance is OUTSIDE numerical process."""
    expected = compile_cycle_solve_plan(plan.initializer, prepared, plan.metadata)
    _same(plan, expected, "cycle-saved-solve-plan-binding")
    if result.plan_digest != plan.digest:
        raise ValueError("cycle-result-plan-binding")
    _same(result.environment, cycle_environment(), "cycle-result-runtime-binding")
    audits = []
    eligible = []
    for branch in result.branches:
        _audit_trace(plan, branch)
        pair = None
        if branch.raw and branch.rounded:
            pair = measure_cycle_pair(plan.initializer.problem, prepared, plan.metadata, plan.initializer,
                                       branch.index, branch.raw, branch.rounded)
        selectable = (result.aborted_reason is None and branch.status == "coordinates-returned"
                      and branch.stats is not None and branch.stats.success and pair is not None and pair.passed)
        if selectable:
            assert branch.final_residuals is not None and pair is not None and pair.raw_maximum_movement_angstrom is not None
            eligible.append((squared_residual(branch.final_residuals), pair.raw_maximum_movement_angstrom, branch.index))
        audits.append(CycleBranchAudit(branch.index, expected.initial_audits[branch.index], pair, selectable,
            "independent-declared-geometry-only" if selectable else "retained-not-selectable:"+branch.status))
    selected = min(eligible)[2] if eligible else None
    status: Literal["selected-declared-geometry-only", "rejected", "unavailable", "aborted"]
    if result.aborted_reason is not None:
        status = "aborted"
    elif selected is not None:
        status = "selected-declared-geometry-only"
    elif any(b.status == "unavailable" for b in result.branches):
        status = "unavailable"
    else:
        status = "rejected"
    return CycleConstructionAudit(plan.digest, result.digest, tuple(audits), selected, status)


def audit_cycle_construction(plan: CycleSolvePlan, prepared: bytes, result: CycleConstructionResult,
                             report: CycleConstructionAudit) -> None:
    _same(report, measure_cycle_construction(plan, prepared, result), "cycle-saved-construction-audit-binding")
