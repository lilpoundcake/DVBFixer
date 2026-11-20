"""Independent projection transport/geometry/movement audit, never re-solves."""
from __future__ import annotations

import math

from dvbfixer.domain.structure_identity import ExactAtomRef

from .candidate_audit import audit_candidate
from .conformer_search import _rounded
from .contracts import ReconstructionRequest, digest_record
from .projection_contracts import ProjectionAudit, ProjectionPlan, ProjectionResult
from .projection_geometry import compile_projection_plan, residual_values
from .proposals import Candidate, GeometryPlan, ProposalResult


def measure_projection(request: ReconstructionRequest, plan: ProjectionPlan, raw: Candidate, rounded: Candidate,
                         solver_success: bool) -> ProjectionAudit:
    errors = set()
    raw_geometry, rounded_geometry = audit_candidate(request, plan, raw), audit_candidate(request, plan, rounded)
    if not raw_geometry.passed or not rounded_geometry.passed:
        errors.add("projection-raw-or-rounded-geometry")
    if not solver_success:
        errors.add("projection-solver-not-converged")
    try:
        if digest_record(_rounded(raw)) != digest_record(rounded):
            errors.add("projection-rounding-transport")
    except ValueError:
        errors.add("projection-rounding-transport")
    assert plan.baseline.rounded is not None
    original_points = {o.identity: o.position for o in plan.baseline.rounded.additions}
    movements: list[float | None] = []
    for candidate in (raw, rounded):
        if ({a.identity for a in candidate.additions} != original_points.keys()
                or len(candidate.additions) != len(original_points)):
            errors.add("projection-movement-mask")
            movements.append(None)
            continue
        maximum = max(math.dist(a.position, original_points[a.identity]) for a in candidate.additions)
        movements.append(maximum)
        if maximum > plan.per_atom_movement_cap_angstrom:
            errors.add("projection-per-atom-movement-cap")
    return ProjectionAudit(not errors, tuple(sorted(errors)), raw_geometry, rounded_geometry, movements[0], movements[1])


def audit_projection_result(original: ReconstructionRequest, request: ReconstructionRequest, baseline_plan: GeometryPlan,
                             baseline: ProposalResult, anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef],
                             result: ProjectionResult) -> tuple[str, ...]:
    """Saved-result integrity is separate from acceptance of a rejected projection."""
    try:
        expected = compile_projection_plan(original, request, baseline_plan, baseline, anchors)
    except ValueError as exc:
        if (result.status != "refused" or result.reason != str(exc)
                or any(x is not None for x in (result.plan, result.raw, result.rounded, result.audit,
                                                result.solver_success, result.solver_status, result.solver_message, result.nfev,
                                                result.initial_squared_residual, result.final_squared_residual))):
            return ("projection-refusal-record-mismatch",)
        return ()
    if result.plan is None or digest_record(result.plan) != digest_record(expected):
        return ("projection-plan-authority-or-source-binding",)
    if result.raw is None or result.rounded is None:
        if (result.status != "solver-failed" or not result.reason or result.audit is not None
                or result.raw is not None or result.rounded is not None
                or any(x is not None for x in (result.solver_success, result.solver_status, result.solver_message,
                                                result.nfev, result.final_squared_residual))):
            return ("projection-incomplete-solver-failure-record",)
        if result.initial_squared_residual is not None:
            assert result.plan.baseline.rounded is not None
            try:
                initial = math.fsum(x*x for x in residual_values(request, result.plan, result.plan.baseline.rounded))
            except (ValueError, KeyError):
                return ("projection-initial-objective-unmeasurable",)
            if initial != result.initial_squared_residual:
                return ("projection-initial-objective-metadata",)
        return ()
    errors = set()
    if (type(result.solver_success) is not bool or type(result.solver_status) is not int
            or result.solver_status not in {0, 1, 2, 3, 4} or result.solver_success != (result.solver_status > 0)
            or type(result.nfev) is not int or not 1 <= result.nfev <= result.plan.max_nfev or not result.solver_message):
        errors.add("projection-solver-metadata")
    audit = measure_projection(request, result.plan, result.raw, result.rounded, bool(result.solver_success))
    if digest_record(audit) != digest_record(result.audit):
        errors.add("projection-saved-audit-mismatch")
    status = "projected-geometry-passed" if audit.passed else "rejected"
    if result.status != status or result.reason != ("declared-projection-geometry-only" if audit.passed else "retained-projection-failure"):
        errors.add("projection-status-or-reason")
    assert result.plan.baseline.rounded is not None
    try:
        initial = math.fsum(x*x for x in residual_values(request, result.plan, result.plan.baseline.rounded))
        final = math.fsum(x*x for x in residual_values(request, result.plan, result.raw))
        if (result.initial_squared_residual != initial or result.final_squared_residual != final):
            errors.add("projection-objective-metadata")
    except (ValueError, KeyError):
        errors.add("projection-objective-unmeasurable")
    return tuple(sorted(errors))
