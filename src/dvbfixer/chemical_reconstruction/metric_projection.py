"""One optional private generated-heavy least-squares proposal, never a repair loop."""
from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace

from dvbfixer.domain.structure_identity import ExactAtomRef

from .conformer_search import _rounded
from .contracts import ReconstructionRequest
from .projection_audit import measure_projection
from .projection_contracts import ProjectionResult
from .projection_geometry import compile_projection_plan, observed_frame, residual_values
from .proposals import Candidate, GeometryPlan, ProposalResult, Strategy


@dataclass(frozen=True)
class ProjectionSolution:
    x: tuple[float, ...]
    success: bool
    status: int
    message: str
    nfev: int


def _solve(function: Callable[[Sequence[float]], tuple[float, ...]], variables: int) -> ProjectionSolution:
    from scipy.optimize import least_squares

    result = least_squares(function, [0.] * variables, bounds=(-2., 2.), method="trf", jac="2-point",
                           loss="linear", x_scale=1., ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=500)
    return ProjectionSolution(tuple(float(x) for x in result.x), bool(result.success), int(result.status),
                              str(result.message), int(result.nfev))


def propose_projection(original: ReconstructionRequest, request: ReconstructionRequest, baseline_plan: GeometryPlan,
                         baseline: ProposalResult, anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]) -> ProjectionResult:
    try:
        plan = compile_projection_plan(original, request, baseline_plan, baseline, anchors)
    except ValueError as exc:
        return ProjectionResult("refused", str(exc))
    assert baseline.rounded is not None
    base = baseline.rounded.additions
    frame = observed_frame(request, anchors)

    def candidate(x: Sequence[float]) -> Candidate:
        if len(x) != 3*len(base) or not all(math.isfinite(float(value)) for value in x):
            raise ValueError("projection-variable-shape-or-finiteness")
        additions = []
        for i, atom in enumerate(base):
            p = tuple(atom.position[j]+sum(float(x[3*i+k])*frame[k][j] for k in range(3)) for j in range(3))
            additions.append(replace(atom, position=(p[0], p[1], p[2])))
        return Candidate(request.digest, plan.digest, Strategy.PROJECTION, tuple(additions))

    def function(x: Sequence[float]) -> tuple[float, ...]:
        return residual_values(request, plan, candidate(x))

    initial = None
    try:
        initial = math.fsum(x*x for x in function([0.] * (3*len(base))))
        solved = _solve(function, 3*len(base))
        raw = candidate(solved.x)
        rounded = _rounded(raw)
        final = math.fsum(x*x for x in residual_values(request, plan, raw))
    except (ImportError, ValueError, FloatingPointError) as exc:
        return ProjectionResult("solver-failed", str(exc), plan, initial_squared_residual=initial)
    audit = measure_projection(request, plan, raw, rounded, solved.success)
    return ProjectionResult("projected-geometry-passed" if audit.passed else "rejected",
                            "declared-projection-geometry-only" if audit.passed else "retained-projection-failure",
                            plan, raw, rounded, audit, solved.success, solved.status, solved.message, solved.nfev, initial, final)
