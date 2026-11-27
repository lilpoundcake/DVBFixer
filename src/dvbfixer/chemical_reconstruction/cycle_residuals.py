"""Frozen R1 numerical residuals, not independent geometry/acceptance."""
from __future__ import annotations

import math

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import ConstraintKind, ContractError, Observation, Unit, _number, _point
from .cycle_solve_plan import CycleSolvePlan, point_order
from .measurement import Point, _cross, _dot, _sub, _unit, angle_degrees, ordered_stereo_volume


class CycleNumericalError(ValueError):
    """Declared numerical failure; no restart, zero penalty or pose fallback."""


def vector_points(plan: CycleSolvePlan, index: int, x: tuple[float, ...]) -> dict[ExactAtomRef, Point]:
    initializer = plan.initializer
    variables = initializer.problem.coverage.variables
    if (type(index) is not int or not 0 <= index < 8 or not isinstance(x, tuple) or len(x) != 3*len(variables)
            or any(not _number(v) or abs(v) > plan.settings.coordinate_bound_angstrom for v in x)):
        raise CycleNumericalError("cycle-displacement-vector-shape-bounds-or-finiteness")
    start = initializer.starts[index]
    if start.status != "coordinates-initialized" or initializer.system is None:
        raise CycleNumericalError("cycle-no-initialized-start")
    frame = initializer.system.frame
    points = {o.identity: o.angstrom for o in initializer.problem.source.observations}
    for j, addition in enumerate(start.additions):
        try:
            point = tuple(addition.position[k]+math.fsum(x[3*j+t]*frame.axes[t][k] for t in range(3)) for k in range(3))
            _point(point)
        except (ContractError, OverflowError) as exc:
            raise CycleNumericalError("cycle-coordinate-representation-or-finiteness") from exc
        points[addition.identity] = point[0], point[1], point[2]
    return points


def vector_additions(plan: CycleSolvePlan, index: int, x: tuple[float, ...]) -> tuple[Observation, ...]:
    points = vector_points(plan, index, x)
    return tuple(Observation(o.identity, o.element, points[o.identity], Unit.ANGSTROM) for o in plan.initializer.starts[index].additions)


def residual_values(plan: CycleSolvePlan, index: int, x: tuple[float, ...]) -> tuple[float, ...]:
    branch = plan.branches[index]
    if not branch.ready:
        raise CycleNumericalError("cycle-branch-not-numerically-admitted")
    points = vector_points(plan, index, x)
    problem = plan.initializer.problem
    values = []
    try:
        for term in problem.coverage.metrics:
            p = tuple(points[a] for a in term.atoms)
            value = math.dist(p[0], p[1]) if term.kind == "bond" else angle_degrees(*p)
            values.append((value-term.target)/term.tolerance)
        for q in (*problem.source.constraints, *problem.extra_constraints):
            p = tuple(points[a] for a in q.atoms)
            if q.kind in {ConstraintKind.STEREO_FRAME, ConstraintKind.STEREO_LIGANDS}:
                assert q.minimum_volume is not None and q.orientation is not None
                values.append(max(0., q.minimum_volume-q.orientation*ordered_stereo_volume(q.kind, p))/q.minimum_volume)
            elif q.kind == ConstraintKind.PLANARITY:
                normal = _unit(_cross(_sub(p[1], p[0]), _sub(p[2], p[0])))
                assert q.upper is not None
                scale = q.upper*(10. if q.unit == Unit.NANOMETER else 1.)
                values.extend(_dot(_sub(point, p[0]), normal)/scale for point in p[3:])
            else:
                value = angle_degrees(*p) if q.kind == ConstraintKind.ANGLE else math.dist(p[0], p[1])
                if q.unit == Unit.NANOMETER:
                    value /= 10.
                assert q.lower is not None and q.upper is not None
                if q.kind == ConstraintKind.EXCLUSION:
                    values.append(max(0., q.lower-value)/q.lower+max(0., value-q.upper)/q.upper)
                else:
                    values.append((value-(q.lower+q.upper)/2)/((q.upper-q.lower)/2))
        order = point_order(plan.initializer)
        separation = problem.policy.minimum_separation_angstrom
        values.extend(max(0., separation-math.dist(points[order[a]], points[order[b]]))/separation for a, b in branch.pairs)
    except (ValueError, ZeroDivisionError, OverflowError) as exc:
        raise CycleNumericalError("cycle-degenerate-or-nonfinite-residual:"+str(exc)) from exc
    if len(values) != branch.residual_count or not all(math.isfinite(v) for v in values):
        raise CycleNumericalError("cycle-residual-size-or-finiteness")
    return tuple(values)


def squared_residual(values: tuple[float, ...]) -> float:
    try:
        result = math.fsum(v*v for v in values)
    except (OverflowError, ValueError) as exc:
        raise CycleNumericalError("cycle-nonfinite-squared-residual") from exc
    if not math.isfinite(result):
        raise CycleNumericalError("cycle-nonfinite-squared-residual")
    return result
