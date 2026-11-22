"""Read-only graph/frame/metric diagnosis; never constructs alternative points."""
from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations

from dvbfixer.domain.structure_identity import ExactAtomRef

from .authority_geometry import geometry_view
from .candidate_audit import measure_additions
from .contracts import (
    Constraint,
    ConstraintKind,
    Observation,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from .measurement import (
    _cross,
    _dot,
    _sub,
    _unit,
    angle_degrees,
    ordered_stereo_volume,
    satisfied,
    torsion_degrees,
)
from .metric_geometry import audit_policy_for, explicit_metrics, metric_key
from .proposals import AuditPolicy, GeometryPlan, PlacementKind


@dataclass(frozen=True)
class MeasuredTerm:
    family: str
    kind: str
    atoms: tuple[ExactAtomRef, ...]
    role: str
    unit: str
    value: float | None
    lower: float | None
    upper: float | None
    passed: bool
    label: str = ""
    error: str | None = None


@dataclass(frozen=True)
class TermTable:
    request_digest: str
    plan_digest: str
    additions_digest: str
    terms: tuple[MeasuredTerm, ...]
    audit_findings: tuple[str, ...]
    explained_findings: tuple[str, ...]
    unexplained_findings: tuple[str, ...]
    scalar_only_findings: tuple[str, ...]


def _positions(request: ReconstructionRequest, additions: tuple[Observation, ...]) -> dict[ExactAtomRef, tuple[float, float, float]]:
    graph = {ExactAtomRef(c.identity, a.name): a for c in request.components if c.resolved is not None for a in c.resolved.graph.atoms}
    if (not isinstance(additions, tuple) or len(additions) != len(request.ownership.generated)
            or {o.identity for o in additions} != request.ownership.generated
            or request.ownership.generated & {o.identity for o in request.observations}
            or any(o.unit != Unit.ANGSTROM or o.occupancy != 1 or o.identity not in graph
                   or graph[o.identity].element in {"H", "D"} or graph[o.identity].element != o.element for o in additions)):
        raise ValueError("diagnosis-needs-complete-exact-heavy-additions")
    points = {o.identity: o.angstrom for o in request.observations}
    points.update({o.identity: o.angstrom for o in additions})
    if any(not all(math.isfinite(v) for v in p) for p in points.values()):
        raise ValueError("diagnosis-nonfinite-coordinates")
    return points


def _constraint_value(constraint: Constraint, points: tuple[tuple[float, float, float], ...]) -> tuple[float, float | None, float | None]:
    kind = constraint.kind
    if kind in {ConstraintKind.STEREO_FRAME, ConstraintKind.STEREO_LIGANDS}:
        assert constraint.orientation is not None
        return constraint.orientation*ordered_stereo_volume(kind, points), constraint.minimum_volume, None
    if kind in {ConstraintKind.DISTANCE, ConstraintKind.RING_CLOSURE, ConstraintKind.EXCLUSION}:
        value = math.dist(points[0], points[1])
    elif kind == ConstraintKind.ANGLE:
        value = angle_degrees(points[0], points[1], points[2])
    elif kind == ConstraintKind.TORSION:
        value = torsion_degrees((points[0], points[1], points[2], points[3]))
    elif kind == ConstraintKind.PLANARITY:
        normal = _unit(_cross(_sub(points[1], points[0]), _sub(points[2], points[0])))
        value = max(abs(_dot(_sub(p, points[0]), normal)) for p in points[3:])
    else:
        raise ValueError("unknown-diagnostic-constraint")
    if constraint.unit == Unit.NANOMETER:
        value /= 10
    return value, constraint.lower, constraint.upper


def measure_terms(request: ReconstructionRequest, plan: GeometryPlan, additions: tuple[Observation, ...]) -> TermTable:
    """Explain existing scalar audit findings without replacing its decision."""
    positions = _positions(request, additions)
    if plan.request_digest != request.digest:
        raise ValueError("diagnostic-plan-request-binding")
    ideals, bonds = geometry_view(request)
    metrics = explicit_metrics(plan)
    terms, explained = [], set()
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in ideals}

    def role(atoms: tuple[ExactAtomRef, ...]) -> str:
        return "generated-involving" if set(atoms) & request.ownership.generated else "fixed-only"

    def metric(family: str, atoms: tuple[ExactAtomRef, ...], target: float, tolerance: float) -> None:
        error = None
        value = None
        try:
            value = math.dist(positions[atoms[0]], positions[atoms[1]]) if len(atoms) == 2 else angle_degrees(*(positions[a] for a in atoms))
        except ValueError as exc:
            error = str(exc)
        passed = value is not None and abs(value-target) <= tolerance
        terms.append(MeasuredTerm(family, "distance" if len(atoms) == 2 else "angle", atoms, role(atoms),
            "angstrom" if len(atoms) == 2 else "degree", value, target-tolerance, target+tolerance, passed, error=error))
        if not passed:
            if error:
                explained.add("degenerate-authority-or-candidate-angle" if family == "authority-angle" else "degenerate-candidate-metric")
            else:
                explained.add("explicit-metric-geometry" if family == "explicit-metric" else "authority-bond-geometry" if len(atoms) == 2 else "authority-angle-geometry")

    for pair, length in sorted(bonds.items(), key=lambda row: tuple(sorted(row[0]))):
        a, b = sorted(pair)
        adjacent[a].add(b)
        adjacent[b].add(a)
        if a.component in request.targets and {a, b} <= positions.keys():
            metric("authority-bond", (a, b), metrics.get(metric_key((a, b)), length), audit_policy_for(plan, (a, b), request.ownership.generated).bond_tolerance_angstrom)
    for centre in sorted(ideals):
        if centre.component not in request.targets or centre not in positions:
            continue
        for a, b in combinations(sorted(adjacent[centre]), 2):
            if {a, b} <= positions.keys():
                key = metric_key((a, centre, b))
                try:
                    target = metrics[key] if key in metrics else angle_degrees(ideals[a], ideals[centre], ideals[b])
                except ValueError as exc:
                    terms.append(MeasuredTerm("authority-angle", "angle", (a, centre, b), role((a, centre, b)), "degree", None, None, None, False, error=str(exc)))
                    explained.add("degenerate-authority-or-candidate-angle")
                else:
                    metric("authority-angle", (a, centre, b), target, audit_policy_for(plan, (a, centre, b), request.ownership.generated).angle_tolerance_degrees)
    for atoms, target in sorted(metrics.items()):
        if set(atoms) <= positions.keys():
            policy = audit_policy_for(plan, atoms, request.ownership.generated)
            metric("explicit-metric", atoms, target, policy.bond_tolerance_angstrom if len(atoms) == 2 else policy.angle_tolerance_degrees)
    for index, q in enumerate(request.constraints):
        value = lower = upper = None
        error = None
        try:
            if not set(q.atoms) <= positions.keys():
                raise ValueError("unmeasured-constraint")
            p = tuple(positions[a] for a in q.atoms)
            value, lower, upper = _constraint_value(q, p)
            passed = satisfied(q, p)
        except ValueError as exc:
            error, passed = str(exc), False
        if not passed:
            explained.add("candidate-constraint-conflict" if error is None else "unmeasured-constraint" if error == "unmeasured-constraint" else "degenerate-candidate-constraint")
        terms.append(MeasuredTerm("explicit-constraint", q.kind.value, q.atoms, role(q.atoms), q.unit.value,
            value, lower, upper, passed, label=f"{index}:{q.label}", error=error))
    audit = measure_additions(request, plan, additions)
    # Admission may stop the independent audit before candidate measurements.
    # Keep additional scalar observations separate, never relabel them as its
    # emitted findings or reinterpret an admission refusal as a candidate pass.
    additional = explained-set(audit.findings)
    explained &= set(audit.findings)
    for term in terms:
        if term.family == "explicit-constraint" and term.role == "fixed-only" and not term.passed:
            code = "degenerate-fixed-constraint" if term.error else "fixed-constraint-conflict"
            if code in audit.findings:
                explained.add(code)
    return TermTable(request.digest, plan.digest, digest_record(additions), tuple(terms), audit.findings,
                     tuple(sorted(explained)), tuple(sorted(set(audit.findings)-explained)), tuple(sorted(additional)))


def is_cycle_edge(edges: tuple[tuple[ExactAtomRef, ExactAtomRef], ...], edge: tuple[ExactAtomRef, ExactAtomRef]) -> bool:
    """An edge lies on a cycle iff its endpoints connect after removing that edge."""
    first, last = edge
    if first == last or frozenset(edge) not in {frozenset(e) for e in edges}:
        raise ValueError("diagnostic-cycle-needs-real-distinct-edge")
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {}
    for a, b in edges:
        if frozenset((a, b)) == frozenset(edge):
            continue
        adjacent.setdefault(a, set()).add(b)
        adjacent.setdefault(b, set()).add(a)
    seen, pending = {first}, [first]
    while pending:
        for neighbour in adjacent.get(pending.pop(), set()) - seen:
            if neighbour == last:
                return True
            seen.add(neighbour)
            pending.append(neighbour)
    return False


@dataclass(frozen=True)
class FrameMetric:
    kind: str
    atoms: tuple[ExactAtomRef, ...]
    target: float
    measured: float
    delta: float


@dataclass(frozen=True)
class FrameTrace:
    atom: ExactAtomRef
    kind: str
    anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]
    anchor_roles: tuple[str, ...]
    metrics: tuple[FrameMetric, ...]


@dataclass(frozen=True)
class UnenforcedLengthEdge:
    atoms: tuple[ExactAtomRef, ExactAtomRef]
    is_cycle: bool
    authority_length_angstrom: float
    measured_length_angstrom: float


def trace_frames(request: ReconstructionRequest, plan: GeometryPlan, additions: tuple[Observation, ...]
                 ) -> tuple[tuple[FrameTrace, ...], tuple[UnenforcedLengthEdge, ...]]:
    positions = _positions(request, additions)
    if (plan.request_digest != request.digest or {p.atom for p in plan.placements} != request.ownership.generated
            or len(plan.placements) != len(request.ownership.generated)):
        raise ValueError("diagnostic-exact-placement-binding")
    ideals, bonds = geometry_view(request)
    traces, length_edges = [], set()
    available = {o.identity for o in request.observations}
    for placement in plan.placements:
        a, b, c = placement.anchors
        d = placement.atom
        if not set(placement.anchors) <= available:
            raise ValueError("diagnostic-placement-anchor-dependency")
        frames: tuple[tuple[str, tuple[ExactAtomRef, ...]], ...]
        if placement.kind == PlacementKind.CHAIN:
            frames = (("distance", (c, d)), ("angle", (b, c, d)), ("torsion", (a, b, c, d)))
            length_edges.add(frozenset((c, d)))
        elif placement.kind == PlacementKind.CONE:
            frames = (("distance", (a, d)), ("angle", (b, a, d)), ("angle", (c, a, d)))
            length_edges.add(frozenset((a, d)))
        else:
            raise ValueError("diagnostic-chain-cone-only-not-an-alternative-constructor")
        values = []
        for kind, atoms in frames:
            if kind == "distance":
                target, measured = bonds[frozenset(atoms)], math.dist(*(positions[x] for x in atoms))
            elif kind == "angle":
                target, measured = angle_degrees(*(ideals[x] for x in atoms)), angle_degrees(*(positions[x] for x in atoms))
            else:
                four = (a, b, c, d)
                target = torsion_degrees(tuple_points(ideals, four))
                measured = torsion_degrees(tuple_points(positions, four))
            delta = (measured-target+180) % 360-180 if kind == "torsion" else measured-target
            values.append(FrameMetric(kind, atoms, target, measured, delta))
        traces.append(FrameTrace(d, placement.kind.value, placement.anchors,
            tuple("generated-conditioning" if x in request.ownership.generated else "fixed-observation" for x in placement.anchors), tuple(values)))
        available.add(d)
    graph = {ExactAtomRef(c.identity, a.name): a for c in request.components if c.resolved is not None for a in c.resolved.graph.atoms}
    edges = tuple((min(pair), max(pair)) for pair in bonds if all(graph[a].element not in {"H", "D"} for a in pair))
    omitted = tuple(UnenforcedLengthEdge((a, b), is_cycle_edge(edges, (a, b)), bonds[frozenset((a, b))], math.dist(positions[a], positions[b]))
                    for a, b in sorted(edges) if {a, b} & request.ownership.generated
                    and frozenset((a, b)) not in length_edges and {a, b} <= positions.keys())
    return tuple(traces), omitted


def tuple_points(points: dict[ExactAtomRef, tuple[float, float, float]], atoms: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef, ExactAtomRef]
                 ) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    return points[atoms[0]], points[atoms[1]], points[atoms[2]], points[atoms[3]]


@dataclass(frozen=True)
class NecessarySupport:
    anchor_distance_angstrom: float
    first_radius_interval: tuple[float, float]
    second_radius_interval: tuple[float, float]
    anchor_distance_interval: tuple[float, float]
    conservative_cosine_interval: tuple[float, float]
    allowed_cosine_interval: tuple[float, float]
    distance_shell_possible: bool
    angle_interval_possible: bool
    scope: str = "necessary-only; not a jointly feasible pose or permission to move fixed observations"


def necessary_support(distance: float, first_length: float, second_length: float, bond_tolerance: float,
                         angle: float, angle_tolerance: float) -> NecessarySupport:
    values = (distance, first_length, second_length, bond_tolerance, angle, angle_tolerance)
    if (any(type(v) not in {int, float} or not math.isfinite(v) for v in values)
            or distance <= 1e-12 or bond_tolerance <= 0 or min(first_length, second_length) <= bond_tolerance
            or not 0 < angle < 180 or not 0 < angle_tolerance < 180):
        raise ValueError("diagnostic-needs-finite-nondegenerate-positive-intervals")
    # Outward rounding prevents cancellation at a valid boundary from becoming
    # a spurious necessary-test rejection. This does NOT widen any audit bound.
    def down(v: float) -> float:
        return math.nextafter(v, -math.inf)
    def up(v: float) -> float:
        return math.nextafter(v, math.inf)
    l1, u1 = down(first_length-bond_tolerance), up(first_length+bond_tolerance)
    l2, u2 = down(second_length-bond_tolerance), up(second_length+bond_tolerance)
    chord = max(0., down(l1-u2), down(l2-u1)), up(u1+u2)
    numerator = (down(down(down(l1*l1)+down(l2*l2))-up(distance*distance)),
                 up(up(up(u1*u1)+up(u2*u2))-down(distance*distance)))
    denominator = down(2*down(l1*l2)), up(2*up(u1*u2))
    if min(denominator) <= 0 or not all(math.isfinite(v) for v in (*numerator, *denominator)):
        raise ValueError("diagnostic-interval-numeric-overflow-or-underflow")
    quotients = tuple(a/b for a in numerator for b in denominator)
    cosine = down(min(quotients)), up(max(quotients))
    allowed = down(math.cos(math.radians(min(180., up(angle+angle_tolerance))))), up(math.cos(math.radians(max(0., down(angle-angle_tolerance)))))
    return NecessarySupport(distance, (l1, u1), (l2, u2), chord, cosine, allowed,
        chord[0] <= distance <= chord[1], max(-1., cosine[0], allowed[0]) <= min(1., cosine[1], allowed[1]))


def fixed_neighbour_support(request: ReconstructionRequest, policy: AuditPolicy) -> tuple[tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef, NecessarySupport | str], ...]:
    ideals, bonds = geometry_view(request)
    observed = {o.identity: o.angstrom for o in request.observations}
    output = []
    for missing in sorted(request.ownership.generated):
        neighbours = sorted(next(iter(pair-{missing})) for pair in bonds if missing in pair and pair-{missing} <= observed.keys())
        for first, second in combinations(neighbours, 2):
            try:
                support: NecessarySupport | str = necessary_support(math.dist(observed[first], observed[second]),
                    bonds[frozenset((first, missing))], bonds[frozenset((second, missing))], policy.bond_tolerance_angstrom,
                    angle_degrees(ideals[first], ideals[missing], ideals[second]), policy.angle_tolerance_degrees)
            except ValueError as exc:
                support = str(exc)
            output.append((missing, first, second, support))
    return tuple(output)
