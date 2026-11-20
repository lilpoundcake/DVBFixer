"""Pure graph-metric compilation and residual measurements, no optimizer import."""
from __future__ import annotations

import math
from itertools import combinations

from dvbfixer.domain.structure_identity import ExactAtomRef

from .authority_geometry import geometry_view
from .candidate_audit import audit_candidate
from .construction import proposal_admission
from .contracts import (
    ConstraintKind,
    Operation,
    PreparedContextRequest,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from .measurement import Point, _cross, _dot, _sub, _unit, angle_degrees, ordered_stereo_volume
from .projection_contracts import (
    POLICY,
    ProjectionMetric,
    ProjectionPlan,
    append_projection_constraints,
)
from .proposals import Candidate, GeometryPlan, ProposalResult, Strategy


def observed_frame(request: ReconstructionRequest, anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]
                   ) -> tuple[Point, Point, Point]:
    observed = {o.identity: o for o in request.observations}
    if any(a not in observed or observed[a].occupancy != 1 for a in anchors):
        raise ValueError("projection-frame-needs-exact-fixed-unit-observations")
    a, b, c = (observed[x].angstrom for x in anchors)
    u = _unit(_sub(b, a))
    w = _unit(_cross(u, _sub(c, a)))
    return u, _cross(w, u), w


def compile_projection_plan(original: ReconstructionRequest, request: ReconstructionRequest, baseline_plan: GeometryPlan,
                              baseline: ProposalResult, anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]) -> ProjectionPlan:
    if (type(original) not in {ReconstructionRequest, PreparedContextRequest} or original.operation != Operation.HEAVY
            or original.ownership.mobile or original.ownership.temporary or not 1 <= len(original.ownership.generated) <= 16
            or len(original.targets) != 1 or baseline_plan.policy != POLICY):
        raise ValueError("projection-needs-bounded-fixed-heavy-input-and-policy")
    prefix = len(original.constraints)
    if (digest_record(request) != digest_record(append_projection_constraints(original, request.constraints[prefix:]))
            or baseline.status != "declared-geometry-passed" or baseline.plan_digest != baseline_plan.digest
            or baseline.admission != proposal_admission(original, baseline_plan)
            or baseline.raw is None or baseline.rounded is None):
        raise ValueError("projection-original-baseline-or-new-request-binding")
    for candidate, saved in ((baseline.raw, baseline.raw_audit), (baseline.rounded, baseline.rounded_audit)):
        audit = audit_candidate(original, baseline_plan, candidate)
        if not audit.passed or audit != saved:
            raise ValueError("projection-original-baseline-geometry-not-passed")
    from .conformer_search import _rounded
    if digest_record(_rounded(baseline.raw)) != digest_record(baseline.rounded):
        raise ValueError("projection-baseline-rounding-transport")
    observed_frame(request, anchors)
    ideals, bonds = geometry_view(request)
    elements = {ExactAtomRef(c.identity, a.name): a.element for c in request.components if c.resolved is not None
                for a in c.resolved.graph.atoms}
    target = next(iter(request.targets))
    heavy = {a for a in ideals if a.component == target and elements[a] not in {"H", "D"}}
    selected = {o.identity for o in request.observations if o.identity.component == target}
    if (selected | request.ownership.generated != heavy or selected & request.ownership.generated
            or not request.ownership.generated <= heavy or any(o.occupancy != 1 for o in request.observations if o.identity in selected)
            or {o.identity for o in baseline.rounded.additions} != request.ownership.generated):
        raise ValueError("projection-exact-complete-heavy-only-inventory")
    terms = []
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in heavy}
    for pair, length in sorted(bonds.items(), key=lambda row: tuple(sorted(row[0]))):
        if not pair <= heavy:
            continue
        a, b = sorted(pair)
        adjacent[a].add(b)
        adjacent[b].add(a)
        terms.append(ProjectionMetric((a, b), length, POLICY.bond_tolerance_angstrom))
    for centre in sorted(heavy):
        for first, last in combinations(sorted(adjacent[centre]), 2):
            terms.append(ProjectionMetric((first, centre, last), angle_degrees(ideals[first], ideals[centre], ideals[last]),
                                          POLICY.angle_tolerance_degrees))
    count = len(terms)
    present = {o.identity for o in request.observations} | request.ownership.generated
    for constraint in request.constraints:
        if constraint.kind == ConstraintKind.TORSION or not set(constraint.atoms) <= present:
            raise ValueError("projection-unqualified-constraint-or-incomplete-frame")
        count += len(constraint.atoms)-3 if constraint.kind == ConstraintKind.PLANARITY else 1
        if constraint.kind in {ConstraintKind.DISTANCE, ConstraintKind.RING_CLOSURE, ConstraintKind.ANGLE}:
            if constraint.lower is None or constraint.upper == constraint.lower:
                raise ValueError("projection-unqualified-unbounded-or-zero-width-interval")
        if constraint.kind == ConstraintKind.EXCLUSION and (constraint.lower is None or constraint.lower <= 0):
            raise ValueError("projection-exclusion-needs-positive-lower-bound")
    if count > 512:
        raise ValueError("projection-residual-budget")
    return ProjectionPlan(request.digest, Strategy.PROJECTION, POLICY, original=original, baseline_plan=baseline_plan,
                          baseline=baseline, anchors=anchors, terms=tuple(terms))


def residual_values(request: ReconstructionRequest, plan: ProjectionPlan, candidate: Candidate) -> tuple[float, ...]:
    """Measure the fixed residual definition on supplied coordinates only."""
    points = {o.identity: o.angstrom for o in request.observations}
    points.update({o.identity: o.position for o in candidate.additions})
    residuals = []
    for term in plan.terms:
        p = tuple(points[a] for a in term.atoms)
        measured = math.dist(p[0], p[1]) if len(p) == 2 else angle_degrees(p[0], p[1], p[2])
        residuals.append((measured-term.value)/term.scale)
    for constraint in request.constraints:
        p = tuple(points[a] for a in constraint.atoms)
        kind = constraint.kind
        if kind in {ConstraintKind.STEREO_FRAME, ConstraintKind.STEREO_LIGANDS}:
            assert constraint.orientation is not None and constraint.minimum_volume is not None
            residuals.append(max(0., constraint.minimum_volume-constraint.orientation*ordered_stereo_volume(kind, p))/constraint.minimum_volume)
        elif kind == ConstraintKind.PLANARITY:
            normal = _unit(_cross(_sub(p[1], p[0]), _sub(p[2], p[0])))
            assert constraint.upper is not None
            scale = constraint.upper * (10 if constraint.unit == Unit.NANOMETER else 1)
            residuals.extend(_dot(_sub(point, p[0]), normal)/scale for point in p[3:])
        else:
            measured = angle_degrees(p[0], p[1], p[2]) if kind == ConstraintKind.ANGLE else math.dist(p[0], p[1])
            if constraint.unit == Unit.NANOMETER:
                measured /= 10
            assert constraint.lower is not None and constraint.upper is not None
            if kind == ConstraintKind.EXCLUSION:
                residuals.append(max(0., constraint.lower-measured)/constraint.lower)
            else:
                residuals.append((measured-(constraint.lower+constraint.upper)/2)/((constraint.upper-constraint.lower)/2))
    if not all(math.isfinite(value) for value in residuals):
        raise ValueError("nonfinite-projection-residual")
    return tuple(residuals)
