"""Measurement-only candidate audit; never calls a proposal constructor."""
from __future__ import annotations

import math
from collections import defaultdict
from itertools import combinations, product

from dvbfixer.domain.structure_identity import ExactAtomRef

from .admission import admit
from .authority_geometry import geometry_view
from .contracts import LinkKind, ReconstructionRequest, Unit
from .joint_hydrogen_geometry import joint_group_findings
from .measurement import Point, angle_degrees, satisfied
from .metric_geometry import audit_policy_for, explicit_metrics, metric_findings, metric_key
from .proposals import Candidate, CandidateAudit, GeometryPlan


def audit_candidate(request: ReconstructionRequest, plan: GeometryPlan, candidate: Candidate) -> CandidateAudit:
    planning = admit(request)
    errors = {f.code for f in planning.findings if f.code != "no-declared-capability-matches"}
    errors.update(metric_findings(request, plan))
    errors.update(joint_group_findings(request, plan))
    if candidate.request_digest != request.digest or plan.request_digest != request.digest:
        errors.add("request-binding")
    if candidate.plan_digest != plan.digest:
        errors.add("plan-binding")
    if candidate.strategy != plan.strategy:
        errors.add("strategy-binding")
    additions = {o.identity: o for o in candidate.additions}
    if len(additions) != len(candidate.additions) or set(additions) != request.ownership.generated:
        errors.add("generated-mask")
    if request.ownership.temporary or request.ownership.mobile:
        errors.add("unsupported-ownership")
    chemical = {ExactAtomRef(c.identity, a.name): a for c in request.components if c.resolved is not None
                for a in c.resolved.graph.atoms}
    for identity, addition in additions.items():
        if identity not in chemical or addition.element.upper() != chemical[identity].element.upper():
            errors.add("generated-element")
        if addition.unit != Unit.ANGSTROM or addition.occupancy != 1:
            errors.add("generated-unit-or-occupancy")
    if errors:
        return CandidateAudit(False, tuple(sorted(errors)), (), planning.unassessed_context)

    positions = {o.identity: o.angstrom for o in request.observations}
    positions.update({a: o.position for a, o in additions.items()})
    ideals, bonds = geometry_view(request)
    metrics = explicit_metrics(plan)
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in ideals}
    for pair, length in bonds.items():
        a, b = sorted(pair)
        adjacent[a].add(b)
        adjacent[b].add(a)
        if a.component in request.targets and a in positions and b in positions:
            target = metrics.get(metric_key((a, b)), length)
            policy = audit_policy_for(plan, (a, b), request.ownership.generated)
            if abs(math.dist(positions[a], positions[b])-target) > policy.bond_tolerance_angstrom:
                errors.add("authority-bond-geometry")
    for centre in sorted(ideals):
        if centre.component not in request.targets or centre not in positions:
            continue
        for a, b in combinations(sorted(adjacent[centre]), 2):
            if a not in positions or b not in positions:
                continue
            try:
                key = metric_key((a, centre, b))
                ideal = metrics[key] if key in metrics else angle_degrees(ideals[a], ideals[centre], ideals[b])
                measured = angle_degrees(positions[a], positions[centre], positions[b])
            except ValueError:
                errors.add("degenerate-authority-or-candidate-angle")
            else:
                policy = audit_policy_for(plan, (a, centre, b), request.ownership.generated)
                if abs(measured-ideal) > policy.angle_tolerance_degrees:
                    errors.add("authority-angle-geometry")
    # Explicit cross-component metrics have no shared ideal-coordinate frame.
    # Measure every declared metric, independently of constructor completion.
    for atoms, target in metrics.items():
        if not set(atoms) <= positions.keys():
            errors.add("unmeasured-metric")
            continue
        try:
            measured = math.dist(positions[atoms[0]], positions[atoms[1]]) if len(atoms) == 2 else angle_degrees(
                positions[atoms[0]], positions[atoms[1]], positions[atoms[2]])
        except ValueError:
            errors.add("degenerate-candidate-metric")
        else:
            policy = audit_policy_for(plan, atoms, request.ownership.generated)
            tolerance = policy.bond_tolerance_angstrom if len(atoms) == 2 else policy.angle_tolerance_degrees
            if abs(measured-target) > tolerance:
                errors.add("explicit-metric-geometry")
    checked = []
    for index, constraint in enumerate(request.constraints):
        if not set(constraint.atoms) <= positions.keys():
            errors.add("unmeasured-constraint")
            continue
        checked.append(index)
        try:
            passed = satisfied(constraint, tuple(positions[a] for a in constraint.atoms))
        except ValueError:
            errors.add("degenerate-candidate-constraint")
        else:
            if not passed:
                errors.add("candidate-constraint-conflict")

    # Links extend graph-distance exclusions only; missing link geometry authority
    # is never manufactured from the two separately framed component templates.
    for connection in request.connections:
        if connection.kind == LinkKind.COVALENT:
            a, b = connection.link.component_atom, connection.link.other_atom
            adjacent.setdefault(a, set()).add(b)
            adjacent.setdefault(b, set()).add(a)
    separation = plan.policy.minimum_separation_angstrom
    buckets: dict[tuple[int, int, int], list[ExactAtomRef]] = defaultdict(list)

    def cell(position: Point) -> tuple[int, int, int]:
        return tuple_floor(position[0]/separation, position[1]/separation, position[2]/separation)

    for a, position in positions.items():
        buckets[cell(position)].append(a)
    for a in sorted(additions):
        near = adjacent.get(a, set())
        excluded = {a} | near | {n for b in near for n in adjacent.get(b, set())}
        x, y, z = cell(positions[a])
        for dx, dy, dz in product((-1, 0, 1), repeat=3):
            for b in buckets.get((x+dx, y+dy, z+dz), ()):
                if b not in excluded and math.dist(positions[a], positions[b]) < separation:
                    errors.add("generated-context-minimum-separation")
    return CandidateAudit(not errors, tuple(sorted(errors)), tuple(checked), planning.unassessed_context)


def tuple_floor(x: float, y: float, z: float) -> tuple[int, int, int]:
    return math.floor(x), math.floor(y), math.floor(z)
