"""Private complete-H group planning; no measured-H replacement or phase search."""
from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations

from dvbfixer.domain.structure_identity import ExactAtomRef

from .amide_authority import AmideMetricAuthority
from .ccd_authority import verify_cip
from .contracts import (
    Constraint,
    ConstraintKind,
    Operation,
    Ownership,
    ReconstructionRequest,
    Unit,
    canonical_graph_digest,
    digest_record,
)
from .joint_hydrogen_geometry import (
    DirectionGroup,
    JointMetricGeometryPlan,
    JointSourceMetricGeometryPlan,
)
from .measurement import _cross, _dot, _sub, _unit, angle_degrees, torsion_degrees
from .metric_geometry import GeneratedMetricGeometryPlan, GeometryMetric, SourceMetricGeometryPlan
from .peptide_state import LocalPeptideAuthority, PeptideStateError
from .proposals import AuditPolicy, Candidate, Placement, PlacementKind, Strategy


@dataclass(frozen=True)
class HydrogenGroup:
    parent: ExactAtomRef
    kind: str
    hydrogens: tuple[ExactAtomRef, ...]
    heavy_neighbours: tuple[ExactAtomRef, ...]
    anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]


@dataclass(frozen=True)
class CompleteHydrogenView:
    request: ReconstructionRequest
    plan: GeneratedMetricGeometryPlan
    groups: tuple[HydrogenGroup, ...]


def compile_complete_hydrogens(base: ReconstructionRequest, authority: LocalPeptideAuthority,
                                incoming_carbon: ExactAtomRef, parameters: AmideMetricAuthority, *,
                                backbone_types: tuple[str, str, str], fixed_policy: AuditPolicy,
                                generated_policy: AuditPolicy) -> CompleteHydrogenView:
    return _compile_hydrogens(base, authority, incoming_carbon, parameters, backbone_types=backbone_types,
                              fixed_policy=fixed_policy, generated_policy=generated_policy, extended=False)


def compile_extended_hydrogens(base: ReconstructionRequest, authority: LocalPeptideAuthority,
                                incoming_carbon: ExactAtomRef, parameters: AmideMetricAuthority, *,
                                backbone_types: tuple[str, str, str] | None, fixed_policy: AuditPolicy,
                                generated_policy: AuditPolicy) -> CompleteHydrogenView:
    """Separate opt-in: authority-proven planar N and neutral thiol, not new states."""
    return _compile_hydrogens(base, authority, incoming_carbon, parameters, backbone_types=backbone_types,
                              fixed_policy=fixed_policy, generated_policy=generated_policy, extended=True)


def _compile_hydrogens(base: ReconstructionRequest, authority: LocalPeptideAuthority,
                       incoming_carbon: ExactAtomRef, parameters: AmideMetricAuthority, *,
                       backbone_types: tuple[str, str, str] | None, fixed_policy: AuditPolicy,
                       generated_policy: AuditPolicy, extended: bool) -> CompleteHydrogenView:
    target = authority.identity
    component = next((c for c in base.components if c.identity == target), None)
    if (base.targets != {target} or component is None or component.resolved != authority.state
            or base.operation != Operation.HEAVY or base.ownership.generated or base.ownership.mobile
            or base.ownership.temporary):
        raise PeptideStateError("complete-H-needs-exact-authority-and-fixed-complete-heavy-input")
    graph = authority.state.graph
    atoms = {a.name: a for a in graph.atoms}
    observed = {o.identity: o for o in base.observations}
    selected = {a.name: o for a, o in observed.items() if a.component == target}
    if any(o.element.upper() in {"H", "D"} for o in selected.values()):
        raise PeptideStateError("measured-target-H-transport-unqualified")
    heavy = {a.name for a in graph.atoms if a.element not in {"H", "D"}}
    if set(selected) != heavy or any(o.occupancy != 1 for o in selected.values()):
        raise PeptideStateError("complete-H-needs-exact-unit-occupancy-heavy-inventory")
    if incoming_carbon not in observed or observed[incoming_carbon].element != "C":
        raise PeptideStateError("complete-H-needs-real-observed-incoming-C")
    adjacent: dict[str, set[str]] = {name: set() for name in atoms}
    bonds = {}
    for bond in graph.bonds:
        adjacent[bond.atom1].add(bond.atom2)
        adjacent[bond.atom2].add(bond.atom1)
        bonds[frozenset((bond.atom1, bond.atom2))] = bond
    hydrogen_names = {a.name for a in graph.atoms if a.element == "H"}
    if not hydrogen_names or any(atoms[n].element == "D" for n in atoms):
        raise PeptideStateError("explicit-protium-inventory-required")
    if any(len(adjacent[h]) != 1 or not adjacent[h] <= heavy for h in hydrogen_names):
        raise PeptideStateError("hydrogen-parent-connectivity-unqualified")
    def ref(name: str) -> ExactAtomRef:
        return ExactAtomRef(target, name)

    positions = {a.name: a.ideal_position for a in graph.atoms}
    constraints = list(base.constraints)
    for q in authority.constraints(minimum_volume=.05, planar_tolerance_angstrom=.05, include_hydrogens=True):
        if q.kind == ConstraintKind.PLANARITY and any(atoms[a.name].element == "H" for a in q.atoms):
            q = replace(q, upper=.005)
        constraints.append(q)
    placements: list[Placement] = []
    metrics: list[GeometryMetric] = []
    groups: list[HydrogenGroup] = []
    for parent in sorted(heavy):
        hs = tuple(sorted(adjacent[parent] & hydrogen_names))
        if not hs:
            continue
        neighbours = tuple(sorted(adjacent[parent] & heavy))
        element = atoms[parent].element
        single = all(bonds[frozenset((parent, n))].order == 1 for n in adjacent[parent])
        planar_n = (extended and element == "N" and atoms[parent].formal_charge in {0, 1}
                    and len(adjacent[parent]) == 3
                    and any(set(p) == {parent, *adjacent[parent]} for p in authority.parent.planar_groups)
                    and all(bonds[frozenset((parent, h))].order == 1 for h in hs))
        if parent == "N" and hs == ("H",) and neighbours == ("CA",):
            if backbone_types is None:
                raise PeptideStateError("backbone-H-requires-explicit-parameter-types")
            anchors = (ref("N"), incoming_carbon, ref("CA"))
            metrics.extend(parameters.metrics(ref("N"), ref("H"), ref("CA"), incoming_carbon, atom_types=backbone_types))
            placements.append(Placement(ref("H"), anchors, PlacementKind.PLANAR))
            constraints.append(Constraint(ConstraintKind.PLANARITY, (*anchors, ref("H")), Unit.ANGSTROM, upper=.005))
            kind = "backbone-amide"
        elif ((element == "C" and atoms[parent].aromatic or planar_n)
              and len(neighbours) == 2 and len(hs) == 1):
            anchors = (ref(parent), ref(neighbours[0]), ref(neighbours[1]))
            h = hs[0]
            source = (f"CCD ideal local planar N-H; equal-priority geometric compromise; {authority.parent.digest}" if planar_n else
                      f"CCD ideal local aromatic H; equal-priority geometric compromise; {authority.parent.digest}")
            pin = authority.parent.component.provenance.sha256
            metrics.extend((GeometryMetric((ref(parent), ref(h)), bonds[frozenset((parent, h))].ideal_length, Unit.ANGSTROM, source, pin),
                            *(GeometryMetric((ref(n), ref(parent), ref(h)), angle_degrees(positions[n], positions[parent], positions[h]),
                                             Unit.DEGREE, source, pin) for n in neighbours)))
            placements.append(Placement(ref(h), anchors, PlacementKind.PLANAR))
            constraints.append(Constraint(ConstraintKind.PLANARITY, (*anchors, ref(h)), Unit.ANGSTROM, upper=.005))
            kind = "planar-NH" if planar_n else "aromatic-CH"
        elif element == "C" and single and (len(neighbours), len(hs)) in {(3, 1), (2, 2)}:
            anchors = (ref(parent), ref(neighbours[0]), ref(neighbours[1]))
            for h in hs:
                directions = tuple(_unit(_sub(positions[n], positions[parent])) for n in (*neighbours[:2], h))
                volume = _dot(directions[0], _cross(directions[1], directions[2]))
                if abs(volume) < .05:
                    raise PeptideStateError("degenerate-ideal-H-frame")
                constraints.append(Constraint(ConstraintKind.STEREO_FRAME, (*anchors, ref(h)), Unit.DIMENSIONLESS,
                                              orientation=1 if volume > 0 else -1, minimum_volume=.05,
                                              label=atoms[parent].stereo or "ordered-CCD-H-frame-not-CIP"))
                placements.append(Placement(ref(h), anchors, PlacementKind.CONE))
            kind = "tetrahedral-CH" if len(hs) == 1 else "methylene"
        elif len(neighbours) == 1 and (single and (
                len(hs) == 3 and (element == "C" or element == "N" and atoms[parent].formal_charge == 1)
                or len(hs) == 1 and (element == "O" or extended and element == "S"
                                       and atoms[parent].formal_charge == 0 and not atoms[parent].aromatic))
                or planar_n and len(hs) == 2):
            b = neighbours[0]
            frame = None
            for a in sorted((adjacent[b] & heavy)-{parent}):
                try:
                    for h in hs:
                        torsion_degrees((positions[a], positions[b], positions[parent], positions[h]))
                except ValueError:
                    continue
                frame = (ref(a), ref(b), ref(parent))
                break
            if frame is None:
                raise PeptideStateError("no-ideal-heavy-chain-H-frame")
            anchors = frame
            placements.extend(Placement(ref(h), anchors) for h in hs)
            if planar_n:
                kind = "planar-NH2"
                constraints.append(Constraint(ConstraintKind.PLANARITY, (*anchors, *(ref(h) for h in hs)),
                                              Unit.ANGSTROM, upper=.005))
            else:
                kind = "terminal-three-H" if len(hs) == 3 else "thiol" if element == "S" else "hydroxyl"
        else:
            raise PeptideStateError(f"hydrogen-group-unqualified:{parent}")
        if not set(anchors) <= observed.keys():
            raise PeptideStateError("H-frame-needs-observed-heavy-anchors")
        groups.append(HydrogenGroup(ref(parent), kind, tuple(ref(h) for h in hs), tuple(ref(n) for n in neighbours), anchors))
    if {p.atom for p in placements} != {ref(h) for h in hydrogen_names} or len(placements) != len(hydrogen_names):
        raise PeptideStateError("complete-H-plan-mask-mismatch")
    source_only = extended and not metrics
    plan: GeneratedMetricGeometryPlan
    if source_only:
        # No backbone H means no amide override. Keep the nonempty metric-plan
        # contract by carrying complete explicit CCD H metrics, not a dummy row.
        pin = canonical_graph_digest(graph)
        source = f"Exact derived canonical graph {pin}; CCD {authority.parent.component.provenance.sha256}; parent {authority.parent.digest}"
        for bond in graph.bonds:
            if {bond.atom1, bond.atom2} & hydrogen_names:
                metrics.append(GeometryMetric((ref(bond.atom1), ref(bond.atom2)), bond.ideal_length, Unit.ANGSTROM, source, pin))
        for centre in sorted(adjacent):
            for first, last in combinations(sorted(adjacent[centre]), 2):
                if {first, centre, last} & hydrogen_names:
                    metrics.append(GeometryMetric((ref(first), ref(centre), ref(last)),
                                    angle_degrees(positions[first], positions[centre], positions[last]), Unit.DEGREE, source, pin))
    origin = base.origin
    if origin is not None:
        origin = replace(origin, settings=(*origin.settings,
                          ("complete_H_policy", "complete-H-groups-v2-explicit-sp2-N-thiol" if extended else "complete-H-groups-v1-fixed-heavy"),
                         ("complete_H_authority", digest_record((authority.digest, parameters.digest, backbone_types)))))
    request = replace(base, operation=Operation.HYDROGENS, ownership=Ownership(frozenset(ref(h) for h in hydrogen_names)),
                      constraints=tuple(dict.fromkeys(constraints)), origin=origin)
    if source_only:
        plan = SourceMetricGeometryPlan(request.digest, Strategy.INTERNAL, fixed_policy, tuple(placements), tuple(metrics),
                                         generated_policy=generated_policy, authority_graph_digest=canonical_graph_digest(graph))
    else:
        plan = GeneratedMetricGeometryPlan(request.digest, Strategy.INTERNAL, fixed_policy, tuple(placements), tuple(metrics),
                                           generated_policy=generated_policy)
    return CompleteHydrogenView(request, plan, tuple(groups))


def compile_joint_hydrogens(view: CompleteHydrogenView) -> CompleteHydrogenView:
    """Upgrade only admitted sp3 carbon cone groups; source/constraints stay exact."""
    groups = tuple(DirectionGroup(g.parent, g.heavy_neighbours, g.hydrogens) for g in view.groups
                   if g.kind in {"tetrahedral-CH", "methylene"})
    selected = {h for g in groups for h in g.hydrogens}
    if not groups or any(p.kind != PlacementKind.CONE for p in view.plan.placements if p.atom in selected):
        raise PeptideStateError("joint-H-upgrade-needs-original-complete-carbon-cone-groups")
    placements = tuple(replace(p, kind=PlacementKind.JOINT) if p.atom in selected else p for p in view.plan.placements)
    plan: JointMetricGeometryPlan
    if isinstance(view.plan, SourceMetricGeometryPlan):
        plan = JointSourceMetricGeometryPlan(view.request.digest, view.plan.strategy, view.plan.policy, placements, view.plan.metrics,
                                             generated_policy=view.plan.generated_policy, joint_groups=groups,
                                             authority_graph_digest=view.plan.authority_graph_digest)
    else:
        plan = JointMetricGeometryPlan(view.request.digest, view.plan.strategy, view.plan.policy, placements, view.plan.metrics,
                                       generated_policy=view.plan.generated_policy, joint_groups=groups)
    return replace(view, plan=plan)


@dataclass(frozen=True)
class ConditionalHydrogenCipAudit:
    passed: bool
    findings: tuple[str, ...]
    target_labels: tuple[tuple[str, str], ...]
    scope: str = "conditional target CIP with virtual ALA chemistry, not actual neighbour/full-polymer qualification"


def audit_hydrogen_cip(view: CompleteHydrogenView, authority: LocalPeptideAuthority,
                       candidate: Candidate) -> ConditionalHydrogenCipAudit:
    r = view.request
    component = next((c for c in r.components if c.identity == authority.identity), None)
    if (candidate.request_digest != r.digest or candidate.plan_digest != view.plan.digest
            or view.plan.request_digest != r.digest or component is None or component.resolved != authority.state
            or candidate.strategy != view.plan.strategy or r.targets != {authority.identity}
            or {a.identity for a in candidate.additions} != r.ownership.generated
            or len(candidate.additions) != len(r.ownership.generated)):
        return ConditionalHydrogenCipAudit(False, ("hydrogen-CIP-binding-or-mask",), ())
    points = {o.identity: o.angstrom for o in r.observations}
    if any(o.element != "H" or o.unit != Unit.ANGSTROM or o.occupancy != 1 for o in candidate.additions):
        return ConditionalHydrogenCipAudit(False, ("hydrogen-CIP-element-or-unit",), ())
    points.update({o.identity: o.position for o in candidate.additions})
    flat = authority.witness.flat_fragments
    positions = {a.name: a.ideal_position for a in flat.atoms}
    mapping = dict(authority.witness.node_map)
    for name, atom in mapping.items():
        if atom.component == authority.identity:
            if atom not in points:
                return ConditionalHydrogenCipAudit(False, ("hydrogen-CIP-incomplete-target",), ())
            positions[name] = points[atom]
    try:
        labels = verify_cip(flat, positions, extra_bonds=authority.witness.extra_bonds)
    except ValueError as exc:
        return ConditionalHydrogenCipAudit(False, (str(exc),), ())
    return ConditionalHydrogenCipAudit(True, (), tuple((mapping[name].name, label) for name, label in labels
                                                     if mapping[name].component == authority.identity))
