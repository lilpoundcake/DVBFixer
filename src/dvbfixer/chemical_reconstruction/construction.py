"""Private authority-driven proposals; no registration, repair or optimizer."""
from __future__ import annotations

import math
from dataclasses import replace

from dvbfixer.domain.nonprotein_reconstruction import (
    ComponentClass,
    CoordinateSource,
    GeometryPolicy,
    ObservedAtom,
    UserMappedAuthority,
)
from dvbfixer.domain.nonprotein_reconstruction import (
    ReconstructionRequest as LegacyClassA,
)
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.protein_completion.geometry import extend, two_angle_tip

from .adapters import AdapterError, from_class_a
from .admission import admit
from .authority_geometry import geometry_view
from .candidate_audit import audit_candidate
from .contracts import (
    AdmissionStatus,
    Capability,
    ChemistryKind,
    ConstraintKind,
    LinkKind,
    Observation,
    Operation,
    ReconstructionRequest,
    Unit,
    canonical_graph_digest,
)
from .measurement import angle_degrees, satisfied, torsion_degrees
from .proposals import Candidate, GeometryPlan, PlacementKind, ProposalResult, Strategy


def _plan_reason(request: ReconstructionRequest, plan: GeometryPlan) -> str | None:
    if plan.request_digest != request.digest:
        return "plan-request-binding"
    if request.ownership.temporary or request.ownership.mobile:
        return "temporary-or-mobile-construction-not-implemented"
    generated = request.ownership.generated
    if len(generated) > 4096:
        return "construction-budget"
    if any({c.link.component_atom, c.link.other_atom} & generated for c in request.connections):
        return "generated-external-link-needs-linkage-geometry-strategy"
    if plan.strategy == Strategy.TEMPLATE:
        if plan.placements or len(request.targets) != 1:
            return "rigid-template-needs-one-target-and-no-chain-plan"
        if request.ph is None:
            return "rigid-template-needs-explicit-legacy-ph-metadata"
        target = next(iter(request.targets))
        component = next(c for c in request.components if c.identity == target)
        if component.kind != ChemistryKind.ISOLATED:
            return "rigid-template-isolated-class-a-only"
        if any(target in {c.link.component_atom.component, c.link.other_atom.component} for c in request.connections):
            return "rigid-template-target-must-be-isolated"
        state = component.resolved
        assert state is not None
        present = {o.identity for o in request.observations} | generated
        if any(ExactAtomRef(target, a.name) not in present for a in state.graph.atoms):
            return "rigid-template-would-generate-unowned-atoms"
    else:
        if len({p.atom for p in plan.placements}) != len(plan.placements) or {p.atom for p in plan.placements} != generated:
            return "plan-must-cover-exact-generated-mask-once"
        available = {o.identity for o in request.observations}
        _, bonds = geometry_view(request)
        for placement in plan.placements:
            if any(a.component != placement.atom.component for a in placement.anchors):
                return "cross-component-internal-frame-not-implemented"
            if not set(placement.anchors) <= available:
                return "plan-anchor-dependency"
            a, b, c = placement.anchors
            pairs = ((a, b), (b, c), (c, placement.atom)) if placement.kind == PlacementKind.CHAIN else (
                (a, b), (a, c), (a, placement.atom))
            if any(frozenset(pair) not in bonds for pair in pairs):
                return "construction-anchors-not-authoritative-bonds"
            if placement.kind == PlacementKind.CONE and not any(
                    q.kind == ConstraintKind.STEREO_FRAME and q.atoms == (*placement.anchors, placement.atom)
                    for q in request.constraints):
                return "cone-needs-explicit-ordered-stereo-frame"
            available.add(placement.atom)
    ideals, bonds = geometry_view(request)
    for pair, length in bonds.items():
        a, b = sorted(pair)
        if a.component in request.targets and abs(math.dist(ideals[a], ideals[b])-length) > plan.policy.bond_tolerance_angstrom:
            return "authority-ideal-bond-inconsistent"
    return None


def _internal(request: ReconstructionRequest, plan: GeometryPlan, additions: list[Observation]) -> None:
    ideals, bonds = geometry_view(request)
    positions = {o.identity: o.angstrom for o in request.observations}
    elements = {ExactAtomRef(c.identity, a.name): a.element for c in request.components if c.resolved is not None
                for a in c.resolved.graph.atoms}
    for placement in plan.placements:
        a, b, c = placement.anchors
        d = placement.atom
        if placement.kind == PlacementKind.CHAIN:
            length = bonds[frozenset((c, d))]
            angle = angle_degrees(ideals[b], ideals[c], ideals[d])
            chi = torsion_degrees((ideals[a], ideals[b], ideals[c], ideals[d]))
            point = extend(positions[a], positions[b], positions[c], length, angle, chi)
        else:
            frame = next(q for q in request.constraints if q.kind == ConstraintKind.STEREO_FRAME
                         and q.atoms == (*placement.anchors, d))
            if not satisfied(frame, tuple(ideals[x] for x in frame.atoms)):
                raise ValueError("authority-stereo-frame-conflict")
            assert frame.orientation is not None
            point = two_angle_tip(positions[a], positions[b], positions[c], length=bonds[frozenset((a, d))],
                                  first_angle=angle_degrees(ideals[b], ideals[a], ideals[d]),
                                  second_angle=angle_degrees(ideals[c], ideals[a], ideals[d]),
                                  orientation=frame.orientation)
        additions.append(Observation(d, elements[d], point, Unit.ANGSTROM))
        positions[d] = point


def _template(request: ReconstructionRequest, plan: GeometryPlan, additions: list[Observation]) -> None:
    # Reuse the existing owner, not a second Kabsch implementation or its legacy
    # geometry judgement/authority loader. This import is deliberately lazy.
    from dvbfixer.nonprotein_reconstruction import MAX_ATOMS, _fit_coordinates, _OutcomeError

    target = next(iter(request.targets))
    state = next(c.resolved for c in request.components if c.identity == target)
    assert state is not None
    if len(state.graph.atoms) > MAX_ATOMS:
        raise ValueError("legacy-class-a-authority-bounds")
    ideals, bonds = geometry_view(request)
    graph = replace(state.graph,
                    atoms=tuple(replace(a, ideal_position=ideals[ExactAtomRef(target, a.name)]) for a in state.graph.atoms),
                    bonds=tuple(replace(b, ideal_length=bonds[frozenset((ExactAtomRef(target, b.atom1),
                                           ExactAtomRef(target, b.atom2)))]) for b in state.graph.bonds))
    observations = tuple(ObservedAtom(o.identity, o.element, o.angstrom) for o in request.observations
                         if o.identity.component == target)
    assert request.ph is not None  # Checked transport metadata, never a state-selection heuristic.
    legacy = LegacyClassA(target, observations, UserMappedAuthority(graph, state.authority.source,
                          canonical_graph_digest(graph)), request.ph, ComponentClass.A,
                          microstate_override=graph.state_id,
                          geometry_policy=GeometryPolicy(max_anchor_rmsd=plan.policy.max_anchor_rmsd_angstrom))
    # Compatibility envelope check only: no source/context is deleted from the
    # original shared request used by both independent candidate audits below.
    from_class_a(legacy, source_sha256=request.source_sha256, operation=request.operation, hydrogens=state.hydrogens)
    try:
        coordinates, _, _ = _fit_coordinates(legacy, graph)
    except _OutcomeError as exc:
        raise ValueError(f"legacy-class-a-fit:{exc.finding.code}") from exc
    observed = {o.identity: o for o in observations}
    for coordinate in coordinates:
        if coordinate.identity in observed:
            if coordinate.source != CoordinateSource.OBSERVED or coordinate.position != observed[coordinate.identity].position:
                raise ValueError("legacy-fit-moved-observation")
        else:
            if coordinate.identity not in request.ownership.generated or coordinate.source != CoordinateSource.SOURCE_TEMPLATE:
                raise ValueError("legacy-fit-unowned-output")
            additions.append(Observation(coordinate.identity, coordinate.element, coordinate.position, Unit.ANGSTROM))


def propose(request: ReconstructionRequest, plan: GeometryPlan) -> ProposalResult:
    """One deterministic proposal; retain raw/rounded failures, never rank/repair."""
    # Element-agnostic coordinate math: this declaration is not an element-specific
    # scientific capability/parameter table or global backend registration.
    elements = {o.element.upper() for o in request.observations} | {
        a.element.upper() for c in request.components if c.resolved is not None for a in c.resolved.graph.atoms}
    if not elements:
        return ProposalResult("unsupported", admit(request), plan.digest, "request-not-admitted")
    capability = Capability(plan.strategy.value, frozenset({Operation.HEAVY, Operation.HYDROGENS}),
                            frozenset({ChemistryKind.PROTEIN, ChemistryKind.ISOLATED}), frozenset(elements),
                            frozenset(ConstraintKind), frozenset({LinkKind.COVALENT, LinkKind.COORDINATION}),
                            allow_opaque_context=True, minimum_anchors=3)
    admission = admit(request, (capability,))
    if admission.status != AdmissionStatus.ADMITTED:
        return ProposalResult("unsupported", admission, plan.digest, "request-not-admitted")
    reason = _plan_reason(request, plan)
    if reason is not None:
        return ProposalResult("invalid-plan", admission, plan.digest, reason)
    additions: list[Observation] = []
    try:
        if plan.strategy == Strategy.INTERNAL:
            _internal(request, plan, additions)
        else:
            _template(request, plan, additions)
    except AdapterError as exc:
        return ProposalResult("unsupported", admission, plan.digest, exc.code)
    except (ValueError, ImportError) as exc:
        raw = Candidate(request.digest, plan.digest, plan.strategy, tuple(additions))
        return ProposalResult("construction-failed", admission, plan.digest, str(exc), raw, audit_candidate(request, plan, raw))
    raw = Candidate(request.digest, plan.digest, plan.strategy, tuple(additions))
    raw_audit = audit_candidate(request, plan, raw)
    rounded_additions = []
    for addition in additions:
        fields = tuple(f"{x:8.3f}" for x in addition.position)
        if any(len(field) != 8 for field in fields):
            return ProposalResult("rejected", admission, plan.digest, "pdb-coordinate-range", raw, raw_audit)
        rounded_additions.append(replace(addition, position=(float(fields[0]), float(fields[1]), float(fields[2]))))
    rounded = replace(raw, additions=tuple(rounded_additions))
    rounded_audit = audit_candidate(request, plan, rounded)
    passing = raw_audit.passed and rounded_audit.passed
    return ProposalResult("declared-geometry-passed" if passing else "rejected", admission, plan.digest,
                          "declared-checks-only" if passing else "raw-or-rounded-geometry-failed",
                          raw, raw_audit, rounded, rounded_audit)
