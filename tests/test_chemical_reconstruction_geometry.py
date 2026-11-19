"""WP1 synthetic geometry controls; not native conformation/H qualification."""
from __future__ import annotations

import math
from dataclasses import replace

import pytest

from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemicalState,
    ChemistryKind,
    Component,
    Connection,
    Constraint,
    ConstraintKind,
    ContractError,
    HydrogenInventory,
    LinkKind,
    Observation,
    Operation,
    Ownership,
    ReconstructionRequest,
    Unit,
    canonical_graph_digest,
)
from dvbfixer.chemical_reconstruction.proposals import (
    AuditPolicy,
    GeometryPlan,
    Placement,
    PlacementKind,
    Strategy,
)
from dvbfixer.domain.nonprotein_reconstruction import (
    AuthorityMode,
    AuthorityProvenance,
    ChemicalAtom,
    ChemicalBond,
    ChemicalGraph,
    ExternalLink,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

TARGET = ComponentInstanceRef(1, "A", "82", "B", occurrence=2)
CONTEXT = ComponentInstanceRef(1, "a", "82", "B", occurrence=3)
POLICY = AuditPolicy(.002, .2, .2, .01)  # Synthetic controls, not calibrated chemistry.


def ref(name: str) -> ExactAtomRef:
    return ExactAtomRef(TARGET, name)


def chemical_state(g: ChemicalGraph, unit: Unit = Unit.ANGSTROM, *, h: bool = False) -> ChemicalState:
    return ChemicalState(g, AuthorityProvenance(AuthorityMode.USER_MAPPED, "synthetic geometry control",
                         canonical_graph_digest(g), "literal-coordinate-v1"), unit,
                         HydrogenInventory.EXPLICIT if h else HydrogenInventory.UNSPECIFIED)


def chain(*, h: bool = False, count: int = 4, last_z: float | None = None,
          kind: ChemistryKind = ChemistryKind.PROTEIN) -> tuple[ReconstructionRequest, GeometryPlan]:
    points = [(0., 0., 0.)]
    for i in range(count-1):
        p = list(points[-1])
        p[i % 3] += 1.
        points.append(tuple(p))
    if last_z is not None:
        points[-1] = (points[-1][0], points[-1][1], last_z)
    atoms = tuple(ChemicalAtom(f"X{i}", "H" if h and i == count-1 else "C", 0, point)
                  for i, point in enumerate(points))
    g = ChemicalGraph("synthetic", "explicit", atoms,
                     tuple(ChemicalBond(f"X{i}", f"X{i+1}", 1., math.dist(points[i], points[i+1]))
                           for i in range(count-1)), 0)
    r = ReconstructionRequest("a"*64, Operation.HYDROGENS if h else Operation.HEAVY,
                              (Component(TARGET, kind, (chemical_state(g, h=h),)),),
                              tuple(Observation(ref(a.name), a.element, a.ideal_position, Unit.ANGSTROM) for a in atoms[:3]),
                              frozenset({TARGET}), Ownership(frozenset(ref(a.name) for a in atoms[3:])),
                              ph=7. if kind == ChemistryKind.ISOLATED else None)
    plan = GeometryPlan(r.digest, Strategy.INTERNAL, POLICY,
                        tuple(Placement(ref(f"X{i}"), (ref(f"X{i-3}"), ref(f"X{i-2}"), ref(f"X{i-1}")))
                              for i in range(3, count)))
    return r, plan


def test_protein_internal_coordinate_proposal_preserves_observations() -> None:
    r, plan = chain()
    before = r.observations
    result = propose(r, plan)
    assert result.status == "declared-geometry-passed", result
    assert result.raw.additions[0].position == pytest.approx((1., 1., 1.))
    assert result.raw_audit.passed and result.rounded_audit.passed
    assert result.coordinates_produced and r.observations == before
    assert not result.geometry_approved and not result.md_ready and not result.source_coverage_verified


def test_twelve_missing_atoms_are_constructed_in_explicit_dependency_order() -> None:
    r, plan = chain(count=15)
    result = propose(r, plan)
    assert result.status == "declared-geometry-passed", result
    assert len(result.raw.additions) == 12
    assert result.raw.additions[-1].position == pytest.approx((5., 5., 4.))
    assert {a.identity for a in result.raw.additions} == r.ownership.generated


def test_explicit_h_placement_uses_same_authority_without_state_selection() -> None:
    r, plan = chain(h=True)
    result = propose(r, plan)
    assert result.status == "declared-geometry-passed", result
    assert result.raw.additions[0].element == "H"
    assert r.components[0].resolved.hydrogens == HydrogenInventory.EXPLICIT
    assert not result.md_ready  # No chemical H-count/CIP or physical qualification.


def test_class_a_rigid_fit_reuses_existing_owner(monkeypatch: pytest.MonkeyPatch) -> None:
    import dvbfixer.nonprotein_reconstruction as legacy

    calls = []
    original = legacy._fit_coordinates

    def tracked(*args):
        calls.append(args)
        return original(*args)

    monkeypatch.setattr(legacy, "_fit_coordinates", tracked)
    # Refuse any hidden chemistry resolver, validator or full legacy service.
    def forbidden(*args, **kwargs):
        raise AssertionError("unexpected legacy scientific/authority path")

    monkeypatch.setattr(legacy, "_resolve_authority", forbidden)
    monkeypatch.setattr(legacy, "_validate_geometry", forbidden)
    monkeypatch.setattr(legacy, "reconstruct_nonprotein", forbidden)
    r, _ = chain(kind=ChemistryKind.ISOLATED)
    def transform(p):
        return p[2]+7., p[0]-4., p[1]+2.
    r = replace(r, observations=tuple(replace(o, position=transform(o.position)) for o in r.observations))
    result = propose(r, GeometryPlan(r.digest, Strategy.TEMPLATE, POLICY))
    assert result.status == "declared-geometry-passed", result
    assert calls and result.raw.additions[0].position == pytest.approx(transform((1., 1., 1.)))


@pytest.mark.parametrize("strategy", [Strategy.INTERNAL, Strategy.TEMPLATE])
def test_units_are_converted_as_views_not_source_mutations(strategy: Strategy) -> None:
    r, plan = chain(kind=ChemistryKind.ISOLATED)
    g = r.components[0].resolved.graph
    g = replace(g, atoms=tuple(replace(a, ideal_position=tuple(x/10 for x in a.ideal_position)) for a in g.atoms),
                bonds=tuple(replace(b, ideal_length=b.ideal_length/10) for b in g.bonds))
    r = replace(r, components=(replace(r.components[0], states=(chemical_state(g, Unit.NANOMETER),)),),
                observations=tuple(replace(o, unit=Unit.NANOMETER, position=tuple(x/10 for x in o.position))
                                   for o in r.observations))
    before = r.observations
    plan = GeometryPlan(r.digest, strategy, POLICY, plan.placements if strategy == Strategy.INTERNAL else ())
    result = propose(r, plan)
    assert result.status == "declared-geometry-passed", result
    assert result.raw.additions[0].position == pytest.approx((1., 1., 1.))
    assert r.observations == before and r.components[0].resolved.ideal_unit == Unit.NANOMETER


def cone() -> tuple[ReconstructionRequest, GeometryPlan]:
    atoms = (ChemicalAtom("CA", "C", 0, (0., 0., 0.), stereo="synthetic-ordered-positive"),
             ChemicalAtom("N", "N", 0, (1., 0., 0.)), ChemicalAtom("C", "C", 0, (0., 1., 0.)),
             ChemicalAtom("CB", "C", 0, (0., 0., 1.)))
    g = ChemicalGraph("synthetic-cone", "one", atoms,
                     tuple(ChemicalBond("CA", a.name, 1., 1.) for a in atoms[1:]), 0)
    frame = Constraint(ConstraintKind.STEREO_FRAME, tuple(ref(n) for n in ("CA", "N", "C", "CB")),
                       Unit.DIMENSIONLESS, orientation=1, minimum_volume=.1, label="synthetic-ordered-positive")
    r = ReconstructionRequest("a"*64, Operation.HEAVY, (Component(TARGET, ChemistryKind.PROTEIN, (chemical_state(g),)),),
                              tuple(Observation(ref(a.name), a.element, a.ideal_position, Unit.ANGSTROM) for a in atoms[:3]),
                              frozenset({TARGET}), Ownership(frozenset({ref("CB")})), constraints=(frame,))
    return r, GeometryPlan(r.digest, Strategy.INTERNAL, POLICY,
                           (Placement(ref("CB"), (ref("CA"), ref("N"), ref("C")), PlacementKind.CONE),))


def test_ordered_protein_cone_and_inversion_rejection() -> None:
    r, plan = cone()
    result = propose(r, plan)
    assert result.status == "declared-geometry-passed", result
    bad = replace(result.raw, additions=(replace(result.raw.additions[0], position=(0., 0., -1.)),))
    audit = audit_candidate(r, plan, bad)
    assert not audit.passed and "candidate-constraint-conflict" in audit.findings


@pytest.mark.parametrize("strategy", [Strategy.INTERNAL, Strategy.TEMPLATE])
def test_whole_complex_context_clash_is_retained_and_rejected(strategy: Strategy) -> None:
    r, plan = chain(kind=ChemistryKind.ISOLATED)
    partner = Observation(ExactAtomRef(CONTEXT, "ZN"), "ZN", (1., 1., 1.), Unit.ANGSTROM)
    r = replace(r, components=(*r.components, Component(CONTEXT, ChemistryKind.UNKNOWN)),
                observations=(*r.observations, partner))
    plan = GeometryPlan(r.digest, strategy, POLICY, plan.placements if strategy == Strategy.INTERNAL else ())
    result = propose(r, plan)
    assert result.status == "rejected" and result.raw is not None and result.rounded is not None
    assert "generated-context-minimum-separation" in result.raw_audit.findings
    assert result.raw_audit.unassessed_context == (CONTEXT,)
    assert r.observations[-1] == partner


def test_rounding_is_independent_gate_not_a_repair() -> None:
    r, plan = chain(last_z=1.0004)
    plan = replace(plan, policy=replace(POLICY, bond_tolerance_angstrom=.00001))
    result = propose(r, plan)
    assert result.status == "rejected" and result.raw_audit.passed
    assert not result.rounded_audit.passed
    assert result.raw.additions[0].position[2] == pytest.approx(1.0004)
    assert result.rounded.additions[0].position[2] == 1.


@pytest.mark.parametrize("change", ["binding", "dependency", "duplicate", "bond", "temporary", "linked"])
def test_invalid_plans_do_not_build_coordinates(change: str) -> None:
    r, plan = chain(count=5)
    if change == "binding":
        plan = replace(plan, request_digest="b"*64)
    elif change == "dependency":
        plan = replace(plan, placements=tuple(reversed(plan.placements)))
    elif change == "duplicate":
        plan = replace(plan, placements=(plan.placements[0], plan.placements[0]))
    elif change == "bond":
        plan = replace(plan, placements=(replace(plan.placements[0], anchors=(ref("X0"), ref("X2"), ref("X1"))),
                                        plan.placements[1]))
    elif change == "temporary":
        r = replace(r, ownership=Ownership(frozenset({ref("X3")}), temporary=frozenset({ref("X4")})))
        plan = replace(plan, request_digest=r.digest)
    else:
        partner = Observation(ExactAtomRef(CONTEXT, "C"), "C", (10., 10., 10.), Unit.ANGSTROM)
        r = replace(r, components=(*r.components, Component(CONTEXT, ChemistryKind.UNKNOWN)),
                    observations=(*r.observations, partner),
                    connections=(Connection(ExternalLink(ref("X4"), partner.identity, 1.), LinkKind.COVALENT),))
        plan = replace(plan, request_digest=r.digest)
    result = propose(r, plan)
    assert result.status in {"invalid-plan", "unsupported"} and not result.coordinates_produced


def test_degenerate_observed_anchors_refuse_without_reflecting_source() -> None:
    r, plan = chain()
    r = replace(r, observations=tuple(replace(o, position=(float(i), 0., 0.)) for i, o in enumerate(r.observations)))
    plan = replace(plan, request_digest=r.digest)
    before = r.observations
    result = propose(r, plan)
    assert result.status == "construction-failed" and result.raw is not None
    assert not result.raw_audit.passed and r.observations == before


@pytest.mark.parametrize("change", ["binding", "element", "mask", "units", "occupancy", "duplicate"])
def test_independent_auditor_refuses_forged_candidates(change: str) -> None:
    r, plan = chain()
    raw = propose(r, plan).raw
    if change == "binding":
        raw = replace(raw, request_digest="b"*64)
    elif change == "element":
        raw = replace(raw, additions=(replace(raw.additions[0], element="O"),))
    elif change == "mask":
        raw = replace(raw, additions=())
    elif change == "units":
        raw = replace(raw, additions=(replace(raw.additions[0], unit=Unit.NANOMETER),))
    elif change == "occupancy":
        raw = replace(raw, additions=(replace(raw.additions[0], occupancy=.5),))
    else:
        raw = replace(raw, additions=raw.additions*2)
    assert not audit_candidate(r, plan, raw).passed


def test_validation_does_not_call_construction(monkeypatch: pytest.MonkeyPatch) -> None:
    import dvbfixer.chemical_reconstruction.construction as construction

    r, plan = chain()
    raw = propose(r, plan).raw
    def forbidden(*args, **kwargs):
        raise AssertionError("validator called constructor")
    monkeypatch.setattr(construction, "extend", forbidden)
    monkeypatch.setattr(construction, "two_angle_tip", forbidden)
    monkeypatch.setattr(construction, "_template", forbidden)
    assert audit_candidate(r, plan, raw).passed


def test_policy_and_plan_hashes_bind_geometry_authority_and_bounds() -> None:
    r, plan = chain()
    assert replace(plan, policy=replace(POLICY, angle_tolerance_degrees=.1)).digest != plan.digest
    assert replace(plan, placements=()).digest != plan.digest
    with pytest.raises(ContractError):
        replace(POLICY, minimum_separation_angstrom=1e-310)
    with pytest.raises(ContractError):
        replace(POLICY, bond_tolerance_angstrom=float("nan"))


def test_authority_inconsistency_is_not_fixed_by_rebuilding() -> None:
    r, plan = chain()
    g = r.components[0].resolved.graph
    g = replace(g, bonds=(*g.bonds[:-1], replace(g.bonds[-1], ideal_length=2.)))
    r = replace(r, components=(replace(r.components[0], states=(chemical_state(g),)),))
    result = propose(r, replace(plan, request_digest=r.digest))
    assert result.status == "invalid-plan" and result.reason == "authority-ideal-bond-inconsistent"
    assert not result.coordinates_produced


def test_unknown_empty_context_is_a_typed_refusal() -> None:
    r, _ = chain()
    r = replace(r, components=(Component(TARGET, ChemistryKind.UNKNOWN),), observations=())
    result = propose(r, GeometryPlan(r.digest, Strategy.INTERNAL, POLICY))
    assert result.status == "unsupported" and not result.coordinates_produced


def test_class_a_cannot_materialize_unowned_hydrogen() -> None:
    r, _ = chain(count=5, h=True, kind=ChemistryKind.ISOLATED)
    r = replace(r, operation=Operation.HEAVY, ownership=Ownership(frozenset({ref("X3")})))
    result = propose(r, GeometryPlan(r.digest, Strategy.TEMPLATE, POLICY))
    assert result.reason == "rigid-template-would-generate-unowned-atoms" and not result.coordinates_produced


def test_class_a_does_not_invent_legacy_ph_metadata() -> None:
    r, _ = chain(kind=ChemistryKind.ISOLATED)
    r = replace(r, ph=None)
    result = propose(r, GeometryPlan(r.digest, Strategy.TEMPLATE, POLICY))
    assert result.reason == "rigid-template-needs-explicit-legacy-ph-metadata" and not result.coordinates_produced


def test_class_a_adapter_refusal_is_not_retried_as_internal_geometry() -> None:
    r, _ = chain(kind=ChemistryKind.ISOLATED)
    g = r.components[0].resolved.graph
    g = replace(g, atoms=(replace(g.atoms[0], element="S"), *g.atoms[1:]))
    r = replace(r, components=(replace(r.components[0], states=(chemical_state(g),)),),
                observations=(replace(r.observations[0], element="S"), *r.observations[1:]))
    result = propose(r, GeometryPlan(r.digest, Strategy.TEMPLATE, POLICY))
    assert result.status == "unsupported" and result.reason == "outside-initial-chno-class-a-adapter-envelope"
    assert not result.coordinates_produced


def test_cone_does_not_relabel_or_reflect_inconsistent_authority() -> None:
    r, plan = cone()
    r = replace(r, constraints=(replace(r.constraints[0], orientation=-1),))
    before = r.observations
    result = propose(r, replace(plan, request_digest=r.digest))
    assert result.status == "construction-failed" and result.reason == "authority-stereo-frame-conflict"
    assert r.observations == before


def test_partial_proposal_survives_constructor_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    import dvbfixer.chemical_reconstruction.construction as construction

    r, plan = chain(count=5)
    original = construction.extend
    calls = 0
    def fail_second(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise ValueError("injected second-step failure")
        return original(*args, **kwargs)
    monkeypatch.setattr(construction, "extend", fail_second)
    result = propose(r, plan)
    assert result.status == "construction-failed" and len(result.raw.additions) == 1
    assert not result.raw_audit.passed and "generated-mask" in result.raw_audit.findings


def test_false_constructor_success_is_rejected_without_repair(monkeypatch: pytest.MonkeyPatch) -> None:
    import dvbfixer.chemical_reconstruction.construction as construction

    r, plan = chain()
    monkeypatch.setattr(construction, "extend", lambda *args: (1., 1., 9.))
    result = propose(r, plan)
    assert result.status == "rejected" and result.raw.additions[0].position == (1., 1., 9.)
    assert "authority-bond-geometry" in result.raw_audit.findings
    assert result.rounded.additions[0].position == (1., 1., 9.)


def test_pdb_coordinate_range_refuses_without_normalizing_source() -> None:
    r, plan = chain()
    r = replace(r, observations=tuple(replace(o, position=tuple(x+10000. for x in o.position)) for o in r.observations))
    before = r.observations
    result = propose(r, replace(plan, request_digest=r.digest))
    assert result.status == "rejected" and result.reason == "pdb-coordinate-range"
    assert result.raw_audit.passed and result.rounded is None and r.observations == before
