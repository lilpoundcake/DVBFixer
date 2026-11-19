"""Synthetic WP0 planning controls, not reconstruction/MD qualification."""
from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.chemical_reconstruction.adapters import AdapterError, from_class_a, from_protein
from dvbfixer.chemical_reconstruction.admission import admit
from dvbfixer.chemical_reconstruction.contracts import (
    AdmissionStatus,
    Capability,
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
from dvbfixer.domain.nonprotein_reconstruction import (
    AuthorityMode,
    AuthorityProvenance,
    ChemicalAtom,
    ChemicalBond,
    ChemicalGraph,
    ComponentClass,
    ExternalLink,
    GeometryPolicy,
    LocalPinnedCcdAuthority,
    ObservedAtom,
    OnlineEnrichmentAuthority,
    UserMappedAuthority,
)
from dvbfixer.domain.nonprotein_reconstruction import (
    ReconstructionRequest as LegacyClassA,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.contracts import CompletionRequest
from dvbfixer.protein_completion.contracts import Observation as ProteinObservation

SOURCE = "a"*64
TARGET = ComponentInstanceRef(1, "A", "82", "A", occurrence=2)
PARTNER = ComponentInstanceRef(1, "a", "82", occurrence=3)


def atom(name: str, component: ComponentInstanceRef = TARGET) -> ExactAtomRef:
    return ExactAtomRef(component, name)


def graph(*, state_id: str = "explicit", stereo: str | None = None, hydrogen: bool = False) -> ChemicalGraph:
    atoms = (ChemicalAtom("C", "C", 0, (0., 0., 0.), stereo=stereo),
             ChemicalAtom("N", "N", 0, (1., 0., 0.)),
             ChemicalAtom("O", "O", 0, (0., 1., 0.)),
             ChemicalAtom("Q", "C", 0, (0., 0., 1.)))
    if hydrogen:
        atoms += (ChemicalAtom("H", "H", 0, (-1., -1., -1.)),)
    return ChemicalGraph("synthetic", state_id, atoms,
                         tuple(ChemicalBond("C", a.name, 1., 1.) for a in atoms[1:]), 0)


def state(g: ChemicalGraph | None = None, *, hydrogens: HydrogenInventory = HydrogenInventory.UNSPECIFIED) -> ChemicalState:
    g = graph() if g is None else g
    return ChemicalState(g, AuthorityProvenance(AuthorityMode.USER_MAPPED, "synthetic—not licensed tables",
                         canonical_graph_digest(g), "v1"), Unit.ANGSTROM, hydrogens)


def request() -> ReconstructionRequest:
    g = graph()
    observations = tuple(Observation(atom(a.name), a.element, a.ideal_position, Unit.ANGSTROM) for a in g.atoms[:-1])
    return ReconstructionRequest(SOURCE, Operation.HEAVY,
                                 (Component(TARGET, ChemistryKind.PROTEIN, (state(g),)),), observations,
                                 frozenset({TARGET}), Ownership(frozenset({atom("Q")})))


def capability(**overrides: object) -> Capability:
    default = Capability("synthetic-declaration-only-v1", frozenset(Operation), frozenset(ChemistryKind),
                         frozenset({"C", "N", "O", "H", "ZN", "CL"}), frozenset(ConstraintKind),
                         frozenset({LinkKind.COVALENT, LinkKind.COORDINATION}), allow_opaque_context=True,
                         allow_mobile=True, allow_temporary=True)
    return replace(default, **overrides)


def codes(r: ReconstructionRequest) -> set[str]:
    return {f.code for f in admit(r, (capability(),)).findings}


def test_planning_never_claims_generation_geometry_or_source_coverage() -> None:
    r = request()
    a = admit(r, (capability(),))
    assert a.status == AdmissionStatus.ADMITTED
    assert a.request_digest == r.digest and a.capability_id == capability().identifier
    assert not any((a.coordinates_produced, a.geometry_approved, a.md_ready, a.source_coverage_verified))
    assert admit(r).status == AdmissionStatus.UNSUPPORTED


def test_canonical_digest_matches_existing_class_a() -> None:
    from dvbfixer.nonprotein_reconstruction import chemical_graph_digest

    g = replace(graph(), component_id="α-component")
    assert canonical_graph_digest(g) == chemical_graph_digest(g)


@pytest.mark.parametrize("change", ["charge", "fragments", "duplicate", "bond", "mutable", "nan"])
def test_graph_shapes_fail_closed(change: str) -> None:
    g = graph()
    if change == "charge":
        g = replace(g, net_charge=1)
    elif change == "fragments":
        g = replace(g, component_count=2)
    elif change == "duplicate":
        g = replace(g, atoms=(*g.atoms, g.atoms[0]))
    elif change == "bond":
        g = replace(g, bonds=(*g.bonds, ChemicalBond("N", "missing", 1., 1.)))
    elif change == "mutable":
        g = replace(g, atoms=list(g.atoms))
    else:
        g = replace(g, atoms=(replace(g.atoms[0], ideal_position=(float("nan"), 0., 0.)), *g.atoms[1:]))
    with pytest.raises(ContractError):
        state(g)


def test_forged_graph_authority_refuses() -> None:
    s = state()
    with pytest.raises(ContractError, match="digest mismatch"):
        replace(s, authority=replace(s.authority, source_digest="0"*64))
    with pytest.raises(ContractError, match="digest mismatch"):
        replace(s, graph=replace(s.graph, state_id="different"))


@pytest.mark.parametrize("part", ["source", "units", "occupancy", "context", "origin-state", "identity"])
def test_digest_binds_source_units_context_identity_and_chemistry(part: str) -> None:
    r = request()
    if part == "source":
        changed = replace(r, source_sha256="b"*64)
    elif part == "units":
        changed = replace(r, observations=(replace(r.observations[0], unit=Unit.NANOMETER), *r.observations[1:]))
    elif part == "occupancy":
        changed = replace(r, observations=(replace(r.observations[0], occupancy=.5), *r.observations[1:]))
    elif part == "context":
        changed = replace(r, components=(*r.components, Component(PARTNER, ChemistryKind.UNKNOWN)))
    elif part == "origin-state":
        changed = replace(r, components=(replace(r.components[0], states=(state(graph(state_id="second")),)),))
    else:
        other = replace(TARGET, chain_id="a", insertion_code="B", occurrence=3)
        changed = replace(r, targets=frozenset({other}))
    assert r.digest != changed.digest
    assert request().digest == r.digest


def test_context_and_exact_link_transport_from_protein() -> None:
    r = request()
    partner = ProteinObservation(atom("ZN", PARTNER), "ZN", (8., 9., 10.))
    old = CompletionRequest(SOURCE, TARGET, "ALA", tuple(ProteinObservation(o.identity, o.element, o.position)
                            for o in r.observations)+(partner,), ("Q",), ((atom("Q"), partner.atom),))
    adapted = from_protein(old, state=state())
    assert tuple(o.identity for o in adapted.observations) == tuple(o.atom for o in old.observations)
    assert tuple(o.position for o in adapted.observations) == tuple(o.xyz for o in old.observations)
    assert adapted.components[1].identity == PARTNER
    assert adapted.origin is not None and adapted.origin.request_digest == old.digest
    assert adapted.ownership.generated == frozenset({atom("Q")})
    assert "unresolved-connection-chemistry" in codes(adapted)
    link = Connection(ExternalLink(atom("Q"), partner.atom, 1.), LinkKind.COVALENT)
    resolved = from_protein(old, state=state(), connections=(link,))
    admission = admit(resolved, (capability(),))
    assert admission.status == AdmissionStatus.ADMITTED
    assert admission.unassessed_context == (PARTNER,)
    assert resolved.connections[0].link == link.link
    with pytest.raises(AdapterError, match="legacy-link-omitted"):
        from_protein(old, state=state(), connections=())
    assert admit(resolved, (capability(allow_opaque_context=False),)).status == AdmissionStatus.UNSUPPORTED


def test_protein_name_does_not_manufacture_graph_or_capability() -> None:
    r = request()
    old = CompletionRequest(SOURCE, TARGET, "ALA", tuple(ProteinObservation(o.identity, o.element, o.position)
                            for o in r.observations), ("Q",))
    with pytest.raises(TypeError):
        from_protein(old)
    assert from_protein(old, state=state()).components[0].resolved.graph.component_id == "synthetic"


def test_unknown_and_ambiguous_targets_are_not_completed() -> None:
    r = request()
    unknown = replace(r, components=(Component(TARGET, ChemistryKind.PROTEIN),))
    assert admit(unknown, (capability(),)).status == AdmissionStatus.NOT_ASSESSED
    ambiguous = replace(r, components=(Component(TARGET, ChemistryKind.PROTEIN,
                        (state(), state(graph(state_id="alternative")))),))
    assert admit(ambiguous, (capability(),)).status == AdmissionStatus.AMBIGUOUS
    selected = replace(ambiguous, components=(replace(ambiguous.components[0], selected_state="explicit"),))
    assert admit(selected, (capability(),)).status == AdmissionStatus.ADMITTED
    bad = replace(selected, components=(replace(selected.components[0], selected_state="absent"),))
    assert admit(bad, (capability(),)).status == AdmissionStatus.INVALID_INPUT


@pytest.mark.parametrize("case,expected", [
    ("observed-new", "new-atom-is-observed"), ("missing-mobile", "mobile-atom-not-observed"),
    ("silent-repack", "completion-cannot-repack"), ("missing-heavy", "unowned-missing-heavy-atoms"),
    ("wrong-element", "observed-element-conflicts-with-graph"), ("foreign", "ownership-outside-target"),
    ("wrong-name", "atom-not-in-selected-graph"), ("undeclared", "undeclared-target-component"),
])
def test_ownership_and_semantics(case: str, expected: str) -> None:
    r = request()
    if case == "observed-new":
        r = replace(r, ownership=Ownership(frozenset({atom("C")})))
    elif case == "missing-mobile":
        r = replace(r, operation=Operation.REPACK, ownership=Ownership(mobile=frozenset({atom("Q")})))
    elif case == "silent-repack":
        r = replace(r, ownership=replace(r.ownership, mobile=frozenset({atom("C")})))
    elif case == "missing-heavy":
        r = replace(r, ownership=Ownership(frozenset({atom("absent")})))
    elif case == "wrong-element":
        r = replace(r, observations=(replace(r.observations[0], element="O"), *r.observations[1:]))
    elif case == "foreign":
        r = replace(r, ownership=Ownership(frozenset({atom("Q", PARTNER)})))
    elif case == "wrong-name":
        r = replace(r, observations=(replace(r.observations[0], identity=atom("BOGUS")), *r.observations[1:]))
    else:
        r = replace(r, targets=frozenset({PARTNER}))
    assert expected in codes(r)


def test_roles_immutable_disjoint_and_explicit_repack() -> None:
    with pytest.raises(ContractError, match="roles overlap"):
        Ownership(frozenset({atom("C")}), mobile=frozenset({atom("C")}))
    with pytest.raises(ContractError, match="immutable"):
        Ownership({atom("C")})
    r = request()
    q = Observation(atom("Q"), "C", (0., 0., 1.), Unit.ANGSTROM)
    r = replace(r, observations=(*r.observations, q), operation=Operation.REPACK,
                ownership=Ownership(mobile=frozenset({atom("Q")})))
    assert admit(r, (capability(),)).status == AdmissionStatus.ADMITTED
    assert admit(r, (capability(allow_mobile=False),)).status == AdmissionStatus.UNSUPPORTED


def test_whole_missing_subset_and_nonprotein_categories_fit_same_contract() -> None:
    # Twelve absent atoms exercise planning beyond one/two-atom prototypes;
    # no positions for the absent atoms are generated or scored by admission.
    g = ChemicalGraph("synthetic-chain", "one", tuple(ChemicalAtom(f"C{i}", "C", 0, (float(i), 0., 0.))
                      for i in range(13)), tuple(ChemicalBond(f"C{i}", f"C{i+1}", 1., 1.) for i in range(12)), 0)
    r = ReconstructionRequest(SOURCE, Operation.HEAVY, (Component(TARGET, ChemistryKind.PROTEIN, (state(g),)),),
                              (Observation(atom("C0"), "C", (0., 0., 0.), Unit.ANGSTROM),), frozenset({TARGET}),
                              Ownership(frozenset(atom(f"C{i}") for i in range(1, 13))))
    for kind in (ChemistryKind.PROTEIN, ChemistryKind.NUCLEOTIDE, ChemistryKind.GLYCAN, ChemistryKind.LINKED):
        declared = replace(r, components=(replace(r.components[0], kind=kind),))
        assert admit(declared, (capability(),)).status == AdmissionStatus.ADMITTED
        assert admit(declared).status == AdmissionStatus.UNSUPPORTED


@pytest.mark.parametrize("case,expected", [
    ("duplicate", "duplicate-connection"), ("dangling", "connection-endpoint-not-materialized-or-owned"),
    ("unresolved-order", "unresolved-connection-chemistry"), ("unresolved-coordination", "unresolved-connection-chemistry"),
])
def test_links_refuse_ambiguous_duplicate_and_unmaterialized_endpoints(case: str, expected: str) -> None:
    r = request()
    partner = Observation(atom("ZN", PARTNER), "ZN", (8., 9., 10.), Unit.ANGSTROM)
    link = Connection(ExternalLink(atom("Q"), partner.identity, 1.), LinkKind.COVALENT)
    links = (link,)
    if case == "duplicate":
        links = (link, Connection(ExternalLink(partner.identity, atom("Q"), 1.), LinkKind.COVALENT))
    elif case == "dangling":
        links = (replace(link, link=replace(link.link, other_atom=atom("absent", PARTNER))),)
    elif case == "unresolved-order":
        links = (replace(link, link=replace(link.link, order=None)),)
    else:
        links = (Connection(replace(link.link, order=None), LinkKind.COORDINATION),)
    r = replace(r, components=(*r.components, Component(PARTNER, ChemistryKind.UNKNOWN)),
                observations=(*r.observations, partner), connections=links)
    assert expected in codes(r)


def test_explicit_coordination_retains_opaque_partner_without_reconstructing_it() -> None:
    r = request()
    partner = Observation(atom("ZN", PARTNER), "ZN", (8., 9., 10.), Unit.ANGSTROM)
    link = Connection(ExternalLink(atom("Q"), partner.identity), LinkKind.COORDINATION, "synthetic-model-v1")
    r = replace(r, components=(*r.components, Component(PARTNER, ChemistryKind.UNKNOWN)),
                observations=(*r.observations, partner), connections=(link,))
    a = admit(r, (capability(),))
    assert a.status == AdmissionStatus.ADMITTED and a.unassessed_context == (PARTNER,)
    assert admit(r, (capability(links=frozenset()),)).status == AdmissionStatus.UNSUPPORTED


def test_h_inventory_masks_and_heavy_anchor_ownership() -> None:
    r = request()
    s = state(graph(hydrogen=True), hydrogens=HydrogenInventory.EXPLICIT)
    r = replace(r, components=(replace(r.components[0], states=(s,)),), operation=Operation.HYDROGENS,
                ownership=Ownership(frozenset({atom("H")}), temporary=frozenset({atom("Q")})))
    assert admit(r, (capability(),)).status == AdmissionStatus.ADMITTED
    assert admit(r, (capability(allow_temporary=False),)).status == AdmissionStatus.UNSUPPORTED
    missing_heavy = replace(r, ownership=Ownership(frozenset({atom("H")})))
    assert "unowned-missing-heavy-atoms" in codes(missing_heavy)
    unspecified = replace(r, components=(replace(r.components[0], states=(replace(s, hydrogens=HydrogenInventory.UNSPECIFIED),)),))
    assert "hydrogen-inventory-not-explicit" in codes(unspecified)
    assert "generated-element-conflicts-with-operation" in codes(replace(r, operation=Operation.HEAVY))
    wrong_h = replace(r, ownership=Ownership(frozenset({atom("Q")})))
    assert "generated-element-conflicts-with-operation" in codes(wrong_h)
    temporary_h = replace(r, operation=Operation.HEAVY, ownership=Ownership(
                          frozenset({atom("Q")}), temporary=frozenset({atom("H")})),
                          components=(replace(r.components[0], states=(replace(s, hydrogens=HydrogenInventory.UNSPECIFIED),)),))
    assert "hydrogen-inventory-not-explicit" in codes(temporary_h)


def constraint(kind: ConstraintKind, names: tuple[str, ...], lower: float | None, upper: float,
               unit: Unit = Unit.ANGSTROM) -> Constraint:
    return Constraint(kind, tuple(atom(n) for n in names), unit, lower, upper)


@pytest.mark.parametrize("kind,names,lower,upper,unit", [
    (ConstraintKind.DISTANCE, ("C", "N"), .99, 1.01, Unit.ANGSTROM),
    (ConstraintKind.RING_CLOSURE, ("C", "N"), .099, .101, Unit.NANOMETER),
    (ConstraintKind.EXCLUSION, ("N", "O"), 1.4, 1.5, Unit.ANGSTROM),
    (ConstraintKind.ANGLE, ("N", "C", "O"), 89., 91., Unit.DEGREE),
    (ConstraintKind.TORSION, ("N", "C", "O", "Q"), -91., -89., Unit.DEGREE),
    (ConstraintKind.PLANARITY, ("C", "N", "O", "Q"), None, 1.01, Unit.ANGSTROM),
])
def test_independent_fixed_constraint_measurements(kind: ConstraintKind, names: tuple[str, ...],
                                                   lower: float | None, upper: float, unit: Unit) -> None:
    r = request()
    q = Observation(atom("Q"), "C", (0., 0., 1.), Unit.ANGSTROM)
    # Repack a different atom, leaving the measured frame fixed; mobile C would defer all.
    g = graph(hydrogen=True)
    h = Observation(atom("H"), "H", (-1., -1., -1.), Unit.ANGSTROM)
    r = replace(r, components=(replace(r.components[0], states=(state(g),)),),
                observations=(*r.observations, q, h), operation=Operation.REPACK,
                ownership=Ownership(mobile=frozenset({atom("H")})),
                constraints=(constraint(kind, names, lower, upper, unit),))
    a = admit(r, (capability(),))
    assert a.status == AdmissionStatus.ADMITTED, a.findings
    assert a.checked_constraints == (0,) and not a.deferred_constraints
    bad = replace(r, constraints=(replace(r.constraints[0], upper=.001, lower=None if lower is None else 0.),))
    assert "fixed-constraint-conflict" in codes(bad)


def test_fixed_stereo_conflict_and_deferred_missing_stereo() -> None:
    r = request()
    s = state(graph(stereo="R"))
    stereo = Constraint(ConstraintKind.STEREO_FRAME, tuple(atom(n) for n in ("C", "N", "O", "Q")),
                        Unit.DIMENSIONLESS, orientation=1, minimum_volume=.1, label="R")
    r = replace(r, components=(replace(r.components[0], states=(s,)),), constraints=(stereo,))
    a = admit(r, (capability(),))
    assert a.status == AdmissionStatus.ADMITTED and a.deferred_constraints == (0,)
    # Make the inverted frame fully fixed; generate H instead, not Q.
    s = state(graph(stereo="R", hydrogen=True), hydrogens=HydrogenInventory.EXPLICIT)
    r = replace(r, components=(replace(r.components[0], states=(s,)),), operation=Operation.HYDROGENS,
                observations=(*r.observations, Observation(atom("Q"), "C", (0., 0., -1.), Unit.ANGSTROM)),
                ownership=Ownership(frozenset({atom("H")})))
    before = r.observations
    assert "fixed-constraint-conflict" in codes(r)
    assert r.observations == before
    assert "ordered-stereo-mapping-required" in codes(replace(r, constraints=()))
    assert "stereo-label-conflicts-with-authority" in codes(replace(r, constraints=(replace(stereo, label="S"),)))


def test_deferred_mobile_constraint_is_not_a_fixed_geometry_pass() -> None:
    r = request()
    r = replace(r, ownership=Ownership(r.ownership.generated, mobile=frozenset({atom("C")})),
                constraints=(constraint(ConstraintKind.DISTANCE, ("C", "N"), 20., 21.),))
    # Completion does not silently permit mobility. Explicit repack is a different
    # request; complete its observed inventory before declaring it.
    assert "completion-cannot-repack" in codes(r)
    r = replace(r, operation=Operation.REPACK, ownership=Ownership(mobile=frozenset({atom("C")})),
                observations=(*r.observations, Observation(atom("Q"), "C", (0., 0., 1.), Unit.ANGSTROM)))
    a = admit(r, (capability(),))
    assert a.status == AdmissionStatus.ADMITTED and a.deferred_constraints == (0,)
    assert not a.checked_constraints and not a.geometry_approved


def test_four_ligand_stereo_and_degenerate_constraints() -> None:
    r = request()
    s = state(graph(hydrogen=True, stereo="synthetic-parity"))
    observations = tuple(Observation(atom(a.name), a.element, a.ideal_position, Unit.ANGSTROM) for a in s.graph.atoms)
    stereo = Constraint(ConstraintKind.STEREO_LIGANDS, tuple(atom(n) for n in ("C", "N", "O", "Q", "H")),
                        Unit.DIMENSIONLESS, orientation=1, minimum_volume=.01, label="synthetic-parity")
    # Move an unrelated context? Mobility outside target is forbidden; instead all
    # fixed constraint measurement can be tested independently of repack admission.
    from dvbfixer.chemical_reconstruction.measurement import satisfied

    assert satisfied(stereo, tuple(o.position for o in observations))
    assert not satisfied(replace(stereo, orientation=-1), tuple(o.position for o in observations))
    degenerate = constraint(ConstraintKind.ANGLE, ("C", "N", "O"), 0., 180., Unit.DEGREE)
    r = replace(r, constraints=(degenerate,), observations=(replace(r.observations[0], position=(1., 0., 0.)), *r.observations[1:]))
    assert "degenerate-fixed-constraint" in codes(r)


@pytest.mark.parametrize("kind", [ConstraintKind.DISTANCE, ConstraintKind.ANGLE, ConstraintKind.PLANARITY])
def test_invalid_constraint_units_and_references(kind: ConstraintKind) -> None:
    names = {ConstraintKind.DISTANCE: ("C", "N"), ConstraintKind.ANGLE: ("C", "N", "O"),
             ConstraintKind.PLANARITY: ("C", "N", "O", "Q")}[kind]
    with pytest.raises(ContractError):
        constraint(kind, names, None if kind == ConstraintKind.PLANARITY else 0., 1., Unit.DIMENSIONLESS)
    r = replace(request(), constraints=(constraint(ConstraintKind.DISTANCE, ("C", "ABSENT"), .9, 1.1),))
    assert "constraint-reference-not-materialized-or-owned" in codes(r)


def legacy_class_a() -> LegacyClassA:
    r = request()
    g = graph()
    return LegacyClassA(TARGET, tuple(ObservedAtom(o.identity, o.element, o.position) for o in r.observations),
                        UserMappedAuthority(g, "synthetic", canonical_graph_digest(g)), 7., ComponentClass.A)


def test_class_a_transport_preserves_context_policy_and_explicit_state() -> None:
    old = legacy_class_a()
    extra = Observation(atom("ZN", PARTNER), "ZN", (8., 9., 10.), Unit.ANGSTROM)
    originals = tuple(Observation(o.identity, o.element, o.position, Unit.ANGSTROM) for o in old.observed_atoms)
    r = from_class_a(old, source_sha256=SOURCE, context_observations=(*originals, extra))
    assert r.observations == (*originals, extra) and r.ownership.generated == frozenset({atom("Q")})
    assert r.ph == old.ph and r.components[0].resolved.graph == old.authority.graph
    assert r.origin is not None
    assert dict(r.origin.settings) == {"requested_class": "A", "max_anchor_rmsd": "0.35",
                                      "max_bond_deviation": "0.25", "min_nonbonded_distance": "0.7",
                                      "max_observed_movement": "0.0"}
    assert admit(r, (capability(),)).status == AdmissionStatus.ADMITTED
    changed = from_class_a(replace(old, geometry_policy=GeometryPolicy(max_anchor_rmsd=.4)), source_sha256=SOURCE)
    assert changed.digest != from_class_a(old, source_sha256=SOURCE).digest
    bad_state = from_class_a(replace(old, microstate_override="absent"), source_sha256=SOURCE)
    assert "unknown-selected-state" in codes(bad_state)
    with pytest.raises(AdapterError, match="omitted-or-changed"):
        from_class_a(old, source_sha256=SOURCE, context_observations=(extra,))


@pytest.mark.parametrize("case", ["local", "online", "link", "wrong-class", "movement", "digest", "coordination", "non-chno"])
def test_class_a_refuses_unresolved_and_out_of_scope_transport(case: str) -> None:
    old = legacy_class_a()
    if case == "local":
        old = replace(old, authority=LocalPinnedCcdAuthority(Path("/must-not-open"), SOURCE, "synthetic"))
    elif case == "online":
        old = replace(old, authority=OnlineEnrichmentAuthority("must-not-fetch"))
    elif case == "link":
        old = replace(old, external_links=(ExternalLink(atom("Q"), atom("ZN", PARTNER), 1.),))
    elif case == "wrong-class":
        old = replace(old, requested_class=ComponentClass.B)
    elif case == "movement":
        old = replace(old, geometry_policy=GeometryPolicy(max_observed_movement=.01))
    elif case == "digest":
        old = replace(old, authority=replace(old.authority, source_digest="0"*64))
    else:
        g = replace(graph(), coordination_model="declared") if case == "coordination" else replace(
            graph(), atoms=(replace(graph().atoms[0], element="S"), *graph().atoms[1:]))
        old = replace(old, authority=UserMappedAuthority(g, "synthetic", canonical_graph_digest(g)))
    with pytest.raises(ContractError):
        from_class_a(old, source_sha256=SOURCE)


@pytest.mark.parametrize("change", ["operation", "chemistry", "elements", "constraints", "anchors"])
def test_capabilities_are_explicit_not_implicit_engine_registration(change: str) -> None:
    r, cap = request(), capability()
    if change == "operation":
        cap = replace(cap, operations=frozenset({Operation.HYDROGENS}))
    elif change == "chemistry":
        cap = replace(cap, chemistry=frozenset({ChemistryKind.GLYCAN}))
    elif change == "elements":
        cap = replace(cap, elements=frozenset({"C"}))
    elif change == "constraints":
        r = replace(r, constraints=(constraint(ConstraintKind.DISTANCE, ("C", "N"), .9, 1.1),))
        cap = replace(cap, constraints=frozenset())
    else:
        cap = replace(cap, minimum_anchors=100)
    assert admit(r, (cap,)).status == AdmissionStatus.UNSUPPORTED


def test_pure_boundary_import_and_execution_do_not_load_science_or_perform_io() -> None:
    # Fresh interpreter prevents already-imported toolkits from masking a leak.
    program = '''
import sys
forbidden = {"numpy", "scipy", "openmm", "pdbfixer", "rdkit", "Bio", "gemmi", "pyrosetta"}
importing = True
def audit(event, args):
    if event.startswith("socket.") or event in {"subprocess.Popen", "os.system"}:
        raise AssertionError(event)
    if not importing and event in {"os.listdir", "os.scandir"}:
        raise AssertionError(event)
    if event == "open":
        path = str(args[0])
        if not importing or not path.endswith((".py", ".pyc")):
            raise AssertionError((event, path))
    if event == "import" and args[0].split(".")[0] in forbidden:
        raise AssertionError(args[0])
sys.addaudithook(audit)
from dvbfixer.chemical_reconstruction.admission import admit
from dvbfixer.chemical_reconstruction.adapters import from_class_a, from_protein
from dvbfixer.chemical_reconstruction.contracts import *
from dvbfixer.domain.nonprotein_reconstruction import *
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation as OldObservation
from pathlib import Path
from dataclasses import replace
importing = False
target = ComponentInstanceRef(1, "a", "1")
g = ChemicalGraph("synthetic", "one", (ChemicalAtom("C", "C", 0, (0.,0.,0.)),
    ChemicalAtom("Q", "C", 0, (1.,0.,0.))), (ChemicalBond("C","Q",1.,1.),), 0)
old = ReconstructionRequest(target, (ObservedAtom(ExactAtomRef(target,"C"),"C",(0.,0.,0.)),),
    UserMappedAuthority(g,"synthetic",canonical_graph_digest(g)), 7.)
r = from_class_a(old, source_sha256="a"*64)
assert admit(r).status == AdmissionStatus.UNSUPPORTED
bad = replace(r, ownership=Ownership(frozenset({ExactAtomRef(target,"Z"), ExactAtomRef(target,"Y")})))
invalid = admit(bad)
assert invalid.status == AdmissionStatus.INVALID_INPUT
print(invalid)
protein = CompletionRequest("a"*64, target, "ALA", (OldObservation(ExactAtomRef(target,"C"),"C",(0.,0.,0.)),), ("Q",))
assert admit(from_protein(protein, state=r.components[0].resolved)).status == AdmissionStatus.UNSUPPORTED
for authority in (LocalPinnedCcdAuthority(Path("/must-not-read"), "a"*64, "synthetic"),
                  OnlineEnrichmentAuthority("must-not-connect")):
    unresolved = ReconstructionRequest(target, old.observed_atoms, authority, 7.)
    try:
        from_class_a(unresolved, source_sha256="a"*64)
    except ValueError:
        pass
    else:
        raise AssertionError("unresolved authority admitted")
assert not forbidden.intersection(sys.modules)
'''
    outputs = []
    for seed in ("1", "2"):
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED=seed)
        result = subprocess.run([sys.executable, "-B", "-c", program], capture_output=True, text=True, env=env, timeout=30)
        assert result.returncode == 0, result.stderr
        outputs.append(result.stdout)
    assert outputs[0] == outputs[1], "admission order/digests must not depend on set hashing"
