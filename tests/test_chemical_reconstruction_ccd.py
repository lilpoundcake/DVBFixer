"""Pinned real CCD geometry transport, not experimental/polymer qualification."""
from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import gemmi
import pytest
from scripts.chemical_reconstruction_ccd_controls import (
    COMPONENTS,
    FIXTURES,
    HEAVY_MASKS,
    POLICY,
    run_controls,
)

from dvbfixer.ccd import CcdError, ccd_to_all_atom_graph, ccd_to_heavy_atom_graph
from dvbfixer.chemical_reconstruction.ccd_authority import (
    CcdAuthorityError,
    load_geometry_authority,
    verify_cip,
)
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    ConstraintKind,
    ContractError,
    HydrogenInventory,
    Observation,
    Operation,
    Ownership,
    ReconstructionRequest,
    Unit,
)
from dvbfixer.chemical_reconstruction.plan_compiler import PlanCompilationError, compile_chain_plan
from dvbfixer.domain.nonprotein_reconstruction import ChemicalAtom, ChemicalBond, ChemicalGraph
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

PINS = json.loads((FIXTURES / "manifest.json").read_text())["components"]
IDENTITY = ComponentInstanceRef(1, "a", "82", "B", occurrence=2)


def authority(name: str):
    return load_geometry_authority(FIXTURES / f"{name}.cif", name, expected_sha256=PINS[name]["sha256"])


def request(name: str, *, hydrogens: bool = False) -> ReconstructionRequest:
    a = authority(name)
    missing = {x.name for x in a.state.graph.atoms if x.element == "H"} if hydrogens else set(HEAVY_MASKS[name])
    return ReconstructionRequest(
        PINS[name]["sha256"], Operation.HYDROGENS if hydrogens else Operation.HEAVY,
        (Component(IDENTITY, ChemistryKind.ISOLATED, (a.state,)),),
        tuple(Observation(ExactAtomRef(IDENTITY, x.name), x.element, x.ideal_position, Unit.ANGSTROM)
              for x in a.state.graph.atoms if x.element != "H" and x.name not in missing),
        frozenset({IDENTITY}), Ownership(frozenset(ExactAtomRef(IDENTITY, n) for n in missing)),
        constraints=a.constraints(IDENTITY, include_hydrogens=hydrogens, minimum_volume=.05, planar_tolerance_angstrom=.02),
    )


@pytest.mark.parametrize("name", COMPONENTS)
def test_pinned_all_atom_authority_retains_source_state_and_old_heavy_projection(name: str) -> None:
    a = authority(name)
    assert a.component.provenance.sha256 == PINS[name]["sha256"]
    assert a.component.modified_date == PINS[name]["modified_date"]
    assert len(a.state.graph.atoms) == PINS[name]["atom_count"]
    assert len(a.state.graph.bonds) == PINS[name]["bond_count"]
    assert a.state.hydrogens == HydrogenInventory.EXPLICIT
    assert {x.name for x in a.state.graph.atoms} == {x.name for x in a.component.atoms}
    old = ccd_to_heavy_atom_graph(a.component)
    assert tuple(x for x in a.state.graph.atoms if x.element != "H") == old.atoms
    assert tuple(b for b in a.state.graph.bonds if b.atom1 in {x.name for x in old.atoms}
                 and b.atom2 in {x.name for x in old.atoms}) == old.bonds
    assert a.state.graph.net_charge == a.component.formal_charge == old.net_charge
    assert old.state_id != a.state.graph.state_id


@pytest.mark.parametrize("name", ["ALA", "THR", "ILE", "LYS", "PHE"])
def test_free_amino_acid_leaving_groups_are_not_silently_polymerized(name: str) -> None:
    a = authority(name)
    assert set(a.leaving_atoms) == {"OXT", "H2", "HXT"}
    assert set(a.leaving_atoms) <= {x.name for x in a.state.graph.atoms}


def test_independent_cip_matches_threonine_and_isoleucine_source_labels() -> None:
    assert dict(authority("THR").verified_cip) == {"CA": "S", "CB": "R"}
    assert dict(authority("ILE").verified_cip) == {"CA": "S", "CB": "S"}
    a = authority("THR")
    mirrored = {x.name: (-x.ideal_position[0], x.ideal_position[1], x.ideal_position[2]) for x in a.state.graph.atoms}
    with pytest.raises(CcdAuthorityError, match="CIP mismatch"):
        verify_cip(a.state.graph, mirrored)


def test_missing_explicit_h_is_not_accepted_as_closed_shell_authority() -> None:
    a = authority("ALA")
    removed = next(x.name for x in a.state.graph.atoms if x.element == "H")
    g = replace(a.state.graph, atoms=tuple(x for x in a.state.graph.atoms if x.name != removed),
                bonds=tuple(b for b in a.state.graph.bonds if removed not in {b.atom1, b.atom2}))
    with pytest.raises(CcdAuthorityError, match="radical or incomplete"):
        verify_cip(g, {x.name: x.ideal_position for x in g.atoms})


def test_rdkit_cleanup_cannot_silently_repair_wrong_source_charges() -> None:
    # RDKit normalizes neutral five-valent nitro N to [N+]/[O-]. That is useful
    # for permissive chemistry import, but forbidden for this authoritative state.
    atoms = (ChemicalAtom("C", "C", 0, (0., 0., 0.)), ChemicalAtom("N", "N", 0, (1.5, 0., 0.)),
             ChemicalAtom("O1", "O", 0, (2., 1., 0.)), ChemicalAtom("O2", "O", 0, (2., -1., 0.)),
             ChemicalAtom("H1", "H", 0, (-1., 0., 0.)), ChemicalAtom("H2", "H", 0, (0., 1., 0.)),
             ChemicalAtom("H3", "H", 0, (0., 0., 1.)))
    g = ChemicalGraph("invalid-nitro-source", "one", atoms,
                      (ChemicalBond("C", "N", 1., 1.5), ChemicalBond("N", "O1", 2., 1.),
                       ChemicalBond("N", "O2", 2., 1.), *(ChemicalBond("C", n, 1., 1.) for n in ("H1", "H2", "H3"))), 0)
    with pytest.raises(CcdAuthorityError, match="cleanup changed"):
        verify_cip(g, {a.name: a.ideal_position for a in atoms})


def test_pin_validation_and_no_network_fallback(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import dvbfixer.ccd as ccd

    def forbidden(*args, **kwargs):
        raise AssertionError("authority attempted network lookup")

    monkeypatch.setattr(ccd, "fetch_ccd_component", forbidden)
    with pytest.raises(CcdError, match="required pin"):
        load_geometry_authority(FIXTURES / "ALA.cif", "ALA", expected_sha256="0"*64)
    with pytest.raises(CcdError, match="regular non-symlink"):
        load_geometry_authority(tmp_path / "absent.cif", "ALA", expected_sha256=PINS["ALA"]["sha256"])
    assert authority("ALA").state.graph.atoms


@pytest.mark.parametrize("pin", [None, "", False, 42])
def test_runtime_invalid_pin_cannot_enable_the_older_unpinned_information_route(pin, monkeypatch: pytest.MonkeyPatch) -> None:
    import dvbfixer.chemical_reconstruction.ccd_authority as adapter

    def forbidden(*args, **kwargs):
        raise AssertionError("invalid pin reached the file loader")
    monkeypatch.setattr(adapter, "load_ccd_component", forbidden)
    with pytest.raises(ContractError, match="SHA-256 required"):
        load_geometry_authority(FIXTURES / "ALA.cif", "ALA", expected_sha256=pin)


@pytest.mark.parametrize("mutation", ["stereo", "bond-stereo", "experimental-coordinates", "missing-h-coordinate"])
def test_source_mutations_fail_closed_or_leave_ideal_geometry_unchanged(mutation: str, tmp_path: Path) -> None:
    document = gemmi.cif.read_string((FIXTURES / "ALA.cif").read_text())
    block = document.sole_block()
    if mutation == "stereo":
        for row in block.find(["_chem_comp_atom.atom_id", "_chem_comp_atom.pdbx_stereo_config"]):
            if row[0] == "CA":
                row[1] = "R"
    elif mutation == "bond-stereo":
        block.find(["_chem_comp_bond.pdbx_stereo_config"])[0][0] = "E"
    elif mutation == "experimental-coordinates":
        for row in block.find(["_chem_comp_atom.model_Cartn_x"]):
            row[0] = "9999"
    else:
        for row in block.find(["_chem_comp_atom.type_symbol", "_chem_comp_atom.pdbx_model_Cartn_x_ideal"]):
            if row[0] == "H":
                row[1] = "?"
                break
    content = document.as_string().encode()
    path = tmp_path / "altered.cif"
    path.write_bytes(content)
    if mutation == "experimental-coordinates":
        changed = load_geometry_authority(path, "ALA", expected_sha256=hashlib.sha256(content).hexdigest())
        assert changed.state.graph.atoms == authority("ALA").state.graph.atoms
        assert changed.verified_cip == authority("ALA").verified_cip
    else:
        with pytest.raises((CcdAuthorityError, CcdError)):
            load_geometry_authority(path, "ALA", expected_sha256=hashlib.sha256(content).hexdigest())


def test_full_sidechain_compiler_and_hydrogen_compiler_preserve_exact_ownership() -> None:
    r = request("LYS")
    plan = compile_chain_plan(r, POLICY)
    assert tuple(p.atom.name for p in plan.placements) == ("CB", "CG", "CD", "CE", "NZ")
    assert propose(r, plan).status == "declared-geometry-passed"
    h = request("PHE", hydrogens=True)
    hp = compile_chain_plan(h, POLICY)
    assert {p.atom for p in hp.placements} == h.ownership.generated
    assert any(c.kind == ConstraintKind.PLANARITY for c in h.constraints)
    assert propose(h, hp).status == "declared-geometry-passed"


def test_compiler_does_not_change_frames_based_on_observed_geometry() -> None:
    r = request("LYS")
    plan = compile_chain_plan(r, POLICY)
    # No completely-fixed plane check is present in this controlled invalid-source
    # request. Compilation is topology/authority planning, not a source repair.
    broken = replace(r, observations=tuple(replace(o, position=(0., 0., 0.)) for o in r.observations),
                     constraints=tuple(c for c in r.constraints if c.kind != ConstraintKind.PLANARITY))
    bp = compile_chain_plan(broken, POLICY)
    assert bp.placements == plan.placements
    assert propose(broken, bp).status == "construction-failed"
    assert all(o.position == (0., 0., 0.) for o in broken.observations)


def test_roundtrip_regression_replay_is_not_independent_scientific_evidence(tmp_path: Path) -> None:
    output = tmp_path / "replay"
    summary = run_controls(output)
    assert summary["declared"] == summary["executed_entries"] == 28
    assert summary["geometry_passed"] == 24 and summary["compiler_refused"] == 4
    assert all(c["observations_unchanged"] for c in summary["controls"])
    assert all(c["raw_cip"]["passed"] and c["rounded_cip"]["passed"] for c in summary["controls"]
               if "-hydrogens-" in c["id"])
    before = (output / "summary.json").read_bytes()
    with pytest.raises(FileExistsError):
        run_controls(output)
    assert (output / "summary.json").read_bytes() == before
    assert not summary["held_out_or_native_accuracy_qualification"] and not summary["polymer_state_or_h_qualification"]
    # Independent JSON readback verifies each request/plan binding, not just flags.
    for control in summary["controls"]:
        evidence = json.loads((output / f"{control['id']}.json").read_text())
        payload = json.dumps(["chemical-reconstruction-contract-v1", evidence["request"]],
                             sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
        assert hashlib.sha256(payload.encode()).hexdigest() == control["request_digest"]
        if "plan" in evidence:
            payload = json.dumps(["chemical-geometry-plan-v1", evidence["plan"]],
                                 sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
            assert hashlib.sha256(payload.encode()).hexdigest() == evidence["result"]["plan_digest"]


def test_all_atom_converter_does_not_edit_source_objects() -> None:
    a = authority("ACT")
    before = a.component
    assert ccd_to_all_atom_graph(before) == a.state.graph
    assert a.component == before and a.state.graph.net_charge == -1


def test_empty_or_unknown_target_cannot_compile() -> None:
    r = request("ALA")
    r = replace(r, components=(Component(IDENTITY, ChemistryKind.UNKNOWN),))
    with pytest.raises(PlanCompilationError, match="request-not-admissible"):
        compile_chain_plan(r, POLICY)
