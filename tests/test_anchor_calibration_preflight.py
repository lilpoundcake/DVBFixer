from types import SimpleNamespace

import pytest
from scripts import anchor_calibration_preflight as study

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def atom_line(serial, name, resid, x, *, chain="a", record="ATOM", resname="ALA"):
    return (f"{record:6s}{serial:5d} {name:^4s} {resname:3s} {chain}{resid:4d}    "
            f"{x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{10.:6.2f}          {name[0]:>2s}  \n")


def test_grid_is_complete_and_zero_control_is_not_a_strength_cell():
    cells = study.variants()
    assert len(cells) == len({c["id"] for c in cells}) == 25
    assert cells[0] == {"id": "fixed", "window": 0, "envelope": None, "strength": None}
    assert {c["window"] for c in cells[1:]} == {1, 2}
    assert all(sum(c["window"] == w for c in cells) == 12 for w in (1, 2))
    assert {tuple(c["envelope"]) for c in cells[1:]} == set(study.ENVELOPES)
    assert {c["strength"] for c in cells[1:]} == set(study.STRENGTHS)


def request_for_window(monkeypatch, observed_indices, gaps):
    axis = tuple(ResidueIdentity("a", str(i)) for i in range(1, 13))
    request = SimpleNamespace(
        gaps=tuple(SimpleNamespace(chain="a", left_anchor=axis[left] if left is not None else None,
                                   right_anchor=axis[right] if right is not None else None)
                   for left, right in gaps),
        sequence_placements=(SimpleNamespace(observed_residues=tuple(axis[i] for i in observed_indices)),))
    monkeypatch.setattr(study, "target_axis", lambda _: {"a": axis})
    return request, axis


def test_two_means_two_not_old_plus_minus_two(monkeypatch):
    req, axis = request_for_window(monkeypatch, (0, 1, 2, 6, 7, 8, 9, 10, 11), [(2, 6)])
    assert study.anchor_residues(req, 0) == ()
    assert set(study.anchor_residues(req, 1)) == {axis[2], axis[6]}
    assert set(study.anchor_residues(req, 2)) == {axis[i] for i in (1, 2, 6, 7)}


def test_windows_do_not_cross_gaps_or_release_terminal_anchors(monkeypatch):
    req, axis = request_for_window(monkeypatch, (2, 6, 7, 10, 11), [(None, 2), (2, 6), (7, 10), (11, None)])
    assert set(study.anchor_residues(req, 2)) == {axis[i] for i in (2, 6, 7, 10, 11)}
    assert len(study.anchor_residues(req, 2)) == 5


@pytest.mark.parametrize("width", [-1, 3, True, 1.0, "2"])
def test_undeclared_width_refused(width):
    with pytest.raises(ValueError, match="undeclared"):
        study.anchor_residues(None, width)


def test_case_sensitive_chains_and_insertion_codes(monkeypatch):
    left = ResidueIdentity("a", "82", "A")
    right = ResidueIdentity("a", "84")
    other = ResidueIdentity("A", "82", "A")
    generated = ResidueIdentity("a", "83")
    req = SimpleNamespace(gaps=(SimpleNamespace(chain="a", left_anchor=left, right_anchor=right),),
        sequence_placements=(SimpleNamespace(observed_residues=(left, right, other)),))
    monkeypatch.setattr(study, "target_axis", lambda _: {"a": (left, generated, right), "A": (other,)})
    assert set(study.anchor_residues(req, 2)) == {left, right}


def test_full_context_inventory_and_observations_not_repaired():
    source = atom_line(1, "CA", 1, 0.)
    added = atom_line(3, "CA", 2, 3.)
    water = atom_line(2, "O", 8, 7., record="HETATM", resname="HOH")
    req = SimpleNamespace(fixed_atoms=(AtomIdentity("a", "1", "", "CA"),),
                          generated_atoms=(AtomIdentity("a", "2", "", "CA"),))
    result = study.validate_proposal(req, source, source+water, source+added+water)
    assert result == {"observed_protein_atoms": 1, "context_atoms": 1,
                      "generated_atoms": 1, "base_observations_exact": True}
    with pytest.raises(ValueError, match="source records"):
        study.validate_proposal(req, source, source+water, atom_line(1, "CA", 1, .001)+added+water)
    with pytest.raises(ValueError, match="inventory"):
        study.validate_proposal(req, source, source+water, source+added)
    with pytest.raises(ValueError, match="context observation"):
        study.validate_proposal(req, source, source+water, source+added+water.replace("   7.000", "   8.000"))


def test_only_serial_alias_is_allowed():
    first = atom_line(1, "CA", 1, 0.)
    assert study.same_observation(first, atom_line(99, "CA", 1, 0.))
    assert not study.same_observation(first, first.replace("  1.00", "  0.50"))
    assert not study.same_observation(first, atom_line(1, "CA", 1, .001))


def test_replay_historical_preparation_is_explicit_and_preserves_original_context():
    original_atom = atom_line(1, "CA", 1, 0.)
    prepared_atom = atom_line(1, "CA", 1, .1)
    added = atom_line(3, "CB", 1, 2.)
    water = atom_line(2, "O", 8, 7., record="HETATM", resname="HOH")
    original = original_atom + water + "END\n"
    rendered, evidence = study.prepared_context(prepared_atom+added+"END\n", original)
    assert prepared_atom in rendered and water in rendered
    assert len(evidence["historical_prepared_coordinate_changes"]) == 1
    assert evidence["historical_prepared_coordinate_changes"][0]["displacement_angstrom"] == .1
    assert evidence["historical_added_fixed_atoms"] == [{"chain": "a", "residue_number": "1", "insertion_code": "", "atom_name": "CB"}]
    assert evidence["new_source_repair_executed"] is False
    assert original == original_atom+water+"END\n"
    with pytest.raises(ValueError, match="metadata"):
        study.prepared_context(prepared_atom.replace(" 10.00", " 20.00"), original)
    with pytest.raises(ValueError, match="removed"):
        study.prepared_context(atom_line(3, "CA", 2, 4.), original)


def test_independent_windows_match_declared_outward_not_symmetric_rule():
    from scripts.audit_anchor_calibration_preflight import independent_windows

    axis = [ResidueIdentity("a", str(i)) for i in range(1, 11)]
    indices = (0, 1, 2, 6, 7, 8, 9)
    req = SimpleNamespace(sequence_placements=(SimpleNamespace(chain="a", observed_target_indices=indices,
        observed_residues=tuple(axis[i] for i in indices)),),
        gaps=(SimpleNamespace(chain="a", left_anchor=axis[2], right_anchor=axis[6],
                              target_interval=SimpleNamespace(start=3, stop=6)),))
    assert set(independent_windows(req, 2)) == {axis[i] for i in (1, 2, 6, 7)}
    assert independent_windows(req, 0) == []


def test_independent_snapshot_readback_rejects_wrong_mass_and_nonfinite_force(tmp_path, monkeypatch):
    import json
    from dataclasses import asdict

    import numpy as np
    from openmm import NonbondedForce, System, XmlSerializer
    from scripts import audit_anchor_calibration_preflight as audit

    case, snapshot = tmp_path / "case", tmp_path / "snapshot"
    case.mkdir()
    snapshot.mkdir()
    (case / "request.json").write_text("{}")
    (case / "proposal.pdb").write_text(atom_line(1, "CA", 1, 0.)+atom_line(2, "CA", 2, 3.)
        +atom_line(3, "CA", 3, 6.)+atom_line(4, "O", 9, 12., record="HETATM", resname="HOH"))
    residues = [ResidueIdentity("a", str(i)) for i in (1, 3)]
    req = SimpleNamespace(generated_atoms=(AtomIdentity("a", "2", "", "CA"),),
        sequence_placements=(SimpleNamespace(chain="a", observed_target_indices=(0, 2), observed_residues=residues),),
        gaps=(SimpleNamespace(chain="a", left_anchor=residues[0], right_anchor=residues[1],
                              target_interval=SimpleNamespace(start=1, stop=2)),))
    monkeypatch.setattr(audit, "DiffusionRequest", SimpleNamespace(from_json=lambda _: req))
    identities = [AtomIdentity("a", str(i), "", "CA") for i in (1, 2, 3)]
    identities += [AtomIdentity("a", "9", "", "O"), AtomIdentity("a", "1", "", "H")]
    system = System()
    nonbonded = NonbondedForce()
    rows = []
    for i, (atom, mass) in enumerate(zip(identities, (12., 12., 12., 0., 1.), strict=True)):
        system.addParticle(mass)
        nonbonded.addParticle(0., .3, .1)
        rows.append({"identity": asdict(atom), "mass_dalton": mass, "temporary_addition": i == 4,
                     "element": "H" if i == 4 else "C"})
    system.addForce(nonbonded)
    (snapshot / "temporary-axis.json").write_text(json.dumps({"atoms": rows}))
    (snapshot / "temporary-system.xml").write_text(XmlSerializer.serialize(system))
    positions = np.array([[0, 0, 0], [.3, 0, 0], [.6, 0, 0], [1.2, 0, 0], [.1, 0, 0]])
    np.savez(snapshot / "initial-state.npz", positions_nm=positions, full_force_kj_mol_nm=np.zeros((5, 3)))
    report = audit.audit_snapshot(case, snapshot, 2)
    assert report["generated_mobile_heavy_atoms"] == 1
    assert report["observed_mobile_heavy_atoms"] == 2
    assert report["fixed_heavy_atoms"] == 1
    assert report["temporary_heavy_additions"] == []
    system.setParticleMass(1, 0.)
    (snapshot / "temporary-system.xml").write_text(XmlSerializer.serialize(system))
    with pytest.raises(ValueError, match="masks"):
        audit.audit_snapshot(case, snapshot, 2)
    (snapshot / "temporary-system.xml").write_text(XmlSerializer.serialize(System()))
    np.savez(snapshot / "initial-state.npz", positions_nm=positions, full_force_kj_mol_nm=np.full((5, 3), np.nan))
    with pytest.raises(ValueError, match="nonfinite"):
        audit.audit_snapshot(case, snapshot, 2)


def test_declaration_binds_immutable_protocol_not_living_handoff(tmp_path, monkeypatch):
    import json

    monkeypatch.setattr(study, "ENTRIES", ())
    output = tmp_path / "declaration"
    study.declare(output)
    plan = json.loads((output / "plan.json").read_text())
    bindings = plan["source_code_and_input_sha256"]
    assert str(output / "protocol.md") in bindings
    assert str(study.PROJECT / "docs/plans/experimental-structure-completion-handoff.md") not in bindings
    assert (output / "protocol.md").read_bytes() == (study.PROJECT / "docs/plans/experimental-structure-completion-handoff.md").read_bytes()
    assert plan["study_optimizer_calls"] == 0 and plan["grid_execution_enabled"] is False


def ff_terminus_fixture():
    from openmm.app import Topology, element

    topology = Topology()
    chain = topology.addChain("a")
    residue = topology.addResidue("ALA", chain, "82", insertionCode="A")
    axis = {name: topology.addAtom(name, elem, residue) for name, elem in
            (("N", element.nitrogen), ("CA", element.carbon), ("C", element.carbon), ("O", element.oxygen))}
    water = topology.addResidue("HOH", chain, "99")
    topology.addAtom("O", element.oxygen, water)
    return SimpleNamespace(topology=topology, missingTerminals={}), residue, axis


def test_explicit_ff_terminus_is_found_before_water_and_not_added_to_source():
    fixer, residue, axis = ff_terminus_fixture()
    before = tuple(fixer.topology.atoms())
    events = study.declare_ff_c_termini(fixer, (ResidueIdentity("a", "82", "A"),))
    assert fixer.missingTerminals == {residue: ["OXT"]}
    assert tuple(fixer.topology.atoms()) == before  # Only plan; seeded helper owns coordinates.
    assert len(events) == 1 and events[0]["output_ownership_extended"] is False
    assert events[0]["already_planned_by_pdbfixer"] is False
    assert study.declare_ff_c_termini(fixer, (ResidueIdentity("a", "82", "A"),))[0]["already_planned_by_pdbfixer"]
    assert fixer.missingTerminals == {residue: ["OXT"]}


def test_ff_terminus_covalent_partner_is_not_silently_dropped():
    from openmm.app import element

    fixer, residue, axis = ff_terminus_fixture()
    cap = fixer.topology.addResidue("NME", residue.chain, "83")
    cap_n = fixer.topology.addAtom("N", element.nitrogen, cap)
    fixer.topology.addBond(axis["C"], cap_n)
    with pytest.raises(ValueError, match="external covalent"):
        study.declare_ff_c_termini(fixer, (ResidueIdentity("a", "82", "A"),))
    assert len(list(fixer.topology.bonds())) == 1 and fixer.missingTerminals == {}


def test_ff_terminus_case_and_insertion_code_are_exact():
    fixer, _, _ = ff_terminus_fixture()
    for identity in (ResidueIdentity("A", "82", "A"), ResidueIdentity("a", "82")):
        with pytest.raises(ValueError, match="absent"):
            study.declare_ff_c_termini(fixer, (identity,))
