"""Check the measurements used by the bounded dependency experiments."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("openmm")


def load_script(name):
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_frozen_drift_distinguishes_translation_rotation_and_deformation():
    audit = load_script("audit_xtb_controls")
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    moved = points + [2, 0, 0]
    metrics = audit.drift_metrics(points, moved)
    assert metrics["max_angstrom"] == pytest.approx(2)
    assert metrics["rigid_fit_rms_angstrom"] < 1e-12
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert audit.drift_metrics(points, moved @ rotation)["rigid_fit_rms_angstrom"] < 1e-12
    moved[0, 0] += 0.1
    assert audit.drift_metrics(points, moved)["rigid_fit_rms_angstrom"] > 0.01


@pytest.mark.parametrize("contents", ["2\n\nC 0 0 0\n", "1\n\nC nan 0 0\n"])
def test_research_xyz_reader_rejects_invalid_results(tmp_path, contents):
    path = tmp_path / "invalid.xyz"
    path.write_text(contents)
    with pytest.raises(ValueError):
        load_script("audit_xtb_controls").read_xyz(path)


def test_probe_pair_identity_preserves_case_and_insertion_code():
    identity = load_script("audit_probe_contacts").probe_identity
    assert identity(" H  82ACYX  HB3 ") == ("H", "82", "A", "HB3")
    assert identity(" h  82ACYX  HB3 ") == ("h", "82", "A", "HB3")


def test_partial_trajectory_is_not_mistaken_for_a_complete_last_frame(tmp_path):
    path = tmp_path / "xtbopt.log"
    first = "2\nfirst complete frame\nC 0 0 0\nO 1.4 0 0\n"
    path.write_text(first + "2\ntruncated frame\nC 1 2 3\n")
    assert load_script("audit_xtb_controls").last_complete_xyz_frame(path) == first


def test_pair_explanations_resolve_blank_insertion_codes():
    from types import SimpleNamespace

    from openmm import Vec3
    from openmm.app import Element, Topology
    from openmm.unit import Quantity, nanometer

    top = Topology()
    residue = top.addResidue("LIG", top.addChain("a"), "1")
    top.addAtom("C1", Element.getBySymbol("C"), residue)
    top.addAtom("C2", Element.getBySymbol("C"), residue)
    pdb = SimpleNamespace(topology=top, positions=Quantity([Vec3(0, 0, 0), Vec3(.2, 0, 0)], nanometer))
    pair = (("a", "1", "", "C1"), ("a", "1", "", "C2"))
    rows = load_script("audit_probe_contacts").explain_discordant_pairs(pdb, set(), {pair}, [])
    assert rows[0]["identity_status"] == "resolved"
    assert "same_residue" in rows[0]["python_exclusion_filters"]


def test_pocket_contact_loss_is_measured_without_claiming_affinity():
    from openmm.app import Element, Topology

    top = Topology()
    chain = top.addChain("a")
    protein = top.addResidue("ALA", chain, "1")
    ligand = top.addResidue("LIG", chain, "2")
    top.addAtom("CA", Element.getBySymbol("C"), protein)
    top.addAtom("C1", Element.getBySymbol("C"), ligand)
    before = np.array([[0., 0., 0.], [3.8, 0., 0.]])
    after = np.array([[0., 0., 0.], [6., 0., 0.]])
    rows = load_script("audit_xtb_controls").pocket_snapshot(top, before, after)
    assert len(rows) == 1
    assert rows[0]["lost_contacts"] == [(1, 0)]
    assert rows[0]["heavy_displacement_max_angstrom"] == pytest.approx(2.2)


def test_qm_controls_do_not_silently_replace_gfnff():
    arguments = load_script("audit_xtb_controls").method_arguments
    assert arguments("gfnff", None, None, "none") == ["--gfnff"]
    assert arguments("gfn2", -1, 0, "water") == [
        "--gfn", "2", "--chrg", "-1", "--uhf", "0", "--alpb", "water"
    ]
    with pytest.raises(ValueError, match="explicit"):
        arguments("gfn2", -1, None, "none")
    with pytest.raises(ValueError, match="explicit"):
        arguments("gfn2", None, 0, "none")
    with pytest.raises(ValueError, match="nonnegative"):
        arguments("gfn2", -1, -1, "none")


def test_directional_control_preserves_distance_and_changes_donor_angle(tmp_path):
    from openmm.app import PDBFile
    from openmm.unit import nanometer

    controls = load_script("audit_probe_contacts").build_directional_controls(tmp_path)
    for control in controls.values():
        pdb = PDBFile(str(control["path"]))
        lookup = {(a.residue.chain.id, a.name): a for a in pdb.topology.atoms()}
        identities = control["target_pair"]
        acceptor = lookup[(identities[0][0], identities[0][3])]
        hydrogen = lookup[(identities[1][0], identities[1][3])]
        parent = next(a for a in pdb.topology.atoms() if a.element.symbol == "N")
        coords = np.asarray(pdb.positions.value_in_unit(nanometer)) * 10
        left, right = coords[parent.index] - coords[hydrogen.index], coords[acceptor.index] - coords[hydrogen.index]
        angle = np.degrees(np.arccos(np.clip(np.dot(left, right) / np.linalg.norm(left) / np.linalg.norm(right), -1, 1)))
        assert np.linalg.norm(right) == pytest.approx(1.8, abs=.002)
        assert angle == pytest.approx(control["donor_hydrogen_acceptor_angle_degrees"], abs=.1)
        assert pdb.topology.getNumAtoms() == 17  # methylamine CH5N + acetone C3H6O


def test_native_failure_logs_preserve_invalid_utf8(tmp_path):
    path = tmp_path / "stdout.raw"
    raw = b"native failure: \xff\xf0\n"
    text = load_script("audit_xtb_controls").save_stream(path, raw)
    assert path.read_bytes() == raw
    assert "native failure" in text
    assert path.with_suffix(".txt").read_text() == text


def test_frozen_ranges_match_production_convention():
    spec = load_script("audit_xtb_controls").frozen_atom_spec
    assert spec([0, 1, 2, 5, 7, 8]) == "1-3,6,8-9"
    assert spec(range(180)) == "1-180"
    assert spec([]) == ""


def test_bonded_controls_have_declared_geometry_and_graph_separation(tmp_path):
    from openmm.app import PDBFile
    from openmm.unit import nanometer

    controls = load_script("audit_probe_contacts").build_bonded_controls(tmp_path)
    for name, control in controls.items():
        pdb = PDBFile(str(control["path"]))
        atoms = {(a.residue.chain.id, a.residue.id, "", a.name): a for a in pdb.topology.atoms()}
        a, b = (atoms[identity] for identity in control["target_pair"])
        coords = np.asarray(pdb.positions.value_in_unit(nanometer)) * 10
        distance = np.linalg.norm(coords[a.index] - coords[b.index])
        assert distance == pytest.approx(control["target_distance_angstrom"], abs=.002)
        neighbors = {}
        for left, right in pdb.topology.bonds():
            neighbors.setdefault(left.index, set()).add(right.index)
            neighbors.setdefault(right.index, set()).add(left.index)
        frontier, seen, depth = {a.index}, set(), 0
        while b.index not in frontier:
            seen |= frontier
            frontier = {j for i in frontier for j in neighbors.get(i, ())} - seen
            assert frontier
            depth += 1
        assert depth == control["target_bond_graph_distance"]
        assert (depth >= 4 and 2.4 < distance < 2.9) if "octane" in name else depth == 1


def test_portable_runner_checks_sources_and_binary_labels(tmp_path):
    audit = load_script("run_probe_xtb_audit")
    assert len(audit.verified_sources()) == 4
    executable = tmp_path / "xtb"
    executable.write_text("#!/bin/sh\nexit 0\n")
    executable.chmod(0o700)
    assert audit.binary_bindings([f"fixed={executable}"]) == {"fixed": str(executable)}
    with pytest.raises(ValueError):
        audit.binary_bindings([f"../unsafe={executable}"])
    with pytest.raises(ValueError):
        audit.binary_bindings([f"fixed={executable}", f"fixed={executable}"])


def test_geometry_deltas_detect_reflection_without_bond_distortion():
    from openmm.app import Element, Topology

    top = Topology()
    residue = top.addResidue("LIG", top.addChain("a"), "82", insertionCode="A")
    atoms = [top.addAtom(str(i), Element.getBySymbol(s), residue)
             for i, s in enumerate(["C", "N", "O", "C", "H"])]
    for atom in atoms[1:]:
        top.addBond(atoms[0], atom)
    before = np.array([[0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1.]])
    delta = load_script("audit_xtb_geometry_deltas").geometry_deltas
    reflected = before * [-1, 1, 1]
    report = delta(top, before, reflected)
    assert report["bond_count"] == 4
    assert report["angle_count"] == 6
    assert len(report["potential_stereovolume_flips"]) == 1
    assert report["bond_delta_rms_angstrom"] == 0
    assert report["angle_delta_rms_degrees"] == 0
    moved = before.copy()
    moved[1] *= 2
    assert len(delta(top, before, moved)["bonds_changed_over_0_1_angstrom"]) == 1
    with pytest.raises(ValueError, match="shape"):
        delta(top, before, moved[:2])
    moved[0, 0] = np.nan
    with pytest.raises(ValueError, match="Nonfinite"):
        delta(top, before, moved)


def test_subsystem_components_preserve_atom_order_and_isolates():
    from openmm.app import Element, Topology

    top = Topology()
    residue = top.addResidue("LIG", top.addChain("a"), "1")
    atoms = [top.addAtom(str(i), Element.getBySymbol("C"), residue) for i in range(5)]
    top.addBond(atoms[0], atoms[2])
    top.addBond(atoms[1], atoms[3])
    assert load_script("audit_xtb_controls").subsystem_components(top) == [[0, 2], [1, 3], [4]]
