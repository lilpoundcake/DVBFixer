"""Collapsed source geometry must not be hidden by bonded-pair exclusions."""

import pytest
from openmm import Vec3
from openmm.app import Topology, element
from openmm.unit import nanometer

from dvbfixer.diagnose.structural import check_coincident_heavy_atoms
from dvbfixer.ffutils.geometry import detect_coincident_heavy_atoms


def test_bonded_cross_chain_pairs_keep_case_and_insertion_code_without_mutation():
    topology = Topology()
    atoms = []
    for chain_id in ("D", "d"):
        residue = topology.addResidue("ALA", topology.addChain(chain_id), "82", "A")
        atoms.append(topology.addAtom("CA", element.carbon, residue))
    topology.addBond(*atoms)
    positions = [Vec3(-0.01, 0, 0), Vec3(0.01, 0, 0)] * nanometer
    before = tuple(tuple(v) for v in positions.value_in_unit(nanometer))
    pair, = detect_coincident_heavy_atoms(topology, positions)
    assert pair.atom1.residue.chain.id == "D"
    assert pair.atom2.residue.chain.id == "d"
    assert pair.distance_nm == pytest.approx(0.02)
    finding, = check_coincident_heavy_atoms(topology, positions)
    assert [atom["chain"] for atom in finding.extra["atoms"]] == ["D", "d"]
    assert all(atom["icode"] == "A" for atom in finding.extra["atoms"])
    assert tuple(tuple(v) for v in positions.value_in_unit(nanometer)) == before
    assert list(topology.bonds()) == [(atoms[0], atoms[1])]


def test_malformed_input_pair_report_is_bounded():
    topology = Topology()
    residue = topology.addResidue("ALA", topology.addChain("A"), "1")
    for index in range(50):
        topology.addAtom(f"C{index}", element.carbon, residue)
    positions = [Vec3(0, 0, 0)] * 50 * nanometer
    assert len(detect_coincident_heavy_atoms(topology, positions, max_pairs=7)) == 7
    findings = check_coincident_heavy_atoms(topology, positions)
    assert len(findings) == 1001
    assert "truncated" in findings[-1].message


def test_grid_matches_brute_force_without_hydrogens_or_duplicate_pairs():
    import random

    topology = Topology()
    residue = topology.addResidue("UNK", topology.addChain("A"), "1")
    generator = random.Random(42)
    points = [Vec3(*(generator.uniform(-0.15, 0.15) for _ in range(3))) for _ in range(150)]
    for index in range(len(points)):
        topology.addAtom(str(index), element.hydrogen if index % 7 == 0 else element.carbon, residue)
    atoms = list(topology.atoms())
    expected = {(i, j) for j in range(len(points)) for i in range(j)
                if atoms[i].element is not element.hydrogen and atoms[j].element is not element.hydrogen
                and sum((points[i][k] - points[j][k]) ** 2 for k in range(3)) < 0.05 ** 2}
    actual = detect_coincident_heavy_atoms(topology, points * nanometer)
    assert {(pair.atom1.index, pair.atom2.index) for pair in actual} == expected
    assert len(actual) == len(expected)


def test_ordinary_diagnose_dispatches_shared_heavy_atom_check(tmp_path):
    from openmm.app import PDBFile

    from dvbfixer.diagnose.structural import run_all

    topology = Topology()
    residue = topology.addResidue("ALA", topology.addChain("A"), "1")
    topology.addAtom("CA", element.carbon, residue)
    topology.addAtom("CB", element.carbon, residue)
    positions = [Vec3(0, 0, 0), Vec3(0, 0, 0)] * nanometer
    path = tmp_path / "collapsed.pdb"
    with path.open("w") as handle:
        PDBFile.writeFile(topology, positions, handle, keepIds=True)
    findings = run_all(topology, positions, path)
    assert any(f.category == "coincident_atoms" and f.extra.get("metric") == "heavy_atom_coincidence.v1"
               for f in findings)


@pytest.mark.parametrize("tolerance", [0, -0.1, float("nan"), float("inf")])
def test_invalid_tolerance_is_not_silently_accepted(tolerance):
    with pytest.raises(ValueError):
        detect_coincident_heavy_atoms(Topology(), [] * nanometer, tol_nm=tolerance)


def test_invalid_pair_limit_is_not_silently_accepted():
    with pytest.raises(ValueError):
        detect_coincident_heavy_atoms(Topology(), [] * nanometer, max_pairs=0)
