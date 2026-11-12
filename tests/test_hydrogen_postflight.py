"""Oxygen-specific distance contract; no force field or stochastic sampling."""
import math

import pytest

openmm = pytest.importorskip('openmm')
from openmm import Vec3, unit  # noqa: E402
from openmm.app import Topology, element  # noqa: E402

from dvbfixer.ffutils.geometry import (  # noqa: E402
    detect_misplaced_hydrogens,
    repair_misplaced_hydrogens,
)

CHEMISTRIES = ('SER', 'THR', 'TYR', 'ASH', 'GLH', 'HOH', 'LIG', 'NAG')
DISTANCES = (0.97, 1.199, 1.2, 1.201, 1.2458133889150464, 1.7)


def oh_graph(resname='SER', distance=1.2458133889150464, chain='B', icode=''):
    topology = Topology()
    residue = topology.addResidue(resname, topology.addChain(chain), '126', icode)
    oxygen_name, hydrogen_name, carbon_name = {
        'SER': ('OG', 'HG', 'CB'), 'THR': ('OG1', 'HG1', 'CB'),
        'TYR': ('OH', 'HH', 'CZ'), 'ASH': ('OD2', 'HD2', 'CG'),
        'GLH': ('OE2', 'HE2', 'CD'), 'HOH': ('O', 'H1', None),
        'LIG': ('O1', 'H1', 'C1'), 'NAG': ('O6', 'HO6', 'C6'),
    }[resname]
    oxygen = topology.addAtom(oxygen_name, element.oxygen, residue)
    hydrogen = topology.addAtom(hydrogen_name, element.hydrogen, residue)
    topology.addBond(oxygen, hydrogen)
    positions = [Vec3(0., 0., 0.), Vec3(distance/10, 0., 0.)]
    if resname == 'HOH':
        sibling = topology.addAtom('H2', element.hydrogen, residue)
        topology.addBond(oxygen, sibling)
        theta = math.radians(104.5)
        positions.append(Vec3(0.097*math.cos(theta), 0.097*math.sin(theta), 0.))
    else:
        carbon = topology.addAtom(carbon_name, element.carbon, residue)
        topology.addBond(oxygen, carbon)
        # 109.5-degree C-O-H, a non-clashing direction whose torsion can
        # be preserved even for the deliberately over-long H bond.
        theta = math.radians(109.5)
        positions.append(Vec3(0.14*math.cos(theta), 0., 0.14*math.sin(theta)))
        tail = topology.addAtom('C_tail', element.carbon, residue)
        topology.addBond(carbon, tail)
        positions.append(Vec3(-0.18, 0., 0.20))
        if resname in ('ASH', 'GLH'):
            other_oxygen = topology.addAtom('O_other', element.oxygen, residue)
            topology.addBond(carbon, other_oxygen)
            positions.append(Vec3(-0.14, -0.09, 0.12))
        elif resname in ('THR', 'TYR'):
            branch = topology.addAtom('C_branch', element.carbon, residue)
            topology.addBond(carbon, branch)
            positions.append(Vec3(-0.04, -0.15, 0.14))
    return topology, positions*unit.nanometer


def xyz(positions):
    return [list(p) for p in positions.value_in_unit(unit.nanometer)]


def test_moderate_oh_is_detected_and_radially_repaired():
    topology, positions = oh_graph()
    before = xyz(positions)
    findings = detect_misplaced_hydrogens(topology, positions)
    assert len(findings) == 1
    assert findings[0].reason == 'parent_distance'
    assert repair_misplaced_hydrogens(topology, positions) == 1
    after = xyz(positions)
    assert after[0] == before[0] and after[2] == before[2]
    assert after[1] == pytest.approx([0.097, 0., 0.])
    assert detect_misplaced_hydrogens(topology, positions) == []
    assert repair_misplaced_hydrogens(topology, positions) == 0


@pytest.mark.parametrize('resname', CHEMISTRIES)
@pytest.mark.parametrize('distance', DISTANCES)
def test_complete_oh_local_graph_grid(resname, distance):
    topology, positions = oh_graph(resname, distance)
    before = xyz(positions)
    count = repair_misplaced_hydrogens(topology, positions)
    after = xyz(positions)
    if distance < 1.2:
        assert count == 0
        assert after == before
    else:
        assert count == 1
        assert 10*math.dist(after[0], after[1]) == pytest.approx(0.97)
        assert after[0] == before[0]
        assert after[2:] == before[2:]
        if distance <= 1.552 or resname == 'HOH':
            assert after[1] == pytest.approx([0.097, 0., 0.])
    assert detect_misplaced_hydrogens(topology, positions) == []


@pytest.mark.parametrize('symbol,distance', [('C', 1.2458133889150464), ('N', 1.2458133889150464), ('S', 1.34), ('P', 1.42)])
def test_nonoxygen_distance_policy_is_not_tightened(symbol, distance):
    topology = Topology()
    residue = topology.addResidue('LIG', topology.addChain('a'), '126', 'A')
    parent = topology.addAtom('X', element.Element.getBySymbol(symbol), residue)
    hydrogen = topology.addAtom('H', element.hydrogen, residue)
    topology.addBond(parent, hydrogen)
    positions = [Vec3(0., 0., 0.), Vec3(distance/10, 0., 0.)]*unit.nanometer
    before = xyz(positions)
    assert detect_misplaced_hydrogens(topology, positions) == []
    assert repair_misplaced_hydrogens(topology, positions) == 0
    assert xyz(positions) == before


def test_radial_collision_uses_existing_tetrahedral_fallback():
    topology = Topology()
    residue = topology.addResidue('SER', topology.addChain('b'), '126', 'A')
    oxygen = topology.addAtom('OG', element.oxygen, residue)
    hydrogen = topology.addAtom('HG', element.hydrogen, residue)
    carbon = topology.addAtom('CB', element.carbon, residue)
    topology.addBond(oxygen, hydrogen)
    topology.addBond(oxygen, carbon)
    positions = [Vec3(0., 0., 0.), Vec3(0.12458133889150464, 0., 0.), Vec3(0.14, 0., 0.)]*unit.nanometer
    before = xyz(positions)
    assert repair_misplaced_hydrogens(topology, positions) == 1
    after = xyz(positions)
    assert after[1] != pytest.approx([0.097, 0., 0.])
    assert after[0] == before[0] and after[2] == before[2]
    assert detect_misplaced_hydrogens(topology, positions) == []


def test_unrepairable_long_oh_fails_closed_without_moving_heavy_atoms():
    from dvbfixer.ffutils.geometry import HydrogenGeometryError

    topology = Topology()
    residue = topology.addResidue('HOH', topology.addChain('b'), '126', 'A')
    oxygen = topology.addAtom('O', element.oxygen, residue)
    hydrogen = topology.addAtom('H1', element.hydrogen, residue)
    topology.addAtom('blocker', element.carbon, residue)  # Deliberately unbonded.
    topology.addBond(oxygen, hydrogen)
    positions = [Vec3(0., 0., 0.), Vec3(0.12458133889150464, 0., 0.), Vec3(0.097, 0., 0.)]*unit.nanometer
    before = xyz(positions)
    with pytest.raises(HydrogenGeometryError, match='b/HOH126A/H1'):
        repair_misplaced_hydrogens(topology, positions)
    assert xyz(positions) == before


def test_chain_case_and_insertion_code_do_not_cross_contaminate():
    topology = Topology()
    positions = []
    targets = []
    for chain, icode, distance, offset in [('B', '', 0.97, 0.), ('B', 'A', 0.97, 1.), ('b', '', 1.2458133889150464, 2.)]:
        residue = topology.addResidue('HOH', topology.addChain(chain), '126', icode)
        oxygen = topology.addAtom('O', element.oxygen, residue)
        hydrogen = topology.addAtom('H1', element.hydrogen, residue)
        topology.addBond(oxygen, hydrogen)
        positions.extend([Vec3(offset, 0., 0.), Vec3(offset+distance/10, 0., 0.)])
        targets.append(hydrogen.index)
    quantity = positions*unit.nanometer
    before = xyz(quantity)
    findings = detect_misplaced_hydrogens(topology, quantity)
    assert len(findings) == 1
    assert findings[0].h_atom.index == targets[-1]
    assert repair_misplaced_hydrogens(topology, quantity) == 1
    after = xyz(quantity)
    assert after[:5] == before[:5]
    assert after[5] == pytest.approx([2.097, 0., 0.])
