"""Nominal, moderate and gross nonoxygen bonds retain their existing policy."""
import math

import pytest


@pytest.mark.parametrize('symbol,target,moderate,gross', [
    ('C', 1.09, 1.2458133889150464, 1.9),
    ('N', 1.01, 1.2458133889150464, 1.8),
    ('S', 1.34, 1.6, 2.3),
    ('P', 1.42, 1.8, 2.5),
])
@pytest.mark.parametrize('kind', ['nominal', 'moderate', 'gross'])
def test_nonoxygen_nominal_moderate_and_gross(symbol, target, moderate, gross, kind):
    pytest.importorskip('openmm')
    from openmm import Vec3, unit
    from openmm.app import Topology, element

    from dvbfixer.ffutils.geometry import (
        detect_misplaced_hydrogens,
        repair_misplaced_hydrogens,
    )

    distance = dict(nominal=target, moderate=moderate, gross=gross)[kind]
    topology = Topology()
    residue = topology.addResidue('LIG', topology.addChain('b'), '126', 'A')
    parent = topology.addAtom('X', element.Element.getBySymbol(symbol), residue)
    hydrogen = topology.addAtom('H', element.hydrogen, residue)
    carbon = topology.addAtom('C', element.carbon, residue)
    topology.addBond(parent, hydrogen)
    topology.addBond(parent, carbon)
    positions = [Vec3(0., 0., 0.), Vec3(distance/10, 0., 0.), Vec3(-0.15, 0., 0.)]*unit.nanometer
    before = [list(p) for p in positions.value_in_unit(unit.nanometer)]
    count = repair_misplaced_hydrogens(topology, positions)
    after = [list(p) for p in positions.value_in_unit(unit.nanometer)]
    if kind == 'gross':
        assert count == 1
        assert 10*math.dist(after[0], after[1]) == pytest.approx(target)
        assert after[0] == before[0] and after[2] == before[2]
    else:
        assert count == 0
        assert after == before
    assert detect_misplaced_hydrogens(topology, positions) == []
    assert repair_misplaced_hydrogens(topology, positions) == 0
