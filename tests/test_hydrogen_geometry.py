"""Deterministic coverage of stretched H bonds after hydrogen placement."""

import numpy as np
import pytest

pytest.importorskip("openmm")

from openmm import Vec3  # noqa: E402
from openmm.app import Topology, element  # noqa: E402
from openmm.unit import nanometer  # noqa: E402

from dvbfixer.ffutils.geometry import (  # noqa: E402
    detect_misplaced_hydrogens,
    repair_misplaced_hydrogens,
)


@pytest.mark.parametrize(
    ("symbol", "bond_nm"),
    [("O", 0.097), ("N", 0.101), ("C", 0.109), ("S", 0.134), ("P", 0.142)],
)
@pytest.mark.parametrize("stretch", [1.0, 1.1, 1.29])
def test_hydrogen_parent_distance(symbol: str, bond_nm: float, stretch: float) -> None:
    """Catch the CI's ~1.25 Å O-H bond without depending on OpenMM's RNG.

    Normal bonds (including a 10% tolerance) must remain untouched. Bonds
    stretched by 29% used to escape the 60% repair threshold entirely.
    """
    topology = Topology()
    residue = topology.addResidue("UNK", topology.addChain("A"), "1")
    neighbor = topology.addAtom("C", element.carbon, residue)
    parent = topology.addAtom(symbol, element.get_by_symbol(symbol), residue)
    hydrogen = topology.addAtom("H", element.hydrogen, residue)
    topology.addBond(neighbor, parent)
    topology.addBond(parent, hydrogen)
    positions = [Vec3(-0.15, 0, 0), Vec3(0, 0, 0), Vec3(bond_nm * stretch, 0, 0)] * nanometer
    original = np.array(positions.value_in_unit(nanometer))

    findings = detect_misplaced_hydrogens(topology, positions)
    if stretch <= 1.1:
        assert findings == []
        assert repair_misplaced_hydrogens(topology, positions) == 0
        np.testing.assert_array_equal(positions.value_in_unit(nanometer), original)
        return

    assert len(findings) == 1
    assert findings[0].reason == "parent_distance"
    assert findings[0].h_atom is hydrogen
    assert repair_misplaced_hydrogens(topology, positions) == 1
    repaired = np.array(positions.value_in_unit(nanometer))
    np.testing.assert_array_equal(repaired[:2], original[:2])
    assert np.linalg.norm(repaired[2] - repaired[1]) == pytest.approx(bond_nm)
    assert detect_misplaced_hydrogens(topology, positions) == []
