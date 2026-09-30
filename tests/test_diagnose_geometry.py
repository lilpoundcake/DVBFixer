"""Pure geometry and exact-scope contracts for diagnose."""

from __future__ import annotations

import math

import pytest

from dvbfixer.diagnose.cli import _parse_boundary_selector, _parse_residue_selector
from dvbfixer.diagnose.geometry import (
    AtomIdentity,
    BoundaryIdentity,
    GeometryMeasurement,
    ResidueIdentity,
    angle_degrees,
    boundary_measurements,
    distance_angstrom,
    signed_dihedral_degrees,
)


def test_distance_uses_explicit_angstrom_coordinates() -> None:
    assert distance_angstrom((0, 0, 0), (3, 4, 0)) == 5.0
    assert distance_angstrom((math.nan, 0, 0), (0, 0, 0)) is None


def test_angle_clamps_and_rejects_degenerate_geometry() -> None:
    assert angle_degrees((1, 0, 0), (0, 0, 0), (0, 1, 0)) == pytest.approx(90.0)
    assert angle_degrees((0, 0, 0), (0, 0, 0), (1, 0, 0)) is None
    assert angle_degrees((math.inf, 0, 0), (0, 0, 0), (1, 0, 0)) is None


def test_signed_dihedral_distinguishes_sign_and_degeneracy() -> None:
    positive = signed_dihedral_degrees((1, 0, 0), (0, 0, 0), (0, 1, 0), (0, 1, 1))
    negative = signed_dihedral_degrees((1, 0, 0), (0, 0, 0), (0, 1, 0), (0, 1, -1))
    assert positive == pytest.approx(-90.0)
    assert negative == pytest.approx(90.0)
    assert signed_dihedral_degrees((0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)) is None


def test_full_identities_are_case_and_insertion_code_sensitive() -> None:
    assert ResidueIdentity("D", "82", "") != ResidueIdentity("d", "82", "")
    assert ResidueIdentity("D", "82", "") != ResidueIdentity("D", "82", "A")
    assert AtomIdentity("D", "82", "A", "CA") != AtomIdentity("D", "82", "", "CA")
    assert ResidueIdentity("", "1").to_dict()["chain"] == ""


def test_scope_selector_is_exact_and_boundary_is_ordered() -> None:
    assert _parse_residue_selector("d:82:A") == ResidueIdentity("d", "82", "A")
    boundary = _parse_boundary_selector("D:82,d:82:A")
    assert boundary.left == ResidueIdentity("D", "82", "")
    assert boundary.right == ResidueIdentity("d", "82", "A")
    with pytest.raises(Exception):
        _parse_residue_selector("D")


def test_measurement_rejects_non_finite_values() -> None:
    atom = AtomIdentity("A", "1", "", "CA")
    with pytest.raises(ValueError):
        GeometryMeasurement("distance.v1", (atom,), math.inf, "angstrom")


def test_measurement_requires_version_unit_atoms_and_ordered_range() -> None:
    atom = AtomIdentity("A", "1", "", "CA")
    with pytest.raises(ValueError):
        GeometryMeasurement("distance", (atom,), 1.0, "angstrom")
    with pytest.raises(ValueError):
        GeometryMeasurement("distance.v1", (), 1.0, "angstrom")
    with pytest.raises(ValueError):
        GeometryMeasurement("distance.v1", (atom,), 1.0, "")
    with pytest.raises(ValueError):
        GeometryMeasurement("distance.v1", (atom,), 1.0, "angstrom", expected_min=2, expected_max=1)


def test_boundary_measurement_reports_bond_and_exact_insertion_codes() -> None:
    from openmm.app import Element, Topology
    from openmm.unit import Quantity, nanometer

    topology = Topology()
    chain = topology.addChain("D")
    left = topology.addResidue("ALA", chain, id="82")
    right = topology.addResidue("GLY", chain, id="82", insertionCode="A")
    c_atom = topology.addAtom("C", Element.getBySymbol("C"), left)
    n_atom = topology.addAtom("N", Element.getBySymbol("N"), right)
    topology.addBond(c_atom, n_atom)
    positions = Quantity([(0, 0, 0), (0.133, 0, 0)], nanometer)
    boundary = BoundaryIdentity(ResidueIdentity("D", "82", ""), ResidueIdentity("D", "82", "A"))

    result = boundary_measurements(topology, positions, [boundary])[0]

    assert result.status == "measured"
    assert result.topology_bonded is True
    assert result.c_n_distance is not None
    assert result.c_n_distance.value == pytest.approx(1.33)
    assert result.c_n_distance.atoms[1].icode == "A"


def test_boundary_measurement_does_not_guess_missing_or_duplicate_atoms() -> None:
    from openmm.app import Element, Topology
    from openmm.unit import Quantity, nanometer

    topology = Topology()
    chain = topology.addChain("d")
    left = topology.addResidue("ALA", chain, id="82")
    right = topology.addResidue("GLY", chain, id="83")
    topology.addAtom("C", Element.getBySymbol("C"), left)
    topology.addAtom("N", Element.getBySymbol("N"), right)
    topology.addAtom("N", Element.getBySymbol("N"), right)
    positions = Quantity([(0, 0, 0), (0.133, 0, 0), (0.14, 0, 0)], nanometer)
    boundary = BoundaryIdentity(ResidueIdentity("d", "82"), ResidueIdentity("d", "83"))

    result = boundary_measurements(topology, positions, [boundary])[0]

    assert result.status == "ambiguous_atom"
    assert result.c_n_distance is None
    assert result.atom_match_counts == (1, 2)
