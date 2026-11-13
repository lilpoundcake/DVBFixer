"""Retained-input adapters only; no native sampling/optimization."""
from dataclasses import asdict

import pytest
from scripts.hydrogen_initialization_guard import Atom, Input
from scripts.research_post_h_after_repair import distance_findings, input_from_table, reconstruct

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def example():
    component = ComponentInstanceRef(1, 'A', '0')
    atoms = tuple(Atom(ExactAtomRef(component, n), 'SER', e) for n, e in
                  [('CA', 'C'), ('N', 'N'), ('C', 'C'), ('CB', 'C'), ('HA', 'H')])
    positions = ((0., 0., 0.), (0.14, 0., 0.01), (-0.07, 0.121, 0.01),
                 (-0.07, -0.121, 0.01), (0., 0., 0.21))
    bonds = ((0, 1), (0, 2), (0, 3), (0, 4))
    return Input(atoms, positions, bonds, frozenset([atoms[4].identity]))


def test_all_unsupported_distance_findings_are_accounted_without_changing_mask():
    value = example()
    findings = distance_findings(value)
    assert len(findings) == 1
    assert findings[0]['radius_nm'] == pytest.approx(0.21)
    assert value.new_hydrogens == frozenset([value.atoms[4].identity])


def test_raw_table_roundtrip_rebinds_graph_and_preserves_coordinates():
    from openmm.unit import nanometer

    value = example()
    table = dict(atoms=[asdict(a) for a in value.atoms], positions_nm=value.positions_nm, bonds=value.bonds)
    metadata = [(a, b, None, None) for a, b in value.bonds]
    parsed = input_from_table(table)
    assert parsed.atoms == value.atoms and parsed.positions_nm == value.positions_nm
    topology, positions = reconstruct(table, metadata)
    assert tuple(tuple(p) for p in positions.value_in_unit(nanometer)) == value.positions_nm
    actual_atoms = set(topology.atoms())
    assert all(a in actual_atoms and b in actual_atoms for a, b in topology.bonds())
