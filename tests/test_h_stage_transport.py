"""Exact identity/bond transport and native post-H repair borrowing controls."""
import copy

import pytest
from scripts.research_h_stage_transport import post_h_repair, topology_from_rows


def source():
    rows = [dict(key=['A', '24', '', 'THR', n], element=e)
            for n, e in [('CG2', 'C'), ('CB', 'C'), ('HG21', 'H'), ('HG22', 'H'), ('HG23', 'H')]]
    bonds = [[0, 1], [0, 2], [0, 3], [0, 4]]
    points = [[0., 0., 0.], [0., 0., .15], [0., 0., -.1], [0., 0., -.1], [0., 0., -.1]]
    return rows, bonds, points


def test_existing_helper_repairs_hydrogens_without_mutating_borrowed_heavy_input():
    pytest.importorskip('openmm')
    rows, bonds, points = source()
    original = copy.deepcopy(points)
    after, repaired = post_h_repair(rows, bonds, points)
    assert repaired >= 2
    assert after[:2] == points[:2] == original[:2]
    assert len({tuple(p) for p in after[2:]}) == 3
    assert points == original


def test_topology_preserves_case_sensitive_chains_residue_numbers_and_icodes():
    pytest.importorskip('openmm')
    rows = [dict(key=[chain, '82', icode, 'GLN', 'CA'], element='C')
            for chain, icode in [('D', ''), ('d', ''), ('d', 'A')]]
    top = topology_from_rows(rows, [])
    assert [c.id for c in top.chains()] == ['D', 'd']
    assert [(r.chain.id, r.id, r.insertionCode) for r in top.residues()] == [('D', '82', ''), ('d', '82', ''), ('d', '82', 'A')]
    assert [a.index for a in top.atoms()] == [0, 1, 2]


@pytest.mark.parametrize('fault', ['duplicate_atom', 'duplicate_bond', 'dangling_bond', 'self_bond'])
def test_ambiguous_transport_is_refused(fault):
    pytest.importorskip('openmm')
    rows, bonds, _ = source()
    if fault == 'duplicate_atom':
        rows[-1] = copy.deepcopy(rows[0])
    elif fault == 'duplicate_bond':
        bonds.append([1, 0])
    elif fault == 'dangling_bond':
        bonds.append([0, 100])
    else:
        bonds.append([0, 0])
    with pytest.raises(ValueError):
        topology_from_rows(rows, bonds)
