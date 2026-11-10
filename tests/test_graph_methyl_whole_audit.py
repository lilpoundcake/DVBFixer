"""Native-volume shape descriptors preserve thresholds, signs and identities."""
import pytest
from scripts.audit_graph_methyl_whole import global_shape


def source(volume):
    return [dict(key=['d', '82', 'A', 'THR', name], element='C', xyz=point) for name, point in
        [('CA', [0., 0., 0.]), ('N', [1., 0., 0.]), ('C', [0., 1., 0.]), ('CB', [0., 0., volume])]]


def test_shape_threshold_is_the_existing_half_A3_not_a_new_complete_gate():
    ordinary = global_shape(source(1.))
    assert not ordinary['CA_at_or_below_0_5']
    assert ordinary['CA_signed_volumes'][0]['volume_angstrom3'] == pytest.approx(1.)
    flat = global_shape(source(.4))
    assert len(flat['CA_at_or_below_0_5']) == 1 and not flat['negative_CA']
    negative = global_shape(source(-.4))
    assert negative['negative_CA'][0]['residue'] == ('d', '82', 'A', 'THR')
    assert negative['negative_CA'][0]['volume_angstrom3'] == pytest.approx(-.4)
    assert 'accepted' not in negative


def test_rounded_descriptor_uses_PDB_precision_without_mutating_observations():
    rows = source(.5004)
    assert not global_shape(rows)['CA_at_or_below_0_5']
    assert len(global_shape(rows, rounded=True)['CA_at_or_below_0_5']) == 1
    assert rows[-1]['xyz'][-1] == .5004
    assert global_shape(rows)['CA_signed_volumes'][0]['volume_angstrom3'] == pytest.approx(.5004)
