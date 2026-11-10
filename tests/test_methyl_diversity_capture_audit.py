"""Captured graph/lineage summaries fail closed on index/ownership changes."""
import copy

import pytest
from scripts.audit_methyl_diversity_capture import bijection, validate_graph


def fixture():
    rows = [dict(key=['a', '82', 'A', 'THR', n], element=e) for n, e in
        [('CG2', 'C'), ('CB', 'C'), ('HG21', 'H'), ('HG22', 'H'), ('HG23', 'H'), ('CA', 'C')]]
    bonds = [[0, 1], [0, 2], [0, 3], [0, 4], [1, 5]]
    group = dict(parent=0, axis=1, frame=5, hydrogens=[2, 3, 4], parent_key=rows[0]['key'], frame_key=rows[5]['key'])
    return rows, bonds, dict(new_h=[2, 3, 4], mobile=[2, 3, 4], graph=dict(groups=[group]))


def test_independent_enumeration_and_frame_check():
    rows, bonds, ownership = fixture()
    before = copy.deepcopy(ownership)
    assert validate_graph(rows, bonds, ownership) == 1
    assert ownership == before
    bijection([2, 0, 1], 3)


@pytest.mark.parametrize('fault', ['omit', 'frame', 'ownership', 'axis'])
def test_incomplete_graphs_are_not_accepted(fault):
    rows, bonds, ownership = fixture()
    if fault == 'omit':
        ownership['graph']['groups'] = []
    elif fault == 'frame':
        ownership['graph']['groups'][0]['frame'] = 0
    elif fault == 'ownership':
        ownership['new_h'] = [2, 3]
    else:
        ownership['graph']['groups'][0]['axis'] = 5
    with pytest.raises(ValueError):
        validate_graph(rows, bonds, ownership)


@pytest.mark.parametrize('mapping', [[0, 0], [0], [0, 2], [False, 1]])
def test_index_lineage_must_be_complete_bijection(mapping):
    with pytest.raises(ValueError):
        bijection(mapping, 2)
