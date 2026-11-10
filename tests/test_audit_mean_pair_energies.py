"""Independent saved-pair deltas preserve absent edges and context accounting."""
import pytest
from scripts.audit_mean_pair_energies import difference, pair_map

TARGET=['A',1,' ','VAL']
OTHER=['A',2,' ','ALA']


def edge(value):
    return dict(first=TARGET,second=OTHER,weighted_terms={'fa_rep':value})


def test_added_removed_and_changed_target_edge_accounting():
    assert difference([], [edge(2.)],['fa_rep'],TARGET)[0]['weighted_deltas']=={'fa_rep':2.}
    assert difference([edge(2.)], [],['fa_rep'],TARGET)[0]['weighted_deltas']=={'fa_rep':-2.}
    result=difference([edge(2.)],[edge(1.)],['fa_rep'],TARGET)
    assert result[0]['contains_target'] and result[0]['weighted_deltas']=={'fa_rep':-1.}


def test_unchanged_pair_is_not_a_changed_context_edge():
    assert difference([edge(2.)],[edge(2.)],['fa_rep'],TARGET)==[]


def test_duplicate_or_incomplete_term_inventory_refuses():
    with pytest.raises(ValueError,match='duplicate pair'):
        pair_map([edge(1.),edge(2.)],['fa_rep'])
    with pytest.raises(ValueError,match='term inventory'):
        pair_map([edge(1.)],['fa_atr'])
