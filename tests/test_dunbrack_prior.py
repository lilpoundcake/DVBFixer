"""Synthetic family-only ranking, no installed-library data in tests."""
import copy

import pytest
from scripts.research_dunbrack_prior import select_prior


def candidates():
    return [dict(index=0,well=3,common_accepted=True),dict(index=1,well=1,common_accepted=False),
            dict(index=2,well=2,common_accepted=True)]


def test_geometry_precedes_probability_and_native_error_never_enters_ranking():
    c=candidates()
    before=copy.deepcopy(c)
    assert select_prior(c,{1:.8,2:.15,3:.05})==2
    for row in c:
        row['native_rmsd']=100. if row['index']==2 else 0.
    assert select_prior(c,{1:.8,2:.15,3:.05})==2
    assert [{k:v for k,v in r.items() if k!='native_rmsd'} for r in c]==before


def test_stable_tie_and_empty_feasible_pool():
    c=candidates()
    assert select_prior(c,{1:0.,2:.5,3:.5})==0
    for row in c:
        row['common_accepted']=False
    assert select_prior(c,{1:0.,2:.5,3:.5}) is None


@pytest.mark.parametrize('probabilities',[{1:.5,2:.5},{1:.5,2:.5,3:.5},{1:float('nan'),2:0.,3:0.}])
def test_incomplete_nonfinite_or_unnormalized_library_result_refuses(probabilities):
    with pytest.raises(ValueError):
        select_prior(candidates(),probabilities)


def test_duplicate_family_does_not_silently_select_a_subset():
    c=candidates()
    c[2]['well']=3
    with pytest.raises(ValueError,match='one saved candidate per family'):
        select_prior(c,{1:.2,2:.3,3:.5})
