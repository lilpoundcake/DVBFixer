"""Saved total-score arithmetic only; no PyRosetta or native data."""
import pytest
from scripts.analyze_dunbrack_mean_terms import deltas


def score(a,b):
    return dict(score_terms={'a':dict(raw=a,weight=2.),'b':dict(raw=b,weight=.5)},total_score=2*a+.5*b)


def test_weighted_deltas_reconstruct_saved_total_and_sort_both_signs():
    result=deltas(score(3.,1.),score(1.,5.))
    assert result['total_delta']==-2.
    assert [r['weighted_delta'] for r in result['terms']]==[-4.,2.]


def test_changed_weights_and_inconsistent_totals_refuse():
    a,b=score(1.,2.),score(3.,4.)
    b['total_score']+=1.
    with pytest.raises(ValueError,match='delta sum'):
        deltas(a,b)
    b=score(3.,4.)
    b['score_terms']['a']['weight']=3.
    with pytest.raises(ValueError,match='changed score weights'):
        deltas(a,b)
