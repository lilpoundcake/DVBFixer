"""Synthetic all-disagreement enrollment and explicit accounting closure policy."""
import pytest
from scripts.research_mean_pair_energies import closure, disagreements


def test_all_disagreements_not_native_error_or_margin_filter():
    geometry=dict(results={n:dict(status='generated',prior_selected_index=a,native_error=e)
        for n,a,e in [('first',0,0.),('second',1,999.),('same',0,1.)]})
    scores=dict(results={n+'-means':dict(selection=dict(selected_index=b)) for n,b in [('first',1),('second',0),('same',0)]})
    assert [r['id'] for r in disagreements(geometry,scores)]==['first','second']


def test_closed_and_incomplete_accounting_are_distinct_not_silently_adjusted():
    assert closure([1.,2.],3.)['closed']
    result=closure([1.,2.],4.)
    assert not result['closed'] and result['pair_sum']==3. and result['recorded_total']==4.


def test_nonfinite_accounting_is_not_a_pass():
    with pytest.raises(ValueError,match='nonfinite'):
        closure([float('nan')],0.)
