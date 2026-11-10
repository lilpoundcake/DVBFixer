"""JSON container normalization must not weaken exact scalar comparisons."""
import pytest
from scripts.research_mean_pair_replay import transport


def test_tuple_inventory_matches_json_without_rounding_a_score():
    value=dict(import_audit=dict(extras=[('A',63,' ','THR','OXT')]),score=1.2345678901234567)
    result=transport(value)
    assert result['import_audit']['extras']==[['A',63,' ','THR','OXT']]
    assert result['score']==value['score']
    assert result['score']!=value['score']+1e-12


def test_nonfinite_metadata_is_not_hidden_by_serialization():
    with pytest.raises(ValueError):
        transport(dict(score=float('nan')))
