"""Raw-weight amendment: preserve API outputs, no native/geometry tuning."""
import pytest
from scripts.research_dunbrack_prior_weights import select_weights


def rows():
    return [dict(index=i,well=i+1,common_accepted=True) for i in range(3)]


def test_values_above_one_and_unnormalized_sums_are_weights_not_probabilities():
    weights={1:1.0358,2:.005,3:.0006}
    assert select_weights(rows(),weights)==0
    assert weights=={1:1.0358,2:.005,3:.0006}
    assert select_weights(rows(),{k:v/sum(weights.values()) for k,v in weights.items()})==0


def test_feasibility_precedes_weight_and_stable_index_breaks_ties():
    c=rows()
    c[0]['common_accepted']=False
    assert select_weights(c,{1:100.,2:1.,3:1.})==1
    for r in c:
        r['common_accepted']=False
    assert select_weights(c,{1:100.,2:1.,3:1.}) is None


@pytest.mark.parametrize('weights',[{1:0.,2:0.,3:0.},{1:-1.,2:2.,3:1.},{1:float('inf'),2:1.,3:1.},{1:1.,2:1.}])
def test_invalid_or_missing_raw_weights_still_fail_closed(weights):
    with pytest.raises(ValueError):
        select_weights(rows(),weights)
