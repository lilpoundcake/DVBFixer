"""Independent raw-weight ranking controls without PyRosetta/table dependencies."""
import pytest
from scripts.audit_dunbrack_prior import ranking


def pool():
    return [dict(index=0,well=3,common_accepted=True),dict(index=1,well=1,common_accepted=False),
            dict(index=2,well=2,common_accepted=True)]


def test_independent_ranking_preserves_geometry_and_stable_index():
    assert ranking(pool(),[dict(well=1,probability=1.1),dict(well=2,probability=.3),dict(well=3,probability=.2)])==2
    assert ranking(pool(),[dict(well=1,probability=1.1),dict(well=2,probability=.2),dict(well=3,probability=.2)])==0


def test_empty_feasible_pool_is_not_an_arbitrary_choice():
    c=pool()
    for row in c:
        row['common_accepted']=False
    assert ranking(c,[dict(well=i,probability=.5) for i in (1,2,3)]) is None


@pytest.mark.parametrize('samples',[[dict(well=1,probability=1.)],
    [dict(well=1,probability=1.),dict(well=2,probability=-1.),dict(well=3,probability=1.)]])
def test_partial_or_negative_raw_library_result_remains_invalid(samples):
    with pytest.raises(ValueError):
        ranking(pool(),samples)
