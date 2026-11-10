"""Pre-download complete-pool contract cannot cherry-pick guard successes."""
import pytest
from scripts.prospective_chi_guard_cohort import SYSTEMS, choose_complete


def rows():
    return [dict(index=0,status='evaluated',common_accepted=True,score=2.,guard_decision='endpoint_gradient_met'),
            dict(index=1,status='evaluated',common_accepted=True,score=1.,guard_decision='endpoint_gradient_met')]


def test_guard_exhaustion_blocks_pool_instead_of_falling_back_to_success():
    values=rows()
    values[1]['guard_decision']='budget_exhausted'
    assert choose_complete(values,guarded=True)==dict(status='incomplete_guard_pool',selected_index=None)


def test_completed_geometry_failure_remains_distinct_from_incomplete_execution():
    values=rows()
    values[1]['common_accepted']=False
    assert choose_complete(values)['selected_index']==0
    values[1]['status']='execution_blocker'
    assert choose_complete(values)['status']=='incomplete_pool'


def test_native_error_field_cannot_change_frozen_score_choice():
    values=rows()
    values[1]['native_error']=999.
    values[0]['native_error']=0.
    assert choose_complete(values)['selected_index']==1
    assert choose_complete(values,guarded=True)['selected_index']==1


def test_nonfinite_scores_refuse_and_no_extra_structure_substitution():
    values=rows()
    values[0]['score']=float('nan')
    with pytest.raises(ValueError,match='finite'):
        choose_complete(values)
    assert list(SYSTEMS)==['1ENH','1CRN','1QYS']
