"""Native-free complete-pool selection cannot reward incomplete optimization."""
from scripts.audit_chi_stopping_controls import select_complete


def test_all_common_feasible_refined_candidates_rank_by_score_not_native_error():
    rows=[dict(index=0,status='audited',common_audit=dict(accepted=True),final_raw_score=2.,native_error=0.),
          dict(index=1,status='audited',common_audit=dict(accepted=True),final_raw_score=1.,native_error=999.)]
    assert select_complete(rows)==1
    rows[1]['common_audit']['accepted']=False
    assert select_complete(rows)==0


def test_runtime_failure_blocks_incomplete_pool_and_empty_geometry_is_none():
    assert select_complete([dict(index=0,status='execution_or_mobility_blocker')]) is None
    assert select_complete([dict(index=0,status='audited',common_audit=dict(accepted=False),final_raw_score=1.)]) is None
