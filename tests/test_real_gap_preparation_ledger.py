import pytest
from scripts.summarize_real_gap_preparation import account


def test_completed_refusals_blocks_and_pending_are_not_imputed_passes():
    rows = [{"entry_id": "A", "status": "private_source_projection_prepared"},
            {"entry_id": "B", "status": "source_preparation_refused"},
            {"entry_id": "C", "status": "source_time_blocked"}]
    counts = account(["A", "B", "C", "D"], rows, ["D"])
    assert counts == {"private_source_projection_prepared": 1, "source_preparation_refused": 1, "source_time_blocked": 1}


@pytest.mark.parametrize("pending", [["A", "D"], [], ["X"]])
def test_duplicate_retries_missing_tasks_or_unknown_tasks_are_rejected(pending):
    with pytest.raises(ValueError, match="partition"):
        account(["A", "D"], [{"entry_id": "A", "status": "source_time_blocked"}], pending)


def test_final_ledger_requires_no_pending_tasks_when_every_outcome_is_settled():
    assert account(["A"], [{"entry_id": "A", "status": "source_time_blocked"}], []) == {"source_time_blocked": 1}
