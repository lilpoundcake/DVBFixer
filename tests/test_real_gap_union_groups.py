import pytest
from scripts.merge_real_gap_candidate_groups import regroup


def test_new_bridge_merges_old_groups_before_split_assignment():
    old = [{"entry_id": "old1", "leakage_keys": ["cluster:a"], "conservative_entry_component": 0},
           {"entry_id": "old2", "leakage_keys": ["cluster:b"], "conservative_entry_component": 1}]
    new = [{"entry_id": "new", "leakage_keys": ["cluster:a", "cluster:b"], "conservative_entry_component": 0}]
    rows, groups = regroup([old, new])
    assert groups == [["new", "old1", "old2"]]
    assert {r["conservative_entry_component"] for r in rows} == {0}
    assert all(r["split_assignment"] is None for r in rows)
    assert old[1]["conservative_entry_component"] == 1  # Source evidence stays unchanged.


def test_duplicate_entries_are_not_silently_adopted():
    row = {"entry_id": "same", "leakage_keys": []}
    with pytest.raises(ValueError, match="overlap"):
        regroup([[row], [row]])
