import pytest
from scripts.acquire_real_gap_hard_reserves import select


def test_hard_recruitment_preserves_unbounded_spans_and_excludes_previous_sources():
    rows = [{"entry_id": name, "longest_gap_metadata": length}
            for name, length in (("old", 812), ("short", 49), ("new", 50), ("huge", 1200))]
    selected, strata = select(rows, {"old"})
    assert [r["entry_id"] for r in selected] == ["huge", "new"]
    assert strata == [{"longest_span": [50, 99], "selected": 1}, {"longest_span": [100, None], "selected": 1}]


def test_resource_budget_is_not_a_silent_gap_length_filter():
    with pytest.raises(ValueError, match="budget"):
        select([{"entry_id": str(i), "longest_gap_metadata": 50} for i in range(201)], set())
