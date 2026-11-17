import pytest
from scripts.acquire_real_gap_medium_reserves import select


def test_medium_recruitment_excludes_prior_and_does_not_consume_long_cases():
    rows = [{"entry_id": name, "longest_gap_metadata": length}
            for name, length in (("old", 26), ("short", 24), ("new", 25), ("long", 50), ("edge", 49))]
    selected, strata = select(rows, {"old"})
    assert [r["entry_id"] for r in selected] == ["edge", "new"]
    assert strata == [{"longest_span": [25, 49], "selected": 2}]


def test_medium_resource_budget_fails_explicitly():
    with pytest.raises(ValueError, match="budget"):
        select([{"entry_id": str(i), "longest_gap_metadata": 30} for i in range(201)], set())
