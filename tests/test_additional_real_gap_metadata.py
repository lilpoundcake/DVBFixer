from scripts.collect_additional_long_gap_candidates import STRATA, rank_new


def test_additional_metadata_rank_excludes_prior_entries_and_is_stable():
    assert rank_new(["NEW", "OLD", "OTHER", "NEW"], {"OLD"}, 0) == rank_new(["OTHER", "NEW"], set(), 0)
    assert len(rank_new(["NEW", "OLD", "OTHER", "NEW"], {"OLD"}, 0)) == 2


def test_declared_budget_is_metadata_only_and_high_stratum_has_no_ceiling():
    assert sum(n for _, _, n in STRATA) == 500
    assert STRATA[-1][1] is None
