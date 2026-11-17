from scripts import expand_real_gap_source_reserves as expansion


def row(name, span, component):
    return {"entry_id": name, "longest_gap_metadata": span,
            "natural_gap_count_metadata": 1, "conservative_entry_component": component}


def test_recruitment_never_repeats_sources_and_reports_absent_hard_strata():
    rows = [row("OLD", 4, 1), row("NEW", 4, 1), row("OTHER", 4, 2), row("LONG", 812, 3)]
    selected, strata = expansion.select(rows, {"OLD", "LONG"})
    assert [r["entry_id"] for r in selected] == ["OTHER", "NEW"]
    assert strata[-1]["available"] == 0
    assert strata[-1]["selected"] == 0
    assert strata[-1]["longest_span"] == [100, None]
    again, _ = expansion.select(list(reversed(rows)), {"OLD", "LONG"})
    assert selected == again
