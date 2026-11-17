import pytest
from scripts.apply_real_gap_loop_scope import scoped_rows


def test_huge_explicit_exclusion_does_not_erase_readiness_or_mutate_request():
    original = {"entry_id": "7N8X", "small_assembly_control": False, "longest_gap": 812,
                "request_path": "retained/request.json", "request_sha256": "unchanged"}
    row, = scoped_rows([original], {"7N8X"})
    assert row["current_scope_status"] == "excluded_by_user"
    assert row["request_sha256"] == "unchanged"
    assert not row["current_test_admitted"]
    assert "current_scope_status" not in original


def test_no_implicit_cutoff_or_short_gap_admission_or_control_promotion():
    rows = [{"entry_id": str(i), "small_assembly_control": control, "longest_gap": length}
            for i, (control, length) in enumerate(((False, 3), (False, 92), (False, 535), (True, 8)))]
    scoped = scoped_rows(rows, set())
    assert [r["current_scope_status"] for r in scoped] == ["loop_scope_and_context_review_pending"]*3 + ["small_control_only"]
    assert not any(r["new_sampling_authorized"] or r["current_test_admitted"] for r in scoped)


def test_duplicate_source_readiness_cannot_be_counted_twice():
    row = {"entry_id": "same", "small_assembly_control": False}
    with pytest.raises(ValueError, match="duplicate"):
        scoped_rows([row, row], set())
