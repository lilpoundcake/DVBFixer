from scripts.retained_return_review_library import candidate_id, filter_rows


def test_library_filter_is_case_sensitive_and_does_not_select_or_approve():
    rows = [{"entry_id": entry, "native_backend": "MODELLER", "native_seed": 42,
             "stage": "generated-only", "manual_approval": None} for entry in ("A", "a")]
    result = filter_rows(rows, entry="a", backend="MODELLER", seed=42, stage="generated-only")
    assert result == [rows[1]]
    assert result[0]["manual_approval"] is None
    assert not filter_rows(rows, entry="a", backend="Protenix")


def test_library_lineage_id_does_not_double_count_a_second_export_directory():
    row = {"entry_id": "1NON", "native_backend": "MODELLER", "native_protocol": "baseline",
           "native_seed": 42, "optimizer_seed": 11, "stage": "generated-only", "sha256": "a"*64,
           "scientific_task_sha256": "b"*64, "path": "/first/export/file.pdb"}
    assert candidate_id(row) == candidate_id({**row, "path": "/second/export/file.pdb"})
