from scripts.audit_retained_construct_differences import gap_evidence, map_difference


def test_difference_seq_num_is_label_not_author_number_and_chain_case_is_preserved():
    candidate = {"source_author_chains": ["a"], "source_label_chain": "X", "entity_id": "1"}
    row = {"seq_num": "5", "pdbx_auth_seq_num": "82", "pdbx_pdb_strand_id": "a", "pdbx_pdb_ins_code": "A", "mon_id": "HIS", "details": "expression tag"}
    scheme = [{"asym_id": "X", "entity_id": "1", "pdb_strand_id": "a", "pdb_seq_num": "82",
               "pdb_ins_code": "A", "mon_id": "HIS", "seq_id": "5"}]
    assert map_difference(row, candidate, scheme)["label_index_zero_based"] == 4
    assert map_difference({**row, "pdbx_pdb_strand_id": "A"}, candidate, scheme)["label_index_zero_based"] is None
    assert map_difference({**row, "seq_num": "82"}, candidate, scheme)["label_index_zero_based"] is None
    assert map_difference(row, candidate, [*scheme, *scheme])["label_index_zero_based"] is None


def test_complete_artifact_evidence_requires_all_source_candidate_agreement():
    mapped = [{"label_index_zero_based": i, "is_explicit_artifact_description": True} for i in (0, 1)]
    gap = {"start": 0, "stop": 2, "terminal": True}
    assert gap_evidence(gap, [mapped, mapped])["coverage"] == "complete-explicit-artifact-annotation"
    assert gap_evidence(gap, [mapped, []])["coverage"] == "candidate-dependent-artifact-annotation"
    assert gap_evidence(gap, [mapped, mapped[:1]])["coverage"] == "partial-explicit-artifact-annotation"
    assert not gap_evidence(gap, [mapped])["request_mask_or_admission_changed"]


def test_missing_mapping_never_becomes_no_tag_claim():
    evidence = gap_evidence({"start": 0, "stop": 1}, [[{"label_index_zero_based": None,
                                                        "is_explicit_artifact_description": True}]])
    assert evidence["coverage"] == "no-mapped-artifact-annotation"
    assert evidence["no_annotation_is_not_proof_of_native_sequence"]


def test_explicit_database_deletion_does_not_invent_a_deposited_coordinate_target():
    candidate = {"source_author_chains": ["a"], "source_label_chain": "X", "entity_id": "1"}
    row = {"details": "deletion", "seq_num": None, "mon_id": None,
           "db_mon_id": "ALA", "pdbx_seq_db_seq_num": "83", "pdbx_pdb_strand_id": "a"}
    result = map_difference(row, candidate, [])
    assert result["label_index_zero_based"] is None
    assert result["mapping_status"] == "not-applicable-database-deletion"
    incomplete = map_difference({**row, "details": "conflict"}, candidate, [])
    assert incomplete["mapping_status"] == "unresolved"
