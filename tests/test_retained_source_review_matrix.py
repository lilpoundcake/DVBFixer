from scripts.build_retained_source_review_matrix import review_topics


def test_review_topics_do_not_admit_a_short_internal_case_or_invent_a_gap_limit():
    row = {"current_scope_status": "loop_scope_and_context_review_pending", "small_assembly_control": False,
           "terminal_construct_review_required": False, "prepared_source_derivative": False,
           "group_membership_complete": True, "longest_gap": 1}
    topics = review_topics(row, {"HOH": 10}, False)
    assert "loop-versus-domain-and-construct-adjudication" in topics
    assert "omitted-context-other-than-HOH-WAT" not in topics
    assert topics == review_topics({**row, "longest_gap": 1000}, {"HOH": 10}, False)


def test_scope_exclusion_and_context_exposure_remain_explicit():
    excluded = {"current_scope_status": "excluded_by_user"}
    assert review_topics(excluded, {}, False) == ["explicitly-excluded-retain-research-only"]
    row = {"current_scope_status": "loop_scope_and_context_review_pending", "small_assembly_control": False,
           "terminal_construct_review_required": True, "prepared_source_derivative": True,
           "group_membership_complete": False}
    topics = review_topics(row, {"ATP": 31}, True)
    assert {"omitted-context-other-than-HOH-WAT", "incomplete-family-membership",
            "known-development-component-no-final-assignment", "terminal-construct-identity"} <= set(topics)
