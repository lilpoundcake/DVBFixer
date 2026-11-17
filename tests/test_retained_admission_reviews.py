from scripts.build_retained_admission_reviews import review_questions


def questions(**overrides):
    values = dict(terminal_gaps=0, other_context=False, water_only=False, artifact_overlap=False,
                  complete_artifact_gap=False, unresolved_annotations=0, prepared=False,
                  group_complete=True, exposed=False)
    values.update(overrides)
    return review_questions("loop_scope_and_context_review_pending", **values)


def test_verified_facts_do_not_auto_admit_or_create_length_cutoff():
    result = questions()
    assert result == ["adjudicate-loop-versus-domain-in-observed-scaffold", "backend-and-whole-context-admission-pending"]


def test_artifact_evidence_requires_task_decision_not_silent_trimming():
    result = questions(artifact_overlap=True, complete_artifact_gap=True, terminal_gaps=1)
    assert "whole-gap-artifact-annotation-no-silent-mask-trimming" in result
    assert "adjudicate-terminal-construct-not-just-gap-length" in result


def test_controls_exclusions_and_ambiguity_are_preserved():
    result = questions(unresolved_annotations=1, exposed=True, group_complete=False)
    assert "resolve-unmapped-difference-annotations" in result
    assert "complete-family-grouping-before-split" in result
    assert "known-development-component-not-untouched-final-evaluation" in result
    assert review_questions("excluded_by_user", terminal_gaps=0, other_context=False, water_only=False,
                            artifact_overlap=False, complete_artifact_gap=False, unresolved_annotations=0,
                            prepared=False, group_complete=True, exposed=False) == ["retain-explicit-exclusion-no-active-task"]


def test_database_deletion_is_not_an_unmapped_coordinate_target():
    result = questions(database_deletions=2)
    assert "retain-database-deletions-as-outside-deposited-coordinate-task" in result
    assert "resolve-unmapped-difference-annotations" not in result
