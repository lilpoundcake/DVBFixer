from types import SimpleNamespace

from scripts.audit_real_gap_cohort_capacity import capacity, task_inventory


def test_terminal_classification_uses_target_axis_not_author_numbering():
    request = SimpleNamespace(target_sequences=[SimpleNamespace(chain="a", sequence="A"*500)],
        sequence_placements=[SimpleNamespace(observed_residues=tuple(range(350)))],
        generated_atoms=tuple(range(100)), fixed_atoms=tuple(range(1000)), gaps=[
            SimpleNamespace(chain="a", target_interval=SimpleNamespace(start=0, stop=100)),
            SimpleNamespace(chain="a", target_interval=SimpleNamespace(start=200, stop=250))])
    row = task_inventory(request)
    assert [g["kind"] for g in row["gaps"]] == ["N-terminal", "internal"]
    assert row["terminal_construct_review_required"]
    assert row["internal_gaps_present"]
    assert not row["small_assembly_control"]


def test_small_controls_do_not_fill_quota_and_shared_groups_are_not_dropped():
    rows = [{"longest_gap": 60, "small_assembly_control": small, "union_component": 1,
             "group_membership_complete": True, "terminal_construct_review_required": False}
            for small in (True, False, False)]
    row = capacity(rows)[3]
    assert row["ready_projections"] == 3
    assert row["not_small_controls"] == 2
    assert row["union_components_represented"] == 1
    assert row["minimum_additional_projections_before_other_admission"] == 38


def test_longest_gap_bins_partition_ready_rows_at_every_boundary():
    rows = [{"longest_gap": length, "small_assembly_control": False, "union_component": i,
             "group_membership_complete": True, "terminal_construct_review_required": False}
            for i, length in enumerate((1, 8, 9, 24, 25, 49, 50, 99, 100, 812))]
    quotas = capacity(rows)
    assert [q["ready_projections"] for q in quotas] == [2, 2, 2, 2, 2]
    assert sum(q["ready_projections"] for q in quotas) == len(rows)
    assert sum(q["required_entries"] for q in quotas) == 200
