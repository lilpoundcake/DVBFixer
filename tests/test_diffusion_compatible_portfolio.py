from dataclasses import replace

import pytest

from dvbfixer.model.diffusion.compatible_portfolio import DonorAlternative, select_compatible
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment, RepairIssue


def alternative(name, group, gaps, family=""):
    return DonorAlternative(name, group, gaps, "a"*64, family)


def test_joint_groups_are_atomic_and_every_complete_plan_is_assessed():
    alternatives = (alternative("joint7", "joint", (0, 1), "seed7"),
                    alternative("joint42", "joint", (0, 1), "seed42"),
                    alternative("solo7", "solo", (2,), "seed7"), alternative("solo42", "solo", (2,), "seed42"))
    seen = []
    def assess(keys):
        seen.append(keys)
        return RepairAssessment((), () if keys == ("joint42", "solo7") else ("whole-assembly-failure",))
    result = select_compatible(alternatives, assess_complete=assess)
    assert result.selected.alternatives == ("joint42", "solo7")
    assert len(seen) == 4 and all(len(keys) == 2 for keys in seen)
    assert result.optimality_established and not result.publication_qualified


def test_manual_allowed_sets_locks_and_conflicts_use_same_search():
    alternatives = tuple(alternative(f"{g}{i}", g, (j,)) for j, g in enumerate(("a", "b")) for i in range(2))
    conflict = frozenset((frozenset(("a1", "b0")),))
    result = select_compatible(alternatives, assess_complete=lambda _: RepairAssessment(), conflicts=conflict,
                               locks={"a": "a1"}, allowed={"b": frozenset(("b0", "b1"))})
    assert result.selected.alternatives == ("a1", "b1")
    assert result.pairwise_pruned == 1 and result.constrained_space_size == 2
    with pytest.raises(ValueError, match="contradicts"):
        select_compatible(alternatives, assess_complete=lambda _: RepairAssessment(), locks={"a": "a1"},
                          allowed={"a": frozenset(("a0",))})


def test_pairwise_pass_is_not_whole_assembly_acceptance():
    alternatives = (alternative("a", "one", (0,)), alternative("b", "two", (1,)))
    result = select_compatible(alternatives, assess_complete=lambda _: RepairAssessment((RepairIssue("clash", (), "collective"),)))
    assert not result.selected.assessment.geometry_passed


def test_unknown_conflicts_and_overlapping_joint_ownership_are_refused():
    a, b = alternative("a", "one", (0, 1)), alternative("b", "two", (1, 2))
    with pytest.raises(ValueError, match="overlap"):
        select_compatible((a, b), assess_complete=lambda _: RepairAssessment())
    b = replace(b, gap_indices=(2,))
    with pytest.raises(ValueError, match="conflicts"):
        select_compatible((a, b), assess_complete=lambda _: RepairAssessment(), conflicts=frozenset((frozenset(("a", "unknown")),)))


def test_bounded_large_search_is_deterministic_and_not_claimed_optimal():
    alternatives = tuple(alternative(f"{group}-{seed}", str(group), (group,), str(seed))
                         for group in range(4) for seed in range(5))
    def call(values):
        return select_compatible(values, assess_complete=lambda _: RepairAssessment(), max_evaluations=8, beam_width=2)
    first, second = call(alternatives), call(tuple(reversed(alternatives)))
    assert first.selected == second.selected
    assert len(first.trials) <= 8 and not first.optimality_established
    assert first.constrained_space_size == 625


def test_locally_inadmissible_locked_donor_does_not_silently_unlock():
    a = replace(alternative("a", "one", (0,)), local_admission_failures=("source-motion",))
    result = select_compatible((a,), assess_complete=lambda _: pytest.fail("inadmissible donor was evaluated"), locks={"one": "a"})
    assert result.selected is None and result.constrained_space_size == 0
