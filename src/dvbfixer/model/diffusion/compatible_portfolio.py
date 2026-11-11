"""Private atomic-group selection plans; scientific assembly stays in the caller.

Manual locks and allowed alternatives constrain the same deterministic search.
Pairwise compatibility is only pruning: every scored complete plan must pass
the supplied whole-assembly assessor. This is not a schema-4 candidate, a sampler
trace, coordinate materializer or publication authorization.
"""
from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass

from dvbfixer.model.diffusion.refinement_stages import RepairAssessment


@dataclass(frozen=True)
class DonorAlternative:
    alternative_id: str
    group_id: str
    gap_indices: tuple[int, ...]
    lineage_sha256: str
    control_family: str
    local_defects: int = 0
    local_admission_failures: tuple[str, ...] = ()

    def __post_init__(self):
        if (not self.alternative_id or not self.group_id or not self.gap_indices or
                len(set(self.gap_indices)) != len(self.gap_indices) or any(type(i) is not int or i < 0 for i in self.gap_indices)):
            raise ValueError("alternative requires a unique atomic group of gap indices")
        if len(self.lineage_sha256) != 64 or any(c not in "0123456789abcdef" for c in self.lineage_sha256):
            raise ValueError("alternative requires content-addressed private lineage")
        if type(self.local_defects) is not int or self.local_defects < 0:
            raise ValueError("invalid local defect count")


@dataclass(frozen=True)
class CompatibleTrial:
    alternatives: tuple[str, ...]
    assessment: RepairAssessment

    @property
    def score(self):
        return len(self.assessment.admission_failures), sum(issue.actionable for issue in self.assessment.issues)


@dataclass(frozen=True)
class CompatibleSelection:
    groups: tuple[str, ...]
    selected: CompatibleTrial | None
    trials: tuple[CompatibleTrial, ...]
    constrained_space_size: int
    pairwise_pruned: int
    optimality_established: bool
    budget_exhausted: bool
    method: str = "experimental-atomic-group-compatible-selection-v1"
    publication_qualified: bool = False


def select_compatible(alternatives: tuple[DonorAlternative, ...], *,
                      assess_complete: Callable[[tuple[str, ...]], RepairAssessment],
                      conflicts: frozenset[frozenset[str]] = frozenset(),
                      allowed: Mapping[str, frozenset[str]] | None = None,
                      locks: Mapping[str, str] | None = None,
                      max_evaluations: int = 64, beam_width: int = 4) -> CompatibleSelection:
    """Enumerate small spaces; otherwise deterministic local-defect beam.

    Ranking is lexicographic: whole-assembly admission failures, actionable
    geometry findings, then alternative IDs. No confidence/reference reward or
    sum of incomparable scientific scores. Multi-gap groups are indivisible.
    Homogeneous families are whole-assembly controls under the SAME manual mask.
    """
    if (type(max_evaluations) is not int or type(beam_width) is not int or
            not 1 <= max_evaluations <= 4096 or not 1 <= beam_width <= 32):
        raise ValueError("invalid compatibility search budget")
    by_id = {a.alternative_id: a for a in alternatives}
    if not alternatives or len(by_id) != len(alternatives):
        raise ValueError("alternative IDs must be unique and nonempty")
    groups = tuple(sorted({a.group_id for a in alternatives}))
    allowed, locks = dict(allowed or {}), dict(locks or {})
    if (allowed.keys() | locks.keys()) - set(groups):
        raise ValueError("manual selection refers to an unknown group")
    owned = set()
    options = {}
    for group in groups:
        candidates = [a for a in alternatives if a.group_id == group]
        group_gaps = set(candidates[0].gap_indices)
        if any(set(a.gap_indices) != group_gaps for a in candidates) or owned & group_gaps:
            raise ValueError("atomic groups overlap or alternatives change group membership")
        owned.update(group_gaps)
        permitted = allowed.get(group, frozenset(a.alternative_id for a in candidates))
        if not permitted <= {a.alternative_id for a in candidates}:
            raise ValueError("allowed alternative belongs to an unknown/different group")
        if group in locks:
            if locks[group] not in permitted:
                raise ValueError("manual lock contradicts allowed alternatives")
            permitted = frozenset((locks[group],))
        options[group] = tuple(sorted((a.alternative_id for a in candidates if a.alternative_id in permitted
                                      and not a.local_admission_failures), key=lambda k: (by_id[k].local_defects, k)))
    if any(len(pair) != 2 or not pair <= by_id.keys() or len({by_id[k].group_id for k in pair}) != 2 for pair in conflicts):
        raise ValueError("conflicts must reference two alternatives in different groups")
    space = math.prod(len(options[g]) for g in groups)
    pruned = 0

    def compatible(indices):
        return not any(frozenset(pair) in conflicts for pair in itertools.combinations(indices, 2))

    trials = {}

    def evaluate(indices):
        if indices not in trials and len(trials) < max_evaluations:
            assessment = assess_complete(indices)
            if not isinstance(assessment, RepairAssessment):
                raise ValueError("whole-assembly assessor returned an invalid assessment")
            trials[indices] = CompatibleTrial(indices, assessment)

    # Atomic homogeneous controls, not fabricated mixtures of a joint group.
    families = sorted({a.control_family for a in alternatives if a.control_family})
    for family in families:
        selected = [tuple(k for k in options[g] if by_id[k].control_family == family) for g in groups]
        if all(len(keys) == 1 for keys in selected):
            indices = tuple(keys[0] for keys in selected)
            if compatible(indices):
                evaluate(indices)
    exhaustive = space <= max_evaluations
    if exhaustive:
        for indices in itertools.product(*(options[g] for g in groups)):
            if compatible(indices):
                evaluate(indices)
            else:
                pruned += 1
    else:
        beam = [()]
        for group in groups:
            expanded = []
            for partial in beam:
                for key in options[group]:
                    proposal = (*partial, key)
                    if compatible(proposal):
                        expanded.append(proposal)
                    else:
                        pruned += 1
            beam = sorted(expanded, key=lambda keys: (sum(by_id[k].local_defects for k in keys), keys))[:beam_width]
        for indices in beam:
            evaluate(indices)
        # Improve complete compatible recipes using whole-assembly findings,
        # rather than selecting solely on independent local quality. The
        # partial beam above supplies valid starting points under pair pruning.
        beam = sorted(trials, key=lambda keys: (*trials[keys].score, keys))[:beam_width]
        for group_index, group in enumerate(groups):
            expanded = set(beam)
            for current in beam:
                for key in options[group]:
                    proposal = current[:group_index] + (key,) + current[group_index+1:]
                    if compatible(proposal):
                        evaluate(proposal)
                        if proposal in trials:
                            expanded.add(proposal)
                    else:
                        pruned += 1
            beam = sorted(expanded, key=lambda keys: (*trials[keys].score, keys))[:beam_width]
            if len(trials) >= max_evaluations:
                break
    ranked = tuple(sorted(trials.values(), key=lambda t: (*t.score, t.alternatives)))
    return CompatibleSelection(groups, ranked[0] if ranked else None, tuple(trials.values()), space, pruned,
                               exhaustive, not exhaustive and len(trials) >= max_evaluations)
