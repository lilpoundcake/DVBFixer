"""Exact-identity heavy holdout measurements; no fitting, ranking or proposal input."""
from __future__ import annotations

import math
from dataclasses import dataclass

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import (
    ContractError,
    Observation,
    Operation,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from .proposals import Candidate, GeometryPlan


class HoldoutScoringError(ContractError):
    """Incomplete/incompatible reference or candidate; never invent a partial score."""


@dataclass(frozen=True)
class HeavyHoldoutScore:
    atom_count: int
    rmsd_angstrom: float
    maximum_angstrom: float
    request_digest: str
    plan_digest: str
    candidate_digest: str
    reference_digest: str
    scope: str = "exact-name-fixed-frame-heavy-only-no-symmetry-fit-or-acceptance-gate"


def score_heavy_holdout(request: ReconstructionRequest, plan: GeometryPlan, candidate: Candidate,
                        reference: tuple[Observation, ...]) -> HeavyHoldoutScore:
    """Measure complete generated-only heavy coordinates, including rejected geometry.

    Call AFTER proposal completion. The reference is not accepted by any compiler
    or constructor and must never select a phase, candidate, mask or chemical state.
    """
    if (request.operation != Operation.HEAVY or request.ownership.temporary or request.ownership.mobile
            or candidate.request_digest != request.digest or plan.request_digest != request.digest
            or candidate.plan_digest != plan.digest or candidate.strategy != plan.strategy):
        raise HoldoutScoringError("heavy-holdout-request-plan-candidate-binding")
    count, rmsd, maximum = measure_heavy_holdout(request, candidate.additions, reference)
    return HeavyHoldoutScore(count, rmsd, maximum, request.digest, plan.digest, digest_record(candidate),
                             digest_record(tuple(sorted(reference, key=lambda o: o.identity))))


def measure_heavy_holdout(request: ReconstructionRequest, additions: tuple[Observation, ...],
                            reference: tuple[Observation, ...]) -> tuple[int, float, float]:
    """Shared exact heavy coordinates measurement, without constructor relabelling."""
    if request.operation != Operation.HEAVY or request.ownership.temporary or request.ownership.mobile:
        raise HoldoutScoringError("heavy-holdout-request-plan-candidate-binding")
    generated = request.ownership.generated
    if not generated or len(generated) > 4096:
        raise HoldoutScoringError("heavy-holdout-needs-nonempty-bounded-mask")
    if generated & {o.identity for o in request.observations}:
        raise HoldoutScoringError("heavy-holdout-mask-overlaps-observation")
    if not isinstance(reference, tuple) or len(reference) != len(generated) or len(additions) != len(generated):
        raise HoldoutScoringError("heavy-holdout-needs-complete-exact-mask")
    if any(not isinstance(o, Observation) for o in reference):
        raise HoldoutScoringError("heavy-holdout-needs-typed-reference")
    chemical = {ExactAtomRef(c.identity, a.name): a for c in request.components if c.resolved is not None for a in c.resolved.graph.atoms}
    for observations in (reference, additions):
        if {o.identity for o in observations} != generated:
            raise HoldoutScoringError("heavy-holdout-identity-or-duplicate")
        for o in observations:
            atom = chemical.get(o.identity)
            if (atom is None or o.identity.component not in request.targets or atom.element in {"H", "D"}
                    or o.element != atom.element or o.unit != Unit.ANGSTROM or o.occupancy != 1):
                raise HoldoutScoringError("heavy-holdout-heavy-element-unit-or-occupancy")
    truth = {o.identity: o.position for o in reference}
    distances = [math.dist(o.position, truth[o.identity]) for o in additions]
    if not all(math.isfinite(d) for d in distances):
        raise HoldoutScoringError("heavy-holdout-nonfinite-displacement")
    scale = math.sqrt(len(distances))
    return len(distances), math.hypot(*(d/scale for d in distances)), max(distances)
