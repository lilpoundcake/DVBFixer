"""Distinct generated-conditioning transport for an already committed projection."""
from __future__ import annotations

from dataclasses import dataclass

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import (
    ContractError,
    Operation,
    Origin,
    Ownership,
    PreparedContextRequest,
    ReconstructionRequest,
    digest_record,
)
from .heavy_conditioning import conditioning_fractional_atoms, validate_heavy_conditioning
from .projection_contracts import ProjectionResult, append_projection_constraints
from .proposals import Candidate, Strategy


def projected_snapshot_digest(original: ReconstructionRequest, checkpoint: Candidate, projection: ProjectionResult) -> str:
    return digest_record(("rounded-projected-heavy-conditioning-snapshot-v1", original.digest, checkpoint, projection))


@dataclass(frozen=True, kw_only=True)
class ProjectedConditionedRequest(ReconstructionRequest):
    original: ReconstructionRequest
    checkpoint: Candidate
    projection: ProjectionResult

    def __post_init__(self) -> None:
        super().__post_init__()
        p = self.projection
        if type(self.original) not in {ReconstructionRequest, PreparedContextRequest}:
            raise ContractError("conditioning-needs-original-source-request")
        if (not isinstance(p, ProjectionResult) or p.plan is None or p.status != "projected-geometry-passed"
                or p.audit is None or not p.audit.passed or p.rounded is None
                or p.plan.request_digest != self.original.digest or digest_record(p.rounded) != digest_record(self.checkpoint)):
            raise ContractError("projected-conditioning-needs-bound-passing-rounded-parent")
        source = p.plan.original
        expected = append_projection_constraints(source, self.original.constraints[len(source.constraints):])
        if expected.digest != self.original.digest:
            raise ContractError("projected-conditioning-original-port-request-binding")
        validate_heavy_conditioning(self, self.original, self.checkpoint, strategy=Strategy.PROJECTION,
            adapter="generated-projected-heavy-conditioning-v1",
            snapshot_digest=lambda original, checkpoint: projected_snapshot_digest(original, checkpoint, p))

    @property
    def allowed_fractional_atoms(self) -> frozenset[ExactAtomRef]:
        return conditioning_fractional_atoms(self.original)

    @property
    def digest(self) -> str:
        return digest_record(("chemical-projected-conditioned-request-v1", self))


def condition_on_projected(original: ReconstructionRequest, projection: ProjectionResult) -> ProjectedConditionedRequest:
    """Pure transport only; scientific parent/declaration audit remains mandatory."""
    if projection.rounded is None:
        raise ContractError("projected-conditioning-needs-rounded-checkpoint")
    snapshot = projected_snapshot_digest(original, projection.rounded, projection)
    return ProjectedConditionedRequest(
        snapshot, Operation.HEAVY, original.components, original.observations+projection.rounded.additions,
        original.targets, Ownership(frozenset()), original.connections, original.constraints, original.ph,
        Origin("generated-projected-heavy-conditioning-v1", snapshot), original=original,
        checkpoint=projection.rounded, projection=projection)
