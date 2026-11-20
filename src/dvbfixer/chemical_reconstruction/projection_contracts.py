"""Private bounded generated-heavy projection contracts, not a torsion parent."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import Constraint, Origin, ReconstructionRequest, _tuple, digest_record
from .proposals import (
    AuditPolicy,
    Candidate,
    CandidateAudit,
    GeometryPlan,
    ProposalResult,
    Strategy,
)

POLICY = AuditPolicy(.10, 12., .7, .02)


def append_projection_constraints(original: ReconstructionRequest, extra: tuple[Constraint, ...]) -> ReconstructionRequest:
    _tuple(extra, Constraint, "explicit projection constraints", maximum=512)
    if not extra:
        raise ValueError("projection-needs-new-declared-constraints")
    settings = original.origin.settings if original.origin else ()
    if any(key.startswith("metric_projection_") for key, _ in settings):
        raise ValueError("nested-projection-request-refused")
    settings = (*settings, ("metric_projection_parent", original.digest), ("metric_projection_constraints", digest_record(extra)))
    origin = (replace(original.origin, settings=settings) if original.origin else
              Origin("generated-heavy-projection-request-v1", original.digest, settings))
    return replace(original, constraints=original.constraints+extra, origin=origin)


@dataclass(frozen=True)
class ProjectionMetric:
    atoms: tuple[ExactAtomRef, ...]
    value: float
    scale: float


@dataclass(frozen=True, kw_only=True)
class ProjectionPlan(GeometryPlan):
    original: ReconstructionRequest
    baseline_plan: GeometryPlan
    baseline: ProposalResult
    anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]
    terms: tuple[ProjectionMetric, ...]
    max_nfev: int = 500
    local_component_bound_angstrom: float = 2.
    per_atom_movement_cap_angstrom: float = 2.

    def __post_init__(self) -> None:
        super().__post_init__()
        _tuple(self.anchors, ExactAtomRef, "observed projection frame", maximum=3)
        _tuple(self.terms, ProjectionMetric, "projection terms", maximum=512)
        if (self.strategy != Strategy.PROJECTION or self.policy != POLICY or self.placements
                or len(self.anchors) != 3 or len(set(self.anchors)) != 3 or not self.terms
                or type(self.max_nfev) is not int or self.max_nfev != 500
                or type(self.local_component_bound_angstrom) is not float or self.local_component_bound_angstrom != 2.
                or type(self.per_atom_movement_cap_angstrom) is not float or self.per_atom_movement_cap_angstrom != 2.):
            raise ValueError("projection-frozen-plan-shape-or-policy")

    @property
    def digest(self) -> str:
        return digest_record(("generated-only-heavy-projection-plan-v1", self))


@dataclass(frozen=True)
class ProjectionAudit:
    passed: bool
    findings: tuple[str, ...]
    raw_geometry: CandidateAudit | None = None
    rounded_geometry: CandidateAudit | None = None
    maximum_raw_movement_angstrom: float | None = None
    maximum_rounded_movement_angstrom: float | None = None


@dataclass(frozen=True)
class ProjectionResult:
    status: Literal["refused", "solver-failed", "rejected", "projected-geometry-passed"]
    reason: str
    plan: ProjectionPlan | None = None
    raw: Candidate | None = None
    rounded: Candidate | None = None
    audit: ProjectionAudit | None = None
    solver_success: bool | None = None
    solver_status: int | None = None
    solver_message: str | None = None
    nfev: int | None = None
    initial_squared_residual: float | None = None
    final_squared_residual: float | None = None

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def heavy_H_or_native_quality_qualified(self) -> Literal[False]:
        return False
