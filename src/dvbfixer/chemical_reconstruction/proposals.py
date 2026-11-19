"""Private geometry proposal records, separate from WP0 planning success."""
from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Literal

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

from .contracts import (
    Admission,
    ContractError,
    Observation,
    _identity,
    _number,
    _sha,
    _tuple,
    digest_record,
)


class Strategy(StrEnum):
    INTERNAL = "graph-internal-authority-v1"
    TEMPLATE = "class-a-rigid-authority-v1"


class PlacementKind(StrEnum):
    CHAIN = "chain"
    CONE = "two-angle-cone"


@dataclass(frozen=True)
class AuditPolicy:
    bond_tolerance_angstrom: float
    angle_tolerance_degrees: float
    minimum_separation_angstrom: float
    max_anchor_rmsd_angstrom: float

    def __post_init__(self) -> None:
        values = (self.bond_tolerance_angstrom, self.angle_tolerance_degrees,
                  self.minimum_separation_angstrom, self.max_anchor_rmsd_angstrom)
        if (any(not _number(v) or v < 0 for v in values) or min(values[:3]) <= 0
                or self.angle_tolerance_degrees >= 180 or self.minimum_separation_angstrom < 1e-6):
            raise ContractError("explicit finite positive geometry bounds and nonnegative anchor bound required")


@dataclass(frozen=True)
class Placement:
    atom: ExactAtomRef
    anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]
    kind: PlacementKind = PlacementKind.CHAIN

    def __post_init__(self) -> None:
        _identity(self.atom)
        _tuple(self.anchors, ExactAtomRef, "ordered construction anchors", maximum=3)
        if len(self.anchors) != 3 or len(set((self.atom, *self.anchors))) != 4:
            raise ContractError("three distinct ordered anchors and a distinct generated atom required")
        for a in self.anchors:
            _identity(a)
        if not isinstance(self.kind, PlacementKind):
            raise ContractError("typed placement strategy required")


@dataclass(frozen=True)
class GeometryPlan:
    request_digest: str
    strategy: Strategy
    policy: AuditPolicy
    placements: tuple[Placement, ...] = ()

    def __post_init__(self) -> None:
        _sha(self.request_digest)
        if not isinstance(self.strategy, Strategy) or not isinstance(self.policy, AuditPolicy):
            raise ContractError("typed strategy and explicit audit policy required")
        _tuple(self.placements, Placement, "construction plan", maximum=4096)

    @property
    def digest(self) -> str:
        return digest_record(("chemical-geometry-plan-v1", self))


@dataclass(frozen=True)
class Candidate:
    request_digest: str
    plan_digest: str
    strategy: Strategy
    additions: tuple[Observation, ...]

    def __post_init__(self) -> None:
        _sha(self.request_digest)
        _sha(self.plan_digest)
        if not isinstance(self.strategy, Strategy):
            raise ContractError("typed candidate strategy required")
        _tuple(self.additions, Observation, "generated-only coordinates", maximum=4096)


@dataclass(frozen=True)
class CandidateAudit:
    passed: bool
    findings: tuple[str, ...]
    checked_constraints: tuple[int, ...]
    unassessed_context: tuple[ComponentInstanceRef, ...]


@dataclass(frozen=True)
class ProposalResult:
    status: Literal["declared-geometry-passed", "unsupported", "invalid-plan", "construction-failed", "rejected"]
    admission: Admission
    plan_digest: str
    reason: str
    raw: Candidate | None = None
    raw_audit: CandidateAudit | None = None
    rounded: Candidate | None = None
    rounded_audit: CandidateAudit | None = None

    @property
    def coordinates_produced(self) -> bool:
        return self.raw is not None and bool(self.raw.additions)

    @property
    def geometry_approved(self) -> Literal[False]:
        return False

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False
