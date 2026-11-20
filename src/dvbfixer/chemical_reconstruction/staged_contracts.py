"""Private role-separated heavy/H evidence, never a mixed generic operation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .contracts import Observation, _sha, _tuple, digest_record
from .hydrogen_completion import CompleteHydrogenView, ConditionalHydrogenCipAudit
from .proposals import AuditPolicy, CandidateAudit, ProposalResult

# Frozen WP2L profile; not caller-adjustable acceptance knobs.
FIXED_POLICY = AuditPolicy(.10, 12., .7, .02)
HYDROGEN_POLICY = AuditPolicy(.02, 5., .7, .02)


@dataclass(frozen=True)
class StagedEnvelope:
    original_request_digest: str
    parent_digest: str
    conditioning_digest: str
    hydrogen_candidate_digest: str
    phase: Literal["raw-H-after-rounded-heavy", "rounded-H-after-rounded-heavy"]
    generated_heavy: tuple[Observation, ...]
    generated_hydrogens: tuple[Observation, ...]

    def __post_init__(self) -> None:
        for value in (self.original_request_digest, self.parent_digest, self.conditioning_digest,
                      self.hydrogen_candidate_digest):
            _sha(value)
        _tuple(self.generated_heavy, Observation, "staged generated heavy", maximum=4096)
        _tuple(self.generated_hydrogens, Observation, "staged generated H", maximum=4096)
        if self.phase not in {"raw-H-after-rounded-heavy", "rounded-H-after-rounded-heavy"}:
            raise ValueError("staged-envelope-phase")


@dataclass(frozen=True)
class PairCoverage:
    heavy_source: int
    heavy_heavy: int
    hydrogen_source: int
    hydrogen_heavy: int
    hydrogen_hydrogen: int
    eligible_pairs_digest: str
    scope: str = "authoritative 1-2/1-3 excluded; parent covers heavy pairs, H audit covers every H pair"


@dataclass(frozen=True)
class StagedAudit:
    passed: bool
    findings: tuple[str, ...]
    raw_geometry: CandidateAudit | None = None
    rounded_geometry: CandidateAudit | None = None
    raw_cip: ConditionalHydrogenCipAudit | None = None
    rounded_cip: ConditionalHydrogenCipAudit | None = None
    pair_coverage: PairCoverage | None = None


@dataclass(frozen=True)
class StagedResult:
    original_request_digest: str
    parent_digest: str
    status: Literal["parent-refused", "H-refused", "H-rejected", "staged-geometry-passed"]
    reason: str
    view: CompleteHydrogenView | None = None
    hydrogen_result: ProposalResult | None = None
    raw: StagedEnvelope | None = None
    rounded: StagedEnvelope | None = None
    audit: StagedAudit | None = None

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def native_quality_qualified(self) -> Literal[False]:
        return False


def parent_digest(value: object) -> str:
    return digest_record(("staged-frozen-heavy-parent-v1", value))
