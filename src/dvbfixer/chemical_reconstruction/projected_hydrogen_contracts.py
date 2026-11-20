"""Separate projected-parent/H evidence kinds; no torsion-result relabelling."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .contracts import PreparedContextRequest, ReconstructionRequest, _sha, digest_record
from .declared_peptide import DeclaredPeptideAuthority
from .hydrogen_completion import CompleteHydrogenView
from .planar_port import PlanarPortAuthority, PlanarPortDeclaration
from .projection_contracts import ProjectionResult
from .proposals import ProposalResult
from .staged_contracts import StagedAudit, StagedEnvelope


@dataclass(frozen=True)
class ProjectedHeavyParent:
    original: ReconstructionRequest
    request: ReconstructionRequest
    projection: ProjectionResult
    peptide_authority: DeclaredPeptideAuthority
    port_declaration: PlanarPortDeclaration
    port_parameters: PlanarPortAuthority

    def __post_init__(self) -> None:
        if (type(self.original) not in {ReconstructionRequest, PreparedContextRequest}
                or type(self.request) not in {ReconstructionRequest, PreparedContextRequest}
                or not isinstance(self.projection, ProjectionResult)
                or not isinstance(self.peptide_authority, DeclaredPeptideAuthority)
                or not isinstance(self.port_declaration, PlanarPortDeclaration)
                or not isinstance(self.port_parameters, PlanarPortAuthority)):
            raise ValueError("typed-projected-heavy-parent-required")

    @property
    def digest(self) -> str:
        return digest_record(("real-port-projected-heavy-parent-v1", self))


@dataclass(frozen=True)
class ProjectedHydrogenResult:
    original_request_digest: str
    port_request_digest: str
    parent_digest: str
    status: Literal["parent-refused", "H-refused", "H-rejected", "projected-H-geometry-passed"]
    reason: str
    view: CompleteHydrogenView | None = None
    hydrogen_result: ProposalResult | None = None
    raw: StagedEnvelope | None = None
    rounded: StagedEnvelope | None = None
    audit: StagedAudit | None = None

    def __post_init__(self) -> None:
        for pin in (self.original_request_digest, self.port_request_digest, self.parent_digest):
            _sha(pin)

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def native_quality_qualified(self) -> Literal[False]:
        return False
