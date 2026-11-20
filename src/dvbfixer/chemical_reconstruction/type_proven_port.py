"""Explicit v2 port declaration, carrying a source/family-invariant type proof."""
from __future__ import annotations

from dataclasses import dataclass

from .contracts import ReconstructionRequest, digest_record
from .incoming_type import IncomingCarbonylTypeProof
from .planar_port import PlanarPortDeclaration


@dataclass(frozen=True, kw_only=True)
class TypeProvenPlanarPortDeclaration(PlanarPortDeclaration):
    incoming_proof: IncomingCarbonylTypeProof

    def __post_init__(self) -> None:
        super().__post_init__()
        self.incoming_proof.validate()
        if (self.incoming != self.incoming_proof.source.incoming or self.nitrogen != self.incoming_proof.source.recipient
                or self.atom_types[0] != self.incoming_proof.atom_type):
            raise ValueError("type-proven-port-incoming-identity-or-type-binding")

    def validate_incoming_type(self, request: ReconstructionRequest) -> None:
        self.incoming_proof.validate()
        if self.incoming_proof.request.digest != request.digest:
            raise ValueError("type-proven-port-source-request-binding")

    @property
    def digest(self) -> str:
        return digest_record(("type-proven-real-tertiary-peptide-port-v2", self))
