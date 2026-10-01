"""Backend-neutral contracts for fail-closed nonprotein reconstruction."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

from .parameterization import ParameterizationRoute
from .structure_identity import ComponentInstanceRef, ExactAtomRef

Coordinate = tuple[float, float, float]


class ReconstructionStatus(StrEnum):
    SUCCEEDED = "succeeded"
    UNSUPPORTED = "unsupported"
    AMBIGUOUS = "ambiguous"
    INVALID_INPUT = "invalid-input"
    FAILED = "failed"


class ComponentClass(StrEnum):
    AUTO = "auto"
    A = "A"
    B = "B"
    C = "C"
    D = "D"


class AuthorityMode(StrEnum):
    LOCAL_PINNED_CCD = "local-pinned-ccd"
    USER_MAPPED = "user-mapped"
    ONLINE_ENRICHMENT = "online-enrichment"


class CoordinateSource(StrEnum):
    OBSERVED = "observed"
    SOURCE_TEMPLATE = "source-template"


@dataclass(frozen=True)
class ChemicalAtom:
    name: str
    element: str
    formal_charge: int
    ideal_position: Coordinate
    aromatic: bool = False
    stereo: str | None = None


@dataclass(frozen=True)
class ChemicalBond:
    atom1: str
    atom2: str
    order: float
    ideal_length: float
    aromatic: bool = False


@dataclass(frozen=True)
class ChemicalGraph:
    component_id: str
    state_id: str
    atoms: tuple[ChemicalAtom, ...]
    bonds: tuple[ChemicalBond, ...]
    net_charge: int
    ph_min: float | None = None
    ph_max: float | None = None
    priority: int = 0
    redox_state: str | None = None
    coordination_model: str | None = None
    component_count: int = 1


@dataclass(frozen=True)
class LocalPinnedCcdAuthority:
    path: Path
    sha256: str
    component_id: str
    mode: AuthorityMode = field(default=AuthorityMode.LOCAL_PINNED_CCD, init=False)


@dataclass(frozen=True)
class UserMappedAuthority:
    graph: ChemicalGraph
    source_label: str
    source_digest: str
    mode: AuthorityMode = field(default=AuthorityMode.USER_MAPPED, init=False)


@dataclass(frozen=True)
class OnlineEnrichmentAuthority:
    provider: str
    mode: AuthorityMode = field(default=AuthorityMode.ONLINE_ENRICHMENT, init=False)


ChemicalAuthority = LocalPinnedCcdAuthority | UserMappedAuthority | OnlineEnrichmentAuthority


@dataclass(frozen=True)
class ObservedAtom:
    identity: ExactAtomRef
    element: str
    position: Coordinate


@dataclass(frozen=True)
class ExternalLink:
    component_atom: ExactAtomRef
    other_atom: ExactAtomRef
    order: float | None = None


@dataclass(frozen=True)
class GeometryPolicy:
    max_anchor_rmsd: float = 0.35
    max_bond_deviation: float = 0.25
    min_nonbonded_distance: float = 0.70
    max_observed_movement: float = 0.0


@dataclass(frozen=True)
class ReconstructionRequest:
    component: ComponentInstanceRef
    observed_atoms: tuple[ObservedAtom, ...]
    authority: ChemicalAuthority
    ph: float
    requested_class: ComponentClass = ComponentClass.AUTO
    external_links: tuple[ExternalLink, ...] = ()
    microstate_override: str | None = None
    geometry_policy: GeometryPolicy = field(default_factory=GeometryPolicy)


@dataclass(frozen=True)
class AtomCoordinate:
    identity: ExactAtomRef
    element: str
    position: Coordinate
    source: CoordinateSource
    authority_atom: str


@dataclass(frozen=True)
class AtomMovement:
    identity: ExactAtomRef
    distance: float


@dataclass(frozen=True)
class ValidationFinding:
    code: str
    message: str
    atom_names: tuple[str, ...] = ()


@dataclass(frozen=True)
class MicrostateProvenance:
    policy: str
    ph: float
    selected_state: str | None
    ranked_states: tuple[str, ...]
    candidate_scores: tuple[tuple[str, int | None], ...]
    override: bool
    reason: str


@dataclass(frozen=True)
class AuthorityProvenance:
    mode: AuthorityMode
    source: str
    source_digest: str | None
    source_version: str | None


@dataclass(frozen=True)
class ReconstructionResult:
    status: ReconstructionStatus
    component_class: ComponentClass
    component: ComponentInstanceRef
    resolved_graph: ChemicalGraph | None = None
    external_links: tuple[ExternalLink, ...] = ()
    coordinates: tuple[AtomCoordinate, ...] = ()
    movements: tuple[AtomMovement, ...] = ()
    findings: tuple[ValidationFinding, ...] = ()
    authority: AuthorityProvenance | None = None
    microstate: MicrostateProvenance | None = None
    input_digest: str | None = None
    output_digest: str | None = None
    geometry_approved: bool = False
    md_ready: bool = False


@dataclass(frozen=True)
class ParameterizationRequest:
    geometry: ReconstructionResult
    geometry_approved: bool
    force_field_family: str
    policy: str
    exact_native_match: bool = False
    validated_user_template: bool = False
    allow_isolated_ligand_candidate: bool = False
    explicit_strip: bool = False


@dataclass(frozen=True)
class ParameterizationDecision:
    status: ReconstructionStatus
    route: ParameterizationRoute | None
    reason: str
    geometry_digest: str | None
    force_field_family: str
    policy: str
    md_ready: bool = False
