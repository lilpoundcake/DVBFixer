"""Versioned immutable planning records reusing existing domain chemistry types."""
from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass, fields, is_dataclass
from enum import StrEnum
from typing import Literal, TypeGuard

from dvbfixer.domain.nonprotein_reconstruction import (
    AuthorityMode,
    AuthorityProvenance,
    ChemicalAtom,
    ChemicalBond,
    ChemicalGraph,
    ExternalLink,
    ObservedAtom,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

REVISION = "chemical-reconstruction-contract-v1"
MAX_ITEMS = 50_000


class ContractError(ValueError):
    """Malformed immutable DTO, not a request for a scientific fallback."""


class Unit(StrEnum):
    ANGSTROM = "angstrom"
    NANOMETER = "nanometer"
    DEGREE = "degree"
    DIMENSIONLESS = "dimensionless"


class Operation(StrEnum):
    HEAVY = "complete-heavy"
    HYDROGENS = "place-hydrogens"
    REPACK = "repack-observed"


class ChemistryKind(StrEnum):
    PROTEIN = "protein"
    NUCLEOTIDE = "nucleotide"
    GLYCAN = "glycan"
    ISOLATED = "isolated-component"
    LINKED = "linked-component"
    COORDINATION = "coordination-complex"
    UNKNOWN = "unknown"


class HydrogenInventory(StrEnum):
    UNSPECIFIED = "unspecified"
    EXPLICIT = "explicit-complete"


class LinkKind(StrEnum):
    COVALENT = "covalent"
    COORDINATION = "coordination"
    UNKNOWN = "unknown"


class ConstraintKind(StrEnum):
    DISTANCE = "distance"
    ANGLE = "angle"
    PLANARITY = "planarity"
    TORSION = "torsion"
    STEREO_FRAME = "stereo-center-frame"
    STEREO_LIGANDS = "stereo-four-ligands"
    RING_CLOSURE = "ring-closure"
    EXCLUSION = "exclusion"


def _text(value: object, label: str, *, empty: bool = False) -> None:
    if not isinstance(value, str) or (not empty and not value) or len(value) > 4096:
        raise ContractError(f"bounded {label} required")


def _number(value: object) -> TypeGuard[int | float]:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _sha(value: object) -> None:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ContractError("SHA-256 required")


def _component(value: ComponentInstanceRef) -> None:
    if not isinstance(value, ComponentInstanceRef):
        raise ContractError("exact component identity required")
    if type(value.model) is not int or type(value.occurrence) is not int or min(value.model, value.occurrence) < 1:
        raise ContractError("positive integer model/occurrence required")
    for field in ("chain_id", "sequence_number", "insertion_code", "alternate_location"):
        _text(getattr(value, field), field, empty=field in {"insertion_code", "alternate_location"})


def _identity(value: ExactAtomRef) -> None:
    if not isinstance(value, ExactAtomRef):
        raise ContractError("exact atom identity required")
    _component(value.component)
    _text(value.name, "atom name")


def _point(value: object) -> None:
    if not isinstance(value, tuple) or len(value) != 3 or any(not _number(x) or abs(x) > 1e5 for x in value):
        raise ContractError("bounded immutable finite position required")


def _element(value: object) -> None:
    # Syntax only, not an independent periodic-table/valence authority.
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z]{1,2}", value):
        raise ContractError("explicit element token required")


def _tuple(value: object, item_type: type, label: str, *, maximum: int = MAX_ITEMS) -> None:
    if not isinstance(value, tuple) or len(value) > maximum or any(not isinstance(x, item_type) for x in value):
        raise ContractError(f"bounded immutable {label} required")


def _enum_set(value: object, kind: type, label: str) -> None:
    if not isinstance(value, frozenset) or len(value) > MAX_ITEMS or any(not isinstance(x, kind) for x in value):
        raise ContractError(f"immutable typed {label} required")


def validate_graph(graph: ChemicalGraph) -> None:
    """Validate immutable graph shape, not chemical valence or ideal geometry."""
    if not isinstance(graph, ChemicalGraph):
        raise ContractError("existing domain ChemicalGraph required")
    _text(graph.component_id, "chemical component id")
    _text(graph.state_id, "state id")
    _tuple(graph.atoms, ChemicalAtom, "graph atoms", maximum=4096)
    _tuple(graph.bonds, ChemicalBond, "graph bonds", maximum=16384)
    if not graph.atoms:
        raise ContractError("empty intended graph")
    names: set[str] = set()
    for atom in graph.atoms:
        _text(atom.name, "graph atom name")
        if atom.name in names:
            raise ContractError("duplicate graph atom identity")
        names.add(atom.name)
        _element(atom.element)
        _point(atom.ideal_position)
        if type(atom.formal_charge) is not int or abs(atom.formal_charge) > 128 or type(atom.aromatic) is not bool:
            raise ContractError("explicit integer charge and boolean aromaticity required")
        if atom.stereo is not None:
            _text(atom.stereo, "stereo label")
    if type(graph.net_charge) is not int or sum(a.formal_charge for a in graph.atoms) != graph.net_charge:
        raise ContractError("graph charge mismatch")
    seen: set[frozenset[str]] = set()
    adjacent: dict[str, set[str]] = {name: set() for name in names}
    for bond in graph.bonds:
        _text(bond.atom1, "bond atom name")
        _text(bond.atom2, "bond atom name")
        pair = frozenset((bond.atom1, bond.atom2))
        if len(pair) != 2 or not pair <= names or pair in seen:
            raise ContractError("duplicate/self/unknown graph bond")
        if (not _number(bond.order) or not 0 < bond.order <= 3
                or not _number(bond.ideal_length) or bond.ideal_length <= 0 or type(bond.aromatic) is not bool):
            raise ContractError("explicit positive bond order/length required")
        seen.add(pair)
        adjacent[bond.atom1].add(bond.atom2)
        adjacent[bond.atom2].add(bond.atom1)
    unvisited, fragments = set(names), 0
    while unvisited:
        fragments += 1
        pending = [unvisited.pop()]
        while pending:
            for other in adjacent[pending.pop()] & unvisited:
                unvisited.remove(other)
                pending.append(other)
    if type(graph.component_count) is not int or graph.component_count != fragments:
        raise ContractError("declared covalent fragment count differs from graph")
    if type(graph.priority) is not int:
        raise ContractError("integer metadata priority required")
    if ((graph.ph_min is None) != (graph.ph_max is None)
            or graph.ph_min is not None and (not _number(graph.ph_min) or not _number(graph.ph_max)
                or not 0 <= graph.ph_min <= graph.ph_max <= 14)):
        raise ContractError("explicit ordered pH range required")
    for value in (graph.redox_state, graph.coordination_model):
        if value is not None:
            _text(value, "chemical model metadata")


def canonical_graph_digest(graph: ChemicalGraph) -> str:
    """Same graph payload/hash as existing Class A, without importing its solver."""
    validate_graph(graph)
    raw = json.dumps(asdict(graph), sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    return hashlib.sha256(raw.encode()).hexdigest()


@dataclass(frozen=True)
class ChemicalState:
    graph: ChemicalGraph
    authority: AuthorityProvenance
    ideal_unit: Unit
    hydrogens: HydrogenInventory = HydrogenInventory.UNSPECIFIED

    def __post_init__(self) -> None:
        validate_graph(self.graph)
        if not isinstance(self.authority, AuthorityProvenance) or not isinstance(self.authority.mode, AuthorityMode):
            raise ContractError("explicit authority provenance required")
        _text(self.authority.source, "authority source")
        _sha(self.authority.source_digest)
        if self.authority.source_version is not None:
            _text(self.authority.source_version, "authority version")
        if self.authority.mode == AuthorityMode.USER_MAPPED and self.authority.source_digest != canonical_graph_digest(self.graph):
            raise ContractError("user graph authority digest mismatch")
        if not isinstance(self.ideal_unit, Unit) or self.ideal_unit not in {Unit.ANGSTROM, Unit.NANOMETER}:
            raise ContractError("ideal-coordinate length unit required")
        if not isinstance(self.hydrogens, HydrogenInventory):
            raise ContractError("explicit H inventory semantics required")


@dataclass(frozen=True)
class Component:
    identity: ComponentInstanceRef
    kind: ChemistryKind
    states: tuple[ChemicalState, ...] = ()
    selected_state: str | None = None

    def __post_init__(self) -> None:
        _component(self.identity)
        if not isinstance(self.kind, ChemistryKind):
            raise ContractError("declared chemistry category required")
        _tuple(self.states, ChemicalState, "graph states", maximum=64)
        if len({s.graph.state_id for s in self.states}) != len(self.states):
            raise ContractError("duplicate state ids; do not merge conflicting authorities")
        if self.selected_state is not None:
            _text(self.selected_state, "explicit state choice")

    @property
    def resolved(self) -> ChemicalState | None:
        if self.selected_state is not None:
            return next((s for s in self.states if s.graph.state_id == self.selected_state), None)
        return self.states[0] if len(self.states) == 1 else None


@dataclass(frozen=True)
class Observation(ObservedAtom):
    unit: Unit
    occupancy: float = 1.0

    def __post_init__(self) -> None:
        _identity(self.identity)
        _element(self.element)
        _point(self.position)
        if not isinstance(self.unit, Unit) or self.unit not in {Unit.ANGSTROM, Unit.NANOMETER}:
            raise ContractError("observed-coordinate length unit required")
        if not _number(self.occupancy) or not 0 < self.occupancy <= 1:
            raise ContractError("occupancy in (0,1] required")

    @property
    def angstrom(self) -> tuple[float, float, float]:
        factor = 10.0 if self.unit == Unit.NANOMETER else 1.0
        x, y, z = self.position
        return x*factor, y*factor, z*factor


@dataclass(frozen=True)
class Ownership:
    generated: frozenset[ExactAtomRef] = frozenset()
    temporary: frozenset[ExactAtomRef] = frozenset()
    mobile: frozenset[ExactAtomRef] = frozenset()

    def __post_init__(self) -> None:
        for value in (self.generated, self.temporary, self.mobile):
            _enum_set(value, ExactAtomRef, "atom ownership")
            for atom in value:
                _identity(atom)
        if self.generated & self.temporary or self.mobile & (self.generated | self.temporary):
            raise ContractError("ownership roles overlap")


@dataclass(frozen=True)
class Connection:
    link: ExternalLink
    kind: LinkKind
    model: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.link, ExternalLink) or not isinstance(self.kind, LinkKind):
            raise ContractError("typed explicit connection required")
        _identity(self.link.component_atom)
        _identity(self.link.other_atom)
        if self.link.component_atom == self.link.other_atom:
            raise ContractError("self connection")
        if self.link.order is not None and (not _number(self.link.order) or not 0 < self.link.order <= 3):
            raise ContractError("invalid explicit connection order")
        if self.kind != LinkKind.COVALENT and self.link.order is not None:
            raise ContractError("do not substitute covalent order for unknown/coordination semantics")
        if self.model is not None:
            _text(self.model, "connection model")


@dataclass(frozen=True)
class Constraint:
    kind: ConstraintKind
    atoms: tuple[ExactAtomRef, ...]
    unit: Unit
    lower: float | None = None
    upper: float | None = None
    orientation: int | None = None
    minimum_volume: float | None = None
    label: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.kind, ConstraintKind) or not isinstance(self.unit, Unit):
            raise ContractError("typed constraint kind/unit required")
        _tuple(self.atoms, ExactAtomRef, "constraint references", maximum=64)
        for atom in self.atoms:
            _identity(atom)
        arity = {ConstraintKind.DISTANCE: 2, ConstraintKind.EXCLUSION: 2, ConstraintKind.RING_CLOSURE: 2,
                 ConstraintKind.ANGLE: 3, ConstraintKind.TORSION: 4, ConstraintKind.STEREO_FRAME: 4,
                 ConstraintKind.STEREO_LIGANDS: 5}
        if len(set(self.atoms)) != len(self.atoms) or (self.kind in arity and len(self.atoms) != arity[self.kind]):
            raise ContractError("distinct exact constraint arity required")
        if self.kind == ConstraintKind.PLANARITY and len(self.atoms) < 4:
            raise ContractError("planarity needs plane anchors and at least one tested atom")
        stereo = self.kind in {ConstraintKind.STEREO_FRAME, ConstraintKind.STEREO_LIGANDS}
        if stereo:
            if (self.unit != Unit.DIMENSIONLESS or type(self.orientation) is not int or self.orientation not in {-1, 1}
                    or not _number(self.minimum_volume) or not 0 < self.minimum_volume <= 1
                    or self.lower is not None or self.upper is not None or not self.label):
                raise ContractError("explicit ordered stereo parity/label/normalized bound required")
        else:
            if self.orientation is not None or self.minimum_volume is not None or not _number(self.upper):
                raise ContractError("bounded non-stereo interval required")
            if self.kind == ConstraintKind.PLANARITY:
                if self.lower is not None or self.upper <= 0:
                    raise ContractError("positive plane distance tolerance required")
            elif not _number(self.lower) or self.lower > self.upper:
                raise ContractError("ordered constraint interval required")
            angular = self.kind in {ConstraintKind.ANGLE, ConstraintKind.TORSION}
            if angular:
                assert self.lower is not None
                if self.unit != Unit.DEGREE or not -180 <= self.lower <= self.upper <= 180:
                    raise ContractError("angular interval in degrees required")
                if self.kind == ConstraintKind.ANGLE and self.lower < 0:
                    raise ContractError("negative bond angle")
            elif self.unit not in {Unit.ANGSTROM, Unit.NANOMETER} or (self.lower is not None and self.lower < 0):
                raise ContractError("nonnegative explicit length interval required")
        _text(self.label, "constraint label", empty=not stereo)


@dataclass(frozen=True)
class Origin:
    adapter: str
    request_digest: str
    settings: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        _text(self.adapter, "adapter revision")
        _sha(self.request_digest)
        _tuple(self.settings, tuple, "origin settings", maximum=64)
        for row in self.settings:
            if len(row) != 2:
                raise ContractError("origin key/value pairs required")
            for item in row:
                _text(item, "origin setting", empty=True)
        if len({row[0] for row in self.settings}) != len(self.settings):
            raise ContractError("duplicate origin setting")


@dataclass(frozen=True)
class ReconstructionRequest:
    source_sha256: str
    operation: Operation
    components: tuple[Component, ...]
    observations: tuple[Observation, ...]
    targets: frozenset[ComponentInstanceRef]
    ownership: Ownership
    connections: tuple[Connection, ...] = ()
    constraints: tuple[Constraint, ...] = ()
    ph: float | None = None
    origin: Origin | None = None

    def __post_init__(self) -> None:
        _sha(self.source_sha256)
        if not isinstance(self.operation, Operation) or not isinstance(self.ownership, Ownership):
            raise ContractError("typed operation/ownership required")
        _tuple(self.components, Component, "components", maximum=4096)
        _tuple(self.observations, Observation, "observations")
        _tuple(self.connections, Connection, "connections")
        _tuple(self.constraints, Constraint, "constraints")
        _enum_set(self.targets, ComponentInstanceRef, "target components")
        for target in self.targets:
            _component(target)
        if not self.components or not self.targets:
            raise ContractError("explicit components/targets required")
        if len({c.identity for c in self.components}) != len(self.components):
            raise ContractError("duplicate component identity")
        if len({o.identity for o in self.observations}) != len(self.observations):
            raise ContractError("duplicate observation")
        if self.ph is not None and (not _number(self.ph) or not 0 <= self.ph <= 14):
            raise ContractError("explicit bounded pH metadata required")
        if self.origin is not None and not isinstance(self.origin, Origin):
            raise ContractError("typed adapter origin required")
        if sum(sum(len(s.graph.atoms) for s in c.states) for c in self.components) > MAX_ITEMS:
            raise ContractError("aggregate authority inventory exceeds bound")

    @property
    def digest(self) -> str:
        return digest_record((REVISION, self))


@dataclass(frozen=True)
class OpaqueContextPolicy:
    prepared_source_sha256: str
    preparation_digest: str
    fractional_atoms: frozenset[ExactAtomRef]
    observations_digest: str

    def __post_init__(self) -> None:
        _sha(self.prepared_source_sha256)
        _sha(self.preparation_digest)
        _sha(self.observations_digest)
        _enum_set(self.fractional_atoms, ExactAtomRef, "exact fractional opaque context")
        for atom in self.fractional_atoms:
            _identity(atom)


@dataclass(frozen=True, kw_only=True)
class PreparedContextRequest(ReconstructionRequest):
    """Explicit conservative fractional context; old requests keep old semantics."""
    context_policy: OpaqueContextPolicy

    def __post_init__(self) -> None:
        super().__post_init__()
        policy = self.context_policy
        if not isinstance(policy, OpaqueContextPolicy) or policy.prepared_source_sha256 != self.source_sha256:
            raise ContractError("prepared-context-source-binding-mismatch")
        if (self.origin is None or self.origin.adapter != "prepared-observed-peptide-v1"
                or self.origin.request_digest != policy.preparation_digest):
            raise ContractError("prepared-context-evidence-binding-mismatch")
        if policy.observations_digest != digest_record(self.observations):
            raise ContractError("prepared-context-observation-binding-mismatch")
        fractional = frozenset(o.identity for o in self.observations if o.occupancy != 1)
        if fractional != policy.fractional_atoms:
            raise ContractError("fractional-context-exact-scope-required")
        components = {c.identity: c for c in self.components}
        forbidden = self.ownership.generated | self.ownership.temporary | self.ownership.mobile
        forbidden |= frozenset(a for c in self.constraints for a in c.atoms)
        forbidden |= frozenset(a for c in self.connections for a in (c.link.component_atom, c.link.other_atom))
        if fractional & forbidden:
            raise ContractError("fractional-context-cannot-be-owned-constrained-or-linked")
        for atom in fractional:
            component = components.get(atom.component)
            if (component is None or atom.component in self.targets or component.kind != ChemistryKind.UNKNOWN
                    or component.states or atom.component.alternate_location):
                raise ContractError("fractional-context-must-be-unselected-opaque")

    @property
    def digest(self) -> str:
        return digest_record(("chemical-prepared-context-request-v1", self))


def digest_record(value: object) -> str:
    def canonical(item: object) -> object:
        if isinstance(item, StrEnum):
            return item.value
        if is_dataclass(item) and not isinstance(item, type):
            return {f.name: canonical(getattr(item, f.name)) for f in fields(item)}
        if isinstance(item, tuple):
            return [canonical(x) for x in item]
        if isinstance(item, frozenset):
            result = [canonical(x) for x in item]
            return sorted(result, key=lambda x: json.dumps(x, sort_keys=True))
        if item is None or type(item) in (str, int, float, bool):
            return item
        raise ContractError("unsupported mutable/foreign digest value")

    raw = json.dumps(canonical(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    return hashlib.sha256(raw.encode()).hexdigest()


@dataclass(frozen=True)
class Capability:
    """Caller-supplied strategy declaration, not registered/qualified execution."""
    identifier: str
    operations: frozenset[Operation]
    chemistry: frozenset[ChemistryKind]
    elements: frozenset[str]
    constraints: frozenset[ConstraintKind]
    links: frozenset[LinkKind] = frozenset()
    allow_opaque_context: bool = False
    allow_mobile: bool = False
    allow_temporary: bool = False
    minimum_anchors: int = 1

    def __post_init__(self) -> None:
        _text(self.identifier, "strategy id")
        for value, kind, label in ((self.operations, Operation, "operations"), (self.chemistry, ChemistryKind, "chemistry"),
                                   (self.constraints, ConstraintKind, "constraint kinds"), (self.links, LinkKind, "link kinds")):
            _enum_set(value, kind, label)
        _enum_set(self.elements, str, "elements")
        if not self.operations or not self.chemistry or not self.elements:
            raise ContractError("explicit nonempty strategy scope required")
        for element in self.elements:
            _element(element)
            if element != element.upper():
                raise ContractError("capability element tokens must be uppercase")
        if any(type(x) is not bool for x in (self.allow_opaque_context, self.allow_mobile, self.allow_temporary)):
            raise ContractError("explicit boolean capability switches required")
        if type(self.minimum_anchors) is not int or not 0 <= self.minimum_anchors <= MAX_ITEMS:
            raise ContractError("bounded anchor count required")


class AdmissionStatus(StrEnum):
    ADMITTED = "admitted-contract-only"
    INVALID_INPUT = "invalid-input"
    AMBIGUOUS = "ambiguous"
    NOT_ASSESSED = "not-assessed"
    UNSUPPORTED = "unsupported"


@dataclass(frozen=True)
class Finding:
    status: AdmissionStatus
    code: str
    atoms: tuple[ExactAtomRef, ...] = ()


@dataclass(frozen=True)
class Admission:
    status: AdmissionStatus
    request_digest: str
    findings: tuple[Finding, ...]
    capability_id: str | None
    unassessed_context: tuple[ComponentInstanceRef, ...]
    checked_constraints: tuple[int, ...]
    deferred_constraints: tuple[int, ...]

    @property
    def coordinates_produced(self) -> Literal[False]:
        return False

    @property
    def geometry_approved(self) -> Literal[False]:
        return False

    @property
    def md_ready(self) -> Literal[False]:
        return False

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False
