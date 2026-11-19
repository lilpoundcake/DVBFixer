"""Explicit, in-memory compatibility adapters; no authority resolution or solving."""
from __future__ import annotations

from dvbfixer.domain.nonprotein_reconstruction import (
    AuthorityProvenance,
    ComponentClass,
    ExternalLink,
    GeometryPolicy,
    UserMappedAuthority,
)
from dvbfixer.domain.nonprotein_reconstruction import (
    ReconstructionRequest as ClassARequest,
)
from dvbfixer.domain.parameterization import COMPLEX_COFACTORS
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.protein_completion.contracts import CompletionRequest

from .contracts import (
    ChemicalState,
    ChemistryKind,
    Component,
    Connection,
    Constraint,
    ContractError,
    HydrogenInventory,
    LinkKind,
    Observation,
    Operation,
    Origin,
    Ownership,
    ReconstructionRequest,
    Unit,
    digest_record,
)


class AdapterError(ContractError):
    """Explicit refused transport; not permission to silently use another engine."""

    def __init__(self, code: str) -> None:
        self.code = code
        super().__init__(code)


def _components(target: Component, context: tuple[Component, ...], observations: tuple[Observation, ...],
                connections: tuple[Connection, ...]) -> tuple[Component, ...]:
    if not isinstance(context, tuple) or any(not isinstance(c, Component) for c in context):
        raise AdapterError("typed-immutable-context-required")
    declared = {target.identity, *(c.identity for c in context)}
    if len(declared) != len(context)+1:
        raise AdapterError("duplicate-context-or-target-chemistry")
    identities = {o.identity.component for o in observations}
    for connection in connections:
        identities.update((connection.link.component_atom.component, connection.link.other_atom.component))
    opaque = tuple(Component(i, ChemistryKind.UNKNOWN) for i in sorted(identities-declared))
    return (target, *context, *opaque)


def from_protein(request: CompletionRequest, *, state: ChemicalState,
                 context: tuple[Component, ...] = (), connections: tuple[Connection, ...] | None = None,
                 constraints: tuple[Constraint, ...] = ()) -> ReconstructionRequest:
    """Graph/state authority is mandatory; residue names supply no chemistry."""
    if not isinstance(request, CompletionRequest) or not isinstance(state, ChemicalState):
        raise AdapterError("typed-protein-request-and-explicit-state-required")
    observations = tuple(Observation(o.atom, o.element, o.xyz, Unit.ANGSTROM, o.occupancy)
                         for o in request.observations)
    if connections is None:
        connections = tuple(Connection(ExternalLink(a, b), LinkKind.UNKNOWN) for a, b in request.links)
    if not isinstance(connections, tuple) or any(not isinstance(c, Connection) for c in connections):
        raise AdapterError("typed-immutable-connections-required")
    pairs = {frozenset((c.link.component_atom, c.link.other_atom)) for c in connections}
    if any(frozenset(link) not in pairs for link in request.links):
        raise AdapterError("legacy-link-omitted")
    target = Component(request.residue, ChemistryKind.PROTEIN, (state,))
    return ReconstructionRequest(
        request.source_sha256, Operation.HEAVY, _components(target, context, observations, connections),
        observations, frozenset({request.residue}),
        Ownership(frozenset(ExactAtomRef(request.residue, name) for name in request.missing_atoms)),
        connections, constraints,
        origin=Origin("protein-contract-adapter-v1", request.digest, (("residue_name", request.residue_name),)),
    )


def from_class_a(request: ClassARequest, *, source_sha256: str,
                 context_observations: tuple[Observation, ...] | None = None,
                 context: tuple[Component, ...] = (), context_connections: tuple[Connection, ...] = (),
                 constraints: tuple[Constraint, ...] = (), operation: Operation = Operation.HEAVY,
                 hydrogens: HydrogenInventory = HydrogenInventory.UNSPECIFIED) -> ReconstructionRequest:
    """Transport an isolated, in-memory CHNO Class A subset, never read a CCD.

    This deliberately narrow adapter envelope is not the existing solver's full
    Class A classification policy. Other organic elements remain unsupported at
    this compatibility boundary until a shared classification policy is extracted.
    """
    if not isinstance(request, ClassARequest):
        raise AdapterError("typed-class-a-request-required")
    if not isinstance(request.authority, UserMappedAuthority):
        raise AdapterError("unresolved-authority-requires-explicit-resolution-outside-adapter")
    if (not isinstance(request.requested_class, ComponentClass)
            or request.requested_class not in {ComponentClass.AUTO, ComponentClass.A} or request.external_links):
        raise AdapterError("isolated-class-a-only")
    graph = request.authority.graph
    state = ChemicalState(graph, AuthorityProvenance(request.authority.mode, request.authority.source_label,
                          request.authority.source_digest, None), Unit.ANGSTROM, hydrogens)
    if (graph.component_count != 1 or graph.redox_state is not None or graph.coordination_model is not None
            or graph.component_id.upper() in COMPLEX_COFACTORS
            or not {a.element.upper() for a in graph.atoms} <= {"C", "H", "N", "O"}
            or not any(a.element.upper() == "C" for a in graph.atoms)):
        raise AdapterError("outside-initial-chno-class-a-adapter-envelope")
    if not isinstance(request.geometry_policy, GeometryPolicy):
        raise AdapterError("typed-legacy-geometry-policy-required")
    if request.geometry_policy.max_observed_movement != 0:
        raise AdapterError("legacy-movement-policy-requires-explicit-repack-contract")
    if operation not in {Operation.HEAVY, Operation.HYDROGENS}:
        raise AdapterError("class-a-adapter-completion-only")
    original = tuple(Observation(o.identity, o.element, o.position, Unit.ANGSTROM) for o in request.observed_atoms)
    if any(o.identity.component != request.component for o in original):
        raise AdapterError("legacy-observation-outside-target")
    observations = original if context_observations is None else context_observations
    if not isinstance(observations, tuple) or any(not isinstance(o, Observation) for o in observations):
        raise AdapterError("typed-immutable-observations-required")
    supplied = {o.identity: o for o in observations}
    if any(supplied.get(o.identity) != o for o in original):
        raise AdapterError("legacy-observations-omitted-or-changed")
    if not isinstance(context_connections, tuple) or any(not isinstance(c, Connection) for c in context_connections):
        raise AdapterError("typed-immutable-connections-required")
    if any(request.component in {c.link.component_atom.component, c.link.other_atom.component}
           for c in context_connections):
        raise AdapterError("isolated-target-cannot-have-context-link")
    missing = frozenset(ExactAtomRef(request.component, a.name) for a in graph.atoms
                        if ExactAtomRef(request.component, a.name) not in supplied
                        and (a.element.upper() == "H") == (operation == Operation.HYDROGENS))
    target = Component(request.component, ChemistryKind.ISOLATED, (state,), request.microstate_override)
    policy = request.geometry_policy
    settings = (("requested_class", request.requested_class.value),
                ("max_anchor_rmsd", repr(policy.max_anchor_rmsd)),
                ("max_bond_deviation", repr(policy.max_bond_deviation)),
                ("min_nonbonded_distance", repr(policy.min_nonbonded_distance)),
                ("max_observed_movement", repr(policy.max_observed_movement)))
    return ReconstructionRequest(
        source_sha256, operation, _components(target, context, observations, context_connections),
        observations, frozenset({request.component}), Ownership(missing), context_connections, constraints,
        request.ph, Origin("class-a-memory-adapter-v1", digest_record(request), settings),
    )
