"""Pure fail-closed request semantics and caller-declared capability matching.

No solver, toolkit, authority lookup, source reader or production registration.
Admission is deliberately weaker than scientific or full-source qualification.
"""
from __future__ import annotations

from collections import Counter

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import (
    Admission,
    AdmissionStatus,
    Capability,
    ChemistryKind,
    ConstraintKind,
    Finding,
    HydrogenInventory,
    LinkKind,
    Operation,
    PreparedContextRequest,
    ReconstructionRequest,
)
from .heavy_conditioning import HeavyConditionedRequest
from .measurement import satisfied


def admit(request: ReconstructionRequest, capabilities: tuple[Capability, ...] = ()) -> Admission:
    """Check supplied records without generating coordinates or reading source."""
    if (not isinstance(request, ReconstructionRequest) or not isinstance(capabilities, tuple)
            or len(capabilities) > 256 or any(not isinstance(c, Capability) for c in capabilities)
            or len({c.identifier for c in capabilities}) != len(capabilities)):
        raise ValueError("typed request and bounded unique capability declarations required")
    findings: list[Finding] = []

    def fail(code: str, atoms: tuple[ExactAtomRef, ...] = (),
             status: AdmissionStatus = AdmissionStatus.INVALID_INPUT) -> None:
        findings.append(Finding(status, code, atoms))

    components = {c.identity: c for c in request.components}
    states = {c.identity: c.resolved for c in request.components}
    observed = {o.identity: o for o in request.observations}
    ownership = request.ownership
    new = ownership.generated | ownership.temporary
    active = new | ownership.mobile
    present = set(observed) | new
    opaque = tuple(c.identity for c in request.components if c.resolved is None and c.identity not in request.targets)
    checked: list[int] = []
    deferred: list[int] = []

    if not request.targets <= components.keys():
        fail("undeclared-target-component")
    if any(a.component not in components for a in present | ownership.mobile):
        fail("undeclared-atom-component")
    if any(a.component not in request.targets for a in active):
        fail("ownership-outside-target")
    if new & observed.keys():
        fail("new-atom-is-observed", tuple(sorted(new & observed.keys())))
    if not ownership.mobile <= observed.keys():
        fail("mobile-atom-not-observed")
    if request.operation != Operation.REPACK and ownership.mobile:
        fail("completion-cannot-repack")
    if request.operation == Operation.REPACK and (new or not ownership.mobile):
        fail("repack-requires-only-observed-mobile-atoms")
    if request.operation != Operation.REPACK and not ownership.generated:
        fail("completion-requires-generated-mask")
    allowed_fractional = request.context_policy.fractional_atoms if isinstance(request, PreparedContextRequest) else frozenset()
    if isinstance(request, HeavyConditionedRequest):
        allowed_fractional = request.allowed_fractional_atoms
    if (any(o.occupancy != 1 and o.identity not in allowed_fractional for o in request.observations)
            or any(c.identity.alternate_location for c in request.components)):
        fail("partial-occupancy-or-altloc-requires-explicit-selection", status=AdmissionStatus.UNSUPPORTED)
    if len({c.identity.model for c in request.components}) != 1:
        fail("multiple-models-require-selection", status=AdmissionStatus.UNSUPPORTED)

    for component in request.components:
        if component.selected_state is not None and component.resolved is None:
            fail("unknown-selected-state")
        if component.identity in request.targets:
            if not component.states:
                fail("target-chemistry-not-assessed", status=AdmissionStatus.NOT_ASSESSED)
            elif component.resolved is None:
                fail("target-state-ambiguous", status=AdmissionStatus.AMBIGUOUS)
            if component.kind == ChemistryKind.UNKNOWN:
                fail("target-chemistry-category-not-assessed", status=AdmissionStatus.NOT_ASSESSED)

    graph_atoms = {ExactAtomRef(identity, atom.name): atom
                   for identity, state in states.items() if state is not None for atom in state.graph.atoms}
    intended_by_component = {identity: set[ExactAtomRef]() for identity in components}
    for atom in graph_atoms:
        intended_by_component[atom.component].add(atom)
    new_h_components = {a.component for a in new if a in graph_atoms and graph_atoms[a].element.upper() in {"H", "D"}}
    for atom in sorted(present | ownership.mobile):
        if states.get(atom.component) is not None:
            if atom not in graph_atoms:
                fail("atom-not-in-selected-graph", (atom,))
            elif atom in observed and observed[atom].element.upper() != graph_atoms[atom].element.upper():
                fail("observed-element-conflicts-with-graph", (atom,))
    for atom in sorted(active):
        if atom not in graph_atoms:
            fail("owned-atom-requires-resolved-graph", (atom,), AdmissionStatus.NOT_ASSESSED)
    for atom in sorted(ownership.generated):
        if atom in graph_atoms:
            hydrogen = graph_atoms[atom].element.upper() in {"H", "D"}
            if request.operation == Operation.HEAVY and hydrogen or request.operation == Operation.HYDROGENS and not hydrogen:
                fail("generated-element-conflicts-with-operation", (atom,))
    for identity in sorted(request.targets):
        state = states.get(identity)
        if state is None:
            continue
        intended = intended_by_component[identity]
        heavy = {a for a in intended if graph_atoms[a].element.upper() not in {"H", "D"}}
        if not heavy <= present:
            fail("unowned-missing-heavy-atoms", tuple(sorted(heavy-present)))
        places_h = request.operation == Operation.HYDROGENS or identity in new_h_components
        if places_h and state.hydrogens != HydrogenInventory.EXPLICIT:
            fail("hydrogen-inventory-not-explicit", status=AdmissionStatus.NOT_ASSESSED)
        if request.operation == Operation.HYDROGENS:
            if not intended <= present:
                fail("unowned-missing-hydrogens", tuple(sorted(intended-present)))

    adjacency: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in graph_atoms}
    for identity, state in states.items():
        if state is not None:
            for bond in state.graph.bonds:
                a, b = ExactAtomRef(identity, bond.atom1), ExactAtomRef(identity, bond.atom2)
                adjacency[a].add(b)
                adjacency[b].add(a)
    pairs: set[frozenset[ExactAtomRef]] = set()
    for connection in request.connections:
        a, b = connection.link.component_atom, connection.link.other_atom
        pair = frozenset((a, b))
        if pair in pairs:
            fail("duplicate-connection", (a, b))
        pairs.add(pair)
        if a not in present or b not in present:
            fail("connection-endpoint-not-materialized-or-owned", (a, b))
        if a.component == b.component:
            fail("external-connection-within-component", (a, b))
        if (connection.kind == LinkKind.UNKNOWN
                or connection.kind == LinkKind.COVALENT and connection.link.order is None
                or connection.kind == LinkKind.COORDINATION and connection.model is None):
            fail("unresolved-connection-chemistry", (a, b), AdmissionStatus.AMBIGUOUS)
        if connection.kind == LinkKind.COVALENT:
            adjacency.setdefault(a, set()).add(b)
            adjacency.setdefault(b, set()).add(a)

    stereo_centres: set[ExactAtomRef] = set()
    for index, constraint in enumerate(request.constraints):
        if not set(constraint.atoms) <= present:
            fail("constraint-reference-not-materialized-or-owned", constraint.atoms)
            continue
        if constraint.kind in {ConstraintKind.STEREO_FRAME, ConstraintKind.STEREO_LIGANDS}:
            centre, *ligands = constraint.atoms
            if not set(ligands) <= adjacency.get(centre, set()):
                fail("stereo-ligands-not-authoritative-neighbours", constraint.atoms)
            if centre in graph_atoms and graph_atoms[centre].stereo is not None:
                if constraint.label != graph_atoms[centre].stereo:
                    fail("stereo-label-conflicts-with-authority", constraint.atoms)
            stereo_centres.add(centre)
        if any(a not in observed or a in ownership.mobile for a in constraint.atoms):
            deferred.append(index)
        else:
            checked.append(index)
            try:
                valid = satisfied(constraint, tuple(observed[a].angstrom for a in constraint.atoms))
            except ValueError:
                fail("degenerate-fixed-constraint", constraint.atoms)
            else:
                if not valid:
                    fail("fixed-constraint-conflict", constraint.atoms)
    for atom, chemical in graph_atoms.items():
        if atom.component in request.targets and chemical.stereo is not None and atom not in stereo_centres:
            fail("ordered-stereo-mapping-required", (atom,), AdmissionStatus.NOT_ASSESSED)

    selected: str | None = None
    if not findings:
        elements = {o.element.upper() for o in request.observations} | {a.element.upper() for a in graph_atoms.values()}
        anchors = Counter(a.component for a in observed if a not in ownership.mobile)
        for capability in capabilities:
            anchors_ok = all(anchors[target] >= capability.minimum_anchors for target in request.targets)
            if (request.operation in capability.operations
                    and all(components[t].kind in capability.chemistry for t in request.targets)
                    and elements <= capability.elements
                    and {c.kind for c in request.constraints} <= capability.constraints
                    and {c.kind for c in request.connections} <= capability.links
                    and (not opaque or capability.allow_opaque_context)
                    and (not ownership.mobile or capability.allow_mobile)
                    and (not ownership.temporary or capability.allow_temporary) and anchors_ok):
                selected = capability.identifier
                break
        if selected is None:
            fail("no-declared-capability-matches", status=AdmissionStatus.UNSUPPORTED)
    priority = (AdmissionStatus.INVALID_INPUT, AdmissionStatus.AMBIGUOUS,
                AdmissionStatus.NOT_ASSESSED, AdmissionStatus.UNSUPPORTED)
    status = next((s for s in priority if any(f.status == s for f in findings)), AdmissionStatus.ADMITTED)
    return Admission(status, request.digest, tuple(findings), selected, opaque, tuple(checked), tuple(deferred))
