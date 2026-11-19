"""Explicit leaving-atom/state declarations; no protonation or graph-edit inference."""
from __future__ import annotations

from dataclasses import dataclass, replace

from dvbfixer.domain.nonprotein_reconstruction import (
    AuthorityMode,
    AuthorityProvenance,
    ChemicalGraph,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils.variants import VARIANT_TO_PARENT

from .amide_authority import AmideMetricAuthority
from .ccd_authority import CcdGeometryAuthority
from .contracts import (
    ChemicalState,
    Constraint,
    ConstraintKind,
    HydrogenInventory,
    LinkKind,
    Origin,
    ReconstructionRequest,
    Unit,
    _identity,
    _sha,
    _text,
    canonical_graph_digest,
    digest_record,
)
from .hydrogen_completion import (
    CompleteHydrogenView,
    compile_extended_hydrogens,
    compile_joint_hydrogens,
)
from .peptide_state import (
    ChemicalWitness,
    LocalPeptideAuthority,
    PeptideStateError,
    PeptideSubstitution,
    TemplateInventory,
    assemble_witness,
)
from .proposals import AuditPolicy


@dataclass(frozen=True)
class PeptideStateDeclaration:
    target: ComponentInstanceRef
    requested_state: str
    parent_digest: str
    alanine_digest: str
    inventory_digest: str
    substitutions: tuple[PeptideSubstitution, ...]

    def __post_init__(self) -> None:
        _identity(ExactAtomRef(self.target, "N"))
        _text(self.requested_state, "declared peptide state")
        for pin in (self.parent_digest, self.alanine_digest, self.inventory_digest):
            _sha(pin)
        if (not isinstance(self.substitutions, tuple) or len(self.substitutions) != 2
                or any(not isinstance(s, PeptideSubstitution) for s in self.substitutions)):
            raise PeptideStateError("exact-two-peptide-substitutions-required")

    @property
    def digest(self) -> str:
        return digest_record(("explicit-peptide-state-declaration-v1", self))


def _values(parent: CcdGeometryAuthority, alanine: CcdGeometryAuthority, identity: ComponentInstanceRef,
             inventory: TemplateInventory, declaration: PeptideStateDeclaration
             ) -> tuple[ChemicalWitness, ChemicalState, frozenset[str]]:
    if (declaration.target != identity or declaration.parent_digest != parent.digest
            or declaration.alanine_digest != alanine.digest or declaration.inventory_digest != digest_record(inventory)
            or declaration.requested_state != inventory.state_name):
        raise PeptideStateError("declared-peptide-source-state-inventory-binding")
    if VARIANT_TO_PARENT.get(inventory.state_name, inventory.state_name) != parent.component.component_id:
        raise PeptideStateError("declared-state-parent-name-mismatch")
    if alanine.component.component_id != "ALA":
        raise PeptideStateError("this version declares neutral ALA witnesses only")
    left = ComponentInstanceRef(identity.model, "__left_chemical_witness__", "1")
    right = ComponentInstanceRef(identity.model, "__right_chemical_witness__", "1")
    if identity in {left, right}:
        raise PeptideStateError("source identity collides with chemical witness")
    edges = tuple((s.carbon.component, s.nitrogen.component) for s in declaration.substitutions)
    if edges != ((left, identity), (identity, right)):
        raise PeptideStateError("declared-peptide-directed-role-mismatch")
    witness = assemble_witness(((left, alanine), (identity, parent), (right, alanine)), declaration.substitutions)
    graph: ChemicalGraph = dict(witness.projected)[identity]
    if inventory.atoms != {a.name for a in graph.atoms} or inventory.bonds != {
            frozenset((b.atom1, b.atom2)) for b in graph.bonds}:
        raise PeptideStateError("declared-state-unchanged-graph-inventory-mismatch")
    state = ChemicalState(graph, AuthorityProvenance(AuthorityMode.USER_MAPPED,
                          f"declared conditional peptide; declaration {declaration.digest}; witness {witness.digest}",
                          canonical_graph_digest(graph), "declared-peptide-state-v1"), Unit.ANGSTROM, HydrogenInventory.EXPLICIT)
    removed = frozenset(a.name for a in parent.state.graph.atoms if a.name not in inventory.atoms)
    return witness, state, removed


@dataclass(frozen=True)
class DeclaredPeptideAuthority(LocalPeptideAuthority):
    declaration: PeptideStateDeclaration

    def __post_init__(self) -> None:
        if len(self.witness.source_units) != 3:
            raise PeptideStateError("declared-peptide-needs-three-source-units")
        alanine = self.witness.source_units[0][1]
        witness, state, removed = _values(self.parent, alanine, self.identity, self.inventory, self.declaration)
        if (self.witness, self.state, self.removed_atoms) != (witness, state, removed):
            raise PeptideStateError("declared-peptide-authority-projection-binding")

    @property
    def digest(self) -> str:
        return digest_record(("declared-local-peptide-authority-v1", self))


def derive_declared_peptide(parent: CcdGeometryAuthority, alanine_witness: CcdGeometryAuthority,
                             identity: ComponentInstanceRef, inventory: TemplateInventory,
                             declaration: PeptideStateDeclaration) -> DeclaredPeptideAuthority:
    """Prove the caller's exact substitutions; never choose a state or edit its graph."""
    witness, state, removed = _values(parent, alanine_witness, identity, inventory, declaration)
    return DeclaredPeptideAuthority(parent, inventory, witness, identity, state, removed, declaration)


def compile_declared_hydrogens(base: ReconstructionRequest, authority: DeclaredPeptideAuthority,
                                incoming_carbon: ExactAtomRef, parameters: AmideMetricAuthority, *,
                                backbone_types: tuple[str, str, str] | None, fixed_policy: AuditPolicy,
                                generated_policy: AuditPolicy) -> CompleteHydrogenView:
    if not isinstance(authority, DeclaredPeptideAuthority):
        raise PeptideStateError("explicit-declared-peptide-authority-required")
    nitrogen = ExactAtomRef(authority.identity, "N")
    pair = frozenset((incoming_carbon, nitrogen))
    if (incoming_carbon.component == nitrogen.component or incoming_carbon.name != "C"
            or not any(c.kind == LinkKind.COVALENT and c.link.order == 1
                       and frozenset((c.link.component_atom, c.link.other_atom)) == pair for c in base.connections)):
        raise PeptideStateError("declared-H-needs-bound-real-incoming-peptide-bond")
    incoming_observation = next((o for o in base.observations if o.identity == incoming_carbon), None)
    active = base.ownership.generated | base.ownership.temporary | base.ownership.mobile
    if (incoming_observation is None or incoming_observation.element != "C" or incoming_observation.occupancy != 1
            or incoming_carbon in active):
        raise PeptideStateError("declared-H-needs-fixed-unit-occupancy-incoming-carbon")
    graph = authority.state.graph
    atoms = {a.name: a for a in graph.atoms}
    adjacent = sorted(b.atom2 if b.atom1 == "N" else b.atom1 for b in graph.bonds if "N" in {b.atom1, b.atom2})
    if not adjacent or "N" not in atoms:
        raise PeptideStateError("declared-H-needs-explicit-backbone-nitrogen")
    has_H = any(atoms[n].element in {"H", "D"} for n in adjacent)
    constraints = base.constraints
    if not has_H:
        if backbone_types is not None or len(adjacent) != 2 or "CA" not in adjacent:
            raise PeptideStateError("no-backbone-H-needs-two-heavy-neighbours-and-no-H-typing")
        plane = (nitrogen, incoming_carbon, *(ExactAtomRef(authority.identity, n) for n in adjacent))
        if not set(plane) <= {o.identity for o in base.observations}:
            raise PeptideStateError("declared-backbone-plane-needs-real-fixed-observations")
        constraints = (*constraints, Constraint(ConstraintKind.PLANARITY, plane, Unit.ANGSTROM, upper=.05))
    elif backbone_types is None:
        raise PeptideStateError("declared-backbone-H-requires-explicit-parameter-types")
    settings = (("peptide_state_declaration", authority.declaration.digest), ("declared_peptide_authority", authority.digest))
    origin = (replace(base.origin, settings=(*base.origin.settings, *settings)) if base.origin is not None else
              Origin("declared-peptide-H-v1", authority.digest, settings))
    bound = replace(base, constraints=constraints, origin=origin)
    return compile_joint_hydrogens(compile_extended_hydrogens(bound, authority, incoming_carbon, parameters,
                                  backbone_types=backbone_types, fixed_policy=fixed_policy, generated_policy=generated_policy))
