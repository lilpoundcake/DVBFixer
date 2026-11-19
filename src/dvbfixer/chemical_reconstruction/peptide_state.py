"""Explicit conditional peptide chemistry; never creates witness coordinates for output."""
from __future__ import annotations

from dataclasses import dataclass, replace

from dvbfixer.domain.nonprotein_reconstruction import (
    AuthorityMode,
    AuthorityProvenance,
    ChemicalGraph,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

from .ccd_authority import CcdGeometryAuthority, ordered_constraints, verify_cip
from .contracts import (
    ChemicalState,
    Constraint,
    ConstraintKind,
    ContractError,
    HydrogenInventory,
    Unit,
    _identity,
    _sha,
    canonical_graph_digest,
    digest_record,
    validate_graph,
)


class PeptideStateError(ContractError):
    """A declared state/substitution is unqualified; do not guess or repair it."""


@dataclass(frozen=True)
class TemplateInventory:
    source: str
    source_sha256: str
    state_name: str
    atoms: frozenset[str]
    bonds: frozenset[frozenset[str]]
    aliases: tuple[tuple[str, str], ...] = ()
    external_tokens: tuple[tuple[str, ...], ...] = ()

    def __post_init__(self) -> None:
        _sha(self.source_sha256)
        if (not self.source or not self.state_name or not isinstance(self.atoms, frozenset)
                or not self.atoms or any(not isinstance(a, str) or not a for a in self.atoms)
                or not isinstance(self.bonds, frozenset)
                or any(not isinstance(b, frozenset) or len(b) != 2 or not b <= self.atoms for b in self.bonds)
                or not isinstance(self.aliases, tuple)):
            raise PeptideStateError("immutable exact template inventory required")


@dataclass(frozen=True)
class PeptideSubstitution:
    carbon: ExactAtomRef
    nitrogen: ExactAtomRef
    removed: frozenset[ExactAtomRef]

    def __post_init__(self) -> None:
        _identity(self.carbon)
        _identity(self.nitrogen)
        if (self.carbon.name != "C" or self.nitrogen.name != "N" or self.carbon.component == self.nitrogen.component
                or not isinstance(self.removed, frozenset) or not self.removed):
            raise PeptideStateError("explicit distinct peptide C/N ports and removals required")
        for atom in self.removed:
            _identity(atom)


@dataclass(frozen=True)
class ChemicalWitness:
    source_units: tuple[tuple[ComponentInstanceRef, CcdGeometryAuthority], ...]
    substitutions: tuple[PeptideSubstitution, ...]
    projected: tuple[tuple[ComponentInstanceRef, ChemicalGraph], ...]
    flat_fragments: ChemicalGraph
    node_map: tuple[tuple[str, ExactAtomRef], ...]
    extra_bonds: tuple[tuple[str, str, float], ...]
    verified_cip: tuple[tuple[str, str], ...]

    @property
    def digest(self) -> str:
        return digest_record(("conditional-peptide-chemical-witness-v1", self))


def assemble_witness(units: tuple[tuple[ComponentInstanceRef, CcdGeometryAuthority], ...],
                     substitutions: tuple[PeptideSubstitution, ...]) -> ChemicalWitness:
    """Prove the joined chemistry; independent source frames are NOT a 3D peptide."""
    if not isinstance(units, tuple) or not 2 <= len(units) <= 16 or not isinstance(substitutions, tuple):
        raise PeptideStateError("bounded explicit units and substitutions required")
    by_id = dict(units)
    if len(by_id) != len(units):
        raise PeptideStateError("duplicate source component identity")
    atoms = {ExactAtomRef(identity, a.name): a for identity, source in units for a in source.state.graph.atoms}
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in atoms}
    for identity, source in units:
        if "PEPTIDE" not in source.component.component_type.upper():
            raise PeptideStateError("peptide source component type required")
        for bond in source.state.graph.bonds:
            a, b = ExactAtomRef(identity, bond.atom1), ExactAtomRef(identity, bond.atom2)
            adjacent[a].add(b)
            adjacent[b].add(a)
    removed: set[ExactAtomRef] = set()
    ports: set[ExactAtomRef] = set()
    for substitution in substitutions:
        c, n = substitution.carbon, substitution.nitrogen
        if c not in atoms or n not in atoms or atoms[c].element != "C" or atoms[n].element != "N":
            raise PeptideStateError("unknown or wrong-element peptide port")
        if {c, n} & ports:
            raise PeptideStateError("peptide port consumed twice")
        ports.update((c, n))
        if removed & substitution.removed or {c, n} & substitution.removed:
            raise PeptideStateError("overlapping or deleted substitution ports")
        for atom in substitution.removed:
            if atom.component not in {c.component, n.component} or atom not in atoms:
                raise PeptideStateError("removal outside named source components")
            if atom.name not in by_id[atom.component].leaving_atoms:
                raise PeptideStateError("removed atom is not a declared source leaving atom")
            endpoint = c if atom.component == c.component else n
            if not adjacent[atom] <= (set(substitution.removed) | {endpoint}):
                raise PeptideStateError("leaving branch crosses another retained atom")
        removed.update(substitution.removed)
    if ports & removed:
        raise PeptideStateError("a substitution removes another consumed port")
    reachable = {units[0][0]}
    while True:
        expanded = reachable | {s.nitrogen.component for s in substitutions if s.carbon.component in reachable} | {
            s.carbon.component for s in substitutions if s.nitrogen.component in reachable}
        if expanded == reachable:
            break
        reachable = expanded
    if reachable != set(by_id):
        raise PeptideStateError("chemical witness is not one explicitly connected assembly")

    projected = []
    flat_atoms, flat_bonds, mapping = [], [], []
    to_name: dict[ExactAtomRef, str] = {}
    plan_digest = digest_record((units, substitutions))
    for index, (identity, source) in enumerate(units):
        graph = source.state.graph
        names = {a.name for a in graph.atoms if ExactAtomRef(identity, a.name) not in removed}
        graph = replace(graph, state_id=f"peptide-{plan_digest[:16]}",
                        atoms=tuple(a for a in graph.atoms if a.name in names),
                        bonds=tuple(b for b in graph.bonds if {b.atom1, b.atom2} <= names),
                        net_charge=sum(a.formal_charge for a in graph.atoms if a.name in names))
        validate_graph(graph)
        projected.append((identity, graph))
        for chemical_atom in graph.atoms:
            ref = ExactAtomRef(identity, chemical_atom.name)
            name = f"unit{index}:{chemical_atom.name}"
            to_name[ref] = name
            mapping.append((name, ref))
            flat_atoms.append(replace(chemical_atom, name=name))
        for bond in graph.bonds:
            flat_bonds.append(replace(bond, atom1=to_name[ExactAtomRef(identity, bond.atom1)],
                                      atom2=to_name[ExactAtomRef(identity, bond.atom2)]))
    flat = ChemicalGraph("conditional-peptide-witness", plan_digest, tuple(flat_atoms), tuple(flat_bonds),
                         sum(a.formal_charge for a in flat_atoms), component_count=len(units))
    validate_graph(flat)
    extra = tuple((to_name[s.carbon], to_name[s.nitrogen], 1.) for s in substitutions)
    labels = verify_cip(flat, {a.name: a.ideal_position for a in flat.atoms}, extra_bonds=extra)
    return ChemicalWitness(units, substitutions, tuple(projected), flat, tuple(mapping), extra, labels)


@dataclass(frozen=True)
class LocalPeptideAuthority:
    parent: CcdGeometryAuthority
    inventory: TemplateInventory
    witness: ChemicalWitness
    identity: ComponentInstanceRef
    state: ChemicalState
    removed_atoms: frozenset[str]

    @property
    def digest(self) -> str:
        return digest_record(("local-peptide-state-v1", self))

    def constraints(self, *, minimum_volume: float, planar_tolerance_angstrom: float) -> tuple[Constraint, ...]:
        mapping = dict(self.witness.node_map)
        labels = tuple((mapping[name].name, label) for name, label in self.witness.verified_cip
                       if mapping[name].component == self.identity)
        present = {a.name for a in self.state.graph.atoms}
        planes = tuple(p for p in self.parent.planar_groups if set(p) <= present)
        constraints = ordered_constraints(self.state, labels, planes, self.identity, include_hydrogens=False,
                                          minimum_volume=minimum_volume, planar_tolerance_angstrom=planar_tolerance_angstrom)
        result = []
        for constraint in constraints:
            if (constraint.kind == ConstraintKind.STEREO_FRAME and constraint.atoms[0].name == "CA"
                    and tuple(a.name for a in constraint.atoms[1:]) == ("C", "CB", "N")):
                # (C,CB,N) -> (N,C,CB) is an even cyclic permutation, not an L/S inference.
                constraint = replace(constraint, atoms=tuple(ExactAtomRef(self.identity, n) for n in ("CA", "N", "C", "CB")))
            result.append(constraint)
        return tuple(result)


def derive_local_peptide(parent: CcdGeometryAuthority, alanine_witness: CcdGeometryAuthority,
                         identity: ComponentInstanceRef, inventory: TemplateInventory) -> LocalPeptideAuthority:
    if inventory.state_name != parent.component.component_id:
        raise PeptideStateError("explicit-variant-projection-unqualified")
    if alanine_witness.component.component_id != "ALA":
        raise PeptideStateError("this version declares neutral ALA witnesses only")
    left = ComponentInstanceRef(identity.model, "__left_chemical_witness__", "1")
    right = ComponentInstanceRef(identity.model, "__right_chemical_witness__", "1")
    if identity in {left, right}:
        raise PeptideStateError("source identity collides with chemical witness")
    units = ((left, alanine_witness), (identity, parent), (right, alanine_witness))
    substitutions = tuple(PeptideSubstitution(ExactAtomRef(a, "C"), ExactAtomRef(b, "N"),
                         frozenset((ExactAtomRef(a, "OXT"), ExactAtomRef(a, "HXT"), ExactAtomRef(b, "H2"))))
                          for a, b in ((left, identity), (identity, right)))
    witness = assemble_witness(units, substitutions)
    graph = dict(witness.projected)[identity]
    if inventory.atoms != {a.name for a in graph.atoms} or inventory.bonds != {
            frozenset((b.atom1, b.atom2)) for b in graph.bonds}:
        raise PeptideStateError("derived-state-template-inventory-mismatch")
    state = ChemicalState(graph, AuthorityProvenance(AuthorityMode.USER_MAPPED,
                          f"conditional local peptide; parent CCD {parent.component.provenance.sha256}; witness {witness.digest}",
                          canonical_graph_digest(graph), "local-peptide-state-v1"), Unit.ANGSTROM, HydrogenInventory.EXPLICIT)
    removed = frozenset(a.name for a in parent.state.graph.atoms if a.name not in inventory.atoms)
    return LocalPeptideAuthority(parent, inventory, witness, identity, state, removed)
