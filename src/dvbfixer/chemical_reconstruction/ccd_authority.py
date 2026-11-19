"""Explicit pinned-local CCD all-atom authority; no fetch or polymer projection."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

from dvbfixer.ccd import CcdComponent, ccd_to_all_atom_graph, load_ccd_component
from dvbfixer.domain.nonprotein_reconstruction import (
    AuthorityMode,
    AuthorityProvenance,
    ChemicalGraph,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

from .contracts import (
    ChemicalState,
    Constraint,
    ConstraintKind,
    ContractError,
    HydrogenInventory,
    Unit,
    _point,
    _sha,
    digest_record,
)
from .measurement import Point, _cross, _dot, _sub, _unit


class CcdAuthorityError(ContractError):
    """Unusable chemical authority, not permission to guess another state."""


def _molecule(graph: ChemicalGraph, positions: dict[str, Point]) -> Any:
    from rdkit import Chem

    builder = Chem.RWMol()
    indices: dict[str, int] = {}
    for atom in graph.atoms:
        if atom.element.upper() == "D":
            raise CcdAuthorityError("isotope authority is not represented by the shared atom DTO")
        item = Chem.Atom(atom.element.title())
        if item.GetAtomicNum() == 0:
            raise CcdAuthorityError("unknown atomic number")
        item.SetFormalCharge(atom.formal_charge)
        item.SetNoImplicit(True)
        item.SetIsAromatic(atom.aromatic)
        indices[atom.name] = builder.AddAtom(item)
    types = {1.: Chem.BondType.SINGLE, 2.: Chem.BondType.DOUBLE, 3.: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}
    for bond in graph.bonds:
        if bond.order not in types:
            raise CcdAuthorityError("unsupported RDKit bond order")
        builder.AddBond(indices[bond.atom1], indices[bond.atom2], types[bond.order])
        builder.GetBondBetweenAtoms(indices[bond.atom1], indices[bond.atom2]).SetIsAromatic(bond.aromatic)
    molecule = builder.GetMol()
    try:
        Chem.SanitizeMol(molecule)
    except (ValueError, RuntimeError) as exc:
        raise CcdAuthorityError(f"CCD valence/aromaticity failed RDKit sanitization: {exc}") from exc
    for index, atom in enumerate(graph.atoms):
        sanitized = molecule.GetAtomWithIdx(index)
        if sanitized.GetFormalCharge() != atom.formal_charge or sanitized.GetNumImplicitHs():
            raise CcdAuthorityError("RDKit cleanup changed the declared charge or implicit-H inventory")
    if any(a.GetNumRadicalElectrons() for a in molecule.GetAtoms()):
        raise CcdAuthorityError("radical or incomplete explicit-H state needs a separate authority policy")
    conformer = Chem.Conformer(len(graph.atoms))
    conformer.Set3D(True)
    for atom in graph.atoms:
        conformer.SetAtomPosition(indices[atom.name], positions[atom.name])
    molecule.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(molecule, confId=0, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return molecule


def verify_cip(graph: ChemicalGraph, positions: dict[str, Point]) -> tuple[tuple[str, str], ...]:
    """Independent of the coordinate constructor; no coordinate repair or H addition."""
    if set(positions) != {a.name for a in graph.atoms}:
        raise CcdAuthorityError("CIP verification requires the complete exact all-atom graph")
    for point in positions.values():
        _point(point)
    molecule = _molecule(graph, positions)
    labels = []
    for index, atom in enumerate(graph.atoms):
        assigned = molecule.GetAtomWithIdx(index)
        label = assigned.GetProp("_CIPCode") if assigned.HasProp("_CIPCode") else None
        if label != atom.stereo:
            raise CcdAuthorityError(f"CCD/RDKit CIP mismatch at {atom.name}: declared={atom.stereo}, measured={label}")
        if label is not None:
            labels.append((atom.name, label))
    return tuple(labels)


@dataclass(frozen=True)
class CcdGeometryAuthority:
    component: CcdComponent
    state: ChemicalState
    rdkit_version: str
    verified_cip: tuple[tuple[str, str], ...]
    planar_groups: tuple[tuple[str, ...], ...]

    @property
    def digest(self) -> str:
        return digest_record(("ccd-all-atom-geometry-v1", self))

    @property
    def leaving_atoms(self) -> tuple[str, ...]:
        return tuple(a.name for a in self.component.atoms if a.leaving)

    def constraints(self, identity: ComponentInstanceRef, *, include_hydrogens: bool,
                    minimum_volume: float, planar_tolerance_angstrom: float) -> tuple[Constraint, ...]:
        """Map CIP-verified ideal geometry to explicit ordered measurement frames."""
        graph = self.state.graph
        atoms = {a.name: a for a in graph.atoms}
        positions = {a.name: a.ideal_position for a in graph.atoms}
        adjacent: dict[str, set[str]] = {a.name: set() for a in graph.atoms}
        for bond in graph.bonds:
            adjacent[bond.atom1].add(bond.atom2)
            adjacent[bond.atom2].add(bond.atom1)
        result = []
        for centre, label in self.verified_cip:
            neighbours = tuple(sorted(n for n in adjacent[centre] if include_hydrogens or atoms[n].element != "H"))
            if len(neighbours) not in {3, 4}:
                raise CcdAuthorityError("stereocentre needs unavailable H or unsupported ligand inventory")
            directions = tuple(_unit(_sub(positions[n], positions[centre])) for n in neighbours)
            kind = ConstraintKind.STEREO_FRAME
            if len(directions) == 4:
                kind = ConstraintKind.STEREO_LIGANDS
                directions = tuple(_unit(_sub(d, directions[3])) for d in directions[:3])
            volume = _dot(directions[0], _cross(directions[1], directions[2]))
            if abs(volume) < minimum_volume:
                raise CcdAuthorityError("ideal stereo frame is degenerate under the declared bound")
            result.append(Constraint(kind, tuple(ExactAtomRef(identity, n) for n in (centre, *neighbours)),
                                     Unit.DIMENSIONLESS, orientation=1 if volume > 0 else -1,
                                     minimum_volume=minimum_volume, label=label))
        for group in self.planar_groups:
            if not include_hydrogens and any(atoms[n].element == "H" for n in group):
                continue
            result.append(Constraint(ConstraintKind.PLANARITY, tuple(ExactAtomRef(identity, n) for n in group),
                                     Unit.ANGSTROM, upper=planar_tolerance_angstrom))
        return tuple(result)


def load_geometry_authority(path: Path, component_id: str, *, expected_sha256: str) -> CcdGeometryAuthority:
    """Explicit local read through the existing pin-checked loader; no network."""
    _sha(expected_sha256)  # The older information loader also supports unpinned reads; this route does not.
    component = load_ccd_component(path, component_id, expected_sha256=expected_sha256)
    from rdkit import rdBase

    if component.provisional_class != "A":
        raise CcdAuthorityError("initial all-atom authority supports unambiguous organic CCD chemistry only")
    if any(b.stereo is not None for b in component.bonds):
        raise CcdAuthorityError("bond stereo is not represented by the current shared graph DTO")
    if any(a.stereo not in {None, "R", "S"} for a in component.atoms):
        raise CcdAuthorityError("unsupported CCD atom stereo label")
    graph = ccd_to_all_atom_graph(component)
    state = ChemicalState(graph, AuthorityProvenance(AuthorityMode.LOCAL_PINNED_CCD, component.provenance.source,
                          component.provenance.sha256, component.modified_date), Unit.ANGSTROM, HydrogenInventory.EXPLICIT)
    positions = {a.name: a.ideal_position for a in graph.atoms}
    labels = verify_cip(graph, positions)
    molecule = _molecule(graph, positions)
    names = tuple(a.name for a in graph.atoms)
    groups: set[tuple[str, ...]] = set()
    for atom in molecule.GetAtoms():
        if str(atom.GetHybridization()) == "SP2" and atom.GetDegree() == 3:
            groups.add(tuple(sorted((names[atom.GetIdx()], *(names[n.GetIdx()] for n in atom.GetNeighbors())))))
    for ring in molecule.GetRingInfo().AtomRings():
        if len(ring) >= 4 and all(molecule.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            groups.add(tuple(sorted(names[i] for i in ring)))
    ordered = []
    for group in sorted(groups):
        for anchors in combinations(group, 3):
            normal = _cross(_sub(positions[anchors[1]], positions[anchors[0]]),
                            _sub(positions[anchors[2]], positions[anchors[0]]))
            if _dot(normal, normal) > 1e-12:
                ordered.append((*anchors, *(n for n in group if n not in anchors)))
                break
        else:
            raise CcdAuthorityError("no noncollinear authority plane anchors")
    return CcdGeometryAuthority(component, state, rdBase.rdkitVersion, labels, tuple(ordered))
