"""Policy-free protein geometry evidence, with explicit identity and connectivity.

Coordinates are in angstroms. Peptide neighbors, disulfides, and terminal states
are caller declarations, never guessed from numbering or spatial proximity.
Biopython supplies canonical chi definitions; OpenMM's standard residue bond
templates supply canonical heavy connectivity. Neither implies a calibrated
geometry reference distribution. This module does not dispatch diagnose checks.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field
from functools import lru_cache
from importlib.resources import files
from itertools import combinations
from types import MappingProxyType
from typing import Any, Literal
from xml.etree import ElementTree

from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB.ic_data import ic_data_sidechains

from dvbfixer.diagnose.geometry import (
    VDW_RADII_ANGSTROM,
    AtomIdentity,
    BoundaryIdentity,
    GeometryMeasurement,
    Point3,
    ResidueIdentity,
    angle_degrees,
    atom_identity,
    distance_angstrom,
    position_angstrom,
    residue_identity,
    signed_dihedral_degrees,
)

RESULT_SCHEMA = "dvbfixer.protein_geometry.v1"
UndefinedReason = Literal[
    "missing_atom",
    "duplicate_atom",
    "missing_residue",
    "duplicate_residue",
    "degenerate_geometry",
    "non_finite_coordinates",
    "unsupported_residue",
    "undefined_terminal_neighbor",
    "ambiguous_peptide_neighbor",
    "unsupported_connectivity",
    "unsupported_reference_data",
    "unsupported_element",
    "undeclared_atom_role",
]


@dataclass(frozen=True)
class MeasurementResult:
    """Exactly one numeric measurement or explicit reason for absent evidence."""

    metric: str
    atoms: tuple[AtomIdentity, ...]
    measurement: GeometryMeasurement | None = None
    reason: UndefinedReason | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if (self.measurement is None) == (self.reason is None):
            raise ValueError("result requires exactly one measurement or undefined reason")
        if self.measurement is not None and (
            self.measurement.metric != self.metric or self.measurement.atoms != self.atoms
        ):
            raise ValueError("result and measurement identities must agree")
        if self.measurement is not None and dict(self.measurement.metadata) != dict(self.metadata):
            raise ValueError("result and measurement metadata must agree")

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": RESULT_SCHEMA,
            "metric": self.metric,
            "atoms": [atom.to_dict() for atom in self.atoms],
            "status": "measured" if self.measurement is not None else "undefined",
            "reason": self.reason,
            "measurement": self.measurement.to_dict() if self.measurement is not None else None,
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class ProteinAtom:
    identity: AtomIdentity
    position: Point3
    element: str | None


@dataclass(frozen=True)
class ProteinResidue:
    identity: ResidueIdentity
    name: str


@dataclass(frozen=True)
class ProteinStructure:
    """Immutable snapshot; tuples deliberately retain duplicate identities.

    Bonds are identity pairs. A bond involving a duplicated identity is never
    accepted as unambiguous evidence by the measurement functions.
    """

    residues: tuple[ProteinResidue, ...]
    atoms: tuple[ProteinAtom, ...]
    bonds: frozenset[frozenset[AtomIdentity]] = frozenset()

    @classmethod
    def from_openmm(cls, topology: Any, positions: Any) -> ProteinStructure:
        return cls(
            tuple(ProteinResidue(residue_identity(r), str(r.name)) for r in topology.residues()),
            tuple(
                ProteinAtom(
                    atom_identity(a),
                    position_angstrom(positions, a),
                    a.element.symbol if a.element is not None else None,
                )
                for a in topology.atoms()
            ),
            frozenset(frozenset((atom_identity(a), atom_identity(b))) for a, b in topology.bonds()),
        )


def _atom(residue: ResidueIdentity, name: str) -> AtomIdentity:
    return AtomIdentity(residue.chain, residue.resid, residue.icode, name)


class GeometryLookup:
    """Duplicate-safe lookup shared across measurements on a snapshot."""

    def __init__(self, structure: ProteinStructure):
        self.structure = structure
        self.atoms: dict[AtomIdentity, list[ProteinAtom]] = defaultdict(list)
        self.residues: dict[ResidueIdentity, list[ProteinResidue]] = defaultdict(list)
        for atom in structure.atoms:
            self.atoms[atom.identity].append(atom)
        for residue in structure.residues:
            self.residues[residue.identity].append(residue)

    def reason(self, atoms: tuple[AtomIdentity, ...]) -> UndefinedReason | None:
        for identity in atoms:
            count = len(self.residues.get(identity.residue, ()))
            if count == 0:
                return "missing_residue"
            if count != 1:
                return "duplicate_residue"
            count = len(self.atoms.get(identity, ()))
            if count == 0:
                return "missing_atom"
            if count != 1:
                return "duplicate_atom"
            if not all(math.isfinite(x) for x in self.atoms[identity][0].position):
                return "non_finite_coordinates"
        return None

    def measure(
        self,
        metric: str,
        atoms: tuple[AtomIdentity, ...],
        kind: Literal["distance", "angle", "dihedral"],
        metadata: Mapping[str, object] | None = None,
    ) -> MeasurementResult:
        details = dict(metadata or {})
        reason = self.reason(atoms)
        if reason is not None:
            return MeasurementResult(metric, atoms, reason=reason, metadata=details)
        points = [self.atoms[a][0].position for a in atoms]
        if kind == "distance":
            value = distance_angstrom(points[0], points[1])
        elif kind == "angle":
            value = angle_degrees(points[0], points[1], points[2])
        else:
            value = signed_dihedral_degrees(points[0], points[1], points[2], points[3])
        if value is None:
            return MeasurementResult(metric, atoms, reason="degenerate_geometry", metadata=details)
        return MeasurementResult(
            metric,
            atoms,
            GeometryMeasurement(
                metric,
                atoms,
                value,
                "angstrom" if kind == "distance" else "degree",
                metadata=details,
            ),
            metadata=details,
        )


def ramachandran_class(resname: str, following_resname: str | None = None) -> str | None:
    """GLY/PRO take precedence over pre-PRO; residue names are canonical PDB names."""
    if resname not in protein_letters_3to1:
        return None
    if resname == "GLY":
        return "gly"
    if resname == "PRO":
        return "pro"
    return "pre-pro" if following_resname == "PRO" else "general"


def _peptide_boundary_reason(
    lookup: GeometryLookup,
    boundary: BoundaryIdentity,
) -> UndefinedReason | None:
    atoms = (_atom(boundary.left, "C"), _atom(boundary.right, "N"))
    reason = lookup.reason(atoms)
    if reason is not None:
        return reason
    residues = [lookup.residues[identity][0] for identity in (boundary.left, boundary.right)]
    if (
        boundary.left == boundary.right
        or boundary.left.chain != boundary.right.chain
        or any(residue.name not in protein_letters_3to1 for residue in residues)
        or frozenset(atoms) not in lookup.structure.bonds
    ):
        return "unsupported_connectivity"
    return None


def backbone_torsions(
    structure: ProteinStructure,
    residue: ResidueIdentity,
    peptides: tuple[BoundaryIdentity, ...],
) -> tuple[MeasurementResult, MeasurementResult]:
    """Phi/psi from declared directed peptide neighbors, including absent endpoints."""
    lookup = GeometryLookup(structure)
    previous = [p.left for p in peptides if p.right == residue]
    following = [p.right for p in peptides if p.left == residue]
    current = lookup.residues.get(residue, [])
    next_residues = lookup.residues.get(following[0], []) if len(following) == 1 else []
    classification_reason: UndefinedReason | None = (
        _peptide_boundary_reason(lookup, BoundaryIdentity(residue, following[0]))
        if len(following) == 1
        else "ambiguous_peptide_neighbor"
        if len(following) > 1
        else None
    )
    cls = (
        ramachandran_class(
            current[0].name,
            next_residues[0].name if len(next_residues) == 1 else None,
        )
        if len(current) == 1 and classification_reason is None
        else None
    )
    metadata: dict[str, object] = {"ramachandran_class": cls}
    results: list[MeasurementResult] = []
    for label, neighbors in (("phi", previous), ("psi", following)):
        core = tuple(_atom(residue, n) for n in ("N", "CA", "C"))
        metric = f"backbone_{label}.v1"
        if len(neighbors) != 1:
            reason: UndefinedReason = (
                "undefined_terminal_neighbor" if not neighbors else "ambiguous_peptide_neighbor"
            )
            # Missing/duplicate current residues take precedence over endpoint absence.
            reason = lookup.reason(core) or reason
            results.append(MeasurementResult(metric, core, reason=reason, metadata=metadata))
            continue
        atoms = (
            (_atom(neighbors[0], "C"), *core)
            if label == "phi"
            else (*core, _atom(neighbors[0], "N"))
        )
        boundary = (
            BoundaryIdentity(neighbors[0], residue)
            if label == "phi"
            else BoundaryIdentity(residue, neighbors[0])
        )
        boundary_reason = _peptide_boundary_reason(lookup, boundary)
        if boundary_reason is not None:
            results.append(
                MeasurementResult(metric, atoms, reason=boundary_reason, metadata=metadata)
            )
        elif classification_reason is not None:
            results.append(
                MeasurementResult(
                    metric,
                    atoms,
                    reason=classification_reason,
                    metadata=metadata,
                )
            )
        elif cls is None and len(current) == 1:
            results.append(
                MeasurementResult(metric, atoms, reason="unsupported_residue", metadata=metadata)
            )
        else:
            results.append(lookup.measure(metric, atoms, "dihedral", metadata))
    return results[0], results[1]


def sidechain_torsions(
    structure: ProteinStructure,
    residue: ResidueIdentity,
) -> tuple[MeasurementResult, ...]:
    """Canonical chi1 through chi5, with no rotamer bins or probability claims."""
    lookup = GeometryLookup(structure)
    matches = lookup.residues.get(residue, [])
    if len(matches) != 1:
        reason: UndefinedReason = "missing_residue" if not matches else "duplicate_residue"
        return (
            MeasurementResult(
                "sidechain_chi.v1",
                (),
                reason=reason,
                metadata={"residue": residue.to_dict()},
            ),
        )
    name = matches[0].name
    if name not in protein_letters_3to1:
        return (
            MeasurementResult(
                "sidechain_chi.v1",
                (),
                reason="unsupported_residue",
                metadata={"residue": residue.to_dict(), "resname": name},
            ),
        )
    definitions = ic_data_sidechains.get(protein_letters_3to1[name], ())
    # ALA and GLY legitimately have no sidechain chi; not a missing measurement.
    return tuple(
        lookup.measure(
            f"sidechain_{definition[4]}.v1",
            tuple(_atom(residue, str(n)) for n in definition[:4]),
            "dihedral",
            {"resname": name, "definition_source": "Bio.PDB.ic_data"},
        )
        for definition in definitions
        if len(definition) == 5 and str(definition[4]).startswith("chi")
    )


def rotamer_reference_result(residue: ResidueIdentity, resname: str) -> MeasurementResult:
    """Pooled Janin occupancy is not a calibrated backbone-dependent probability."""
    return MeasurementResult(
        "backbone_dependent_rotamer_probability.v1",
        (),
        reason="unsupported_reference_data",
        metadata={
            "residue": residue.to_dict(),
            "resname": resname,
            "reference": None,
            "detail": "No calibrated residue/backbone-dependent library",
        },
    )


@lru_cache(maxsize=1)
def canonical_heavy_bonds() -> Mapping[str, tuple[tuple[str, str], ...]]:
    """Reuse OpenMM's canonical connectivity, without hydrogen/variant inference."""
    root = ElementTree.parse(str(files("openmm.app") / "data" / "residues.xml")).getroot()
    result: dict[str, tuple[tuple[str, str], ...]] = {}
    for residue in root.findall("Residue"):
        name = residue.attrib["name"]
        if name not in protein_letters_3to1:
            continue
        pairs = []
        for bond in residue.findall("Bond"):
            a, b = bond.attrib["from"], bond.attrib["to"]
            if any(n.startswith(("-", "+", "H")) for n in (a, b)) or "OXT" in (a, b):
                continue
            pairs.append((a, b))
        result[name] = tuple(sorted(pairs))
    return MappingProxyType(result)


def connectivity_measurements(
    structure: ProteinStructure,
    *,
    peptides: tuple[BoundaryIdentity, ...] = (),
    junctions: tuple[BoundaryIdentity, ...] = (),
    disulfides: tuple[BoundaryIdentity, ...] = (),
    n_termini: frozenset[ResidueIdentity] = frozenset(),
    c_termini: frozenset[ResidueIdentity] = frozenset(),
) -> tuple[MeasurementResult, ...]:
    """Expected heavy bonds + lengths + all angles about expected bond centers.

    Bond presence is numeric 0/1 evidence, not pass/fail. C-terminal OXT is
    required only when declared. Disulfides must be declared, not distance-inferred.
    Unknown intra-residue chemistry and non-CYS disulfides are explicit undefined.
    """
    lookup = GeometryLookup(structure)
    edges: dict[frozenset[AtomIdentity], tuple[str, str]] = {}
    results: list[MeasurementResult] = []
    for terminal_kind, declared_terminals, atom_name in (
        ("n_terminal", n_termini, "N"),
        ("c_terminal", c_termini, "OXT"),
    ):
        for identity in sorted(declared_terminals):
            terminal_atoms = (_atom(identity, atom_name),)
            terminal_metadata: dict[str, object] = {
                "terminal_kind": terminal_kind,
                "residue": identity.to_dict(),
            }
            terminal_reason = lookup.reason(terminal_atoms)
            terminal_matches = lookup.residues.get(identity, [])
            if terminal_reason is None:
                assert len(terminal_matches) == 1
                if terminal_matches[0].name not in protein_letters_3to1:
                    terminal_reason = "unsupported_residue"
            if terminal_reason is not None:
                results.append(
                    MeasurementResult(
                        "terminal_atom_presence.v1",
                        terminal_atoms,
                        reason=terminal_reason,
                        metadata=terminal_metadata,
                    )
                )
            else:
                results.append(
                    MeasurementResult(
                        "terminal_atom_presence.v1",
                        terminal_atoms,
                        GeometryMeasurement(
                            "terminal_atom_presence.v1",
                            terminal_atoms,
                            1.0,
                            "boolean",
                            metadata=terminal_metadata,
                        ),
                        metadata=terminal_metadata,
                    )
                )
    for residue in structure.residues:
        pairs = canonical_heavy_bonds().get(residue.name)
        if pairs is None:
            results.append(
                MeasurementResult(
                    "canonical_connectivity.v1",
                    (),
                    reason="unsupported_residue",
                    metadata={"residue": residue.identity.to_dict(), "resname": residue.name},
                )
            )
            continue
        for a, b in (*pairs, *((("C", "OXT"),) if residue.identity in c_termini else ())):
            edge = frozenset((_atom(residue.identity, a), _atom(residue.identity, b)))
            edges[edge] = (
                "intra_residue",
                "declared_terminal" if "OXT" in (a, b) else "openmm_residue_template",
            )
    for kind, boundaries, names in (
        ("peptide", peptides, ("C", "N")),
        ("junction", junctions, ("C", "N")),
        ("disulfide", disulfides, ("SG", "SG")),
    ):
        for boundary in boundaries:
            boundary_atoms = (_atom(boundary.left, names[0]), _atom(boundary.right, names[1]))
            endpoint_matches = [
                lookup.residues.get(identity, []) for identity in (boundary.left, boundary.right)
            ]
            supported = all(
                len(matches) == 1 and matches[0].name in protein_letters_3to1
                for matches in endpoint_matches
            )
            if kind == "disulfide":
                supported = all(
                    len(matches) == 1 and matches[0].name in {"CYS", "CYM", "CYX", "CSS"}
                    for matches in endpoint_matches
                )
            if (
                not supported
                or boundary.left == boundary.right
                or (kind in {"peptide", "junction"} and boundary.left.chain != boundary.right.chain)
            ):
                results.append(
                    MeasurementResult(
                        "canonical_connectivity.v1",
                        boundary_atoms,
                        reason=lookup.reason(boundary_atoms) or "unsupported_connectivity",
                        metadata={
                            "bond_kind": kind,
                            "definition_source": f"declared_{kind}",
                        },
                    )
                )
                continue
            edges[frozenset(boundary_atoms)] = (kind, f"declared_{kind}")
    neighbors: dict[AtomIdentity, set[AtomIdentity]] = defaultdict(set)
    for edge, (kind, definition_source) in sorted(
        edges.items(), key=lambda item: tuple(sorted(item[0]))
    ):
        a, b = sorted(edge)
        bond_atoms = (a, b)
        residue_names = tuple(
            lookup.residues[identity.residue][0].name
            for identity in bond_atoms
            if len(lookup.residues.get(identity.residue, ())) == 1
        )
        bond_metadata: dict[str, object] = {
            "bond_kind": kind,
            "residue_names": residue_names,
            "definition_source": definition_source,
        }
        bond_reason = lookup.reason(bond_atoms)
        if bond_reason is not None:
            results.append(
                MeasurementResult(
                    "canonical_connectivity.v1",
                    bond_atoms,
                    reason=bond_reason,
                    metadata=bond_metadata,
                )
            )
        else:
            results.append(
                MeasurementResult(
                    "canonical_connectivity.v1",
                    bond_atoms,
                    GeometryMeasurement(
                        "canonical_connectivity.v1",
                        bond_atoms,
                        float(edge in structure.bonds),
                        "boolean",
                        metadata=bond_metadata,
                    ),
                    metadata=bond_metadata,
                )
            )
        results.append(lookup.measure("bond_length.v1", bond_atoms, "distance", bond_metadata))
        neighbors[a].add(b)
        neighbors[b].add(a)
    for center in sorted(neighbors):
        for a, b in combinations(sorted(neighbors[center]), 2):
            angle_atoms = (a, center, b)
            angle_residue_names = tuple(
                lookup.residues[identity.residue][0].name
                for identity in angle_atoms
                if len(lookup.residues.get(identity.residue, ())) == 1
            )
            definition_sources = tuple(
                sorted({edges[frozenset((center, endpoint))][1] for endpoint in (a, b)})
            )
            results.append(
                lookup.measure(
                    "bond_angle.v1",
                    angle_atoms,
                    "angle",
                    {
                        "residue_names": angle_residue_names,
                        "definition_sources": definition_sources,
                    },
                )
            )
    return tuple(results)


def peptide_geometry(
    structure: ProteinStructure,
    boundary: BoundaryIdentity,
    *,
    boundary_kind: Literal["peptide", "junction"] = "peptide",
) -> tuple[MeasurementResult, ...]:
    """Omega and carbonyl improper with policy-free planarity distances."""
    lookup = GeometryLookup(structure)

    def left_atom(name: str) -> AtomIdentity:
        return _atom(boundary.left, name)

    def right_atom(name: str) -> AtomIdentity:
        return _atom(boundary.right, name)

    atom_sets = (
        ("peptide_omega.v1", (left_atom("CA"), left_atom("C"), right_atom("N"), right_atom("CA"))),
        (
            "peptide_carbonyl_improper.v1",
            (left_atom("O"), left_atom("C"), left_atom("CA"), right_atom("N")),
        ),
    )
    boundary_reason = _peptide_boundary_reason(lookup, boundary)
    boundary_metadata: dict[str, object] = {
        "boundary_kind": boundary_kind,
        "boundary": boundary.to_dict(),
    }
    results = [
        MeasurementResult(
            metric,
            atoms,
            reason=boundary_reason,
            metadata=boundary_metadata,
        )
        if boundary_reason is not None
        else lookup.measure(metric, atoms, "dihedral", boundary_metadata)
        for metric, atoms in atom_sets
    ]
    for result in tuple(results):
        if result.measurement is None:
            results.append(
                MeasurementResult(
                    result.metric.replace(".v1", "_planarity.v1"),
                    result.atoms,
                    reason=result.reason,
                    metadata=result.metadata,
                )
            )
            continue
        value = abs(result.measurement.value)
        metric = result.metric.replace(".v1", "_planarity.v1")
        metadata = {
            **result.metadata,
            "distance_to_nearest_plane_degrees": min(value, 180.0 - value),
        }
        if result.metric == "peptide_omega.v1":
            metadata.update(
                distance_to_cis_degrees=value,
                distance_to_trans_degrees=180.0 - value,
            )
        results.append(
            MeasurementResult(
                metric,
                result.atoms,
                GeometryMeasurement(
                    metric, result.atoms, min(value, 180.0 - value), "degree", metadata=metadata
                ),
                metadata=metadata,
            )
        )
    return tuple(results)


def steric_measurements(
    structure: ProteinStructure,
    *,
    generated: frozenset[AtomIdentity],
    fixed: frozenset[AtomIdentity],
    junction: frozenset[AtomIdentity] = frozenset(),
    include_hydrogens: bool = False,
    cutoff_angstrom: float | None = None,
    exclude_bond_hops: int = 3,
) -> tuple[MeasurementResult, ...]:
    """Raw vdW overlap for declared roles; no covalent distance or H-bond inference.

    Junction takes its own disjoint role, and pairs touching it are tagged junction.
    Fixed/fixed pairs are outside this selection. Only declared topology bonds
    define exclusions (1-2/1-3/1-4 by default), so close nonbonded atoms stay visible.
    Unsupported/missing/duplicate atoms yield undefined evidence before search.
    """
    from scipy.spatial import cKDTree

    if (
        cutoff_angstrom is not None and (not math.isfinite(cutoff_angstrom) or cutoff_angstrom <= 0)
    ) or exclude_bond_hops < 0:
        raise ValueError("positive finite cutoff and nonnegative exclusion hops required")
    if generated & fixed or generated & junction or fixed & junction:
        raise ValueError("atom roles must be disjoint")
    lookup = GeometryLookup(structure)
    roles = {
        a: role
        for role, identities in (("generated", generated), ("fixed", fixed), ("junction", junction))
        for a in identities
    }
    results: list[MeasurementResult] = []
    selected: list[ProteinAtom] = []
    for identity in sorted(roles):
        reason = lookup.reason((identity,))
        if reason is None:
            atom = lookup.atoms[identity][0]
            if atom.element is not None and atom.element.upper() == "H" and not include_hydrogens:
                continue
            if atom.element is None or atom.element.upper() not in VDW_RADII_ANGSTROM:
                reason = "unsupported_element"
            else:
                selected.append(atom)
        if reason is not None:
            results.append(
                MeasurementResult(
                    "steric_overlap.v1",
                    (identity,),
                    reason=reason,
                    metadata={"role": roles[identity]},
                )
            )
    neighbors: dict[AtomIdentity, set[AtomIdentity]] = defaultdict(set)
    for edge in structure.bonds:
        if len(edge) != 2:
            continue
        a, b = sorted(edge)
        if lookup.reason((a, b)) is None:
            neighbors[a].add(b)
            neighbors[b].add(a)
    excluded: set[frozenset[AtomIdentity]] = set()
    for atom in selected:
        visited = {atom.identity}
        frontier = {atom.identity}
        for _ in range(exclude_bond_hops):
            frontier = {b for a in frontier for b in neighbors[a]} - visited
            visited.update(frontier)
        excluded.update(frozenset((atom.identity, b)) for b in visited if b != atom.identity)
    if len(selected) < 2:
        return tuple(results)
    search_cutoff = (
        cutoff_angstrom
        if cutoff_angstrom is not None
        else 2.0
        * max(VDW_RADII_ANGSTROM[atom.element.upper()] for atom in selected if atom.element)
    )
    tree = cKDTree([a.position for a in selected])
    for i, j in sorted(tree.query_pairs(search_cutoff)):
        a, b = selected[i], selected[j]
        atoms = (a.identity, b.identity)
        pair_roles = (roles[a.identity], roles[b.identity])
        if pair_roles == ("fixed", "fixed") or frozenset(atoms) in excluded:
            continue
        category = (
            "junction"
            if "junction" in pair_roles
            else "generated/generated"
            if pair_roles == ("generated", "generated")
            else "generated/fixed"
        )
        distance = distance_angstrom(a.position, b.position)
        assert distance is not None and a.element is not None and b.element is not None
        radii = VDW_RADII_ANGSTROM[a.element.upper()] + VDW_RADII_ANGSTROM[b.element.upper()]
        metadata = {
            "pair_class": category,
            "roles": pair_roles,
            "hydrogen_pair": "H" in (a.element, b.element),
            "distance_angstrom": distance,
            "vdw_sum_angstrom": radii,
            "exclude_bond_hops": exclude_bond_hops,
            "cutoff_angstrom": search_cutoff,
            "radius_source": "dvbfixer.diagnose.geometry.VDW_RADII_ANGSTROM",
        }
        results.append(
            MeasurementResult(
                "steric_overlap.v1",
                atoms,
                GeometryMeasurement(
                    "steric_overlap.v1", atoms, radii - distance, "angstrom", metadata=metadata
                ),
                metadata=metadata,
            )
        )
    return tuple(results)
