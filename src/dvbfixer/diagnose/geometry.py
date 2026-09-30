"""Pure, single-structure geometry measurements for diagnostics.

The primitives in this module use explicit ångström coordinates and never attach
severity or pass/fail policy.  OpenMM adaptation stays in the small helpers at
the bottom so callers can reuse the numeric contracts without importing an
experimental reconstruction backend.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

Point3 = tuple[float, float, float]
GEOMETRY_SCHEMA_VERSION = "dvbfixer.geometry.v1"


@dataclass(frozen=True, order=True)
class ResidueIdentity:
    """Exact, case-sensitive PDB residue identity."""

    chain: str
    resid: str
    icode: str = ""

    def __post_init__(self) -> None:
        if not self.resid:
            raise ValueError("residue number must not be empty")
        if len(self.icode) > 1:
            raise ValueError("residue insertion code must be empty or one character")

    @property
    def display_resid(self) -> str:
        return f"{self.resid}{self.icode}"

    def to_dict(self) -> dict[str, str]:
        return {"chain": self.chain, "resid": self.resid, "icode": self.icode}


@dataclass(frozen=True, order=True)
class AtomIdentity:
    """Exact, case-sensitive PDB atom identity."""

    chain: str
    resid: str
    icode: str
    atom: str

    def __post_init__(self) -> None:
        ResidueIdentity(self.chain, self.resid, self.icode)
        if not self.atom:
            raise ValueError("atom name must not be empty")

    @property
    def residue(self) -> ResidueIdentity:
        return ResidueIdentity(self.chain, self.resid, self.icode)

    def to_dict(self) -> dict[str, str]:
        return {
            "chain": self.chain,
            "resid": self.resid,
            "icode": self.icode,
            "atom": self.atom,
        }


@dataclass(frozen=True)
class GeometryMeasurement:
    """Typed numeric evidence with no diagnostic policy attached."""

    metric: str
    atoms: tuple[AtomIdentity, ...]
    value: float
    unit: str
    expected_min: float | None = None
    expected_max: float | None = None
    reference_value: float | None = None
    engine: str = "dvbfixer-python"
    engine_version: str = "1"
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.metric or not self.metric.rpartition(".v")[2].isdigit():
            raise ValueError("geometry metric must have a numeric .vN suffix")
        if not self.atoms:
            raise ValueError("geometry measurement must identify at least one atom")
        if not self.unit:
            raise ValueError("geometry measurement unit must not be empty")
        if not math.isfinite(self.value):
            raise ValueError("geometry measurement value must be finite")
        for value in (self.expected_min, self.expected_max, self.reference_value):
            if value is not None and not math.isfinite(value):
                raise ValueError("geometry measurement expectations must be finite")
        if (
            self.expected_min is not None
            and self.expected_max is not None
            and self.expected_min > self.expected_max
        ):
            raise ValueError("geometry measurement expected range is reversed")

    def to_dict(self) -> dict[str, object]:
        result: dict[str, object] = {
            "schema": GEOMETRY_SCHEMA_VERSION,
            "metric": self.metric,
            "atoms": [atom.to_dict() for atom in self.atoms],
            "value": self.value,
            "unit": self.unit,
            "engine": {"name": self.engine, "version": self.engine_version},
        }
        if self.expected_min is not None or self.expected_max is not None:
            result["expected_range"] = {
                "min": self.expected_min,
                "max": self.expected_max,
            }
        if self.reference_value is not None:
            result["reference_value"] = self.reference_value
        if self.metadata:
            result["metadata"] = dict(self.metadata)
        return result


@dataclass(frozen=True)
class BoundaryIdentity:
    """An explicitly requested ordered left/right residue boundary."""

    left: ResidueIdentity
    right: ResidueIdentity

    def to_dict(self) -> dict[str, object]:
        return {"left": self.left.to_dict(), "right": self.right.to_dict()}


@dataclass(frozen=True)
class DiagnosticScope:
    """Immutable exact residue selection and optional ordered boundaries."""

    residues: frozenset[ResidueIdentity] = frozenset()
    boundaries: tuple[BoundaryIdentity, ...] = ()

    @property
    def selected_residues(self) -> frozenset[ResidueIdentity]:
        boundary_residues = {
            residue for boundary in self.boundaries for residue in (boundary.left, boundary.right)
        }
        return self.residues | frozenset(boundary_residues)


@dataclass(frozen=True)
class BoundaryMeasurement:
    boundary: BoundaryIdentity
    c_n_distance: GeometryMeasurement | None
    topology_bonded: bool
    atom_match_counts: tuple[int, int]

    @property
    def status(self) -> str:
        if 0 in self.atom_match_counts:
            return "missing_atom"
        if self.atom_match_counts != (1, 1):
            return "ambiguous_atom"
        return "measured" if self.c_n_distance is not None else "undefined_geometry"

    def to_dict(self) -> dict[str, object]:
        return {
            **self.boundary.to_dict(),
            "status": self.status,
            "topology_bonded": self.topology_bonded,
            "c_n_distance": (
                self.c_n_distance.to_dict() if self.c_n_distance is not None else None
            ),
            "atoms": [
                {
                    **AtomIdentity(
                        self.boundary.left.chain,
                        self.boundary.left.resid,
                        self.boundary.left.icode,
                        "C",
                    ).to_dict(),
                    "match_count": self.atom_match_counts[0],
                },
                {
                    **AtomIdentity(
                        self.boundary.right.chain,
                        self.boundary.right.resid,
                        self.boundary.right.icode,
                        "N",
                    ).to_dict(),
                    "match_count": self.atom_match_counts[1],
                },
            ],
        }


def _finite_point(point: Sequence[float]) -> Point3 | None:
    if len(point) != 3:
        raise ValueError("geometry coordinates must have exactly three components")
    converted = (float(point[0]), float(point[1]), float(point[2]))
    return converted if all(math.isfinite(value) for value in converted) else None


def distance_angstrom(first: Sequence[float], second: Sequence[float]) -> float | None:
    """Euclidean distance for explicit ångström coordinates."""

    a = _finite_point(first)
    b = _finite_point(second)
    if a is None or b is None:
        return None
    return math.sqrt(sum((left - right) ** 2 for left, right in zip(a, b)))


def angle_degrees(
    first: Sequence[float], vertex: Sequence[float], third: Sequence[float]
) -> float | None:
    """Interior angle in degrees, or ``None`` for invalid/zero-length vectors."""

    a = _finite_point(first)
    b = _finite_point(vertex)
    c = _finite_point(third)
    if a is None or b is None or c is None:
        return None
    left = tuple(x - y for x, y in zip(a, b))
    right = tuple(x - y for x, y in zip(c, b))
    left_norm = math.sqrt(sum(value * value for value in left))
    right_norm = math.sqrt(sum(value * value for value in right))
    if left_norm <= 1e-12 or right_norm <= 1e-12:
        return None
    cosine = sum(x * y for x, y in zip(left, right)) / (left_norm * right_norm)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def signed_dihedral_degrees(
    first: Sequence[float],
    second: Sequence[float],
    third: Sequence[float],
    fourth: Sequence[float],
) -> float | None:
    """Signed four-point dihedral in ``[-180, 180]``, or ``None`` if degenerate."""

    p1 = _finite_point(first)
    p2 = _finite_point(second)
    p3 = _finite_point(third)
    p4 = _finite_point(fourth)
    if p1 is None or p2 is None or p3 is None or p4 is None:
        return None

    def sub(a: Point3, b: Point3) -> Point3:
        return a[0] - b[0], a[1] - b[1], a[2] - b[2]

    def dot(a: Point3, b: Point3) -> float:
        return sum(x * y for x, y in zip(a, b))

    def cross(a: Point3, b: Point3) -> Point3:
        return (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )

    central = sub(p3, p2)
    central_norm = math.sqrt(dot(central, central))
    if central_norm <= 1e-12:
        return None
    unit: Point3 = (
        central[0] / central_norm,
        central[1] / central_norm,
        central[2] / central_norm,
    )
    left = sub(p1, p2)
    right = sub(p4, p3)
    left_projection = dot(left, unit)
    right_projection = dot(right, unit)
    left = sub(left, tuple(left_projection * value for value in unit))
    right = sub(right, tuple(right_projection * value for value in unit))
    if math.sqrt(dot(left, left)) <= 1e-12 or math.sqrt(dot(right, right)) <= 1e-12:
        return None
    return math.degrees(math.atan2(dot(cross(unit, left), right), dot(left, right)))


def residue_identity(residue: Any) -> ResidueIdentity:
    return ResidueIdentity(
        chain=str(residue.chain.id),
        resid=str(residue.id),
        icode=str(getattr(residue, "insertionCode", "") or "").strip(),
    )


def atom_identity(atom: Any) -> AtomIdentity:
    residue = residue_identity(atom.residue)
    return AtomIdentity(residue.chain, residue.resid, residue.icode, str(atom.name))


def position_angstrom(positions: Any, atom: Any) -> Point3:
    """Adapt an OpenMM position to an explicit ångström coordinate."""

    from openmm.unit import angstrom

    point = positions[atom.index].value_in_unit(angstrom)
    return float(point[0]), float(point[1]), float(point[2])


def atom_index_by_identity(topology: Any) -> dict[AtomIdentity, list[Any]]:
    """Index atoms without discarding duplicates; ambiguous identities stay visible."""

    result: dict[AtomIdentity, list[Any]] = {}
    for atom in topology.atoms():
        result.setdefault(atom_identity(atom), []).append(atom)
    return result


def boundary_measurements(
    topology: Any,
    positions: Any,
    boundaries: Iterable[BoundaryIdentity],
) -> list[BoundaryMeasurement]:
    """Measure explicit C(left)-N(right) boundaries without inferring gaps."""

    atoms = atom_index_by_identity(topology)
    bonds = {frozenset((first.index, second.index)) for first, second in topology.bonds()}
    results: list[BoundaryMeasurement] = []
    for boundary in boundaries:
        left_id = AtomIdentity(boundary.left.chain, boundary.left.resid, boundary.left.icode, "C")
        right_id = AtomIdentity(
            boundary.right.chain, boundary.right.resid, boundary.right.icode, "N"
        )
        left_matches = atoms.get(left_id, [])
        right_matches = atoms.get(right_id, [])
        measurement: GeometryMeasurement | None = None
        bonded = False
        if len(left_matches) == 1 and len(right_matches) == 1:
            left_atom, right_atom = left_matches[0], right_matches[0]
            value = distance_angstrom(
                position_angstrom(positions, left_atom),
                position_angstrom(positions, right_atom),
            )
            bonded = frozenset((left_atom.index, right_atom.index)) in bonds
            if value is not None:
                measurement = GeometryMeasurement(
                    metric="peptide_boundary_c_n_distance.v1",
                    atoms=(left_id, right_id),
                    value=value,
                    unit="angstrom",
                    reference_value=1.33,
                    metadata={"topology_bonded": bonded},
                )
        results.append(
            BoundaryMeasurement(
                boundary,
                measurement,
                bonded,
                (len(left_matches), len(right_matches)),
            )
        )
    return results
