"""Independent Cartesian measurements of declared fixed constraints only.

Stereo-frame atoms are (centre, ordered neighbour 1, 2, 3). Four-ligand
atoms are (centre, ordered ligand 1, 2, 3, 4); unit ligand directions are
formed first, then the normalized determinant of (u1-u4, u2-u4, u3-u4).
Neither convention assigns CIP priorities or infers chemical stereochemistry.
"""
from __future__ import annotations

import math

from .contracts import Constraint, ConstraintKind, Unit

Point = tuple[float, float, float]


def _sub(a: Point, b: Point) -> Point:
    return a[0]-b[0], a[1]-b[1], a[2]-b[2]


def _dot(a: Point, b: Point) -> float:
    return sum(x*y for x, y in zip(a, b))


def _cross(a: Point, b: Point) -> Point:
    return a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]


def _unit(a: Point) -> Point:
    length = math.sqrt(_dot(a, a))
    if length < 1e-12:
        raise ValueError("degenerate fixed geometry")
    return a[0]/length, a[1]/length, a[2]/length


def angle_degrees(first: Point, centre: Point, last: Point) -> float:
    a, b = _unit(_sub(first, centre)), _unit(_sub(last, centre))
    return math.degrees(math.acos(max(-1., min(1., _dot(a, b)))))


def torsion_degrees(points: tuple[Point, Point, Point, Point]) -> float:
    b0 = _sub(points[0], points[1])
    b1 = _unit(_sub(points[2], points[1]))
    b2 = _sub(points[3], points[2])
    d0, d2 = _dot(b0, b1), _dot(b2, b1)
    v = _unit((b0[0]-d0*b1[0], b0[1]-d0*b1[1], b0[2]-d0*b1[2]))
    w = _unit((b2[0]-d2*b1[0], b2[1]-d2*b1[1], b2[2]-d2*b1[2]))
    return math.degrees(math.atan2(_dot(_cross(b1, v), w), _dot(v, w)))


def ordered_stereo_volume(kind: ConstraintKind, points: tuple[Point, ...]) -> float:
    """Canonical normalized ordered-frame/ligand determinant, without parity inference."""
    if kind not in {ConstraintKind.STEREO_FRAME, ConstraintKind.STEREO_LIGANDS}:
        raise ValueError("ordered stereo kind required")
    directions = tuple(_unit(_sub(p, points[0])) for p in points[1:])
    if kind == ConstraintKind.STEREO_LIGANDS:
        directions = tuple(_unit(_sub(p, directions[3])) for p in directions[:3])
    return _dot(directions[0], _cross(directions[1], directions[2]))


def satisfied(constraint: Constraint, points: tuple[Point, ...]) -> bool:
    """Positions arrive in angstrom; bounds retain their declared units."""
    kind = constraint.kind
    if kind in {ConstraintKind.STEREO_FRAME, ConstraintKind.STEREO_LIGANDS}:
        volume = ordered_stereo_volume(kind, points)
        assert constraint.orientation is not None and constraint.minimum_volume is not None
        return constraint.orientation*volume >= constraint.minimum_volume
    if kind in {ConstraintKind.DISTANCE, ConstraintKind.RING_CLOSURE, ConstraintKind.EXCLUSION}:
        delta = _sub(points[0], points[1])
        value = math.sqrt(_dot(delta, delta))
    elif kind == ConstraintKind.ANGLE:
        value = angle_degrees(points[0], points[1], points[2])
    elif kind == ConstraintKind.TORSION:
        value = torsion_degrees((points[0], points[1], points[2], points[3]))
    elif kind == ConstraintKind.PLANARITY:
        normal = _unit(_cross(_sub(points[1], points[0]), _sub(points[2], points[0])))
        value = max(abs(_dot(_sub(p, points[0]), normal)) for p in points[3:])
    else:
        raise ValueError("unknown measurement kind")
    if constraint.unit == Unit.NANOMETER:
        value /= 10
    assert constraint.upper is not None
    return (constraint.lower is None or value >= constraint.lower) and value <= constraint.upper
