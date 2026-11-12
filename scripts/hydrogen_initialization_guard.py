"""Private pure proposal only: no production import, OpenMM hook or RNG use."""
from __future__ import annotations

import math
from dataclasses import dataclass

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS

Point = tuple[float, float, float]
REVISION = 'private-selective-carbon-h-v1'


class GuardRefusal(ValueError):
    """Invalid/ambiguous input or no admitted atomic proposal; no partial output."""


@dataclass(frozen=True)
class Atom:
    identity: ExactAtomRef
    residue_name: str
    element: str


@dataclass(frozen=True)
class Input:
    atoms: tuple[Atom, ...]
    positions_nm: tuple[Point, ...]
    bonds: tuple[tuple[int, int], ...]
    new_hydrogens: frozenset[ExactAtomRef]


@dataclass(frozen=True)
class Result:
    positions_nm: tuple[Point, ...]
    moved: tuple[ExactAtomRef, ...]
    supported: tuple[ExactAtomRef, ...]
    unsupported: tuple[ExactAtomRef, ...]


def sub(a: Point, b: Point) -> Point:
    return a[0]-b[0], a[1]-b[1], a[2]-b[2]


def norm(a: Point) -> float:
    return math.sqrt(sum(x*x for x in a))


def normalized_volume(points: tuple[Point, Point, Point], origin: Point) -> float:
    a, b, c = (sub(p, origin) for p in points)
    divisor = norm(a)*norm(b)*norm(c)
    if divisor <= 1e-30:
        raise GuardRefusal('degenerate tetrahedron vectors')
    det = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
           + a[2]*(b[0]*c[1]-b[1]*c[0]))
    return det/divisor


def validate(value: Input) -> list[list[int]]:
    if (not isinstance(value, Input) or not isinstance(value.atoms, tuple)
            or not isinstance(value.positions_nm, tuple) or not isinstance(value.bonds, tuple)
            or not isinstance(value.new_hydrogens, frozenset)
            or not 0 < len(value.atoms) <= 50_000 or len(value.atoms) != len(value.positions_nm)):
        raise GuardRefusal('bounded immutable input required')
    identities = set()
    names = {}
    for atom, point in zip(value.atoms, value.positions_nm, strict=True):
        if (not isinstance(atom, Atom) or not isinstance(atom.identity, ExactAtomRef)
                or not isinstance(atom.residue_name, str) or not atom.residue_name
                or not isinstance(atom.element, str) or not atom.element
                or not isinstance(point, tuple) or len(point) != 3
                or any(isinstance(x, bool) or not isinstance(x, (float, int))
                       or not math.isfinite(x) or abs(x) > 1e4 for x in point)):
            raise GuardRefusal('finite atom/coordinate record required')
        if atom.identity in identities:
            raise GuardRefusal('duplicate exact atom identity')
        identities.add(atom.identity)
        component = atom.identity.component
        if component in names and names[component] != atom.residue_name:
            raise GuardRefusal('inconsistent component residue names')
        names[component] = atom.residue_name
    if not value.new_hydrogens <= identities:
        raise GuardRefusal('new-H mask includes absent identities')
    adjacent: list[list[int]] = [[] for _ in value.atoms]
    seen = set()
    for bond in value.bonds:
        if (not isinstance(bond, tuple) or len(bond) != 2
                or any(type(i) is not int or not 0 <= i < len(value.atoms) for i in bond)):
            raise GuardRefusal('invalid bond endpoints')
        a, b = bond
        key = min(a, b), max(a, b)
        if a == b or key in seen:
            raise GuardRefusal('duplicate/self bond')
        seen.add(key)
        adjacent[a].append(b)
        adjacent[b].append(a)
    for i, atom in enumerate(value.atoms):
        if atom.identity in value.new_hydrogens:
            if (atom.element != 'H' or len(adjacent[i]) != 1
                    or value.atoms[adjacent[i][0]].element in ('H', 'D')):
                raise GuardRefusal('new H must have exactly one heavy parent')
    return adjacent


def propose(value: Input) -> Result:
    """Return all proposals together or refuse, never alter an observation."""
    adjacent = validate(value)
    points = list(value.positions_nm)
    moved, supported, unsupported = [], [], []
    parent_by_moved = {}
    for h, atom in enumerate(value.atoms):
        if atom.identity not in value.new_hydrogens:
            continue
        parent = adjacent[h][0]
        center = value.atoms[parent]
        canonical = VARIANT_TO_PARENT.get(center.residue_name, center.residue_name)
        heavy = [i for i in adjacent[parent] if value.atoms[i].element not in ('H', 'D')]
        hydrogens = [i for i in adjacent[parent] if value.atoms[i].element in ('H', 'D')]
        local = [parent, h]+heavy
        inventory = CANONICAL_HEAVY_ATOMS.get(canonical, ())
        if (center.element != 'C' or len(heavy) != 3 or hydrogens != [h]
                or center.identity.name not in inventory
                or any(value.atoms[i].identity.component != center.identity.component for i in local)
                or any(value.atoms[i].identity.name not in inventory for i in heavy)):
            unsupported.append(atom.identity)
            continue
        supported.append(atom.identity)
        heavy.sort(key=lambda i: value.atoms[i].identity)
        p = value.positions_nm[parent]
        neighbors = tuple(value.positions_nm[i] for i in heavy)
        radius = norm(sub(value.positions_nm[h], p))
        if not 0.05 <= radius <= 0.15 or any(not 0.10 <= norm(sub(q, p)) <= 0.20 for q in neighbors):
            raise GuardRefusal('supported center outside distance admission')
        hv = normalized_volume(neighbors, p)
        fv = normalized_volume(neighbors, value.positions_nm[h])
        if abs(hv) <= 1e-4 or abs(fv) <= 1e-4:
            raise GuardRefusal('degenerate supported center')
        if hv*fv > 0:
            continue
        raw = tuple(sum(p[k]-q[k] for q in neighbors) for k in range(3))
        length = norm(raw)
        if length <= 1e-12:
            raise GuardRefusal('degenerate heavy-neighbor sum')
        candidate = tuple(p[k]+radius*raw[k]/length for k in range(3))
        cv = normalized_volume(neighbors, candidate)
        if abs(cv) <= 1e-4 or hv*cv <= 0:
            raise GuardRefusal('no sign-compatible admitted candidate')
        points[h] = candidate
        moved.append(atom.identity)
        parent_by_moved[h] = parent
    for h, parent in parent_by_moved.items():
        if any(math.dist(points[h], point) < 0.05 for i, point in enumerate(points) if i not in (h, parent)):
            raise GuardRefusal('proposal overlaps a nonparent atom')
    return Result(tuple(points), tuple(moved), tuple(supported), tuple(unsupported))
