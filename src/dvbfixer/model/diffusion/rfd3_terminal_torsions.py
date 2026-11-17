"""Private missing-terminal kinematics; no source motion or public capability.

Use the actual, independently checked topology, not residue-number adjacency.
Only a graph component disconnected by one bond may rotate. A ring or a component
containing any unapproved atom fails closed. This grants no geometry acceptance.
"""
from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, Pair


def terminal_component(
    bonds: tuple[Pair, ...], *, axis: Pair, seed: AtomIdentity,
    allowed: frozenset[AtomIdentity],
) -> frozenset[AtomIdentity]:
    """Return one side of a rotatable bond only if every atom is approved missing."""
    if len(axis) != 2 or axis[0] == axis[1] or not allowed:
        raise CustodyError("invalid terminal torsion axis/mask")
    graph = defaultdict(set)
    seen = set()
    for first, second in bonds:
        pair = tuple(sorted((first, second)))
        if first == second or pair in seen:
            raise CustodyError("ambiguous terminal bond graph")
        seen.add(pair)
        graph[first].add(second)
        graph[second].add(first)
    if tuple(sorted(axis)) not in seen or seed not in graph:
        raise CustodyError("terminal torsion axis/seed absent from native graph")
    graph[axis[0]].remove(axis[1])
    graph[axis[1]].remove(axis[0])
    visited, todo = set(), [seed]
    while todo:
        atom = todo.pop()
        if atom in visited:
            continue
        visited.add(atom)
        todo.extend(graph[atom] - visited)
    if set(axis) <= visited or not set(axis) & visited:
        raise CustodyError("terminal torsion bond is cyclic or disconnected from seed")
    if not visited <= allowed:
        raise CustodyError("terminal torsion would move an unapproved atom")
    return frozenset(visited)


def rotate_terminal_component(
    coordinates: Mapping[AtomIdentity, tuple[float, float, float]],
    *, component: frozenset[AtomIdentity], axis: Pair, degrees: float,
) -> dict[AtomIdentity, tuple[float, float, float]]:
    """Rigid proper rotation; untouched coordinate values remain exactly unchanged."""
    if (isinstance(degrees, bool) or not math.isfinite(degrees) or not component
            or not component < coordinates.keys() or not set(axis) <= coordinates.keys()):
        raise CustodyError("invalid terminal rotation angle/component/axis")
    if any(len(xyz) != 3 or not all(math.isfinite(v) for v in xyz) for xyz in coordinates.values()):
        raise CustodyError("terminal rotation requires finite three-vectors")
    origin, last = (np.asarray(coordinates[atom], dtype=float) for atom in axis)
    direction = last - origin
    norm = float(np.linalg.norm(direction))
    if norm <= 1e-8:
        raise CustodyError("terminal rotation axis is degenerate")
    unit = direction / norm
    angle = math.radians(degrees)
    cosine, sine = math.cos(angle), math.sin(angle)
    result = dict(coordinates)
    for atom in component - set(axis):
        point = np.asarray(coordinates[atom], dtype=float) - origin
        rotated = point * cosine + np.cross(unit, point) * sine + unit * float(np.dot(unit, point)) * (1 - cosine)
        result[atom] = tuple(float(v) for v in rotated + origin)
    return result
