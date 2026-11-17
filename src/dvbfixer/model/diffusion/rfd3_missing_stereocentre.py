"""Private one-point reflection for an explicitly reviewed missing stereocentre.

The canonical diagnostic and target chemistry must establish the problem outside
this generic geometric primitive. Never reflect an observed atom or rename atoms.
No public automatic stereochemistry repair is enabled.
"""
from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError


def reflect_reviewed_missing_atom(
    coordinates: Mapping[AtomIdentity, tuple[float, float, float]],
    *, atom: AtomIdentity, plane: tuple[AtomIdentity, AtomIdentity, AtomIdentity],
    allowed_missing: frozenset[AtomIdentity], observed: frozenset[AtomIdentity],
    max_displacement_angstrom: float,
) -> dict[AtomIdentity, tuple[float, float, float]]:
    """Only one missing point moves; distances to all three plane atoms stay fixed."""
    if (atom not in allowed_missing or atom in observed or allowed_missing & observed
            or atom in plane or len(set(plane)) != 3 or not {atom, *plane} <= coordinates.keys()):
        raise CustodyError("reviewed stereocentre reflection has invalid observed/missing/plane ownership")
    if (isinstance(max_displacement_angstrom, bool) or not math.isfinite(max_displacement_angstrom)
            or max_displacement_angstrom <= 0):
        raise CustodyError("reviewed reflection requires explicit finite positive movement bound")
    if any(len(xyz) != 3 or not all(math.isfinite(v) for v in xyz) for xyz in coordinates.values()):
        raise CustodyError("reviewed reflection coordinates must be finite three-vectors")
    center, first, second = [np.asarray(coordinates[a], dtype=float) for a in plane]
    a, b = first - center, second - center
    normal = np.cross(a, b)
    norm2 = float(np.dot(normal, normal))
    product = float(np.dot(a, a) * np.dot(b, b))
    if product <= 1e-16 or norm2 < .01 * product:
        raise CustodyError("reviewed reflection plane is near-degenerate")
    point = np.asarray(coordinates[atom], dtype=float)
    reflected = point - 2 * float(np.dot(point - center, normal)) / norm2 * normal
    if float(np.linalg.norm(reflected - point)) > max_displacement_angstrom:
        raise CustodyError("reviewed reflection exceeds declared one-point movement scope")
    result = dict(coordinates)
    result[atom] = tuple(float(v) for v in reflected)
    return result
