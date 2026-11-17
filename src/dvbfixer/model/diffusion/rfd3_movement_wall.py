"""Private generated-only inward-buffered wall; final movement audit still required.

This is a finite force, NOT a mathematical hard constraint. It cannot grant an
exception if optimization exceeds the prospectively declared outer movement gate.
"""
from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError


def build_missing_atom_wall(
    topology: Any, references: Mapping[AtomIdentity, tuple[float, float, float]],
    *, mobile: frozenset[AtomIdentity], outer_limit_angstrom: float,
    buffer_angstrom: float, strength_kj_mol_nm2: float,
):
    """Resolve identities on FINAL native topology, never carry pre-rebuild indices."""
    from openmm import CustomExternalForce

    for value in (outer_limit_angstrom, buffer_angstrom, strength_kj_mol_nm2):
        if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise CustodyError("missing-atom wall parameters must be finite positive")
    if buffer_angstrom >= outer_limit_angstrom or not mobile or references.keys() != mobile:
        raise CustodyError("missing-atom wall requires exact mobile references and an inward buffer")
    axis = {}
    for atom in topology.atoms():
        identity = AtomIdentity(atom.residue.chain.id, str(atom.residue.id), atom.residue.insertionCode.strip(), atom.name)
        if identity not in mobile:
            continue
        if identity in axis or atom.element is None or atom.element.symbol in {"H", "D"}:
            raise CustodyError("ambiguous/non-heavy missing-atom wall identity")
        axis[identity] = atom.index
    if axis.keys() != mobile:
        raise CustodyError("missing-atom wall mobile identity absent from final topology")
    force = CustomExternalForce(
        "0.5*k*max(0,sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2+1e-16)-radius)^2")
    force.setName("private-missing-atom-inward-wall-v1")
    for name in ("k", "x0", "y0", "z0", "radius"):
        force.addPerParticleParameter(name)
    radius = (outer_limit_angstrom - buffer_angstrom) / 10
    for identity in sorted(mobile):
        xyz = references[identity]
        if len(xyz) != 3 or not all(math.isfinite(v) for v in xyz):
            raise CustodyError("nonfinite missing-atom wall reference")
        force.addParticle(axis[identity], [strength_kj_mol_nm2, *(v / 10 for v in xyz), radius])
    return force
