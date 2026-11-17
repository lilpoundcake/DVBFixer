"""Private explicit generated/observed S-S torsion geometry restraint.

Only caller-approved actual native bonds with exactly one mobile SG may enter.
Observed-only bonds are never restrained here. This is a declared geometric
restraint, not an unmodified force-field prediction or scientific acceptance.
"""
from __future__ import annotations

import math
from typing import Any

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, Pair


def build_generated_disulfide_torsion_force(
    topology: Any, *, pairs: tuple[Pair, ...], mobile: frozenset[AtomIdentity],
    strength_kj_mol: float,
):
    """Symmetric minima at +/-90deg; preserve target graph and final native indices."""
    from openmm import CustomTorsionForce

    if (isinstance(strength_kj_mol, bool) or not math.isfinite(strength_kj_mol)
            or strength_kj_mol <= 0 or not pairs or not mobile):
        raise CustodyError("generated S-S torsion requires explicit finite positive strength/pairs/mobile mask")
    normalized = tuple(tuple(sorted(pair)) for pair in pairs)
    endpoints = [atom for pair in normalized for atom in pair]
    if (any(len(pair) != 2 or pair[0] == pair[1] or any(a.atom_name != "SG" for a in pair) for pair in normalized)
            or len(set(normalized)) != len(normalized) or len(set(endpoints)) != len(endpoints)
            or any(len(set(pair) & mobile) != 1 for pair in normalized)):
        raise CustodyError("generated S-S torsion requires singly partnered pairs with exactly one mobile SG")
    wanted = set(endpoints) | {AtomIdentity(a.chain, a.residue_number, a.insertion_code, "CB") for a in endpoints}
    axis = {}
    for atom in topology.atoms():
        identity = AtomIdentity(atom.residue.chain.id, str(atom.residue.id), atom.residue.insertionCode.strip(), atom.name)
        if identity not in wanted:
            continue
        if identity in axis or atom.element is None or atom.element.symbol != ("S" if atom.name == "SG" else "C"):
            raise CustodyError("generated S-S torsion identity/element is ambiguous")
        axis[identity] = atom
    if axis.keys() != wanted:
        raise CustodyError("generated S-S torsion atoms absent from final native topology")
    bonds = {frozenset((a.index, b.index)) for a, b in topology.bonds()}
    force = CustomTorsionForce("k*(1+cos(2*theta))")
    force.setName("private-generated-observed-disulfide-torsion-v1")
    force.addPerTorsionParameter("k")
    for first, last in normalized:
        cb_first = AtomIdentity(first.chain, first.residue_number, first.insertion_code, "CB")
        cb_last = AtomIdentity(last.chain, last.residue_number, last.insertion_code, "CB")
        order = [axis[a].index for a in (cb_first, first, last, cb_last)]
        if any(frozenset(pair) not in bonds for pair in zip(order, order[1:], strict=False)):
            raise CustodyError("generated S-S torsion lacks required native CB-SG-SG-CB bonds")
        force.addTorsion(*order, [strength_kj_mol])
    return force
