"""Ligand-agnostic valence/protonation-state correction helpers.

PDB files carry no bond-order information, so RDKit's and OpenBabel's
proximity-based bond perception + naive valence-filling silently
over-protonates resonance-delocalized ionizable groups (carboxylate,
sulfonate/sulfate, phosphate) and misses genuine alkenes whose bond
length looks single-bond-like at typical crystallographic resolution
(e.g. DAN's ring C2=C3 — "2,3-didehydro" sialic acid analog, DANA).

Shared by ``prepare/glycan.py``'s two heterogen H-addition passes
(RDKit + OpenBabel) and ``lig_params.py``'s ``--parametrize-ligands``
SDF extraction, so this knowledge lives in exactly one place instead
of drifting across three copies.
"""

from __future__ import annotations

from typing import Any

# Known ligand double bonds that pure 3D-distance perception can't
# recover — a genuine alkene whose bond length isn't obviously shorter
# than a single bond at typical crystallographic resolution. Not
# generalizable via geometry; this is ligand-specific chemistry.
#
# Used only for `lig_params.py`'s SDF construction, where antechamber
# reads the explicit bond block and genuinely needs the correct order.
# NOT used to control RDKit's `AddHs`/OpenBabel's `addh()` H-count in
# `glycan.py` — verified empirically that both compute "how many H to
# add" from an atom's CONNECTION COUNT (degree), not from the bond-
# order-weighted valence sum, so changing bond order has no effect on
# how many H they add. See `_H_COUNT_OVERRIDES` below for that case.
#
# Applying this ALSO resets every other intra-residue bond to single
# first (see the apply functions below) rather than just adding the
# listed double bond on top of whatever RDKit/OpenBabel's own
# proximity-based perception already guessed. That's necessary, not
# just tidy: conjugation shortens neighboring bonds too (DAN's C1-C2 at
# 1.43 Å and C2-O6 at 1.30 Å are both single bonds in reality, shortened
# by the adjacent C2=C3 alkene), and RDKit's proximityBonding can
# independently misassign double-bond character to those neighbors —
# even cascading further around the ring via its own Kekulization pass.
# Patching individual misassigned bonds one at a time is whack-a-mole;
# resetting the whole residue to single bonds first and then applying
# only the genuinely-known exception is robust regardless of what else
# proximity bonding got wrong nearby.
_KNOWN_DOUBLE_BONDS: dict[str, list[tuple[str, str]]] = {
    'DAN': [
        ('C2', 'C3'),    # ring alkene, "2,3-didehydro" sialic acid analog (DANA)
        ('C10', 'O10'),  # N-acetyl amide carbonyl — OpenBabel's
                          # PerceiveBondOrders also misses this on the
                          # isolated heavy-atom skeleton, leaving both
                          # atoms one bond short of a valid valence
                          # (a "radical") when no H is added to compensate.
    ],
}

# Known ligand atoms whose correct H count differs from what RDKit's/
# OpenBabel's degree-based valence filling computes — an sp2 atom in a
# ring alkene (e.g. DAN's C2=C3) has fewer H slots than the same atom
# would if treated as sp3, but neither library's automatic H-adder
# infers that from geometry alone (see `_KNOWN_DOUBLE_BONDS` above for
# why bond-order overrides don't fix this). Applied as a post-hoc cap
# on however many H the library wanted to add for that atom — same
# proven mechanism as `prepare/glycan.py`'s existing `_NO_H_ATOM_NAMES`
# (which is really just this table's all-zero special case).
_H_COUNT_OVERRIDES: dict[str, dict[str, int]] = {
    'DAN': {'C2': 0, 'C3': 1},
}


def get_h_count_override(resname: str, atom_name: str) -> int | None:
    """Return the correct H count for `resname`'s `atom_name` if known,
    else None (no override — let the library's own count stand)."""
    return _H_COUNT_OVERRIDES.get(resname, {}).get(atom_name)


def find_ionizable_terminal_oxygens_rdkit(mol: Any) -> set[int]:
    """RDKit atom indices of terminal O in a carboxylate/sulfonate-or-
    sulfate/phosphate group.

    Detection is purely connectivity-based: an O with exactly one
    heavy-atom bond, to a C/S/P center that itself has >= 2
    (carboxylate) or >= 3 (sulfonate/sulfate/phosphate) such terminal,
    single-heavy-neighbor O's. No charge or bond-order info needed —
    these groups are ionized at physiological pH regardless of what a
    naive per-atom valence filler would compute. A plain carbonyl/
    amide/ester oxygen never has a second terminal O on the same
    center, so this doesn't false-positive on those.
    """
    ionizable: set[int] = set()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ('C', 'S', 'P'):
            continue
        terminal_o = []
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != 'O':
                continue
            heavy_neighbors = [a for a in nbr.GetNeighbors()
                               if a.GetSymbol() != 'H']
            if len(heavy_neighbors) == 1:
                terminal_o.append(nbr.GetIdx())
        min_count = 2 if symbol == 'C' else 3
        if len(terminal_o) >= min_count:
            ionizable.update(terminal_o)
    return ionizable


def apply_double_bond_overrides_rdkit(mol: Any) -> None:
    """Apply :data:`_KNOWN_DOUBLE_BONDS` to every matching residue
    instance in `mol` (which may contain multiple copies of an
    overridden ligand, e.g. one per crystallographic chain).

    Resets every intra-residue bond to single first (see module
    docstring on `_KNOWN_DOUBLE_BONDS` for why), then applies the known
    double bond(s) on top — robust regardless of what else RDKit's
    proximity-based perception guessed wrong nearby.

    `mol` must have been built with PDB residue info attached to every
    atom (e.g. via ``Chem.MolFromPDBBlock``).
    """
    from rdkit import Chem

    by_residue: dict[tuple, set[int]] = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        resname = info.GetResidueName().strip()
        if resname not in _KNOWN_DOUBLE_BONDS:
            continue
        key = (resname, info.GetChainId(), info.GetResidueNumber())
        by_residue.setdefault(key, set()).add(atom.GetIdx())

    for (resname, _chain, _resnum), atom_indices in by_residue.items():
        name_to_idx = {
            mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetName().strip(): i
            for i in atom_indices
        }
        for bond in mol.GetBonds():
            if (bond.GetBeginAtomIdx() in atom_indices
                    and bond.GetEndAtomIdx() in atom_indices):
                bond.SetBondType(Chem.BondType.SINGLE)
        for a1_name, a2_name in _KNOWN_DOUBLE_BONDS[resname]:
            i1 = name_to_idx.get(a1_name)
            i2 = name_to_idx.get(a2_name)
            if i1 is None or i2 is None:
                continue
            bond = mol.GetBondBetweenAtoms(i1, i2)
            if bond is not None:
                bond.SetBondType(Chem.BondType.DOUBLE)


def find_ionizable_terminal_oxygens_openbabel(obmol: Any) -> set[int]:
    """Same detection as :func:`find_ionizable_terminal_oxygens_rdkit`,
    using the OpenBabel ``OBMol``/``OBAtom`` API (1-based atom
    indices, per OB convention). Used by both `glycan.py`'s OpenBabel
    H-addition pass and `lig_params.py`'s SDF extraction — both call
    this on a single already-isolated residue's OBMol.
    """
    from openbabel import openbabel as ob

    ionizable: set[int] = set()
    for atom in ob.OBMolAtomIter(obmol):
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (6, 16, 15):  # C, S, P
            continue
        terminal_o = []
        for nbr in ob.OBAtomAtomIter(atom):
            if nbr.GetAtomicNum() != 8:
                continue
            heavy_deg = sum(1 for n2 in ob.OBAtomAtomIter(nbr)
                             if n2.GetAtomicNum() != 1)
            if heavy_deg == 1:
                terminal_o.append(nbr.GetIdx())
        min_count = 2 if atomic_num == 6 else 3
        if len(terminal_o) >= min_count:
            ionizable.update(terminal_o)
    return ionizable


def apply_double_bond_override_openbabel(resname: str,
                                         name_to_atom: dict) -> None:
    """Apply :data:`_KNOWN_DOUBLE_BONDS` for `resname` to a single
    already-isolated residue's OBMol, given a ``{atom_name: OBAtom}``
    mapping covering the WHOLE residue (both call sites already have
    the per-residue atom list in hand, in a different shape each).

    Resets every bond between two atoms in `name_to_atom` to single
    first (see module docstring on `_KNOWN_DOUBLE_BONDS` for why), then
    applies the known double bond(s) on top. No-op if `resname` has no
    override entry.
    """
    doubles = _KNOWN_DOUBLE_BONDS.get(resname)
    if not doubles:
        return
    from openbabel import openbabel as ob

    atoms = list(name_to_atom.values())
    atom_idx_set = {a.GetIdx() for a in atoms}
    for atom in atoms:
        for bond in ob.OBAtomBondIter(atom):
            other = bond.GetNbrAtom(atom)
            if other.GetIdx() in atom_idx_set:
                bond.SetBondOrder(1)
    for a1_name, a2_name in doubles:
        a1 = name_to_atom.get(a1_name)
        a2 = name_to_atom.get(a2_name)
        if a1 is None or a2 is None:
            continue
        bond = a1.GetBond(a2)
        if bond is not None:
            bond.SetBondOrder(2)


def assign_ionizable_bond_orders_openbabel(obmol: Any) -> None:
    """Set correct Kekule bond orders + formal charges on every
    carboxylate/sulfonate/sulfate/phosphate group found in `obmol`
    (via the same detection as
    :func:`find_ionizable_terminal_oxygens_openbabel`), so antechamber/
    GAFFTemplateGenerator computes the correct net molecular charge.

    Only relevant for `lig_params.py`'s SDF path — `glycan.py`'s
    H-suppression only needs the atom set, not a specific Kekule
    assignment (both O's end up with 0 H either way there).

    For each detected group: the (n-1) shortest center-O bonds become
    double (order 2, neutral); the single longest becomes a single
    bond (order 1) with formal charge -1 on that O — net -1 per group,
    the mono-anionic default already implicit elsewhere in this
    pipeline. A genuinely dianionic phosphate would need a per-ligand
    override; not attempted here.
    """
    from openbabel import openbabel as ob

    for atom in ob.OBMolAtomIter(obmol):
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (6, 16, 15):
            continue
        group_bonds = []  # (bond, other_atom, length)
        for bond in ob.OBAtomBondIter(atom):
            other = bond.GetNbrAtom(atom)
            if other.GetAtomicNum() != 8:
                continue
            heavy_deg = sum(1 for n2 in ob.OBAtomAtomIter(other)
                             if n2.GetAtomicNum() != 1)
            if heavy_deg != 1:
                continue
            group_bonds.append((bond, other, bond.GetLength()))
        min_count = 2 if atomic_num == 6 else 3
        if len(group_bonds) < min_count:
            continue
        group_bonds.sort(key=lambda t: t[2])  # ascending length
        for bond, other, _len in group_bonds[:-1]:
            bond.SetBondOrder(2)
            other.SetFormalCharge(0)
        last_bond, last_other, _len = group_bonds[-1]
        last_bond.SetBondOrder(1)
        last_other.SetFormalCharge(-1)
