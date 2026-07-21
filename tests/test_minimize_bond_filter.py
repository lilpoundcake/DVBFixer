"""Unit tests for `_drop_spurious_inter_aa_bonds` in
`dvbfixer.minimize.pipeline`. Regression for the 0.6.4 bug where
CONECT inference added a spurious inter-residue bond (HID CG ↔
neighbour's O) that broke `createSystem` with "1 C-O bond too many".
"""

from __future__ import annotations

import pytest

pytest.importorskip("openmm", reason="needs OpenMM Topology")

from openmm.app import Element, Topology  # noqa: E402

from dvbfixer.minimize.pipeline import (  # noqa: E402
    _drop_spurious_inter_aa_bonds,
)


def _make_two_residues() -> tuple[Topology, list]:
    """Two HIS residues on the same chain. Returns (topology, atoms)."""
    top = Topology()
    chain = top.addChain("A")
    r1 = top.addResidue("HIS", chain)
    r2 = top.addResidue("HIS", chain)
    atoms = []
    # Sensible atom set for each residue.
    for r in (r1, r2):
        for name in ("N", "CA", "C", "O", "CB", "CG", "ND1"):
            elem = Element.getBySymbol(name[0])
            atoms.append(top.addAtom(name, elem, r))
    return top, atoms


def test_peptide_bond_kept() -> None:
    """The canonical peptide C(prev) - N(next) bond must survive the filter."""
    top, atoms = _make_two_residues()
    # r1 atoms are 0-6; r2 atoms are 7-13.
    r1_c = atoms[2]         # C of r1
    r2_n = atoms[7]         # N of r2
    top.addBond(r1_c, r2_n)
    n = _drop_spurious_inter_aa_bonds(top)
    assert n == 0
    assert (r1_c, r2_n) in list(top.bonds()) or (r2_n, r1_c) in list(top.bonds())


def test_spurious_cg_to_o_dropped() -> None:
    """A HIS CG ↔ next-residue O bond — the exact failure pattern from
    the 1TM1 report — must be filtered out."""
    top, atoms = _make_two_residues()
    r1_cg = atoms[5]        # CG of r1
    r2_o = atoms[10]        # O of r2
    top.addBond(r1_cg, r2_o)
    # Also add a legit peptide bond as a control.
    r1_c = atoms[2]
    r2_n = atoms[7]
    top.addBond(r1_c, r2_n)

    n = _drop_spurious_inter_aa_bonds(top)
    assert n == 1
    remaining = list(top.bonds())
    # Spurious CG-O gone.
    assert not any(
        (b1 is r1_cg and b2 is r2_o) or (b1 is r2_o and b2 is r1_cg)
        for b1, b2 in remaining
    )
    # Peptide bond retained.
    assert any(
        (b1 is r1_c and b2 is r2_n) or (b1 is r2_n and b2 is r1_c)
        for b1, b2 in remaining
    )


def test_intra_residue_heavy_heavy_bond_kept() -> None:
    """A single legitimate intra-residue heavy-atom bond is kept."""
    top, atoms = _make_two_residues()
    r1_ca = atoms[1]
    r1_cb = atoms[4]
    top.addBond(r1_ca, r1_cb)
    n = _drop_spurious_inter_aa_bonds(top)
    assert n == 0


def test_extra_h_bond_dropped() -> None:
    """An H atom must not have two bonds — the second is dropped
    regardless of residue type. Catches the corrupted CONECT case
    where backbone H got wired to CA/HA/O in addition to N."""
    top = Topology()
    chain = top.addChain("A")
    r = top.addResidue("HIS", chain)
    n = top.addAtom("N", Element.getBySymbol("N"), r)
    ca = top.addAtom("CA", Element.getBySymbol("C"), r)
    h = top.addAtom("H", Element.getBySymbol("H"), r)
    top.addBond(n, h)     # correct: N-H
    top.addBond(ca, h)    # spurious: CA-H (H already has one partner)
    n_dropped = _drop_spurious_inter_aa_bonds(top)
    assert n_dropped == 1
    # First bond (N-H) kept; second (CA-H) dropped.
    remaining = [(b[0], b[1]) for b in top.bonds()]
    assert (n, h) in remaining
    assert (ca, h) not in remaining
