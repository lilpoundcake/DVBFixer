"""Unit tests for `dvbfixer.diagnose.chemistry`.

Covers the previously-untested paths flagged in the 0.6.0 assessment:
- ``_dihedral_deg`` (trans / cis peptide reference geometries)
- ``check_ca_chirality`` (L-Ala positive, D-Ala negative triple)
- ``check_valences`` (spurious 5th bond on carbon)
- ``check_bond_lengths`` (canonical vs stretched)
- ``check_disulfides`` (SS geometry deviations)
"""

from __future__ import annotations

import pytest

pytest.importorskip("openmm", reason="chemistry checks need OpenMM")

from openmm.app import Element, Topology  # noqa: E402
from openmm.unit import Quantity, nanometer  # noqa: E402

from dvbfixer.diagnose.chemistry import (  # noqa: E402
    _dihedral_deg,
    check_bond_lengths,
    check_ca_chirality,
    check_disulfides,
    check_valences,
)
from dvbfixer.diagnose.report import Severity  # noqa: E402

# ---------------------------------------------------------------------------
# _dihedral_deg — trans and cis reference geometries
# ---------------------------------------------------------------------------

def test_dihedral_trans_peptide_returns_180() -> None:
    # Cα(i) - C(i) - N(i+1) - Cα(i+1) coplanar, Cα(i+1) trans.
    omega = _dihedral_deg(
        (0.0, 0.0, 0.0),
        (1.5, 0.0, 0.0),
        (2.4, 1.3, 0.0),
        (3.9, 1.3, 0.0),
    )
    assert abs(abs(omega) - 180.0) < 0.5


def test_dihedral_cis_peptide_returns_0() -> None:
    omega = _dihedral_deg(
        (0.0, 0.0, 0.0),
        (1.5, 0.0, 0.0),
        (2.4, 1.3, 0.0),
        (0.9, 1.3, 0.0),
    )
    assert abs(omega) < 0.5


# ---------------------------------------------------------------------------
# check_ca_chirality — L vs D
# ---------------------------------------------------------------------------

def _make_ca_residue(
    n_pos: tuple[float, float, float],
    ca_pos: tuple[float, float, float],
    c_pos: tuple[float, float, float],
    cb_pos: tuple[float, float, float],
    resname: str = "ALA",
) -> tuple[Topology, Quantity]:
    top = Topology()
    chain = top.addChain("A")
    r = top.addResidue(resname, chain)
    for name, sym in [("N", "N"), ("CA", "C"), ("C", "C"), ("CB", "C")]:
        top.addAtom(name, Element.getBySymbol(sym), r)
    positions_nm = [tuple(x / 10.0 for x in p) for p in
                    (n_pos, ca_pos, c_pos, cb_pos)]
    return top, Quantity(positions_nm, nanometer)


def test_l_alanine_not_flagged_as_d() -> None:
    """Canonical L-Ala geometry (from AMBER14 ALA template)."""
    top, pos = _make_ca_residue(
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        cb_pos=(-0.86, -0.87, -0.87),
    )
    findings = check_ca_chirality(top, pos)
    assert findings == []


def test_d_alanine_flagged_as_error() -> None:
    """D-Ala mirror of L: swap CB and H positions."""
    top, pos = _make_ca_residue(
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        cb_pos=(-0.86, -0.87, 0.87),   # z flipped → D
    )
    findings = check_ca_chirality(top, pos)
    assert len(findings) == 1
    assert findings[0].severity == Severity.ERROR
    assert findings[0].category == "chirality"


# ---------------------------------------------------------------------------
# check_valences
# ---------------------------------------------------------------------------

def test_valence_5_on_carbon_flagged() -> None:
    """A carbon with 5 bonds must be reported as ERROR."""
    top = Topology()
    chain = top.addChain("A")
    r = top.addResidue("LIG", chain)
    c = top.addAtom("C1", Element.getBySymbol("C"), r)
    hs = [top.addAtom(f"H{i}", Element.getBySymbol("H"), r)
          for i in range(5)]
    for h in hs:
        top.addBond(c, h)
    findings = check_valences(top)
    # Expect one valence finding on C1 (5 > 4).
    val_findings = [f for f in findings if f.category == "valence"]
    assert len(val_findings) >= 1
    assert any(f.atom == "C1" for f in val_findings)


def test_valence_4_on_carbon_not_flagged() -> None:
    top = Topology()
    chain = top.addChain("A")
    r = top.addResidue("LIG", chain)
    c = top.addAtom("C1", Element.getBySymbol("C"), r)
    hs = [top.addAtom(f"H{i}", Element.getBySymbol("H"), r)
          for i in range(4)]
    for h in hs:
        top.addBond(c, h)
    findings = check_valences(top)
    assert [f for f in findings if f.category == "valence"] == []


def test_four_coordinate_sulfur_not_flagged() -> None:
    """Sulfonates/sulfates and EPE-like buffers commonly have four S bonds."""
    top = Topology()
    chain = top.addChain("A")
    residue = top.addResidue("EPE", chain)
    sulfur = top.addAtom("S", Element.getBySymbol("S"), residue)
    neighbors = [
        top.addAtom(f"O{i}", Element.getBySymbol("O"), residue)
        for i in range(4)
    ]
    for atom in neighbors:
        top.addBond(sulfur, atom)
    findings = check_valences(top)
    assert [f for f in findings if f.category == "valence"] == []


# ---------------------------------------------------------------------------
# check_bond_lengths
# ---------------------------------------------------------------------------

def _cc_topology_with_bond(distance_a: float) -> tuple[Topology, Quantity]:
    """Two sp3 carbons at ``distance_a`` Å, bonded — for bond-length tests."""
    top = Topology()
    chain = top.addChain("A")
    r = top.addResidue("LIG", chain)
    c1 = top.addAtom("CX", Element.getBySymbol("C"), r)
    c2 = top.addAtom("CY", Element.getBySymbol("C"), r)
    top.addBond(c1, c2)
    positions_nm = [(0.0, 0.0, 0.0), (distance_a / 10.0, 0.0, 0.0)]
    return top, Quantity(positions_nm, nanometer)


def test_canonical_cc_bond_not_flagged() -> None:
    """A canonical sp3 C-C at 1.53 Å must pass silently."""
    top, pos = _cc_topology_with_bond(1.53)
    findings = check_bond_lengths(top, pos)
    assert findings == []


def test_stretched_cc_bond_under_20pct_not_flagged() -> None:
    """A C-C at 1.75 Å (14 % long) sits under the 20 % WARN floor."""
    top, pos = _cc_topology_with_bond(1.75)
    findings = check_bond_lengths(top, pos)
    assert findings == []


def test_broken_cc_bond_over_50pct_flagged_as_error() -> None:
    """A C-C at 2.50 Å (63 % long) exceeds the 50 % ERROR floor."""
    top, pos = _cc_topology_with_bond(2.50)
    findings = check_bond_lengths(top, pos)
    assert len(findings) == 1
    assert findings[0].severity == Severity.ERROR


# ---------------------------------------------------------------------------
# check_disulfides
# ---------------------------------------------------------------------------

def _two_cys_topology(sg_distance_a: float) -> tuple[Topology, Quantity]:
    """Two CYX residues with SG atoms at ``sg_distance_a`` Å along x.
    CA / CB placed to give a reasonable Cα-Cα and χ_ss dihedral.
    """
    top = Topology()
    chain = top.addChain("A")
    positions_a: list[tuple[float, float, float]] = []
    for offset in (0.0, sg_distance_a + 3.0):
        r = top.addResidue("CYX", chain)
        top.addAtom("CA", Element.getBySymbol("C"), r)
        top.addAtom("CB", Element.getBySymbol("C"), r)
        top.addAtom("SG", Element.getBySymbol("S"), r)
    # Positions: two CYX residues placed so SG-SG is along x axis.
    # r1: CA=(0,0,0), CB=(1.5,0,0), SG=(3.0, 0, 0)
    # r2: CA=(3+sg_d+3.0, 0, 0) ... simpler: just place SGs at 0 and sg_d.
    r1_ca = (0.0, 3.0, 0.0)
    r1_cb = (0.5, 1.5, 0.0)
    r1_sg = (0.0, 0.0, 0.0)
    r2_sg = (sg_distance_a, 0.0, 0.0)
    r2_cb = (sg_distance_a - 0.5, 1.5, 0.0)
    r2_ca = (sg_distance_a, 3.0, 0.0)
    positions_a = [r1_ca, r1_cb, r1_sg, r2_ca, r2_cb, r2_sg]
    positions_nm = [tuple(x / 10.0 for x in p) for p in positions_a]
    return top, Quantity(positions_nm, nanometer)


def test_canonical_disulfide_not_flagged() -> None:
    """SG-SG at 2.05 Å — canonical; no findings."""
    top, pos = _two_cys_topology(2.05)
    findings = check_disulfides(top, pos)
    # Canonical bond length; Cα-Cα falls in range; dihedral ~180°
    # (planar test placement) → dihedral WARNING is expected because
    # 180° isn't the ~90° ideal. So we expect only the dihedral flag,
    # not a length flag.
    length_findings = [f for f in findings if "SG-SG to" in f.message]
    assert length_findings == []


def test_stretched_disulfide_flagged() -> None:
    """SG-SG at 2.40 Å is grossly stretched — ERROR."""
    top, pos = _two_cys_topology(2.40)
    findings = check_disulfides(top, pos)
    length_findings = [f for f in findings if "SG-SG to" in f.message]
    assert len(length_findings) == 1
    assert length_findings[0].severity == Severity.ERROR


def test_no_disulfide_no_findings() -> None:
    """Two CYX 5 Å apart aren't bonded — no findings."""
    top, pos = _two_cys_topology(5.00)
    findings = check_disulfides(top, pos)
    assert findings == []
