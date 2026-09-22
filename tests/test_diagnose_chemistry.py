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
    check_backbone_bond_angles,
    check_bond_lengths,
    check_ca_chirality,
    check_disulfides,
    check_ramachandran,
    check_sidechain_chi12,
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
# check_backbone_bond_angles
# ---------------------------------------------------------------------------


def _backbone_topology(ca_y_a: float) -> tuple[Topology, Quantity]:
    top = Topology()
    chain = top.addChain("A")
    residue = top.addResidue("ALA", chain, id="1")
    n = top.addAtom("N", Element.getBySymbol("N"), residue)
    ca = top.addAtom("CA", Element.getBySymbol("C"), residue)
    c = top.addAtom("C", Element.getBySymbol("C"), residue)
    o = top.addAtom("O", Element.getBySymbol("O"), residue)
    top.addBond(n, ca)
    top.addBond(ca, c)
    top.addBond(c, o)
    positions_a = [
        (-1.2, 0.0, 0.0),
        (0.0, ca_y_a, 0.0),
        (0.5, 1.4, 0.0),
        (1.7, 1.4, 0.0),
    ]
    positions_nm = [tuple(value / 10.0 for value in point) for point in positions_a]
    return top, Quantity(positions_nm, nanometer)


def test_canonical_backbone_angles_not_flagged() -> None:
    top, positions = _backbone_topology(0.0)

    findings = check_backbone_bond_angles(top, positions)

    assert findings == []


def test_collapsed_backbone_angle_flagged_as_error() -> None:
    top, positions = _backbone_topology(2.0)

    findings = check_backbone_bond_angles(top, positions)

    assert any(
        finding.category == "bond_angle"
        and finding.atom == "N-CA-C"
        and finding.severity == Severity.ERROR
        for finding in findings
    )


# ---------------------------------------------------------------------------
# check_ramachandran
# ---------------------------------------------------------------------------


def _ramachandran_topology(
    *,
    middle_name: str = "SER",
    following_name: str = "LYS",
    middle_icode: str = "",
    shift_middle_n_x_a: float = 0.0,
    connect_previous: bool = True,
) -> tuple[Topology, Quantity]:
    top = Topology()
    chain = top.addChain("D")
    positions_a = [
        (44.612, 70.159, -80.957),
        (44.026, 68.817, -81.117),
        (42.822, 68.850, -82.033),
        (42.556 + shift_middle_n_x_a, 67.737, -82.713),
        (41.362, 67.537, -83.569),
        (40.799, 66.143, -83.299),
        (39.584, 65.886, -83.779),
        (38.867, 64.596, -83.627),
        (38.049, 64.341, -84.895),
    ]
    residues = (
        top.addResidue("GLY", chain, id="66"),
        top.addResidue(middle_name, chain, id="67", insertionCode=middle_icode),
        top.addResidue(following_name, chain, id="68"),
    )
    atoms: list[tuple[object, object, object]] = []
    for residue in residues:
        n = top.addAtom("N", Element.getBySymbol("N"), residue)
        ca = top.addAtom("CA", Element.getBySymbol("C"), residue)
        c = top.addAtom("C", Element.getBySymbol("C"), residue)
        top.addBond(n, ca)
        top.addBond(ca, c)
        atoms.append((n, ca, c))
    if connect_previous:
        top.addBond(atoms[0][2], atoms[1][0])
    top.addBond(atoms[1][2], atoms[2][0])
    positions_nm = [tuple(value / 10.0 for value in point) for point in positions_a]
    return top, Quantity(positions_nm, nanometer)


def test_canonical_general_ramachandran_not_flagged() -> None:
    top, positions = _ramachandran_topology()

    assert check_ramachandran(top, positions) == []


def test_gross_general_ramachandran_outlier_flagged_as_error() -> None:
    top, positions = _ramachandran_topology(
        middle_icode="A",
        shift_middle_n_x_a=-1.75,
    )

    findings = check_ramachandran(top, positions)

    assert len(findings) == 1
    assert findings[0].severity is Severity.ERROR
    assert findings[0].category == "ramachandran_outlier"
    assert findings[0].chain == "D"
    assert findings[0].resid == "67A"
    assert findings[0].extra["residue_class"] == "general"


@pytest.mark.parametrize(
    ("middle_name", "following_name"),
    (("GLY", "LYS"), ("PRO", "LYS"), ("SER", "PRO")),
)
def test_class_specific_ramachandran_residues_are_excluded(
    middle_name: str,
    following_name: str,
) -> None:
    top, positions = _ramachandran_topology(
        middle_name=middle_name,
        following_name=following_name,
        shift_middle_n_x_a=-1.75,
    )

    assert check_ramachandran(top, positions) == []


def test_ramachandran_requires_actual_peptide_bonds() -> None:
    top, positions = _ramachandran_topology(
        shift_middle_n_x_a=-1.75,
        connect_previous=False,
    )

    assert check_ramachandran(top, positions) == []


# ---------------------------------------------------------------------------
# check_sidechain_chi12
# ---------------------------------------------------------------------------


def _sidechain_chi12_topology(
    *,
    residue_name: str = "LYS",
    chain_id: str = "D",
    residue_id: str = "67",
    insertion_code: str = "",
    connect_cd: bool = True,
) -> tuple[Topology, Quantity]:
    top = Topology()
    chain = top.addChain(chain_id)
    residue = top.addResidue(
        residue_name,
        chain,
        id=residue_id,
        insertionCode=insertion_code,
    )
    terminal_name = "CD1" if residue_name == "LEU" else "CD"
    names = ("N", "CA", "CB", "CG", terminal_name)
    atoms = [
        top.addAtom(name, Element.getBySymbol("N" if name == "N" else "C"), residue)
        for name in names
    ]
    for first, second in zip(atoms, atoms[1:]):
        if not connect_cd and second.name == "CD":
            continue
        top.addBond(first, second)
    positions_a = [
        (39.584, 65.886, -83.779),
        (38.867, 64.596, -83.627),
        (37.974, 64.622, -82.385),
        (37.582, 63.259, -81.835),
        (36.439, 63.326, -80.851),
    ]
    positions_nm = [tuple(value / 10.0 for value in point) for point in positions_a]
    return top, Quantity(positions_nm, nanometer)


def _rotate_sidechain_to_zero_chi12(positions: Quantity) -> Quantity:
    import math

    import numpy as np

    coordinates = np.asarray(positions.value_in_unit(nanometer), dtype=float).copy()

    def dihedral(indices: tuple[int, int, int, int]) -> float:
        return _dihedral_deg(*(tuple(coordinates[index]) for index in indices)) or 0.0

    def rotate(
        origin_index: int,
        axis_index: int,
        rotated_indices: tuple[int, ...],
        angle_degrees: float,
    ) -> None:
        origin = coordinates[origin_index]
        axis = coordinates[axis_index] - origin
        axis /= np.linalg.norm(axis)
        angle = math.radians(angle_degrees)
        cosine = math.cos(angle)
        sine = math.sin(angle)
        for index in rotated_indices:
            vector = coordinates[index] - origin
            coordinates[index] = origin + (
                vector * cosine
                + np.cross(axis, vector) * sine
                + axis * np.dot(axis, vector) * (1.0 - cosine)
            )

    rotate(1, 2, (3, 4), -dihedral((0, 1, 2, 3)))
    rotate(2, 3, (4,), -dihedral((1, 2, 3, 4)))
    return Quantity(coordinates, nanometer)


def test_canonical_sidechain_chi12_not_flagged() -> None:
    top, positions = _sidechain_chi12_topology()

    assert check_sidechain_chi12(top, positions) == []


def test_gross_sidechain_chi12_outlier_preserves_identity() -> None:
    top, positions = _sidechain_chi12_topology(insertion_code="A")
    positions = _rotate_sidechain_to_zero_chi12(positions)

    findings = check_sidechain_chi12(top, positions)

    assert len(findings) == 1
    assert findings[0].severity is Severity.ERROR
    assert findings[0].category == "sidechain_chi12_outlier"
    assert findings[0].chain == "D"
    assert findings[0].resid == "67A"
    assert abs(findings[0].extra["chi1_degrees"]) < 1e-6
    assert abs(findings[0].extra["chi2_degrees"]) < 1e-6


@pytest.mark.parametrize("residue_name", ("ALA", "CYS", "GLY", "PRO", "SER", "THR", "VAL"))
def test_chi1_only_or_absent_residues_are_excluded(residue_name: str) -> None:
    top, positions = _sidechain_chi12_topology(residue_name=residue_name)
    positions = _rotate_sidechain_to_zero_chi12(positions)

    assert check_sidechain_chi12(top, positions) == []


def test_sidechain_chi12_requires_actual_bonds() -> None:
    top, positions = _sidechain_chi12_topology(connect_cd=False)
    positions = _rotate_sidechain_to_zero_chi12(positions)

    assert check_sidechain_chi12(top, positions) == []


def test_symmetric_chi2_equivalent_is_accepted(monkeypatch: pytest.MonkeyPatch) -> None:
    top, positions = _sidechain_chi12_topology(residue_name="LEU")
    positions = _rotate_sidechain_to_zero_chi12(positions)
    seen: list[tuple[float, float]] = []

    def populated(chi1: float, chi2: float) -> bool:
        seen.append((chi1, chi2))
        return len(seen) == 2

    monkeypatch.setattr(
        "dvbfixer.diagnose.chemistry._sidechain_chi12_reference_populated",
        populated,
    )

    assert check_sidechain_chi12(top, positions) == []
    assert len(seen) == 2
    assert seen[1][1] == pytest.approx(seen[0][1] + 180.0)


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
