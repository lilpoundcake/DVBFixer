"""Unit tests for ``dvbfixer.ffutils.ff_names.apply_variants_to_pdb_text``."""

from __future__ import annotations

from pathlib import Path

import pytest

from dvbfixer.ffutils.ff_names import (
    PROTONATION_AMBER_TO_CHARMM,
    apply_variants_to_pdb_text,
)

_CANONICAL_PDB = """\
ATOM      1  N   HIS A   5      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  HIS A   5      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  N   ASP A  12      20.000  10.000  10.000  1.00  0.00           N
ATOM      4  CA  ASP A  12      21.000  10.000  10.000  1.00  0.00           C
ATOM      5  N   GLU A  20      30.000  10.000  10.000  1.00  0.00           N
ATOM      6  N   CYS A  25      40.000  10.000  10.000  1.00  0.00           N
ATOM      7  N   LYN A  30      50.000  10.000  10.000  1.00  0.00           N
ATOM      8  NZ  LYN A  30      52.000  10.000  10.000  1.00  0.00           N
ATOM      9  HZ2 LYN A  30      52.500  10.500  10.500  1.00  0.00           H
ATOM     10  HZ3 LYN A  30      52.500   9.500  10.500  1.00  0.00           H
ATOM     11  N   HIS A  40      60.000  10.000  10.000  1.00  0.00           N
END
"""


def test_amber_target_rewrites_variants(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_CANONICAL_PDB)
    renames = {
        ("A", "5"): "HIE",
        ("A", "12"): "ASH",
        ("A", "20"): "GLH",
        ("A", "25"): "CYX",
        ("A", "40"): "HID",
    }
    n = apply_variants_to_pdb_text(pdb, renames, target_ff="amber",
                                    include_gromacs_lyn=False)
    assert n >= 5
    text = pdb.read_text()
    assert " HIE A   5" in text
    assert " ASH A  12" in text
    assert " GLH A  20" in text
    assert " CYX A  25" in text
    assert " HID A  40" in text


def test_amber_target_gromacs_lyn_rename(tmp_path: Path) -> None:
    """LYN HZ3 → HZ1 is applied by default under target_ff='amber'."""
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_CANONICAL_PDB)
    n = apply_variants_to_pdb_text(pdb, {}, target_ff="amber",
                                    include_gromacs_lyn=True)
    assert n >= 1
    text = pdb.read_text()
    assert " HZ1 LYN A  30" in text
    assert " HZ2 LYN A  30" in text     # untouched
    assert " HZ3 LYN A  30" not in text  # renamed


def test_amber_target_gromacs_lyn_off(tmp_path: Path) -> None:
    """include_gromacs_lyn=False keeps ff14SB HZ2 + HZ3."""
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_CANONICAL_PDB)
    apply_variants_to_pdb_text(pdb, {}, target_ff="amber",
                                include_gromacs_lyn=False)
    text = pdb.read_text()
    assert " HZ2 LYN A  30" in text
    assert " HZ3 LYN A  30" in text


def test_charmm_target_maps_names_and_lyn_atoms(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_CANONICAL_PDB)
    renames = {
        ("A", "5"): "HIE",
        ("A", "12"): "ASH",
        ("A", "20"): "GLH",
        ("A", "25"): "CYX",
        ("A", "30"): "LYN",
    }
    apply_variants_to_pdb_text(pdb, renames, target_ff="charmm")
    text = pdb.read_text()
    # AMBER→CHARMM residue mapping
    assert " HSE A   5" in text                    # HIE → HSE
    assert " ASPP A  12" in text or " ASPPA  12" in text  # ASH → ASPP (4-char)
    assert " GLUP A  20" in text or " GLUPA  20" in text  # GLH → GLUP (4-char)
    assert " CYS A  25" in text                    # CYX → CYS
    assert " LSN A  30" in text                    # LYN → LSN
    # LYN → LSN atom shift: HZ2→HZ1, HZ3→HZ2
    assert " HZ1 LSN A  30" in text or " HZ1 LSN " in text
    assert " HZ2 LSN A  30" in text or " HZ2 LSN " in text


def test_idempotent(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_CANONICAL_PDB)
    renames = {("A", "5"): "HIE", ("A", "12"): "ASH"}
    apply_variants_to_pdb_text(pdb, renames)
    n_second = apply_variants_to_pdb_text(pdb, renames)
    assert n_second == 0


def test_no_renames_no_variants_no_change(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    plain = "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\nEND\n"
    pdb.write_text(plain)
    n = apply_variants_to_pdb_text(pdb, {}, include_gromacs_lyn=False)
    assert n == 0
    assert pdb.read_text() == plain


def test_bad_target_ff_raises(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    pdb.write_text("END\n")
    with pytest.raises(ValueError):
        apply_variants_to_pdb_text(pdb, {}, target_ff="opls")


def test_maps_expose_amber_and_charmm_variants() -> None:
    """Sanity: the exported maps cover the expected AMBER variants."""
    assert PROTONATION_AMBER_TO_CHARMM["HID"] == "HSD"
    assert PROTONATION_AMBER_TO_CHARMM["HIE"] == "HSE"
    assert PROTONATION_AMBER_TO_CHARMM["ASH"] == "ASPP"
    assert PROTONATION_AMBER_TO_CHARMM["LYN"] == "LSN"


# ---------------------------------------------------------------------------
# 0.7.0: comprehensive GROMACS-canonical atom naming
# ---------------------------------------------------------------------------

_LYS_PDB = """\
ATOM      1  N   LYS A   5      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  LYS A   5      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  C   LYS A   5      12.000  10.000  10.000  1.00  0.00           C
ATOM      4  O   LYS A   5      12.500  10.500  10.500  1.00  0.00           O
ATOM      5  CB  LYS A   5      11.000  11.500  10.000  1.00  0.00           C
ATOM      6  HB2 LYS A   5      10.500  11.800  10.500  1.00  0.00           H
ATOM      7  HB3 LYS A   5      11.500  11.800   9.500  1.00  0.00           H
ATOM      8  CG  LYS A   5      11.000  12.500  11.000  1.00  0.00           C
ATOM      9  HG2 LYS A   5      10.500  12.800  11.500  1.00  0.00           H
ATOM     10  HG3 LYS A   5      11.500  12.800  10.500  1.00  0.00           H
ATOM     11  CD  LYS A   5      11.000  13.500  12.000  1.00  0.00           C
ATOM     12  HD2 LYS A   5      10.500  13.800  12.500  1.00  0.00           H
ATOM     13  HD3 LYS A   5      11.500  13.800  11.500  1.00  0.00           H
ATOM     14  CE  LYS A   5      11.000  14.500  13.000  1.00  0.00           C
ATOM     15  HE2 LYS A   5      10.500  14.800  13.500  1.00  0.00           H
ATOM     16  HE3 LYS A   5      11.500  14.800  12.500  1.00  0.00           H
ATOM     17  NZ  LYS A   5      11.000  15.500  14.000  1.00  0.00           N
END
"""


def test_lys_all_methylenes_shifted(tmp_path: Path) -> None:
    """LYS has β/γ/δ/ε methylenes; all should have H3 → H1 rename."""
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_LYS_PDB)
    n = apply_variants_to_pdb_text(pdb, {}, target_ff="amber")
    assert n == 4  # HB3, HG3, HD3, HE3
    text = pdb.read_text()
    # New names (single-source rename: HB3 → HB1 etc.).
    assert " HB1 LYS " in text
    assert " HG1 LYS " in text
    assert " HD1 LYS " in text
    assert " HE1 LYS " in text
    # HB2/HG2/HD2/HE2 stay put (chemically equivalent to HB1 after rename).
    assert " HB2 LYS " in text
    assert " HG2 LYS " in text
    assert " HD2 LYS " in text
    assert " HE2 LYS " in text
    # HB3/HG3/HD3/HE3 gone.
    assert " HB3 LYS " not in text
    assert " HG3 LYS " not in text
    assert " HD3 LYS " not in text
    assert " HE3 LYS " not in text


def test_lys_methylene_shift_idempotent(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_LYS_PDB)
    apply_variants_to_pdb_text(pdb, {}, target_ff="amber")
    n_second = apply_variants_to_pdb_text(pdb, {}, target_ff="amber")
    assert n_second == 0


_GLY_PDB = """\
ATOM      1  N   GLY A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  GLY A   1      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  HA2 GLY A   1      11.500  10.500  10.500  1.00  0.00           H
ATOM      4  HA3 GLY A   1      11.500   9.500  10.500  1.00  0.00           H
ATOM      5  C   GLY A   1      12.000  10.000  10.000  1.00  0.00           C
ATOM      6  O   GLY A   1      12.500  10.500  10.500  1.00  0.00           O
END
"""


def test_gly_alpha_methylene_shift(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_GLY_PDB)
    apply_variants_to_pdb_text(pdb, {}, target_ff="amber")
    text = pdb.read_text()
    assert " HA1 GLY " in text
    assert " HA2 GLY " in text
    assert " HA3 GLY " not in text


_ILE_PDB = """\
ATOM      1  N   ILE A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  ILE A   1      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  CG1 ILE A   1      11.000  11.500  10.000  1.00  0.00           C
ATOM      4 HG12 ILE A   1      10.500  11.800  10.500  1.00  0.00           H
ATOM      5 HG13 ILE A   1      11.500  11.800   9.500  1.00  0.00           H
ATOM      6  CD1 ILE A   1      11.000  12.500  10.000  1.00  0.00           C
END
"""


def test_ile_gamma_methylene_shift(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_ILE_PDB)
    apply_variants_to_pdb_text(pdb, {}, target_ff="amber")
    text = pdb.read_text()
    assert " HG11 ILE " in text
    assert " HG12 ILE " in text
    assert " HG13 ILE " not in text


def _multi_res_pdb(*residues):
    """Build a small PDB with the given (resname, resid, atoms) tuples."""
    lines = []
    serial = 1
    for resname, resid, atoms in residues:
        for name, x, y, z in atoms:
            atomfield = f" {name:<3s}" if len(name) < 4 else name[:4]
            lines.append(f"ATOM  {serial:>5d} {atomfield} {resname:<3s} A{resid:>4d}    "
                         f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           "
                         f"{name[0]}\n")
            serial += 1
    lines.append("END\n")
    return "".join(lines)


def test_nterminal_h_becomes_h1(tmp_path: Path) -> None:
    """First protein residue in a chain: H → H1."""
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_multi_res_pdb(
        ("ALA", 1, [("N", 0, 0, 0), ("H", 0.5, 0.5, 0), ("CA", 1, 0, 0), ("C", 2, 0, 0), ("O", 2.5, 0.5, 0)]),
        ("ALA", 2, [("N", 3, 0, 0), ("H", 3.5, 0.5, 0), ("CA", 4, 0, 0), ("C", 5, 0, 0), ("O", 5.5, 0.5, 0)]),
    ))
    apply_variants_to_pdb_text(pdb, {}, target_ff="amber")
    text = pdb.read_text()
    # First ALA (N-term) has H1; mid-chain ALA has H.
    lines = text.splitlines()
    ala1_h = [l for l in lines if " ALA A   1" in l and l[12:16].strip() in ("H", "H1")]
    ala2_h = [l for l in lines if " ALA A   2" in l and l[12:16].strip() in ("H", "H1")]
    assert len(ala1_h) == 1 and ala1_h[0][12:16].strip() == "H1"
    assert len(ala2_h) == 1 and ala2_h[0][12:16].strip() == "H"


def test_cterminal_o_oxt_becomes_oc1_oc2(tmp_path: Path) -> None:
    """Last protein residue in a chain: O → OC2, OXT → OC1."""
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_multi_res_pdb(
        ("ALA", 1, [("N", 0, 0, 0), ("CA", 1, 0, 0), ("C", 2, 0, 0), ("O", 2.5, 0.5, 0)]),
        ("ALA", 2, [("N", 3, 0, 0), ("CA", 4, 0, 0), ("C", 5, 0, 0), ("O", 5.5, 0.5, 0), ("OXT", 5.5, -0.5, 0)]),
    ))
    apply_variants_to_pdb_text(pdb, {}, target_ff="amber")
    text = pdb.read_text()
    # Last ALA (C-term): O → OC2, OXT → OC1.
    assert " OC1 ALA A   2" in text
    assert " OC2 ALA A   2" in text
    assert " OXT ALA A   2" not in text
    # First ALA's O (backbone carbonyl) untouched.
    assert " O   ALA A   1" in text


def test_charmm_backbone_h_to_hn(tmp_path: Path) -> None:
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_multi_res_pdb(
        ("HIS", 5, [("N", 0, 0, 0), ("H", 0.5, 0.5, 0), ("CA", 1, 0, 0), ("C", 2, 0, 0), ("O", 2.5, 0.5, 0)]),
    ))
    apply_variants_to_pdb_text(pdb, {("A", "5"): "HIE"}, target_ff="charmm")
    text = pdb.read_text()
    assert " HN  HSE " in text
    assert " H   HSE " not in text
    assert " HSE A   5" in text


def test_dna_map_exposed() -> None:
    """DNA/RNA rename map exists — end-to-end test deferred (apostrophe in
    atom names requires column-perfect PDB fixture handling)."""
    from dvbfixer.ffutils.ff_names import GROMACS_AMBER_NA_ATOM_RENAMES
    assert "DA" in GROMACS_AMBER_NA_ATOM_RENAMES
    assert GROMACS_AMBER_NA_ATOM_RENAMES["DA"]["H2'"] == "H2'1"
    assert GROMACS_AMBER_NA_ATOM_RENAMES["DA"]["H2''"] == "H2'2"
    # RNA has an extra HO'2 → HO2' entry.
    assert GROMACS_AMBER_NA_ATOM_RENAMES["A"]["HO'2"] == "HO2'"


def test_two_chain_nterminal_detection_per_chain(tmp_path: Path) -> None:
    """Each chain gets its own N-terminal / C-terminal residue detection."""
    pdb = tmp_path / "in.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\n"
        "ATOM      2  H   ALA A   1      10.500  10.500  10.000  1.00  0.00           H\n"
        "ATOM      3  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00           C\n"
        "TER\n"
        "ATOM      4  N   ALA B   1      20.000  10.000  10.000  1.00  0.00           N\n"
        "ATOM      5  H   ALA B   1      20.500  10.500  10.000  1.00  0.00           H\n"
        "ATOM      6  CA  ALA B   1      21.000  10.000  10.000  1.00  0.00           C\n"
        "END\n"
    )
    apply_variants_to_pdb_text(pdb, {}, target_ff="amber")
    text = pdb.read_text()
    # BOTH chains' first residue should have H → H1.
    assert " H1  ALA A   1" in text
    assert " H1  ALA B   1" in text
