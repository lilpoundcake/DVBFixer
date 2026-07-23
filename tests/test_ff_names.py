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
