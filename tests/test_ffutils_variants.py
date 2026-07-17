"""Unit tests for the shared variant helpers in ``dvbfixer.ffutils.variants``."""

from __future__ import annotations

from pathlib import Path

import pytest

from dvbfixer.ffutils.variants import (
    ALL_VARIANTS,
    AMBER_VARIANT_TO_PARENT,
    CHARMM_VARIANT_TO_PARENT,
    VARIANT_TO_PARENT,
    scan_variant_names,
    text_rename_variants_to_parent,
)


def test_variant_tables_cover_expected_names() -> None:
    for name in ("HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM", "LYN"):
        assert AMBER_VARIANT_TO_PARENT[name]
    for name in ("HSD", "HSE", "HSP", "ASPP", "GLUP", "LSN"):
        assert CHARMM_VARIANT_TO_PARENT[name]
    # Combined lookup is the union.
    for name in ALL_VARIANTS:
        assert VARIANT_TO_PARENT[name] in {"HIS", "ASP", "GLU", "CYS", "LYS"}


def test_variant_tables_do_not_overlap() -> None:
    """No name appears in both AMBER and CHARMM tables (would create ambiguity)."""
    assert not set(AMBER_VARIANT_TO_PARENT) & set(CHARMM_VARIANT_TO_PARENT)


def test_scan_variant_names_finds_cyx(default_pdb: Path) -> None:
    """default.pdb has CYX residues; scan should locate them."""
    saved = scan_variant_names(default_pdb)
    assert saved, "expected at least one variant"
    assert any(name == "CYX" for name in saved.values()), (
        f"expected CYX in scan, got {set(saved.values())}"
    )


def test_scan_variant_names_ignores_standard_residues(small_pdb: Path) -> None:
    """ASN.pdb has one ASN — no variants."""
    saved = scan_variant_names(small_pdb)
    assert saved == {}


def test_text_rename_variants_to_parent_rewrites_cyx(
    tmp_workdir: Path, default_pdb: Path
) -> None:
    """CYX residues in default.pdb → CYS in the rewritten file."""
    out_path, saved = text_rename_variants_to_parent(default_pdb)
    assert saved, "expected renames"
    assert out_path != str(default_pdb), "should have written a temp file"

    # Rewritten file has no CYX left.
    rewritten_names = {
        ln[17:20].strip()
        for ln in Path(out_path).read_text().splitlines()
        if ln.startswith(("ATOM  ", "HETATM"))
    }
    assert "CYX" not in rewritten_names
    assert "CYS" in rewritten_names

    # Clean up temp file.
    Path(out_path).unlink(missing_ok=True)


def test_text_rename_variants_to_parent_no_op_on_standard(
    small_pdb: Path,
) -> None:
    """ASN.pdb has no variants — rename returns the input path unchanged."""
    out_path, saved = text_rename_variants_to_parent(small_pdb)
    assert out_path == str(small_pdb)
    assert saved == {}


def test_text_rename_variants_preserves_columns(
    tmp_workdir: Path, default_pdb: Path
) -> None:
    """Rewritten line must keep atom coord columns byte-identical."""
    orig_lines = {
        ln[6:11].strip(): ln[30:54]
        for ln in default_pdb.read_text().splitlines()
        if ln.startswith(("ATOM  ", "HETATM"))
    }
    out_path, _saved = text_rename_variants_to_parent(default_pdb)
    if out_path == str(default_pdb):
        pytest.skip("no rewrite happened; nothing to compare")

    for ln in Path(out_path).read_text().splitlines():
        if ln.startswith(("ATOM  ", "HETATM")):
            serial = ln[6:11].strip()
            assert ln[30:54] == orig_lines[serial], f"coord drift at serial {serial}"
    Path(out_path).unlink(missing_ok=True)
