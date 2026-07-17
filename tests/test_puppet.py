"""Regression tests for `dvbfixer puppet`."""

from __future__ import annotations

from pathlib import Path

import pytest

from dvbfixer.puppet import BACKBONE_ATOMS, _parse_keep, main


def _iter_atoms(pdb: Path):
    for ln in pdb.read_text().splitlines():
        if ln.startswith("ATOM"):
            yield ln


def test_parse_keep_single() -> None:
    assert _parse_keep(["A:100"]) == {("A", 100)}


def test_parse_keep_range() -> None:
    assert _parse_keep(["A:100-102"]) == {("A", 100), ("A", 101), ("A", 102)}


def test_parse_keep_mixed_single_chain() -> None:
    """One --keep spec covers ONE chain; multiple chains need repeated --keep."""
    kept = _parse_keep(["A:100-101,105"])
    assert kept == {("A", 100), ("A", 101), ("A", 105)}


def test_parse_keep_repeated_for_multiple_chains() -> None:
    """Two chains require two --keep specs."""
    kept = _parse_keep(["A:100", "B:200"])
    assert kept == {("A", 100), ("B", 200)}


def test_puppet_strips_sidechains(tmp_workdir: Path, default_pdb: Path) -> None:
    out = tmp_workdir / "puppet.pdb"
    main([str(default_pdb), "-o", str(out)])
    for ln in _iter_atoms(out):
        atom = ln[12:16].strip()
        assert atom in BACKBONE_ATOMS, f"sidechain atom leaked through: {atom}"


def test_puppet_renames_to_gly(tmp_workdir: Path, default_pdb: Path) -> None:
    out = tmp_workdir / "puppet.pdb"
    main([str(default_pdb), "-o", str(out)])
    for ln in _iter_atoms(out):
        assert ln[17:20] == "GLY", f"non-GLY residue: {ln[17:20]}"


def test_puppet_keep_preserves_original(tmp_workdir: Path, default_pdb: Path) -> None:
    """--keep should retain the residue's original name + all atoms."""
    # default.pdb starts with CYX at H:1
    out = tmp_workdir / "puppet_keep.pdb"
    main([str(default_pdb), "-o", str(out), "--keep", "H:1"])
    kept_lines = [ln for ln in _iter_atoms(out) if ln[21] == "H" and ln[22:26].strip() == "1"]
    assert kept_lines, "expected --keep to preserve chain H residue 1"
    resnames = {ln[17:20] for ln in kept_lines}
    assert resnames == {"CYX"}, f"expected CYX preserved, got {resnames}"


def test_puppet_prints_summary(tmp_workdir: Path, default_pdb: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_workdir / "puppet.pdb"
    main([str(default_pdb), "-o", str(out)])
    captured = capsys.readouterr()
    assert "backbone atoms" in captured.out
