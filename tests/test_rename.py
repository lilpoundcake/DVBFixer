"""Regression tests for `dvbfixer rename`."""

from __future__ import annotations

from pathlib import Path

import pytest

from dvbfixer.rename import CANONICAL_MAP, canonicalize_pdb, main


def _resnames(pdb: Path) -> set[str]:
    names = set()
    for ln in pdb.read_text().splitlines():
        if ln.startswith(("ATOM  ", "HETATM")):
            names.add(ln[17:20].strip())
    return names


def test_canonical_map_is_complete() -> None:
    """AMBER + CHARMM + MSE all map to standard 3-letter codes."""
    for src in ("HIE", "HID", "HIP", "HSD", "HSE", "HSP"):
        assert CANONICAL_MAP[src] == "HIS"
    assert CANONICAL_MAP["ASH"] == "ASP"
    assert CANONICAL_MAP["GLH"] == "GLU"
    assert CANONICAL_MAP["CYX"] == "CYS"
    assert CANONICAL_MAP["CYM"] == "CYS"
    assert CANONICAL_MAP["LYN"] == "LYS"
    assert CANONICAL_MAP["MSE"] == "MET"


def test_canonicalize_pdb_renames_cyx_to_cys(tmp_workdir: Path, default_pdb: Path) -> None:
    """default.pdb has CYX residues; after rename they should be CYS."""
    assert "CYX" in _resnames(default_pdb), "fixture invariant"
    out = tmp_workdir / "renamed.pdb"
    n = canonicalize_pdb(default_pdb, out, verbose=False)
    assert n > 0, "expected at least one rename"
    names = _resnames(out)
    assert "CYX" not in names
    assert "CYS" in names


def test_canonicalize_pdb_is_idempotent(tmp_workdir: Path, small_pdb: Path) -> None:
    """ASN.pdb has only ASN — nothing to rename; second pass is a no-op."""
    out1 = tmp_workdir / "pass1.pdb"
    n1 = canonicalize_pdb(small_pdb, out1, verbose=False)
    assert n1 == 0
    out2 = tmp_workdir / "pass2.pdb"
    n2 = canonicalize_pdb(out1, out2, verbose=False)
    assert n2 == 0
    # ATOM/HETATM lines identical (canonicalize doesn't touch anything else)
    assert out1.read_text() == out2.read_text()


def test_canonicalize_pdb_preserves_coords(tmp_workdir: Path, default_pdb: Path) -> None:
    """Rename must not touch atom coordinates or the atom count."""
    before = [ln for ln in default_pdb.read_text().splitlines() if ln.startswith(("ATOM  ", "HETATM"))]
    out = tmp_workdir / "renamed.pdb"
    canonicalize_pdb(default_pdb, out, verbose=False)
    after = [ln for ln in out.read_text().splitlines() if ln.startswith(("ATOM  ", "HETATM"))]
    assert len(before) == len(after)
    # Coord columns (30-54) are untouched.
    for b, a in zip(before, after):
        assert b[30:54] == a[30:54]


def test_main_cli_writes_output(tmp_workdir: Path, default_pdb: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_workdir / "cli_out.pdb"
    main([str(default_pdb), "-o", str(out)])
    assert out.exists()
    captured = capsys.readouterr()
    assert "Wrote" in captured.out
