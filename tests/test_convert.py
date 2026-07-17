"""Regression tests for `dvbfixer convert` (legacy `glycam`).

Uses `default.pdb` — its CYX residues exercise the AMBER protein-variant
rename path (CYX → CYS on `--to-charmm`) without needing glycan fixtures.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dvbfixer.glycam import main


def _resnames(pdb: Path) -> set[str]:
    names = set()
    for ln in pdb.read_text().splitlines():
        if ln.startswith(("ATOM  ", "HETATM")):
            names.add(ln[17:20].strip())
    return names


def test_convert_to_charmm_renames_cyx(tmp_workdir: Path, default_pdb: Path) -> None:
    """AMBER CYX should become CHARMM CYS under --to-charmm."""
    assert "CYX" in _resnames(default_pdb)
    out = tmp_workdir / "charmm.pdb"
    main([str(default_pdb), "-o", str(out), "--to-charmm"])
    names = _resnames(out)
    assert "CYX" not in names, "CYX should be renamed under --to-charmm"
    assert "CYS" in names


def test_convert_to_amber_is_noop_on_amber_input(tmp_workdir: Path, default_pdb: Path) -> None:
    """default.pdb is already AMBER-named; --to-amber should be idempotent."""
    out = tmp_workdir / "amber.pdb"
    main([str(default_pdb), "-o", str(out), "--to-amber"])
    # Same protein-variant names should still be present.
    before = _resnames(default_pdb)
    after = _resnames(out)
    for name in ("CYX", "HIS", "ASN"):
        if name in before:
            assert name in after, f"{name} should survive --to-amber pass"


def test_convert_round_trip_runs_cleanly(tmp_workdir: Path, default_pdb: Path) -> None:
    """AMBER → CHARMM → AMBER doesn't crash; output has a comparable atom count.

    Exact count preservation isn't guaranteed — variant-aware stale-H drops
    (e.g. HIE's HD1) can trim a few atoms. This test pins that the round-trip
    completes and that the count stays within 1% of the input.
    """
    to_charmm = tmp_workdir / "charmm.pdb"
    to_amber = tmp_workdir / "amber.pdb"
    main([str(default_pdb), "-o", str(to_charmm), "--to-charmm"])
    main([str(to_charmm), "-o", str(to_amber), "--to-amber"])

    def atom_count(pdb: Path) -> int:
        return sum(1 for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM  ", "HETATM")))

    before = atom_count(default_pdb)
    after = atom_count(to_amber)
    assert after > 0
    tolerance = max(5, before // 100)  # 1% or 5 atoms, whichever is larger
    assert abs(before - after) <= tolerance, f"round-trip atom count drifted: {before} -> {after}"
