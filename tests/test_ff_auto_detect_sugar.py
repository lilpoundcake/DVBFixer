"""Unit tests for the 0.7.6 FF auto-detection change: PDB-standard
sugar names should trigger the `amber+glycam` alias (with an
auto-convert reason), not fall back to `amber` with a warning."""
from __future__ import annotations

from pathlib import Path

_TRAST_STUB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
HETATM    3  C1  BGL D1001      50.000  50.000  50.000  1.00  0.00           C
HETATM    4  C1  BMA D1002      52.000  50.000  50.000  1.00  0.00           C
END
"""

_PURE_PROTEIN = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
END
"""


def test_pdb_sugars_trigger_amber_glycam(tmp_workdir: Path) -> None:
    """A structure with PDB-standard sugar names (BGL, BMA) should
    auto-detect as amber+glycam, not amber."""
    from dvbfixer.ffutils import detect_ff_from_pdb

    p = tmp_workdir / "glyco.pdb"
    p.write_text(_TRAST_STUB)
    alias, reason = detect_ff_from_pdb(p)
    assert alias == "amber+glycam", (
        f"expected 'amber+glycam' for PDB-sugar input, got {alias!r}: {reason}"
    )
    assert "auto-convert" in reason.lower() or "glycam" in reason.lower(), (
        f"reason should mention auto-conversion, got: {reason}"
    )


def test_no_sugars_stays_amber(tmp_workdir: Path) -> None:
    """A pure-protein input should stay amber (no glycam upgrade)."""
    from dvbfixer.ffutils import detect_ff_from_pdb

    p = tmp_workdir / "pure.pdb"
    p.write_text(_PURE_PROTEIN)
    alias, _ = detect_ff_from_pdb(p)
    assert alias == "amber"


def test_has_pdb_standard_sugars_helper(tmp_workdir: Path) -> None:
    """The has_pdb_standard_sugars helper should return True on glycos,
    False on pure protein."""
    from dvbfixer.ffutils import has_pdb_standard_sugars

    p_glyco = tmp_workdir / "glyco.pdb"
    p_glyco.write_text(_TRAST_STUB)
    assert has_pdb_standard_sugars(p_glyco) is True

    p_pure = tmp_workdir / "pure.pdb"
    p_pure.write_text(_PURE_PROTEIN)
    assert has_pdb_standard_sugars(p_pure) is False


def test_resolve_ff_auto_on_glycoprotein(tmp_workdir: Path) -> None:
    """resolve_ff('auto', <glycoprotein>) should return the
    amber+glycam alias (not just amber)."""
    from dvbfixer.ffutils import resolve_ff

    p = tmp_workdir / "glyco.pdb"
    p.write_text(_TRAST_STUB)
    xmls, alias, reason = resolve_ff("auto", p, verbose=False)
    assert alias == "amber+glycam"
    assert any("GLYCAM" in x for x in xmls), (
        f"amber+glycam alias should include GLYCAM XML in {xmls}"
    )


def test_resolve_ff_explicit_amber_upgrades_on_glycoprotein(
    tmp_workdir: Path,
) -> None:
    """If user passes --ff amber on a glycoprotein, resolve_ff should
    upgrade to amber+glycam (existing behaviour, now also triggered
    by PDB-standard sugars, not just GLYCAM markers)."""
    from dvbfixer.ffutils import resolve_ff

    p = tmp_workdir / "glyco.pdb"
    p.write_text(_TRAST_STUB)
    xmls, alias, reason = resolve_ff("amber", p, verbose=False)
    assert alias == "amber+glycam", (
        f"amber → amber+glycam upgrade should fire; got {alias!r}"
    )
    assert reason and "upgrade" in reason.lower()
