"""Unit tests for the ``--strip-heterogens`` pre-pass in
``dvbfixer.model.pipeline._strip_hetatm_lines``.

Regression coverage for the 0.6.6 feature: HETATM records + their
orphan CONECT lines are dropped before Modeller runs. Water residues
are optionally preserved.
"""

from __future__ import annotations

from dvbfixer.model.pipeline import _strip_hetatm_lines

_SAMPLE_LINES = [
    "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N  \n",
    "ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00           C  \n",
    "ATOM      3  C   ALA A   1      12.000  10.000  10.000  1.00  0.00           C  \n",
    "HETATM    4  O   HOH A 100      20.000  20.000  20.000  1.00  0.00           O  \n",
    "HETATM    5  O   HOH A 101      21.000  20.000  20.000  1.00  0.00           O  \n",
    "HETATM    6  C1  CIT A 200      30.000  30.000  30.000  1.00  0.00           C  \n",
    "HETATM    7  O1  CIT A 200      31.000  30.000  30.000  1.00  0.00           O  \n",
    "CONECT    6    7\n",     # bond within CIT — should die with it
    "CONECT    3    4\n",     # protein-to-water bond — mixed
    "END\n",
]


def _resnames_and_starts(lines):
    """Return (line_starts, resnames-present) for easier assertions."""
    starts = [ln[:6].rstrip() for ln in lines]
    hetatm_resnames = [ln[17:20].strip() for ln in lines if ln.startswith("HETATM")]
    return starts, hetatm_resnames


def test_strip_without_keep_water_removes_all_hetatm() -> None:
    kept = _strip_hetatm_lines(_SAMPLE_LINES, keep_water=False)
    starts, hets = _resnames_and_starts(kept)
    assert "HETATM" not in starts
    assert hets == []
    # Intra-CIT CONECT (6-7) gone; mixed CONECT (3-4) also gone
    # (references dropped serial 4).
    conect_count = sum(1 for ln in kept if ln.startswith("CONECT"))
    assert conect_count == 0


def test_strip_with_keep_water_preserves_hoh() -> None:
    kept = _strip_hetatm_lines(_SAMPLE_LINES, keep_water=True)
    _starts, hets = _resnames_and_starts(kept)
    assert hets == ["HOH", "HOH"]
    # Mixed CONECT (3-4) references a kept HOH serial and stays.
    conect_lines = [ln for ln in kept if ln.startswith("CONECT")]
    assert len(conect_lines) == 1
    # Intra-CIT CONECT (6-7) gone.
    assert "6" not in conect_lines[0].split() or "7" not in conect_lines[0].split()


def test_no_hetatm_no_change() -> None:
    protein_only = [
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N  \n",
        "ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00           C  \n",
        "END\n",
    ]
    kept = _strip_hetatm_lines(protein_only, keep_water=False)
    assert kept == protein_only


def test_verbose_reports_count(capsys) -> None:
    _strip_hetatm_lines(_SAMPLE_LINES, keep_water=False, verbose=True)
    captured = capsys.readouterr().out
    assert "--strip-heterogens" in captured
    # 4 HETATM + 2 CONECT = 6 dropped.
    assert "6" in captured
