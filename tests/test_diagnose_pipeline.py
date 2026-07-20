"""End-to-end tests for `dvbfixer diagnose` pipeline.

Uses the shipped fixture `test/broken_SER/SER.pdb` (skipped if the
fixture isn't provisioned locally, matching the pattern used by other
tests) plus synthetic in-test inputs so CI can exercise the pipeline
without a checked-in PDB.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("openmm", reason="diagnose needs OpenMM")
pytest.importorskip("pdbfixer", reason="diagnose needs PDBFixer")

from dvbfixer.diagnose import main as diagnose_main  # noqa: E402

# Constructed C-terminal SER with HG placed exactly on OXT — reproduces
# the reported broken_SER/SER.pdb without needing that fixture.
_BROKEN_HG_ON_OXT = """\
ATOM      1  N   SER B 126     -20.246  21.928 -25.860  1.00  0.00           N
ATOM      2  H   SER B 126     -20.762  22.238 -25.051  1.00  0.00           H
ATOM      3  CA  SER B 126     -20.916  22.100 -27.152  1.00  0.00           C
ATOM      4  HA  SER B 126     -20.275  21.751 -27.937  1.00  0.00           H
ATOM      5  C   SER B 126     -21.114  23.596 -27.370  1.00  0.00           C
ATOM      6  O   SER B 126     -20.149  24.343 -27.634  1.00  0.00           O
ATOM      7  CB  SER B 126     -22.268  21.363 -27.174  1.00  0.00           C
ATOM      8  HB2 SER B 126     -22.692  21.265 -28.187  1.00  0.00           H
ATOM      9  HB3 SER B 126     -22.069  20.334 -26.888  1.00  0.00           H
ATOM     10  OG  SER B 126     -23.160  21.914 -26.218  1.00  0.00           O
ATOM     11  HG  SER B 126     -22.356  23.103 -27.085  1.00  0.00           H
ATOM     12  OXT SER B 126     -22.355  23.104 -27.085  1.00  0.00           O
TER
END
"""


# Constructed well-formed single ALA — should be clean.
_CLEAN_ALA = """\
ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  H   ALA A   1      -1.930  -0.891   0.000  1.00  0.00           H
ATOM      3  H2  ALA A   1      -1.930   0.446   0.771  1.00  0.00           H
ATOM      4  H3  ALA A   1      -1.930   0.446  -0.771  1.00  0.00           H
ATOM      5  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      6  HA  ALA A   1       0.360   0.502   0.895  1.00  0.00           H
ATOM      7  CB  ALA A   1       0.508  -1.435   0.000  1.00  0.00           C
ATOM      8  HB1 ALA A   1       0.148  -1.936   0.895  1.00  0.00           H
ATOM      9  HB2 ALA A   1       1.598  -1.435   0.000  1.00  0.00           H
ATOM     10  HB3 ALA A   1       0.148  -1.936  -0.895  1.00  0.00           H
ATOM     11  C   ALA A   1       0.556   0.762  -1.209  1.00  0.00           C
ATOM     12  O   ALA A   1       0.010   0.702  -2.297  1.00  0.00           O
ATOM     13  OXT ALA A   1       1.639   1.435  -1.121  1.00  0.00           O
TER
END
"""


def _run(argv: list[str]) -> int:
    """Run diagnose and return its exit code (converts SystemExit)."""
    try:
        diagnose_main(argv)
        return 0
    except SystemExit as e:
        return int(e.code) if e.code is not None else 0


def test_broken_hg_on_oxt_reports_errors_and_exits_1(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The reported failure mode must be flagged with ERROR severity."""
    in_pdb = tmp_workdir / "broken.pdb"
    in_pdb.write_text(_BROKEN_HG_ON_OXT)
    exit_code = _run([str(in_pdb)])
    out = capsys.readouterr().out
    assert exit_code == 1, out
    assert "ERROR" in out
    assert "coincident_atoms" in out or "coincident with OXT" in out


def test_clean_alanine_zero_findings_exits_0(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    in_pdb = tmp_workdir / "clean.pdb"
    in_pdb.write_text(_CLEAN_ALA)
    exit_code = _run([str(in_pdb)])
    out = capsys.readouterr().out
    # We expect no ERROR-level findings on this trivial single-residue
    # input. WARNING-level clashes on backbone atoms are OK (H2/H3
    # near existing H is a known artefact of the fake terminal geometry).
    assert exit_code == 0, out


def test_only_flag_restricts_categories(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    in_pdb = tmp_workdir / "broken.pdb"
    in_pdb.write_text(_BROKEN_HG_ON_OXT)
    # Only structural checks → no clash / valence findings in output.
    _run([str(in_pdb), "--only", "structural"])
    out = capsys.readouterr().out
    assert "Structural integrity" in out
    assert "Steric analysis" not in out
    assert "Chemistry" not in out


def test_severity_filter_ERROR_drops_lower(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    in_pdb = tmp_workdir / "broken.pdb"
    in_pdb.write_text(_BROKEN_HG_ON_OXT)
    _run([str(in_pdb), "--severity", "ERROR"])
    out = capsys.readouterr().out
    # Every reported finding line should be ERROR (no WARNING/INFO).
    for ln in out.splitlines():
        stripped = ln.strip()
        if stripped.startswith(("WARNING ", "INFO ")):
            pytest.fail(f"WARNING/INFO leaked through --severity ERROR: {ln}")


def test_missing_input_exits_2(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    exit_code = _run([str(tmp_workdir / "does-not-exist.pdb")])
    assert exit_code == 2
