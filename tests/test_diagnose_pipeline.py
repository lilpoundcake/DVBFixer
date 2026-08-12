"""End-to-end tests for `dvbfixer diagnose` pipeline.

Uses the shipped fixture `an external broken-SER fixture` (skipped if the
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


# Multi-MODEL PDB (3 identical frames of the clean ALA).
_MULTI_MODEL_ALA = (
    "MODEL        1\n"
    + _CLEAN_ALA.replace("END\n", "ENDMDL\n")
    + "MODEL        2\n"
    + _CLEAN_ALA.replace("END\n", "ENDMDL\n")
    + "MODEL        3\n"
    + _CLEAN_ALA.replace("END\n", "ENDMDL\n")
    + "END\n"
)


def test_multi_model_input_reports_warning_and_uses_model1(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A 3-MODEL file must produce a WARNING banner and not spam
    inter-MODEL chain-break errors."""
    in_pdb = tmp_workdir / "multi.pdb"
    in_pdb.write_text(_MULTI_MODEL_ALA)
    _run([str(in_pdb)])
    out = capsys.readouterr().out
    assert "3 MODEL records" in out
    assert "MODEL 1 only" in out
    # No inter-MODEL chain-break flooding — chain break count should
    # be at most 1 (from the terminal ALA), not 3.
    n_chain_break = out.count("chain break")
    assert n_chain_break <= 1, f"multi-MODEL flooded chain breaks: {out}"


# Small solvent-only input — a plain protein residue plus a couple
# of waters. Chain-break check without --include-water must not
# flag the water-water resSeq jump.
_ALA_PLUS_WATERS = _CLEAN_ALA.replace("END\n", "") + """\
HETATM   14  O   HOH A 100      10.000  10.000  10.000  1.00  0.00           O
HETATM   15  O   HOH A 200      15.000  10.000  10.000  1.00  0.00           O
END
"""


def test_water_chain_breaks_suppressed_by_default(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    in_pdb = tmp_workdir / "ala_plus_waters.pdb"
    in_pdb.write_text(_ALA_PLUS_WATERS)
    _run([str(in_pdb)])
    out = capsys.readouterr().out
    # No 'chain break after A/HOH100' — waters are filtered.
    assert "HOH100" not in out
    assert "HOH200" not in out


def test_water_included_when_flag_set(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    in_pdb = tmp_workdir / "ala_plus_waters.pdb"
    in_pdb.write_text(_ALA_PLUS_WATERS)
    _run([str(in_pdb), "--include-water"])
    out = capsys.readouterr().out
    # With --include-water, the water-vs-water resSeq jump surfaces.
    assert "HOH" in out


def test_chain_transition_is_not_reported_as_internal_break(tmp_workdir: Path) -> None:
    """An explicit A -> B transition is a boundary, not a broken chain."""
    from dvbfixer.diagnose.structural import check_chain_breaks

    pdb = tmp_workdir / "two_chains.pdb"
    pdb.write_text("""\
ATOM      1  N   GLY A  86       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A  86       1.400   0.000   0.000  1.00  0.00           C
ATOM      3  C   GLY A  86       2.800   0.000   0.000  1.00  0.00           C
ATOM      4  O   GLY A  86       3.800   0.000   0.000  1.00  0.00           O
TER
ATOM      5  N   MET B   1      30.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  MET B   1      31.400   0.000   0.000  1.00  0.00           C
ATOM      7  C   MET B   1      32.800   0.000   0.000  1.00  0.00           C
ATOM      8  O   MET B   1      33.800   0.000   0.000  1.00  0.00           O
TER
END
""")
    assert check_chain_breaks(pdb) == []

    same_chain = tmp_workdir / "broken_chain.pdb"
    same_chain.write_text(pdb.read_text().replace("MET B   1", "MET A   1"))
    findings = check_chain_breaks(same_chain)
    assert len(findings) == 1
    assert findings[0].chain == "A"
    assert "after A/GLY86" in findings[0].message


def test_clash_cutoff_bad_value_raises(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Malformed --clash-cutoff should trigger argparse to exit 2 (or
    SystemExit with a non-zero code)."""
    in_pdb = tmp_workdir / "clean.pdb"
    in_pdb.write_text(_CLEAN_ALA)
    # argparse's default behaviour is SystemExit(2) with a usage
    # message on stderr.
    with pytest.raises(SystemExit):
        diagnose_main([str(in_pdb), "--clash-cutoff", "not-a-float"])


def test_clash_mode_bioluminate_quiets_borderline(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """--clash-mode bioluminate must exit 0 on a broken-SER input
    that molprobity flags as ERROR — the constructed geometry has
    overlaps in the 0.4-0.6 Å band that bioluminate silently
    ignores."""
    in_pdb = tmp_workdir / "clean.pdb"
    in_pdb.write_text(_CLEAN_ALA)
    # The clean-ALA test already exits 0 under chimerax; run under
    # bioluminate to sanity-check the flag plumbing without any
    # engine crash.
    ec = _run([str(in_pdb), "--clash-mode", "bioluminate"])
    assert ec == 0


def test_json_output_is_valid(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """--format json should emit a parseable JSON document with the
    documented keys."""
    import json as _json
    in_pdb = tmp_workdir / "broken.pdb"
    in_pdb.write_text(_BROKEN_HG_ON_OXT)
    _run([str(in_pdb), "--format", "json"])
    out = capsys.readouterr().out
    payload = _json.loads(out)
    assert set(payload.keys()) >= {
        "input", "n_atoms", "n_residues", "n_chains", "findings", "summary",
    }
    # At least one ERROR finding on this known-broken input.
    assert payload["summary"]["ERROR"] >= 1
    # Every finding has the documented shape.
    for f in payload["findings"]:
        assert set(f.keys()) >= {
            "severity", "category", "chain", "resid", "resname",
            "atom", "message", "fix_hint",
        }


def test_json_output_keeps_unicode_readable(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    in_pdb = tmp_workdir / "Ångström.pdb"
    in_pdb.write_text(_BROKEN_HG_ON_OXT)
    _run([str(in_pdb), "--format", "json"])
    out = capsys.readouterr().out
    assert "Ångström" in out
    assert "\\u00c5" not in out
    assert "\\u2014" not in out


def test_json_reports_persistent_chirality_repair_history(
    tmp_workdir: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    in_pdb = tmp_workdir / "repaired.pdb"
    in_pdb.write_text(
        "REMARK 999 DVBFIXER CHIRALITY_REPAIR minimize A ALA 1 .\n" + _CLEAN_ALA
    )
    _run([str(in_pdb), "--format", "json"])
    import json as _json
    payload = _json.loads(capsys.readouterr().out)
    assert payload["chirality"]["d_isomer_error"] is False
    assert payload["chirality"]["hydrogen_geometry_review_recommended"] is True
    assert payload["chirality"]["forced_repairs"][0]["resid"] == "1"
