from pathlib import Path

import pytest

from dvbfixer.zbs import PostflightError, _run_pipeline, parse_args


def _args(source: Path, output: Path, *extra: str):
    return parse_args([
        str(source), "-o", str(output), "--skip-renumber", "--skip-model",
        "--skip-prepare", "--skip-minimize", "--no-align-to-input", *extra,
    ])


def test_postflight_errors_warn_by_default(monkeypatch, tmp_path: Path):
    source = tmp_path / "input.pdb"
    source.write_text("END\n")
    output = tmp_path / "output.pdb"
    monkeypatch.setattr("dvbfixer.diagnose.main", lambda argv: (_ for _ in ()).throw(SystemExit(1)))
    _run_pipeline(_args(source, output), source)
    assert output.exists()


def test_strict_postflight_errors_fail(monkeypatch, tmp_path: Path):
    source = tmp_path / "input.pdb"
    source.write_text("END\n")
    output = tmp_path / "output.pdb"
    monkeypatch.setattr("dvbfixer.diagnose.main", lambda argv: (_ for _ in ()).throw(SystemExit(1)))
    with pytest.raises(PostflightError, match="postflight diagnose failed"):
        _run_pipeline(_args(source, output, "--strict-postflight"), source)
