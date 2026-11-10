from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from dvbfixer import __version__
from dvbfixer.batch import add_runtime_help
from dvbfixer.cli import _diagnose_json_stdout
from dvbfixer.runtime import _severity


def test_log_file_tees_python_and_child_streams(tmp_path: Path) -> None:
    log = tmp_path / "run.log"
    script = (
        "from dvbfixer.runtime import tee_output; import subprocess,sys; "
        f"p={str(log)!r}; "
        "ctx=tee_output(p); ctx.__enter__(); "
        "print('python-out', flush=True); print('python-err', file=sys.stderr, flush=True); "
        "subprocess.run([sys.executable, '-c', "
        "\"import sys; print('child-out'); print('child-err', file=sys.stderr)\"]); "
        "ctx.__exit__(None,None,None)"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).parents[1] / "src")},
    )
    combined_terminal = result.stdout + result.stderr
    combined_log = log.read_text()
    assert "DIAGNOSTIC SUMMARY" not in combined_terminal
    assert "DIAGNOSTIC SUMMARY" not in combined_log
    for message in ("python-out", "python-err", "child-out", "child-err"):
        assert message in combined_terminal
        assert message in combined_log


def test_cli_run_header_goes_to_stderr_and_log(tmp_path: Path) -> None:
    pdb = tmp_path / "one.pdb"
    pdb.write_text("END\n")
    log = tmp_path / "cli.log"
    result = subprocess.run(
        [sys.executable, "-m", "dvbfixer.cli", "rename", str(pdb),
         "--log-file", str(log)],
        capture_output=True, text=True, check=True,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).parents[1] / "src")},
    )
    header = f"=== dvbfixer {__version__}"
    assert header in result.stderr
    assert header in log.read_text()


def test_warnings_are_emphasized_and_summarized(tmp_path: Path) -> None:
    script = (
        "from dvbfixer.runtime import tee_output; "
        "ctx=tee_output(None); ctx.__enter__(); "
        "print('WARNING: suspicious structure', flush=True); "
        "ctx.__exit__(None,None,None)"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).parents[1] / "src")},
    )
    combined = result.stdout + result.stderr
    assert "!!! WARNING !!!" in combined
    assert "DIAGNOSTIC SUMMARY: 0 error(s), 1 warning(s)" in combined
    assert "[WARNING] WARNING: suspicious structure" in combined


def test_runtime_help_states_batch_support_for_every_tool() -> None:
    supported = argparse.ArgumentParser()
    add_runtime_help(supported, batch=True)
    supported_help = supported.format_help()
    assert "Batch mode:" in supported_help
    assert "--input-dir DIR" in supported_help

    unsupported = argparse.ArgumentParser()
    add_runtime_help(unsupported)
    unsupported_help = unsupported.format_help()
    assert "Batch mode:" in unsupported_help
    assert "does not support directory batch input" in unsupported_help
    assert "--input-dir" not in unsupported_help


def test_json_presentation_selection_matches_format_option_forms() -> None:
    from dvbfixer.diagnose.cli import parse_args

    for argv in (
        ["in.pdb", "--format", "json"],
        ["--format=json", "in.pdb"],
        ["--for=json", "in.pdb"],
        ["--format", "text", "--format=json", "in.pdb"],
        ["--format=json", "--", "in.pdb"],
    ):
        assert _diagnose_json_stdout(argv)
        assert parse_args(argv).output_format == "json"
    for argv in (
        ["json.pdb"],
        ["--format=json", "--format", "text", "in.pdb"],
        ["--", "--format=json"],
        ["--output=--format=json", "in.pdb"],
    ):
        assert not _diagnose_json_stdout(argv)
        assert parse_args(argv).output_format == "text"
    assert not _diagnose_json_stdout(["in.pdb", "--format"])


def test_negative_chirality_status_is_not_a_runtime_error() -> None:
    assert _severity(b"D-isomer error: NO\n") is None
    assert _severity(b"D-isomer error: NOT CHECKED\n") is None
    assert _severity(b"D-isomer error: YES\n") == "ERROR"
    assert _severity(b"  ERROR   d/ALA82A:CA  D-amino acid C-alpha\n") == "ERROR"
    assert _severity(b"WARNING: inspect stereochemistry\n") == "WARNING"
    assert _severity(b"  ERROR:   0 findings across 0 residues\n") is None
    assert _severity(b"  WARNING: 2 findings across 1 residues\n") is None


def test_structured_stdout_preserves_data_and_child_stderr_diagnostics(tmp_path: Path) -> None:
    log = tmp_path / "structured.log"
    payload = {"message": "warning inside JSON — not a diagnostic line"}
    script = (
        "from dvbfixer.runtime import tee_output; import json,subprocess,sys; "
        f"ctx=tee_output({str(log)!r}, structured_stdout=True); ctx.__enter__(); "
        f"print(json.dumps({payload!r}, ensure_ascii=False), flush=True); "
        "subprocess.run([sys.executable, '-c', "
        "\"import sys; print('WARNING: child diagnostic', file=sys.stderr)\"]); "
        "ctx.__exit__(None,None,None)"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).parents[1] / "src")},
        timeout=60,
    )
    assert json.loads(result.stdout) == payload
    assert "!!! WARNING !!! WARNING: child diagnostic" in result.stderr
    assert "DIAGNOSTIC SUMMARY: 0 error(s), 1 warning(s)" in result.stderr
    assert "inside JSON" not in result.stderr
    logged = log.read_text()
    assert json.dumps(payload, ensure_ascii=False) in logged
    assert "[WARNING] WARNING: child diagnostic" in logged
    assert "\x1b[" not in logged


def test_structured_stdout_preserves_summary_during_exception_cleanup(tmp_path: Path) -> None:
    log = tmp_path / "failure.log"
    payload = {"literal": "ERROR: report data, not an extra runtime error"}
    script = (
        "from dvbfixer.runtime import tee_output\n"
        "import json, sys\n"
        f"with tee_output({str(log)!r}, structured_stdout=True):\n"
        f"    print(json.dumps({payload!r}), flush=True)\n"
        "    print('ERROR: actual failure', file=sys.stderr, flush=True)\n"
        "    raise SystemExit(7)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=False,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).parents[1] / "src")},
        timeout=60,
    )
    assert result.returncode == 7
    assert json.loads(result.stdout) == payload
    assert "!!! ERROR !!! ERROR: actual failure" in result.stderr
    assert "DIAGNOSTIC SUMMARY: 1 error(s), 0 warning(s)" in result.stderr
    assert "report data" not in result.stderr
    assert "[ERROR] ERROR: actual failure" in log.read_text()


def test_structured_stdout_captures_child_unicode_without_final_newline(tmp_path: Path) -> None:
    log = tmp_path / "child.log"
    child = "import json,sys; sys.stdout.write(json.dumps({'message': 'WARNING: αβ' * 12000}, ensure_ascii=False))"
    script = (
        "from dvbfixer.runtime import tee_output; import subprocess,sys; "
        f"ctx=tee_output({str(log)!r}, structured_stdout=True); ctx.__enter__(); "
        f"subprocess.run([sys.executable, '-c', {child!r}], check=True); "
        "ctx.__exit__(None,None,None)"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).parents[1] / "src")},
        timeout=60,
    )
    assert json.loads(result.stdout) == {"message": "WARNING: αβ" * 12000}
    assert log.read_text() == result.stdout
    assert result.stderr == ""
