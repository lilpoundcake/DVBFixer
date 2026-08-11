from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from dvbfixer.batch import add_runtime_help


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
    assert "=== dvbfixer 0.7.26" in result.stderr
    assert "=== dvbfixer 0.7.26" in log.read_text()


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
