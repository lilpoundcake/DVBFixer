from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


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
    assert "=== dvbfixer 0.7.23" in result.stderr
    assert "=== dvbfixer 0.7.23" in log.read_text()

