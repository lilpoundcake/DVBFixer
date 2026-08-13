"""Process-wide runtime helpers for the unified CLI."""

from __future__ import annotations

import os
import sys
import threading
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import BinaryIO

_RED = b"\x1b[1;97;41m"
_YELLOW = b"\x1b[1;30;43m"
_RESET = b"\x1b[0m"


def run_header(command: str) -> str:
    """Return the stable header printed for every executable CLI run."""
    from dvbfixer import __version__

    timestamp = datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z")
    return f"=== dvbfixer {__version__} | {timestamp} | {command} ==="


def _severity(chunk: bytes) -> str | None:
    upper = chunk.upper()
    if b"ERROR:" in upper or b"ERROR " in upper:
        return "ERROR"
    if b"WARNING:" in upper or b"WARNING " in upper:
        return "WARNING"
    return None


def _copy_stream(
    read_fd: int, console_fd: int, log: BinaryIO | None, lock: threading.Lock,
    findings: list[tuple[str, str]],
) -> None:
    """Copy one redirected file descriptor to its console and the shared log."""
    try:
        pending = b""
        while chunk := os.read(read_fd, 65536):
            pending += chunk
            while b"\n" in pending:
                line, pending = pending.split(b"\n", 1)
                _write_diagnostic_line(line + b"\n", console_fd, log, lock, findings)
        if pending:
            _write_diagnostic_line(pending, console_fd, log, lock, findings)
    finally:
        os.close(read_fd)


def _write_diagnostic_line(
    line: bytes, console_fd: int, log: BinaryIO | None, lock: threading.Lock,
    findings: list[tuple[str, str]],
) -> None:
    severity = _severity(line)
    displayed = line
    logged = line
    if severity:
        message = line.decode("utf-8", errors="replace").strip()
        with lock:
            findings.append((severity, message))
        marker = f"!!! {severity} !!! ".encode()
        logged = marker + line
        if os.isatty(console_fd):
            color = _RED if severity == "ERROR" else _YELLOW
            displayed = color + marker + _RESET + line
        else:
            displayed = marker + line
    os.write(console_fd, displayed)
    if log is not None:
        with lock:
            log.write(logged)  # prominent markers, but no terminal ANSI escapes
            log.flush()


@contextmanager
def tee_output(path: str | Path | None) -> Iterator[None]:
    """Mirror fd-level stdout/stderr to ``path`` while retaining the terminal.

    File-descriptor redirection is intentional: subprocesses that inherit the
    normal stdout/stderr descriptors are captured as well as Python ``print``.
    """
    log_path = Path(path).expanduser() if path is not None else None
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    pipes = [os.pipe(), os.pipe()]
    lock = threading.Lock()
    log_context = open(log_path, "ab", buffering=0) if log_path is not None else None
    try:
        log = log_context
        findings: list[tuple[str, str]] = []
        threads = [
            threading.Thread(
                target=_copy_stream,
                args=(pipes[index][0], saved[index], log, lock, findings),
                daemon=True,
            )
            for index in range(2)
        ]
        try:
            for thread in threads:
                thread.start()
            os.dup2(pipes[0][1], 1)
            os.dup2(pipes[1][1], 2)
            os.close(pipes[0][1])
            os.close(pipes[1][1])
            yield
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            for thread in threads:
                thread.join()
            unique = list(dict.fromkeys(findings))
            errors = sum(level == "ERROR" for level, _ in unique)
            warnings = sum(level == "WARNING" for level, _ in unique)
            summary = [
                "\n" + "!" * 80,
                f"DIAGNOSTIC SUMMARY: {errors} error(s), {warnings} warning(s)",
            ]
            summary.extend(f"  [{level}] {message}" for level, message in unique)
            summary.append("!" * 80)
            summary_bytes = ("\n".join(summary) + "\n").encode("utf-8")
            os.write(saved[1], summary_bytes)
            if log is not None:
                log.write(summary_bytes)
                log.flush()
            os.close(saved[0])
            os.close(saved[1])
    finally:
        if log_context is not None:
            log_context.close()
