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


def run_header(command: str) -> str:
    """Return the stable header printed for every executable CLI run."""
    from dvbfixer import __version__

    timestamp = datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z")
    return f"=== dvbfixer {__version__} | {timestamp} | {command} ==="


def _copy_stream(read_fd: int, console_fd: int, log: BinaryIO, lock: threading.Lock) -> None:
    """Copy one redirected file descriptor to its console and the shared log."""
    try:
        while chunk := os.read(read_fd, 65536):
            os.write(console_fd, chunk)
            with lock:
                log.write(chunk)
                log.flush()
    finally:
        os.close(read_fd)


@contextmanager
def tee_output(path: str | Path | None) -> Iterator[None]:
    """Mirror fd-level stdout/stderr to ``path`` while retaining the terminal.

    File-descriptor redirection is intentional: subprocesses that inherit the
    normal stdout/stderr descriptors are captured as well as Python ``print``.
    """
    if path is None:
        yield
        return

    log_path = Path(path).expanduser()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    pipes = [os.pipe(), os.pipe()]
    lock = threading.Lock()
    with open(log_path, "ab", buffering=0) as log:
        threads = [
            threading.Thread(
                target=_copy_stream,
                args=(pipes[index][0], saved[index], log, lock),
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
            os.close(saved[0])
            os.close(saved[1])
