"""Stdlib-only installation stage reporting; never pass license subprocesses here."""

from __future__ import annotations

import codecs
import contextlib
import subprocess
import threading
import time
from collections.abc import Iterator
from contextvars import ContextVar
from dataclasses import dataclass, field
from pathlib import Path

from dvbfixer import install_result


@dataclass
class _Progress:
    stack: list[object] = field(default_factory=list)
    lock: threading.RLock = field(default_factory=threading.RLock)
    waiting: int = 0
    wait_started: float = 0
    waited: float = 0

    def waiting_time(self) -> float:
        return self.waited + (time.monotonic() - self.wait_started if self.waiting else 0)


_current: ContextVar[_Progress | None] = ContextVar("installation_progress", default=None)


@contextlib.contextmanager
def waiting_for_input() -> Iterator[None]:
    """Pause all enclosing stages before rendering a prompt, not just its input()."""
    progress = _current.get()
    if progress is None:
        yield
        return
    with progress.lock:
        if not progress.waiting:
            progress.wait_started = time.monotonic()
        progress.waiting += 1
    try:
        yield
    finally:
        with progress.lock:
            progress.waiting -= 1
            if not progress.waiting:
                progress.waited += time.monotonic() - progress.wait_started


@contextlib.contextmanager
def stage(label: str, log: Path | None = None, *, interval: float = 10,
          follow: bool = True) -> Iterator[None]:
    progress = _current.get() or _Progress()
    token = _current.set(progress)
    frame = object()
    start = time.monotonic()
    with progress.lock:
        waited_at_start = progress.waiting_time()
        progress.stack.append(frame)
    install_result.stage_started(label, log)
    print(f"\nStarting: {label}", flush=True)
    if log is not None:
        print(f"Log: {log}", flush=True)
    done = threading.Event()
    offset = 0
    decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")

    def drain() -> None:
        nonlocal offset
        if not follow or log is None or not log.exists():
            return
        with log.open("rb") as stream:
            stream.seek(offset)
            while chunk := stream.read(65536):
                offset += len(chunk)
                print(decoder.decode(chunk), end="", flush=True)

    def heartbeat() -> None:
        last_notice = start
        while not done.wait(min(interval, 0.25)):
            with progress.lock:
                if progress.waiting or progress.stack[-1] is not frame:
                    continue
                drain()
                now = time.monotonic() - progress.waiting_time() + waited_at_start
                if now - last_notice >= interval:
                    print(f"  Still running: {label} ({now - start:.0f}s elapsed)", flush=True)
                    last_notice = now

    worker = threading.Thread(target=heartbeat, daemon=True)
    worker.start()
    try:
        yield
    except BaseException:
        install_result.stage_failed(label, log)
        done.set()
        worker.join()
        drain()
        print(f"Failed/interrupted: {label} ({time.monotonic() - start - progress.waiting_time() + waited_at_start:.1f}s active)", flush=True)
        raise
    else:
        done.set()
        worker.join()
        drain()
        print(f"Completed: {label} ({time.monotonic() - start - progress.waiting_time() + waited_at_start:.1f}s active)", flush=True)
    finally:
        done.set()
        worker.join()
        with progress.lock:
            progress.stack.remove(frame)
        _current.reset(token)


def logged_run(command: list[str], *, log: Path, env: dict[str, str],
               cwd: Path | None = None, timeout: float | None = None,
               allowed_codes: tuple[int, ...] = (0,)) -> subprocess.CompletedProcess:
    """Retain full command output; show activity even for silent package solves."""
    with stage(log.stem, log), log.open("w") as stream:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, timeout=timeout, check=False)
        if result.returncode not in allowed_codes:
            raise subprocess.CalledProcessError(result.returncode, command)
        return result
