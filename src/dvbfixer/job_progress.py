"""Best-effort, versioned monitoring events; never scientific acceptance evidence."""

from __future__ import annotations

import json
import os
import stat
from typing import Literal

JobStage = Literal["preparation", "runner", "sampling", "refinement", "validation", "publication"]


def job_progress_descriptor() -> int | None:
    """Return only a live pipe/socket descriptor explicitly owned by the host."""
    try:
        fd = int(os.environ.get("DVBFIXER_JOB_PROGRESS_FD", "-1"))
        if fd >= 3:
            mode = os.fstat(fd).st_mode
            if stat.S_ISFIFO(mode) or stat.S_ISSOCK(mode):
                return fd
    except (OSError, ValueError, OverflowError):
        pass
    return None


def emit_job_progress(
    stage: JobStage, *, completed: int | None = None, total: int | None = None, seed: int | None = None,
) -> None:
    """Write one small NDJSON event to the host-owned pipe, if configured.

    Counts describe this stage only. No counts means indeterminate work. A closed
    monitor must not interrupt a scientific calculation or contaminate stdout.
    """
    if stage not in {"preparation", "runner", "sampling", "refinement", "validation", "publication"}:
        raise ValueError("unknown job stage")
    if (completed is None) != (total is None):
        raise ValueError("completed and total must be supplied together")
    if completed is not None and total is not None:
        if type(completed) is not int or type(total) is not int or not 0 <= completed <= total or total < 1:
            raise ValueError("invalid stage work counts")
        if total > 2**53 - 1:
            raise ValueError("stage work counts exceed the wire integer range")
    if seed is not None and (type(seed) is not int or not 0 <= seed <= 2**53 - 1):
        raise ValueError("seed must be a safe non-negative integer")
    fd = job_progress_descriptor()
    if fd is None:
        return
    event: dict[str, object] = {"version": 1, "stage": stage}
    if completed is not None:
        event.update(completed=completed, total=total)
    if seed is not None:
        event["seed"] = seed
    try:
        os.set_blocking(fd, False)
        os.write(fd, (json.dumps(event, separators=(",", ":")) + "\n").encode("ascii"))
    except (OSError, ValueError, OverflowError):
        # Monitoring is optional, including during host shutdown/cancellation.
        return
