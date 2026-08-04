"""Safe temporary-path helpers used by structure pipelines."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path


def make_temp_path(*, suffix: str = "", prefix: str = "dvbfixer_") -> Path:
    """Reserve and return a unique path without the race in ``mktemp``."""
    descriptor, name = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    os.close(descriptor)
    return Path(name)
