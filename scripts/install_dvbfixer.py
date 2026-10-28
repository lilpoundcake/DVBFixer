#!/usr/bin/env python3
"""Bootstrap the stdlib-only installer from this source checkout."""

import sys
from pathlib import Path

if sys.version_info < (3, 11):  # noqa: UP036 -- bootstrap runs before project installation
    raise SystemExit("The installer bootstrap requires Python 3.11 or newer; target Python is managed separately.")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dvbfixer.install import main

if __name__ == "__main__":
    raise SystemExit(main())
