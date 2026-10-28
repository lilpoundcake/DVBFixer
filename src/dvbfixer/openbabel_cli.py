"""Expose the wheel's bundled obminimize, which has no upstream console entrypoint."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def obminimize() -> None:
    # Import first: the wheel sets BABEL_LIBDIR/BABEL_DATADIR for its plugins.
    import openbabel

    binary = Path(openbabel.__file__).resolve().parent / "bin" / "obminimize"
    if not binary.is_file():
        raise SystemExit("obminimize is missing from openbabel-wheel; reinstall the Standard profile")
    os.execv(str(binary), [str(binary), *sys.argv[1:]])
