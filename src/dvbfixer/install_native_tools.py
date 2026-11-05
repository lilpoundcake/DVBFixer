"""Keep conda-owned Open Babel executables out of pip entry-point ownership."""

from __future__ import annotations

import contextlib
import csv
import json
import os
import shutil
import stat
import tempfile
import uuid
from collections.abc import Iterator
from pathlib import Path


def _native_obminimize(prefix: Path) -> bool:
    for path in (prefix / "conda-meta").glob("*.json"):
        record = json.loads(path.read_text())
        if record.get("name") != "openbabel":
            continue
        files = record.get("files", [])
        paths = record.get("paths_data", {}).get("paths", [])
        return "bin/obminimize" in files or any(item.get("_path") == "bin/obminimize" for item in paths)
    return False


def _release_pip_ownership(prefix: Path, executable: Path) -> None:
    # The package's wheel-specific wrapper must not claim a restored native
    # executable, or a later pip uninstall would delete another package's file.
    for path in prefix.glob("lib/python*/site-packages/dvbfixer-*.dist-info/RECORD"):
        if path.is_symlink() or not path.resolve().is_relative_to(prefix):
            raise ValueError("Refusing to rewrite linked package ownership metadata")
        with path.open(newline="") as stream:
            rows = list(csv.reader(stream))
        retained = [row for row in rows if not row or Path(os.path.abspath(path.parent.parent / row[0])) != executable]
        if retained == rows:
            continue
        temporary = path.with_name(f".RECORD-{uuid.uuid4().hex}")
        try:
            with temporary.open("x", newline="") as stream:
                csv.writer(stream, lineterminator="\n").writerows(retained)
            temporary.chmod(stat.S_IMODE(path.stat().st_mode))
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)


@contextlib.contextmanager
def preserve_native_openbabel_tool(prefix: Path) -> Iterator[None]:
    """Protect the known console-script collision during a package-only install.

    GUI-only/package-only updates intentionally do not install openbabel-wheel.
    Pip nevertheless creates DVBFixer's obminimize wrapper, which would otherwise
    overwrite the conda binary and potentially its hardlinked package-cache copy.
    """
    prefix = prefix.resolve()
    executable = prefix / "bin/obminimize"
    if not _native_obminimize(prefix) or not (executable.exists() or executable.is_symlink()):
        yield
        return
    directory = Path(tempfile.mkdtemp(prefix="dvbfixer-native-tool-"))
    backup = directory / "obminimize"
    linked = executable.is_symlink()
    link = os.readlink(executable) if linked else None
    if linked:
        backup.symlink_to(link)
    else:
        shutil.copy2(executable, backup)
    # Break the native/cache hardlink before pip's script writer can truncate it.
    executable.unlink()
    print("Preserving conda-owned obminimize; its executable will not be replaced by the wheel wrapper.", flush=True)
    retained = False
    try:
        yield
    finally:
        temporary = executable.with_name(f".obminimize-{uuid.uuid4().hex}")
        try:
            if executable.exists() or executable.is_symlink():
                original = executable.is_symlink() and linked and os.readlink(executable) == link
                if not linked and not executable.is_symlink() and executable.is_file():
                    original = executable.read_bytes() == backup.read_bytes()
                if not original:
                    wrapper = (not executable.is_symlink() and executable.is_file()
                               and executable.stat().st_size < 65536
                               and b"dvbfixer.openbabel_cli" in executable.read_bytes())
                    if not wrapper:
                        raise ValueError("obminimize changed unexpectedly during installation; refusing to overwrite it")
            else:
                original = False
            if not original:
                if linked:
                    temporary.symlink_to(link)
                else:
                    shutil.copy2(backup, temporary)
                temporary.replace(executable)
            _release_pip_ownership(prefix, executable)
        except BaseException:
            retained = True
            print(f"Native executable backup retained for recovery: {backup}", flush=True)
            raise
        finally:
            temporary.unlink(missing_ok=True)
            if not retained:
                backup.unlink()
                directory.rmdir()
