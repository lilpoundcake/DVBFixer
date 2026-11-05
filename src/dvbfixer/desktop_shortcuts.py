"""Opt-in, per-user Desktop launchers bound to one installation prefix."""

from __future__ import annotations

import hashlib
import os
import plistlib
import shlex
import sys
from pathlib import Path

from dvbfixer.gui_launcher import GuiLaunchError, installed_gui


def _exclusive_file(path: Path, content: bytes, mode: int = 0o644) -> None:
    """Idempotent exact reuse only; never replace a user's existing launcher."""
    if path.is_symlink():
        raise GuiLaunchError(f"Refusing shortcut symlink: {path}")
    try:
        with path.open("xb") as stream:
            os.fchmod(stream.fileno(), mode)
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError:
        if not path.is_file() or path.read_bytes() != content:
            raise GuiLaunchError(f"Shortcut target already exists with different content: {path}; it was not overwritten") from None


def _desktop_quote(value: str) -> str:
    # Desktop Entry Exec grammar is not shell grammar. Percent is a field code.
    value = value.replace("%", "%%")
    quoted = "".join("\\" + char if char in '\\"`$' else char for char in value)
    # General Desktop Entry string unescaping precedes Exec argument parsing.
    return '"' + quoted.replace("\\", "\\\\") + '"'


def install_desktop_shortcut(prefix: Path, *, home: Path | None = None, platform: str | None = None) -> Path:
    prefix = prefix.resolve()
    installed_gui(prefix)
    python = prefix / "bin/python"
    if not os.access(python, os.X_OK):
        raise GuiLaunchError(f"Desktop Python is unavailable: {python}")
    if any(char in str(prefix) for char in ("\n", "\r", "\0")):
        raise GuiLaunchError("Shortcut paths must not contain newlines or NUL")
    selected = platform or sys.platform
    if selected not in {"darwin", "linux"}:
        raise GuiLaunchError("Application shortcuts currently support macOS and Linux")
    home = home or Path.home()
    key = hashlib.sha256(str(prefix).encode()).hexdigest()[:12]
    name = f"DVBFixer ({prefix.name})"
    script = (
        "#!/bin/sh\n# DVBFIXER_DESKTOP_LAUNCHER_V1\n"
        f"exec {shlex.quote(str(python))} -m dvbfixer.cli gui \"$@\"\n"
    ).encode()
    if selected == "darwin":
        app = home / "Applications" / f"{name}-{key}.app"
        if app.is_symlink():
            raise GuiLaunchError(f"Refusing shortcut symlink: {app}")
        executable = app / "Contents/MacOS/DVBFixer"
        info = plistlib.dumps({
            "CFBundleExecutable": "DVBFixer", "CFBundleIdentifier": f"org.dvbfixer.desktop.{key}",
            "CFBundleName": name, "CFBundleDisplayName": name, "CFBundlePackageType": "APPL",
            "CFBundleVersion": "1", "LSUIElement": True,
        }, sort_keys=True)
        if app.exists():
            # Do not write through user-added symlinks or adopt an unrelated app.
            if (not app.is_dir() or any((app / path).is_symlink() for path in ("Contents", "Contents/MacOS", "Contents/Info.plist"))
                    or executable.is_symlink() or not executable.is_file()
                    or not (app / "Contents/Info.plist").is_file()
                    or executable.read_bytes() != script or (app / "Contents/Info.plist").read_bytes() != info):
                raise GuiLaunchError(f"Application already exists and was not overwritten: {app}")
            return app
        app.parent.mkdir(parents=True, exist_ok=True)
        app.mkdir()  # Exclusive reservation; never rename over an existing directory.
        executable.parent.mkdir(parents=True)
        _exclusive_file(executable, script, 0o755)
        _exclusive_file(app / "Contents/Info.plist", info)
        return app
    script_path = prefix / "share/dvbfixer/desktop-launcher"
    _exclusive_file(script_path, script, 0o755)
    data_home = Path(os.environ.get("XDG_DATA_HOME", str(home / ".local/share"))).expanduser()
    if not data_home.is_absolute():
        raise GuiLaunchError("XDG_DATA_HOME must be absolute")
    target = data_home / "applications" / f"dvbfixer-{key}.desktop"
    target.parent.mkdir(parents=True, exist_ok=True)
    desktop_name = name.replace("\\", "\\\\")
    entry = (
        "[Desktop Entry]\nType=Application\nVersion=1.0\n"
        f"Name={desktop_name}\nComment=Launch the local DVBFixer Desktop GUI\n"
        f"Exec={_desktop_quote(str(script_path))}\nTerminal=false\nCategories=Science;Education;\n"
        "Actions=Stop;\n\n[Desktop Action Stop]\nName=Stop DVBFixer\n"
        f"Exec={_desktop_quote(str(script_path))} --stop\n"
    )
    _exclusive_file(target, entry.encode(), 0o755)
    return target
