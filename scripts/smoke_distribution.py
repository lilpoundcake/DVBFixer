#!/usr/bin/env python
"""Smoke-test an installed dvbfixer wheel without installing runtime dependencies."""

from __future__ import annotations

import argparse
import importlib.metadata
import subprocess
import sys
import sysconfig
import tomllib
from pathlib import Path

REQUIRED_FF_FILES = (
    "aminoacids.rtp",
    "ffbonded.itp",
    "ffnonbonded.itp",
    "forcefield.itp",
)
EXPECTED_FF_DIRS = (
    "amber99sb-ildn-lipid21.ff",
    "charmm36_ljpme-jul2022.ff",
)


def _run_cli(*args: str) -> str:
    scripts_dir = Path(sysconfig.get_path("scripts"))
    executable = scripts_dir / ("dvbfixer.exe" if sys.platform == "win32" else "dvbfixer")
    if not executable.is_file():
        raise RuntimeError(f"installed console entry point is missing: {executable}")
    completed = subprocess.run(
        [str(executable), *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        required=True,
        help="Checkout root containing the authoritative pyproject.toml",
    )
    args = parser.parse_args()

    expected_version = tomllib.loads(
        (args.project_root / "pyproject.toml").read_text()
    )["project"]["version"]

    import dvbfixer
    from dvbfixer.top.cli import FF_CHOICES, bundled_ff_root

    installed_version = importlib.metadata.version("dvbfixer")
    versions = {
        "distribution metadata": installed_version,
        "dvbfixer.__version__": dvbfixer.__version__,
    }
    mismatches = {name: value for name, value in versions.items() if value != expected_version}
    if mismatches:
        details = ", ".join(f"{name}={value!r}" for name, value in mismatches.items())
        raise RuntimeError(f"installed version does not match {expected_version!r}: {details}")

    version_output = _run_cli("--version").strip()
    if version_output != f"dvbfixer {expected_version}":
        raise RuntimeError(f"unexpected --version output: {version_output!r}")
    help_output = _run_cli("--help")
    for expected_text in ("Usage: dvbfixer <command> [options]", "zbs", "top"):
        if expected_text not in help_output:
            raise RuntimeError(f"--help output is missing {expected_text!r}")

    ff_root = bundled_ff_root().resolve()
    expected_root = (Path(sys.prefix) / "share" / "dvbfixer" / "FF").resolve()
    if ff_root != expected_root:
        raise RuntimeError(f"force fields resolved to {ff_root}, expected installed root {expected_root}")
    if set(FF_CHOICES.values()) != set(EXPECTED_FF_DIRS):
        raise RuntimeError(f"unexpected force-field registry: {FF_CHOICES!r}")
    missing = [
        str(ff_root / directory / filename)
        for directory in EXPECTED_FF_DIRS
        for filename in REQUIRED_FF_FILES
        if not (ff_root / directory / filename).is_file()
    ]
    if missing:
        raise RuntimeError("wheel is missing bundled force-field files: " + ", ".join(missing))

    print(f"installed dvbfixer {expected_version} distribution smoke check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
