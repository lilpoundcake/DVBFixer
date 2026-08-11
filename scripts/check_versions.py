#!/usr/bin/env python
"""Ensure release metadata matches the authoritative pyproject version."""

from __future__ import annotations

import json
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    init = (ROOT / "src/dvbfixer/__init__.py").read_text()
    match = re.search(r'^__version__\s*=\s*["\']([^"\']+)', init, re.MULTILINE)
    package_lock = json.loads((ROOT / "gui/package-lock.json").read_text())
    changelog = (ROOT / "CHANGELOG.md").read_text()
    changelog_match = re.search(r"^## \[([^]]+)]", changelog, re.MULTILINE)
    values = {
        "src/dvbfixer/__init__.py": match.group(1) if match else None,
        "gui/package.json": json.loads((ROOT / "gui/package.json").read_text())["version"],
        "gui/package-lock.json": package_lock["version"],
        'gui/package-lock.json packages[""]': package_lock.get("packages", {})
        .get("", {})
        .get("version"),
        "CHANGELOG.md latest release": changelog_match.group(1) if changelog_match else None,
    }
    mismatches = {file: value for file, value in values.items() if value != version}
    if mismatches:
        for file, value in mismatches.items():
            print(f"{file}: {value!r}; expected {version!r}")
        return 1
    print(f"release metadata synchronized at {version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
