#!/usr/bin/env python3
"""Render/check profile specs and the dependency matrix from one catalog."""

from __future__ import annotations

import argparse
from pathlib import Path

from dvbfixer.installation_profiles import dependency_matrix, environment_spec, load_catalog

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    outputs = {
        "environments/standard.yml": environment_spec("standard"),
        "environments/desktop.yml": environment_spec("desktop"),
        "environments/standard-msa-alternatives.yml": environment_spec("standard", ("msa-alternatives",)),
        "environments/standard-probe.yml": environment_spec("standard", ("probe",)),
        "docs/installation-dependency-matrix.md": dependency_matrix(),
    }
    psi4 = load_catalog()["external"]["psi4"]
    outputs["environments/psi4.yml"] = (
        "# Generated; separate runtime, do not install dvbfixer here.\nname: dvbfixer-psi4\nchannels:\n"
        + "".join(f"  - {channel}\n" for channel in psi4["channels"])
        + "  - nodefaults\ndependencies:\n"
        + "".join(f"  - {spec}\n" for spec in psi4["dependencies"])
    )
    failed = False
    for name, content in outputs.items():
        path = ROOT / name
        if args.check:
            if not path.is_file() or path.read_text() != content:
                print(f"Out of date: {name}")
                failed = True
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
