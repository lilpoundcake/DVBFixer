#!/usr/bin/env python3
"""Generate declaration evidence, NOT a resolved SBOM or license clearance."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = "docs/research/license-inventory.json"
CATALOG = "src/dvbfixer/data/installation-profiles.toml"


def tracked_files(root: Path) -> set[str]:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z"],
        check=True, capture_output=True,
    )
    return set(result.stdout.decode().split("\0")) - {""}


def build_inventory(root: Path, tracked: set[str]) -> dict:
    """Read only declared inputs and tracked assets; never scan local prefixes."""
    inputs: dict[str, str] = {}

    def evidence(path: str) -> str:
        data = (root / path).read_bytes()
        inputs[path] = hashlib.sha256(data).hexdigest()
        return data.decode("utf-8")

    project = tomllib.loads(evidence("pyproject.toml"))
    catalog = tomllib.loads(evidence(CATALOG))
    gui = json.loads(evidence("gui/package.json"))
    lock = json.loads(evidence("gui/package-lock.json"))
    project_legal_files = {}
    for name in ("LICENSE", "CONTRIBUTING.md", "CITATION.cff", "gui/LICENSE", "gui/public/LICENSE"):
        if (root / name).is_file():
            evidence(name)
            project_legal_files[name] = inputs[name]
    records = []

    def record(identifier: str, boundary: str, source: str, **fields) -> None:
        records.append({
            "id": identifier, "boundary": boundary, "source": source,
            "license_review": "unknown", "compatibility_review": "unknown",
            **fields,
        })

    for group, requirements in {
        "runtime": project["project"]["dependencies"],
        "build": project["build-system"]["requires"],
        **{f"extra:{key}": value for key, value in
           project["project"].get("optional-dependencies", {}).items()},
    }.items():
        for index, requirement in enumerate(requirements):
            record(f"python:{group}:{index}", f"python:{group}", "pyproject.toml",
                   requirement=requirement, resolution="unresolved")

    def profile_groups(name: str, visiting: tuple[str, ...] = ()) -> set[str]:
        if name in visiting:
            raise ValueError(f"Profile inheritance cycle: {visiting + (name,)}")
        profile = catalog["profiles"][name]
        groups = set(profile.get("groups", []))
        if "inherits" in profile:
            groups |= profile_groups(profile["inherits"], visiting + (name,))
        return groups

    profiles = {name: profile_groups(name) for name in catalog["profiles"]}
    for index, package in enumerate(catalog["packages"]):
        group = package["group"]
        record(f"catalog:{index}", f"catalog:{group}", CATALOG,
               declaration=package, resolution="unresolved",
               profiles=sorted(name for name, groups in profiles.items() if group in groups),
               addons=sorted(name for name, addon in catalog.get("addons", {}).items()
                             if group in addon.get("groups", [])))
    for name, entry in catalog.get("external", {}).items():
        record(f"external:{name}", f"external:{name}", CATALOG,
               declaration=entry, resolution="unresolved")

    for group in ("dependencies", "devDependencies", "optionalDependencies"):
        for name, requirement in gui.get(group, {}).items():
            record(f"gui:{group}:{name}", f"gui:{group}", "gui/package.json",
                   name=name, requirement=requirement, resolution="unresolved")
    for path, package in lock["packages"].items():
        if not path:
            continue
        record(f"npm-lock:{path}", "gui:lockfile-not-shipping-proof",
               "gui/package-lock.json", package_path=path,
               package_metadata={key: package[key] for key in
                                 ("name", "version", "resolved", "integrity", "license",
                                  "dev", "optional", "devOptional", "link", "os", "cpu")
                                 if key in package},
               license_review="declared" if package.get("license") else "unknown",
               distribution_inclusion="unverified")

    # Patterns come from actual setuptools declarations. Only tracked files are
    # approved inputs: never adopt a user's additional local FF tree or notice.
    setuptools = project.get("tool", {}).get("setuptools", {})
    for destination, patterns in setuptools.get("data-files", {}).items():
        for pattern in patterns:
            paths = sorted(path for path in root.glob(pattern)
                           if path.is_file() and path.relative_to(root).as_posix() in tracked)
            if not paths:
                raise ValueError(f"No tracked files for packaged pattern: {pattern}")
            for path in paths:
                relative = path.relative_to(root).as_posix()
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                record(f"data:{destination}:{relative}", "python:shared-data",
                       "pyproject.toml", path=relative, destination=destination,
                       sha256=digest, distribution_inclusion="declared")
    for package, patterns in setuptools.get("package-data", {}).items():
        base = root / "src" / package.replace(".", "/")
        for pattern in patterns:
            paths = sorted(path for path in base.glob(pattern)
                           if path.is_file() and path.relative_to(root).as_posix() in tracked)
            if not paths:
                raise ValueError(f"No tracked files for package data: {package}:{pattern}")
            for path in paths:
                relative = path.relative_to(root).as_posix()
                record(f"package-data:{package}:{relative}", "python:package-data",
                       "pyproject.toml", path=relative,
                       sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                       distribution_inclusion="declared")

    # These files must be reviewed next; hashing YAML does not resolve its deps.
    for path in sorted(tracked):
        candidate = Path(path)
        if (path == "environment.yml" or path == "MANIFEST.in"
                or path.startswith("environments/") and candidate.suffix in {".yml", ".yaml"}
                or path.startswith("deploy/") and (
                    candidate.suffix in {".yml", ".yaml", ".patch"}
                    or candidate.name.startswith("Dockerfile")
                    or candidate.name.endswith("lock.json"))):
            evidence(path)
            record(f"review-input:{path}", "pending-source-review", path,
                   sha256=inputs[path], dependency_enumeration="pending")

    records.sort(key=lambda item: item["id"])
    return {
        "schema_version": 1,
        "kind": "declaration-inventory-not-sbom",
        "project_version": project["project"]["version"],
        "project_license_declaration": project["project"].get("license"),
        "gui_license_declaration": gui.get("license"),
        "project_legal_files_sha256": project_legal_files,
        "limitations": [
            "No installed Python/conda/native/GPU transitive closure or artifact inspection.",
            "npm lock entries are metadata observations, not proof of shipped code or rights.",
            "Untracked files are excluded; actual local build globs can include additional files.",
            "Runtime YAML, image recipes, locks and patches are pending review inputs only.",
            "Copied code, fixtures and downloaded model/CCD/HMM assets need separate inventories.",
            "No upstream legal text is verified and no compatibility decision is made here.",
        ],
        "input_sha256": dict(sorted(inputs.items())),
        "record_count": len(records),
        "records": records,
    }


def render(inventory: dict) -> str:
    return json.dumps(inventory, indent=2, ensure_ascii=False, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check tracked declaration drift only")
    args = parser.parse_args()
    content = render(build_inventory(ROOT, tracked_files(ROOT)))
    output = ROOT / OUTPUT
    if args.check:
        if not output.is_file() or output.read_text(encoding="utf-8") != content:
            print(f"Out of date: {OUTPUT}; run python3 scripts/gen_license_inventory.py")
            return 1
        print("Declaration inventory is current; legal clearance remains unverified")
    else:
        output.write_text(content, encoding="utf-8")
        print(f"Wrote {OUTPUT}; not a resolved SBOM or legal clearance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
