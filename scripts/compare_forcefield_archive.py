#!/usr/bin/env python3
"""Compare tracked parameter files to an archive without extracting or changing them."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path, PurePosixPath

from collect_license_evidence import inspect_archive

ROOT = Path(__file__).resolve().parents[1]


def compare_manifests(local: dict[str, str], upstream: dict[str, str]) -> dict:
    if not local or not upstream:
        raise ValueError("Both comparison manifests must contain files")
    shared = local.keys() & upstream.keys()
    return {
        "matching": sorted(name for name in shared if local[name] == upstream[name]),
        "modified": [{"path": name, "local_sha256": local[name], "upstream_sha256": upstream[name]}
                     for name in sorted(shared) if local[name] != upstream[name]],
        "only_local": sorted(local.keys() - upstream.keys()),
        "only_upstream": sorted(upstream.keys() - local.keys()),
    }


def build_report(root: Path, archive: Path, local_tree: Path, archive_prefix: str,
                 expected_sha256: str) -> dict:
    root = root.resolve()
    local_tree = local_tree.resolve(strict=True)
    if not local_tree.is_relative_to(root) or not local_tree.is_dir():
        raise ValueError("Local tree must be a directory within the repository")
    prefix = PurePosixPath(archive_prefix)
    if prefix.is_absolute() or ".." in prefix.parts or "\\" in archive_prefix or str(prefix) == ".":
        raise ValueError("Expected an explicit safe archive subtree prefix")
    inspection = inspect_archive(archive)
    if inspection["archive_sha256"] != expected_sha256:
        raise ValueError("Archive digest does not match the predeclared identity")
    relative_tree = local_tree.relative_to(root).as_posix()
    tracked = subprocess.run(["git", "-C", str(root), "ls-files", "-z", "--", relative_tree],
                             check=True, capture_output=True).stdout.decode().split("\0")
    local = {}
    for name in sorted(set(tracked) - {""}):
        path = root / name
        if not path.resolve().is_relative_to(local_tree):
            raise ValueError("Tracked parameter path escapes its declared tree")
        # Missing tracked files are errors, not silently skipped observations.
        local[path.relative_to(local_tree).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    archive_root = str(prefix).rstrip("/") + "/"
    upstream = {entry["path"][len(archive_root):]: entry["sha256"]
                for entry in inspection["entries"] if entry["path"].startswith(archive_root)}
    return {
        "schema_version": 1, "kind": "forcefield-archive-content-comparison-not-license-clearance",
        "archive_name": archive.name, "archive_sha256": inspection["archive_sha256"],
        "archive_bytes": archive.stat().st_size, "archive_prefix": str(prefix),
        "local_tree": relative_tree,
        "local_manifest": local, "upstream_manifest": dict(sorted(upstream.items())),
        "outside_archive_prefix": [entry["path"] for entry in inspection["entries"]
                                   if not entry["path"].startswith(archive_root)],
        **compare_manifests(local, upstream),
        "redistribution_clearance": "unknown",
        "limitations": ["Content matches identify a reference, not the owner's original acquisition chain.",
                        "Differences are observations, not errors or authority to replace parameters.",
                        "Untracked local files are excluded; no upstream license is inferred from format."],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--local-tree", type=Path, required=True)
    parser.add_argument("--archive-prefix", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Use a new output path; do not overwrite retained comparisons")
    report = build_report(ROOT, args.archive, args.local_tree, args.archive_prefix, args.expected_sha256)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, sort_keys=True)
        stream.write("\n")
    print(f"Compared {len(report['local_manifest'])} tracked files: "
          f"{len(report['matching'])} matching, {len(report['modified'])} modified; rights remain unverified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
