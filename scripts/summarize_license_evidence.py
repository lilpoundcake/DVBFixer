#!/usr/bin/env python3
"""Build a portable evidence index; full retained texts remain project-local."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from collect_license_evidence import contained_file, digest

ROOT = Path(__file__).resolve().parents[1]
PREFIX_LABELS = ("apple-standard", "apple-desktop", "apple-compat-standard", "apple-compat-desktop",
                 "apple-msa", "apple-psi4", "cc94-mps-sampling", "cc94-cpu-refinement")


def summarize(base: Path, repository: Path) -> dict:
    reports, components = [], {}

    def capture_refs(value, directory: Path) -> None:
        if isinstance(value, dict):
            if "blob" in value:
                blob = contained_file(directory, value["blob"])
                if digest(blob.read_bytes()) != value["sha256"]:
                    raise ValueError("Retained evidence blob digest mismatch")
            for child in value.values():
                capture_refs(child, directory)
        elif isinstance(value, list):
            for child in value:
                capture_refs(child, directory)

    def component(ecosystem: str, entry: dict, label: str) -> None:
        identity = {key: entry.get(key) for key in ("name", "version", "build", "subdir")}
        identity["ecosystem"] = ecosystem
        identity["identity_sha256"] = ((entry.get("sha256") or entry.get("metadata_sha256")) if ecosystem == "conda"
                                       else entry.get("metadata", {}).get("sha256")
                                       or entry.get("metadata_sha256"))
        key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        license_value = entry.get("license_expression") or entry.get("license") or entry.get("license_metadata")
        if not isinstance(license_value, str) or len(license_value) > 200:
            license_value = "see-captured-metadata" if license_value else None
        current = components.setdefault(key, {
            **identity, "declared_license": license_value, "environments": [],
            "legal_text_sha256": [], "license_review": "captured-or-declared-not-reviewed",
            "compatibility_review": "unknown",
        })
        if label not in current["environments"]:
            current["environments"].append(label)
        current["legal_text_sha256"] = sorted(set(current["legal_text_sha256"]) |
                                             {item["sha256"] for item in entry.get("legal_evidence", [])})

    paths = [base / "installed-v4" / label / "inventory.json" for label in PREFIX_LABELS]
    paths += [base / "npm-v3/inventory.json", base / "upstream-v2/inventory.json"]
    paths += [base / "build" / label / "inventory.json" for label in ("inspection-v3-wheel", "inspection-v3-sdist")]
    if not paths or any(not p.is_file() for p in paths):
        raise ValueError("Required retained evidence reports are missing")
    for path in paths:
        data = json.loads(path.read_text())
        capture_refs(data, path.parent)
        relative = path.relative_to(repository).as_posix()
        result = {"path": relative, "sha256": digest(path.read_bytes()), "kind": data["kind"]}
        if "conda" in data:
            label = path.parent.name
            result.update({"label": label, "prefix_label": data["prefix_label"],
                           "site_package_aliases": data.get("site_package_aliases", {}),
                           "conda_records": len(data["conda"]), "python_distributions": len(data["python"]),
                           "conda_with_legal_texts": sum(bool(e["legal_evidence"]) for e in data["conda"]),
                           "python_with_legal_texts": sum(bool(e["legal_evidence"]) for e in data["python"])})
            for entry in data["conda"]:
                component("conda", entry, label)
            for entry in data["python"]:
                component("python-distribution", entry, label)
        elif "packages" in data:
            result.update({"installed_records": len(data["packages"]),
                           "missing_lock_entries": data["missing_lock_entries"],
                           "mismatches": [{key: entry.get(key) for key in
                                           ("name", "version", "locked_version", "package_path")}
                                          for entry in data["packages"] if not entry["matches_locked_version"]]})
            for entry in data["packages"]:
                component("npm-installed", entry, "local-gui-not-lock-clean")
        elif "sources" in data:
            result["sources"] = [{k: v for k, v in entry.items() if k != "blob"}
                                 for entry in data["sources"]]
        elif "entries" in data:
            result.update({key: data[key] for key in
                           ("archive_name", "archive_sha256", "legal_paths", "restricted_or_sensitive_name_findings")})
            result["file_count"] = len(data["entries"])
        reports.append(result)
    return {
        "schema_version": 1, "date": "2026-10-31", "source_baseline": "4a328a6",
        "kind": "captured-evidence-index-not-legal-clearance",
        "limitations": [
            "Conda/Python records overlap; identity count is not a unique-project count.",
            "Installed Python identity uses metadata digest, not the original wheel archive digest.",
            "License strings remain declarations, and all compatibility findings remain unknown.",
            "Linux/CUDA/Protenix environments were not available on this host.",
            "Native dynamic-link closure, unmanaged tools, source delivery and full asset rights remain open.",
            "Full texts/reports are in ignored project-local tmp; index hashes do not replace their retention.",
            "Build snapshot predates the later provenance-only Protpardelle edit; it is not a final release.",
        ],
        "reports": reports,
        "component_identity_count": len(components),
        "components": sorted(components.values(), key=lambda e: (e["ecosystem"], e["name"] or "", e["version"] or "", e["identity_sha256"] or "")),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=ROOT / "tmp/licensing-audit-2026-10-31")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/research/licensing-evidence-2026-10-31.json")
    args = parser.parse_args()
    report = summarize(args.evidence.resolve(), ROOT)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print(f"Indexed {len(report['reports'])} reports and {report['component_identity_count']} component identities; not legal clearance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
