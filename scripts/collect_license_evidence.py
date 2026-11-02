#!/usr/bin/env python3
"""Read installed package metadata/legal files without importing package code.

Outputs are local audit evidence, not permission to redistribute the environment.
Never inspect activation scripts, MODELLER config or registration keys.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import io
import json
import re
import tarfile
import zipfile
from email.parser import Parser
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit, urlunsplit

LEGAL_NAME = re.compile(r"^(licen[sc]e|copying|copyright|notice|authors)([._-].*)?$", re.I)
NONLEGAL_SUFFIXES = {".py", ".pyc", ".pyo", ".js", ".ts", ".c", ".h", ".cpp",
                     ".so", ".dylib", ".dll", ".png", ".jpg", ".ico", ".svg"}


def legal_name(path: Path) -> bool:
    if re.search(r"(^|[._-])(key|token|secret)([._-]|$)", path.name, re.I):
        return False
    named_license = LEGAL_NAME.match(path.name) or re.match(r"^(GNU_)?L?GPL([._-].*)?$", path.name, re.I)
    return bool(named_license) and path.suffix.lower() not in NONLEGAL_SUFFIXES


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def public_url(value: str) -> str:
    """Do not retain URL credentials, query tokens or private conda token paths."""
    parts = urlsplit(value)
    if parts.scheme not in {"http", "https"}:
        return "local-or-non-http-source"
    if re.search(r"/t/[^/]+", parts.path):
        return "redacted-token-channel"
    return urlunsplit((parts.scheme, parts.netloc.rsplit("@", 1)[-1], parts.path, "", ""))


def contained_file(base: Path, relative: str) -> Path:
    path = base / relative
    if not path.resolve().is_relative_to(base.resolve()):
        raise ValueError(f"Evidence path escapes its root: {relative}")
    return path


class EvidenceStore:
    def __init__(self, output: Path):
        self.output = output
        (output / "blobs").mkdir(parents=True, exist_ok=True)

    def capture(self, path: Path, origin: str) -> dict:
        data = path.read_bytes()
        sha = digest(data)
        target = self.output / "blobs" / sha
        if target.exists():
            if target.read_bytes() != data:
                raise ValueError("Existing evidence blob has conflicting content")
        else:
            with target.open("xb") as stream:
                stream.write(data)
        return {"origin": origin, "sha256": sha, "bytes": len(data),
                "blob": f"blobs/{sha}"}


def legal_files(root: Path) -> list[Path]:
    return sorted(p for p in root.rglob("*")
                  if p.is_file() and legal_name(p)
                  and p.resolve().is_relative_to(root.resolve()))


def collect_prefix(prefix: Path, store: EvidenceStore) -> dict:
    prefix = prefix.resolve(strict=True)
    meta = prefix / "conda-meta"
    if not meta.is_dir():
        raise ValueError("Expected an existing conda prefix with conda-meta")
    conda, python = [], []
    for path in sorted(meta.glob("*.json")):
        raw = path.read_bytes()
        package = json.loads(raw)
        entry = {key: package[key] for key in
                 ("name", "version", "build", "build_number", "subdir", "sha256", "md5",
                  "license", "license_family", "depends", "constrains") if key in package}
        entry["metadata_sha256"] = digest(raw)
        entry["url"] = public_url(package.get("url", ""))
        entry["legal_evidence"] = []
        entry["missing_legal_files"] = []
        entry["native_files"] = []
        for relative in package.get("files", []):
            candidate = contained_file(prefix, relative)
            if candidate.suffix in {".so", ".dylib", ".dll", ".a"} or ".so." in candidate.name:
                entry["native_files"].append(relative)
            if legal_name(candidate):
                if candidate.is_file():
                    entry["legal_evidence"].append(store.capture(candidate, f"prefix:{relative}"))
                else:
                    entry["missing_legal_files"].append(relative)
        # Conda often keeps license texts only in package-cache info/licenses.
        cache = package.get("link", {}).get("source")
        if cache:
            base = Path(cache) / "info" / "licenses"
            if base.is_dir():
                for candidate in sorted(p for p in base.rglob("*") if p.is_file()):
                    if not candidate.resolve().is_relative_to(base.resolve()):
                        raise ValueError("Conda license cache symlink escapes evidence root")
                    entry["legal_evidence"].append(store.capture(
                        candidate, f"package-cache:info/licenses/{candidate.relative_to(base)}"))
        entry["license_review"] = "captured-not-reviewed" if entry["legal_evidence"] else "declared-only"
        entry["compatibility_review"] = "unknown"
        conda.append(entry)

    site_aliases: dict[str, list[str]] = {}
    for site in sorted((prefix / "lib").glob("python*/site-packages")):
        resolved = site.resolve(strict=True)
        if not resolved.is_relative_to(prefix):
            raise ValueError("Python site-packages escapes inspected prefix")
        site_aliases.setdefault(resolved.relative_to(prefix).as_posix(), []).append(
            site.relative_to(prefix).as_posix())
    sites = [prefix / relative for relative in sorted(site_aliases)]
    for site in sites:
        for info in sorted(site.glob("*.dist-info")):
            path = info / "METADATA"
            if not path.is_file():
                raise ValueError(f"Missing distribution METADATA: {info.name}")
            metadata = Parser().parsestr(path.read_text(encoding="utf-8"))
            entry = {
                "name": metadata["Name"], "version": metadata["Version"],
                "metadata": store.capture(path, f"prefix:{path.relative_to(prefix)}"),
                "license_expression": metadata["License-Expression"],
                "license": metadata["License"],
                "license_classifiers": [v for v in metadata.get_all("Classifier", [])
                                        if v.startswith("License ::")],
                "requires_dist": metadata.get_all("Requires-Dist", []),
                "legal_evidence": [store.capture(p, f"prefix:{p.relative_to(prefix)}")
                                   for p in legal_files(info)],
                "compatibility_review": "unknown",
            }
            record = info / "RECORD"
            entry["native_file_candidates"] = []
            entry["missing_owned_legal_files"] = []
            if record.is_file():
                entry["record"] = store.capture(record, f"prefix:{record.relative_to(prefix)}")
                captured = {p.resolve() for p in legal_files(info)}
                with record.open(encoding="utf-8", newline="") as stream:
                    for row in csv.reader(stream):
                        if not row:
                            continue
                        candidate = (site / row[0]).resolve()
                        if not candidate.is_relative_to(prefix):
                            raise ValueError("Python RECORD path escapes inspected prefix")
                        relative = candidate.relative_to(prefix).as_posix()
                        if candidate.suffix in {".so", ".dylib", ".dll", ".a"} or ".so." in candidate.name:
                            entry["native_file_candidates"].append(relative)
                        if legal_name(candidate) and candidate not in captured:
                            if candidate.is_file():
                                entry["legal_evidence"].append(store.capture(candidate, f"prefix:{relative}"))
                                captured.add(candidate)
                            else:
                                entry["missing_owned_legal_files"].append(relative)
            direct = info / "direct_url.json"
            if direct.is_file():
                value = json.loads(direct.read_text())
                entry["direct_url"] = public_url(value.get("url", ""))
                entry["direct_url_sha256"] = digest(direct.read_bytes())
            entry["license_review"] = "captured-not-reviewed" if entry["legal_evidence"] else "declared-only"
            python.append(entry)
    return {
        "schema_version": 1, "kind": "installed-metadata-and-legal-evidence",
        "prefix_label": prefix.name,
        "site_package_aliases": site_aliases,
        "limitations": [
            "Snapshot of existing prefix metadata, not a fresh solve or proof of environment usability.",
            "Conda and Python distributions overlap; do not sum them as unique dependencies.",
            "Native filenames are recorded, but dynamic linkage and all installed bytes are not verified.",
            "RECORD hashes are supplied metadata, not independently verified installed-file hashes.",
            "Captured legal texts are not automatically reviewed or redistribution-cleared.",
            "Python egg-info, unmanaged executables and operating-system frameworks need separate review.",
            "RECORD-owned legal files/native candidates are included; data rights and binary source correspondence remain unverified.",
        ],
        "conda": conda, "python": python,
    }


def collect_npm(root: Path, store: EvidenceStore) -> dict:
    lock = json.loads((root / "package-lock.json").read_text())
    entries, missing = [], []
    for relative, locked in sorted(lock["packages"].items()):
        if not relative:
            continue
        base = contained_file(root, relative)
        path = base / "package.json"
        if not path.is_file():
            missing.append(relative)
            continue
        actual = json.loads(path.read_text())
        evidence = []
        # Capture package-owned top-level/legal directories, not nested node_modules.
        for p in sorted(base.iterdir()):
            if p.is_file() and legal_name(p):
                evidence.append(store.capture(contained_file(base, p.name), f"npm:{relative}/{p.name}"))
            elif p.is_dir() and p.name.lower() in {"licenses", "licences"}:
                evidence.extend(store.capture(q, f"npm:{relative}/{q.relative_to(base)}")
                                for q in legal_files(p))
        entries.append({
            "package_path": relative, "name": actual.get("name"), "version": actual.get("version"),
            "locked_version": locked.get("version"),
            "matches_locked_version": actual.get("version") == locked.get("version"),
            "license_metadata": actual.get("license"), "lock_integrity": locked.get("integrity"),
            "metadata_sha256": digest(path.read_bytes()), "legal_evidence": evidence,
            "license_review": "captured-not-reviewed" if evidence else "declared-only",
            "compatibility_review": "unknown", "distribution_inclusion": "unverified",
        })
    return {"schema_version": 1, "kind": "installed-npm-legal-evidence",
            "lock_sha256": digest((root / "package-lock.json").read_bytes()),
            "packages": entries, "missing_lock_entries": missing,
            "limitations": ["Installed lock entries are not a browser or server shipping closure.",
                            "Package version equality does not verify npm archive integrity or all files."]}


def binary_format(data: bytes) -> str | None:
    """Identify candidates by bytes, including extensionless bundled tools."""
    if data.startswith(b"\x7fELF"):
        return "ELF"
    if data[:4] in {b"\xfe\xed\xfa\xce", b"\xce\xfa\xed\xfe",
                    b"\xfe\xed\xfa\xcf", b"\xcf\xfa\xed\xfe"}:
        return "Mach-O"
    if data[:4] in {b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca",
                    b"\xca\xfe\xba\xbf", b"\xbf\xba\xfe\xca"}:
        return "fat-binary-candidate"  # CAFEBABE is also Java's class magic.
    if data.startswith(b"MZ"):
        return "PE-candidate"
    if data.startswith(b"!<arch>\n"):
        return "ar-archive"
    return None


def verify_wheel_record(entries: list[dict], record: bytes) -> dict:
    """Verify internal SHA-256 claims, not authenticity or build provenance."""
    files = {entry["path"]: entry for entry in entries}
    seen, verified, unhashed = set(), [], []
    for row in csv.reader(io.StringIO(record.decode("utf-8"))):
        if len(row) != 3 or row[0] in seen:
            raise ValueError("Malformed or duplicate wheel RECORD row")
        name, claimed_hash, size = row
        seen.add(name)
        if name not in files:
            raise ValueError("Wheel RECORD references a missing archive member")
        item = files[name]
        if size and (not size.isdecimal() or int(size) != item["bytes"]):
            raise ValueError("Wheel RECORD size mismatch")
        if not claimed_hash:
            unhashed.append(name)
            continue
        expected = "sha256=" + base64.urlsafe_b64encode(bytes.fromhex(item["sha256"])).decode().rstrip("=")
        if claimed_hash != expected:
            raise ValueError("Wheel RECORD hash mismatch or unsupported hash algorithm")
        verified.append(name)
    return {"verified_sha256_files": len(verified), "unhashed_files": sorted(unhashed),
            "unlisted_files": sorted(files.keys() - seen),
            "scope": "internal consistency only; not signature, source correspondence or legal clearance"}


def inspect_archive(path: Path, store: EvidenceStore | None = None,
                    expected_sha256: str | None = None) -> dict:
    """Hash small project release archives without extracting or executing them."""
    archive_sha256 = digest(path.read_bytes())
    if expected_sha256 is not None and archive_sha256 != expected_sha256:
        raise ValueError("Archive digest does not match the predeclared identity")
    entries, names, legal_evidence, records = [], set(), [], []

    def entry(name: str, data: bytes) -> None:
        parts = PurePosixPath(name)
        if parts.is_absolute() or ".." in parts.parts or "\\" in name:
            raise ValueError("Unsafe archive member path")
        if name in names:
            raise ValueError("Duplicate archive member path")
        names.add(name)
        entries.append({"path": name, "sha256": digest(data), "bytes": len(data)})
        kind = binary_format(data)
        if kind:
            entries[-1]["binary_format_candidate"] = kind
        if name.endswith(".dist-info/RECORD"):
            records.append(data)
        if store and legal_name(PurePosixPath(name)):
            sha = digest(data)
            target = store.output / "blobs" / sha
            if target.exists():
                if target.read_bytes() != data:
                    raise ValueError("Conflicting retained archive legal text")
            else:
                with target.open("xb") as stream:
                    stream.write(data)
            legal_evidence.append({"origin": f"archive:{name}", "sha256": sha,
                                   "bytes": len(data), "blob": f"blobs/{sha}"})

    total = 0
    if path.suffix == ".whl" or zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            for item in archive.infolist():
                if item.is_dir():
                    continue
                total += item.file_size
                if item.file_size > 64 * 1024 * 1024 or total > 512 * 1024 * 1024:
                    raise ValueError("Archive exceeds bounded project-inspection size")
                entry(item.filename, archive.read(item))
    else:
        with tarfile.open(path) as archive:
            for item in archive:
                if item.issym() or item.islnk():
                    raise ValueError("Archive links require separate review")
                if not item.isfile():
                    continue
                total += item.size
                if item.size > 64 * 1024 * 1024 or total > 512 * 1024 * 1024:
                    raise ValueError("Archive exceeds bounded project-inspection size")
                stream = archive.extractfile(item)
                assert stream is not None
                with stream:
                    entry(item.name, stream.read())
    sensitive = [name for name in sorted(names) if PurePosixPath(name).name in
                 {"modeller_key.txt", ".env", "config.py"}
                 or PurePosixPath(name).suffix.lower() in {".pt", ".pth", ".ckpt", ".onnx", ".safetensors"}]
    wheel_record = None
    if path.suffix == ".whl":
        if len(records) != 1:
            raise ValueError("Expected exactly one wheel RECORD")
        wheel_record = verify_wheel_record(entries, records[0])
    return {"schema_version": 1, "kind": "project-release-archive-inspection",
             "archive_name": path.name, "archive_sha256": archive_sha256,
             "entries": sorted(entries, key=lambda e: e["path"]),
             "native_candidate_paths": sorted(e["path"] for e in entries if "binary_format_candidate" in e),
             "legal_evidence": legal_evidence, "wheel_record_verification": wheel_record,
            "legal_paths": [name for name in sorted(names) if legal_name(PurePosixPath(name))
                            or any(part.lower() in {"notices", "third_party_notices", "licenses"}
                                   for part in PurePosixPath(name).parts)],
            "restricted_or_sensitive_name_findings": sensitive,
            "compatibility_review": "unknown",
             "limitations": ["Filename screening is not a secret-content scan or full legal clearance.",
                             "Binary magic identifies candidates, not versions, linkage, source correspondence or licenses.",
                            "Archive contents do not establish rights to bundled third-party data."]}


def write_report(output: Path, report: dict) -> None:
    target = output / "inventory.json"
    if target.exists():
        raise FileExistsError("Use a new evidence directory; retained reports are never overwritten")
    with target.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, sort_keys=True)
        stream.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--prefix", type=Path)
    source.add_argument("--npm-root", type=Path)
    source.add_argument("--archive", type=Path)
    parser.add_argument("--expected-sha256", help="Require this archive identity before inspecting members")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.expected_sha256 is not None and args.archive is None:
        parser.error("--expected-sha256 requires --archive")
    if args.output.exists():
        raise FileExistsError("Evidence output must be a new directory")
    source_root = args.prefix or args.npm_root
    if source_root and args.output.resolve().is_relative_to(source_root.resolve()):
        raise ValueError("Evidence output must not modify the inspected source/prefix")
    store = EvidenceStore(args.output)
    if args.prefix:
        report = collect_prefix(args.prefix, store)
    elif args.npm_root:
        report = collect_npm(args.npm_root, store)
    else:
        report = inspect_archive(args.archive, store, args.expected_sha256)
    report["collector_sha256"] = digest(Path(__file__).read_bytes())
    report["collector_source"] = store.capture(Path(__file__), "collector-source")
    write_report(args.output, report)
    print(f"Captured {report['kind']}; legal review and compatibility remain pending")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
