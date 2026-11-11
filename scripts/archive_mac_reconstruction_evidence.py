#!/usr/bin/env python3
"""Archive retained research, not runtimes/weights, in bounded Git-sized parts."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import tarfile
from pathlib import Path

ROOTS = (
    "real-complex-cc94-2026-11-06",
    "joint-gap-pilot-1HE8-v1",
    "1HE8-all-gaps-multiseed-local-v2",
    "1HE8-complete-context-v3",
    "1HE8-neighbor-motif-v4",
    "1HE8-full-joint-spatial-v2",
    "1HE8-patchwork-all-generated-anchors-v1",
    "1HE8-motif-seed-ensembles",
    "mac-reconstruction-2026-11-11",
    "mac-reconstruction-comparison-2026-11-11",
)


class Parts:
    def __init__(self, directory, name, limit=32 * 1024**2):
        self.directory, self.name, self.limit = directory, name, limit
        self.stream = None
        self.size = 0
        self.paths = []

    def write(self, data):
        total = len(data)
        while data:
            if self.stream is None or self.size == self.limit:
                if self.stream is not None:
                    self.stream.close()
                path = self.directory / f"{self.name}.tar.gz.part{len(self.paths)+1:03d}"
                self.stream = path.open("xb")
                self.paths.append(path)
                self.size = 0
            count = min(len(data), self.limit - self.size)
            self.stream.write(data[:count])
            data = data[count:]
            self.size += count
        return total

    def flush(self):
        if self.stream is not None:
            self.stream.flush()

    def close(self):
        if self.stream is not None:
            self.stream.close()


class HashedReader:
    def __init__(self, stream):
        self.stream = stream
        self.hash = hashlib.sha256()

    def read(self, size=-1):
        data = self.stream.read(size)
        self.hash.update(data)
        return data


def excluded(path):
    if path.is_symlink():
        return "symlink: never dereference into local runtimes"
    if any(part in {"__pycache__", ".pytest_cache", ".git", "node_modules"} for part in path.parts):
        return "runtime/cache metadata"
    if path.name == ".DS_Store" or path.suffix == ".pyc":
        return "host/cache metadata"
    if "profile" in path.parts and any(part.startswith("chrome-") for part in path.parts):
        return "browser profile: may contain local credentials"
    return None


def archive(project, output):
    output.mkdir(parents=True, exist_ok=False)
    manifest = {"schema": 1, "scope": "private reconstruction research and retained failures",
                "runtimes_weights_and_browser_credentials_included": False,
                "original_paths_relative_to_repository": True, "archives": [], "excluded": []}
    for name in ROOTS:
        root = project / "tmp" / name
        if not root.is_dir():
            raise ValueError(f"required retained study missing: {root}")
        parts = Parts(output, name)
        files = []
        try:
            with gzip.GzipFile(filename="", fileobj=parts, mode="wb", mtime=0, compresslevel=6) as compressed:
                with tarfile.open(fileobj=compressed, mode="w|") as tar:
                    for path in sorted(root.rglob("*")):
                        if not path.is_file() and not path.is_symlink():
                            continue
                        reason = excluded(path.relative_to(project))
                        if path.is_symlink():
                            reason = "symlink: never dereference into local runtimes"
                        if reason:
                            manifest["excluded"].append({"path": str(path.relative_to(project)), "reason": reason})
                            continue
                        before = path.stat()
                        info = tar.gettarinfo(str(path), arcname=str(path.relative_to(project)))
                        info.uid = info.gid = info.mtime = 0
                        info.uname = info.gname = ""
                        info.mode = 0o644
                        with path.open("rb") as stream:
                            reader = HashedReader(stream)
                            tar.addfile(info, reader)
                        after = path.stat()
                        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                            raise ValueError(f"evidence changed during archival: {path}")
                        files.append({"path": info.name, "size": info.size, "sha256": reader.hash.hexdigest()})
        finally:
            parts.close()
        checksums = [{"path": p.name, "size": p.stat().st_size,
                      "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in parts.paths]
        manifest["archives"].append({"study": name, "parts": checksums, "files": files})
        print(json.dumps({"study": name, "files": len(files), "compressed_bytes": sum(p["size"] for p in checksums)}), flush=True)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    archive(Path(__file__).resolve().parents[1], args.output.resolve())
