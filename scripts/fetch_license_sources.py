#!/usr/bin/env python3
"""Archive explicitly registered legal/provenance sources, never model weights."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

from collect_license_evidence import EvidenceStore, public_url, write_report

ROOT = Path(__file__).resolve().parents[1]
REGISTER = ROOT / "docs/research/licensing-upstream-sources.json"
MAX_BYTES = 4 * 1024 * 1024


def fetch_source(source: dict) -> tuple[dict, bytes | None]:
    url = source["url"]
    if public_url(url) != url:
        raise ValueError("Sources must use public HTTP(S) URLs without credentials or tokens")
    request = Request(url, headers={"User-Agent": "dvbfixer-license-evidence/1"})
    try:
        with urlopen(request, timeout=60) as response:
            data = response.read(MAX_BYTES + 1)
            if len(data) > MAX_BYTES:
                raise ValueError("Legal evidence exceeds the bounded download size")
            return ({**source, "status": "captured-not-reviewed",
                     "final_url": public_url(response.url),
                     "content_type": response.headers.get("Content-Type"),
                     "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}, data)
    except (OSError, ValueError) as exc:
        # Do not retain arbitrary HTTP exception messages or redirect credentials.
        return ({**source, "status": "unavailable", "error_type": type(exc).__name__}, None)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Use a new output directory; never overwrite retained evidence")
    sources = json.loads(REGISTER.read_text())["sources"]
    store = EvidenceStore(args.output)
    entries = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for record, data in pool.map(fetch_source, sources):
            if data is not None:
                path = store.output / "blobs" / record["sha256"]
                if not path.exists():
                    with path.open("xb") as stream:
                        stream.write(data)
                elif path.read_bytes() != data:
                    raise ValueError("Conflicting retained evidence blob")
                record["blob"] = f"blobs/{record['sha256']}"
            entries.append(record)
    write_report(args.output, {
        "schema_version": 1, "kind": "upstream-legal-source-capture",
        "captured_at": datetime.now(timezone.utc).isoformat(),
        "register_sha256": hashlib.sha256(REGISTER.read_bytes()).hexdigest(),
        "sources": entries,
        "limitations": ["Mutable pages are timestamped/hash-bound observations, not immutable grants.",
                        "No compatibility or artifact-redistribution decision is made by fetching texts."],
    })
    failures = [entry["id"] for entry in entries if entry["status"] == "unavailable"]
    print(f"Captured {len(entries) - len(failures)}/{len(entries)} sources; unavailable: {failures}")
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
