#!/usr/bin/env python3
"""Freeze a diverse reserve pool, download original CIFs and inventory real gaps."""
from __future__ import annotations

import argparse
import gzip
import io
import json
import shutil
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from scripts.collect_real_gap_candidates import BINS, PINNED, digest, save
from scripts.real_gap_source_inventory import inventory_mmcif

RESERVES = (75, 110, 110, 80, 75)
MAX_COMPRESSED = 16 * 1024**2
MAX_DECOMPRESSED = 128 * 1024**2


def select_reserves(rows):
    eligible = [r for r in rows if r.get("natural_gap_count_metadata", 0) > 0
                and not r.get("small_control_not_primary")]
    by_id = {r["entry_id"]: r for r in eligible}
    selected = [by_id[i] for i in PINNED if i in by_id]
    used = Counter(r["conservative_entry_component"] for r in selected)
    strata = []
    for (lo, hi), target in zip(BINS, RESERVES, strict=True):
        def belongs(row, lo=lo, hi=hi):
            length = row["longest_gap_metadata"]
            return length >= lo and (hi is None or length <= hi)
        pool = [r for r in eligible if belongs(r)]
        chosen_ids = {r["entry_id"] for r in selected}
        remaining = [r for r in pool if r["entry_id"] not in chosen_ids]
        needed = max(0, min(target, len(pool))-sum(belongs(r) for r in selected))
        for _ in range(needed):
            row = min(remaining, key=lambda r: (used[r["conservative_entry_component"]],
                digest(f"real-gap-source-reserves-v1:{r['entry_id']}".encode())))
            selected.append(row)
            used[row["conservative_entry_component"]] += 1
            remaining.remove(row)
        strata.append({"longest_span": [lo, hi], "reserve_target": target, "available": len(pool),
                       "selected": sum(belongs(r) for r in selected)})
    if len(selected) > 450 or len(selected) != len({r["entry_id"] for r in selected}):
        raise ValueError("source reserve count/identity changed")
    return selected, strata


def fetch_source(output, entry_id):
    path = output / "sources" / f"{entry_id}.cif.gz"
    receipt = output / "sources" / f"{entry_id}.receipt.json"
    url = f"https://files.rcsb.org/download/{entry_id}.cif.gz"
    if path.exists() or receipt.exists():
        raise ValueError("source reservation already exists; explicit verified recovery required")
    if shutil.disk_usage(output).free < 2 * 1024**3:
        raise RuntimeError("source acquisition lacks 2 GiB free headroom")
    if sum(p.stat().st_size for p in output.rglob("*") if p.is_file()) >= 2 * 1024**3:
        raise RuntimeError("source acquisition reached 2 GiB evidence ceiling")
    errors = []
    for attempt in range(1, 4):
        try:
            req = Request(url, headers={"User-Agent": "dvbfixer-natural-gap-research/1"})
            with urlopen(req, timeout=60) as response:
                raw = response.read(MAX_COMPRESSED+1)
                status = response.status
            if len(raw) > MAX_COMPRESSED:
                with path.with_suffix(".partial-resource-stop").open("xb") as stream:
                    stream.write(raw)
                raise ValueError("source exceeds declared 16 MiB compressed ceiling; no truncation/admission")
            with path.open("xb") as stream:
                stream.write(raw)
            save(receipt, {"url": url, "sha256": digest(raw), "bytes": len(raw), "status": status,
                "received_utc": datetime.now(UTC).isoformat(), "attempts": attempt, "previous_errors": errors,
                "original_download_unchanged": True})
            return path, raw
        except (HTTPError, URLError, TimeoutError) as exc:
            error = {"attempt": attempt, "type": type(exc).__name__, "message": str(exc)}
            save(output / "sources" / f"{entry_id}.attempt{attempt}.failure.json", error)
            errors.append(error)
    raise RuntimeError(f"source download failed after three retained attempts: {entry_id}")


def run(grouped, output):
    data = (grouped / "grouped-census.json").read_bytes()
    complete = json.loads((grouped / "completed.json").read_text())
    if digest(data) != complete["grouped_census_sha256"]:
        raise ValueError("source-reserve metadata changed")
    selected, strata = select_reserves(json.loads(data)["rows"])
    output.mkdir(parents=True, exist_ok=False)
    (output / "sources").mkdir()
    (output / "inventory").mkdir()
    save(output / "protocol.json", {"protocol": "natural-gap-original-source-reserve-v1",
        "grouped_census_sha256": digest(data), "source_entries_declared": len(selected), "selection_strata": strata,
        "target_primary_entries": 200, "final_cohort_selected": False, "split_assignment": None,
        "source_code_sha256": {str(p): digest(p.read_bytes()) for p in (Path(__file__), Path(inventory_mmcif.__code__.co_filename))},
        "source_compressed_max_bytes": MAX_COMPRESSED, "source_decompressed_max_bytes": MAX_DECOMPRESSED,
        "disk_evidence_cap_bytes": 2*1024**3, "minimum_free_disk_bytes": 2*1024**3,
        "network_concurrency": 1, "request_timeout_seconds": 60, "attempts_per_source": 3,
        "coordinate_edits": False, "model_sampling": False,
        "purpose": "source verification before final selection; biological assembly rendering/chemistry admission still required"})
    save(output / "selection.json", {"rows": selected, "strata": strata})
    rows = []
    for metadata in selected:
        entry_id = metadata["entry_id"]
        path, raw = fetch_source(output, entry_id)
        report = {"entry_id": entry_id, "source_path": str(path), "source_sha256": digest(raw),
                  "conservative_entry_component": metadata["conservative_entry_component"]}
        try:
            with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
                text = stream.read(MAX_DECOMPRESSED+1)
            if len(text) > MAX_DECOMPRESSED:
                raise ValueError("source exceeds 128 MiB decompressed inventory ceiling")
            inventory = inventory_mmcif(text.decode("utf-8"), metadata)
            report.update(status="source_inventory_complete", inventory=inventory)
        except (ValueError, RuntimeError, UnicodeError, OSError) as exc:
            report.update(status="source_inventory_refused", error_type=type(exc).__name__, error=str(exc),
                          scientific_admission_complete=False)
        save(output / "inventory" / f"{entry_id}.json", report)
        rows.append(report)
        if len(rows) % 20 == 0 or len(rows) == len(selected):
            print(json.dumps({"original_sources_downloaded": len(rows), "declared": len(selected),
                "inventory_complete": sum(r["status"] == "source_inventory_complete" for r in rows)}), flush=True)
    save(output / "source-results.json", {"rows": rows, "final_cohort_selected": False, "scientific_admission_complete": False})
    save(output / "completed.json", {"downloaded": len(rows), "inventory_complete": sum(r["status"] == "source_inventory_complete" for r in rows),
        "refusals": sum(r["status"] == "source_inventory_refused" for r in rows),
        "source_results_sha256": digest((output / "source-results.json").read_bytes()),
        "final_200_admitted": False, "checkpoint_sampling_started": False})
    print((output / "completed.json").read_text(), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grouped", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.grouped.resolve(), args.output.resolve())
