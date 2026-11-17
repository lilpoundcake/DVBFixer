#!/usr/bin/env python3
"""Recruit remaining25-49 originals; do not substitute them for long-gap quotas."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts import acquire_real_gap_sources as acquisition
from scripts.collect_real_gap_candidates import digest, save


def select(rows, excluded):
    selected = sorted((r for r in rows if r["entry_id"] not in excluded
                       and 25 <= r.get("longest_gap_metadata", 0) <= 49), key=lambda r: r["entry_id"])
    if len(selected) > 200:
        raise ValueError("medium-source budget exceeds200; explicit new declaration required")
    return selected, [{"longest_span": [25, 49], "selected": len(selected)}]


def declare(grouped, previous, output):
    data_path = grouped / "grouped-census.json"
    complete = json.loads((grouped / "completed.json").read_text())
    if digest(data_path.read_bytes()) != complete["grouped_census_sha256"]:
        raise ValueError("union grouped census changed")
    bindings = {str(data_path): digest(data_path.read_bytes())}
    excluded = set()
    for root in previous:
        path = root / "source-results.json"
        completed = json.loads((root / "completed.json").read_text())
        if digest(path.read_bytes()) != completed["source_results_sha256"]:
            raise ValueError("prior source ledger changed")
        ids = {r["entry_id"] for r in json.loads(path.read_text())["rows"]}
        if ids & excluded:
            raise ValueError("prior source acquisitions overlap")
        excluded.update(ids)
        bindings[str(path)] = digest(path.read_bytes())
    selected, strata = select(json.loads(data_path.read_text())["rows"], excluded)
    for path in (Path(__file__), Path(acquisition.__file__), Path(acquisition.inventory_mmcif.__code__.co_filename)):
        bindings[str(path)] = digest(path.read_bytes())
    output.mkdir(parents=True, exist_ok=False)
    save(output / "declaration.json", {"grouped_directory": str(grouped), "bindings": bindings,
        "selected_rows": selected, "strata": strata, "count": len(selected), "excluded_prior_source_ids": sorted(excluded),
        "reason": "existing source program cannot reach200;25-49 ready quota remains underfilled",
        "max_sources": 200, "max_evidence_bytes": 2*1024**3, "timeout_seconds": 60,
        "max_HTTP_attempts": 3, "serial_network": True,
        "max_source_compressed_bytes": acquisition.MAX_COMPRESSED,
        "max_source_decompressed_bytes": acquisition.MAX_DECOMPRESSED,
        "band_is_recruitment_stratum_not_method_length_cutoff": True,
        "long_quotas_unchanged": True, "new_neural_tasks": 0, "scientific_admission": False})
    print(json.dumps({"additional_medium_originals_declared": len(selected), "strata": strata}, indent=2))


def run(output):
    plan = json.loads((output / "declaration.json").read_text())
    for name, expected in plan["bindings"].items():
        if digest(Path(name).read_bytes()) != expected:
            raise ValueError(f"medium-source declaration changed: {name}")
    rows = json.loads((Path(plan["grouped_directory"]) / "grouped-census.json").read_text())["rows"]
    selected, strata = select(rows, set(plan["excluded_prior_source_ids"]))
    if selected != plan["selected_rows"] or strata != plan["strata"]:
        raise ValueError("medium-source selection changed")
    original = acquisition.select_reserves
    try:
        acquisition.select_reserves = lambda rows: (selected, strata)
        acquisition.run(Path(plan["grouped_directory"]), output / "acquisition")
    finally:
        acquisition.select_reserves = original


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--grouped", type=Path)
    parser.add_argument("--previous", type=Path, action="append")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.grouped.resolve(), [p.resolve() for p in args.previous], args.output.resolve())
    else:
        run(args.output.resolve())
