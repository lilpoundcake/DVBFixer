#!/usr/bin/env python3
"""Freeze and acquire unused50+ natural-gap originals after union regrouping."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts import acquire_real_gap_sources as acquisition
from scripts.collect_real_gap_candidates import digest, save


def select(rows, excluded):
    selected = sorted((r for r in rows if r["entry_id"] not in excluded
                       and r.get("longest_gap_metadata", 0) >= 50), key=lambda r: r["entry_id"])
    if len(selected) > 200:
        raise ValueError("hard-source budget exceeds200 originals; a new declaration is required")
    strata = [{"longest_span": [lo, hi], "selected": sum(lo <= r["longest_gap_metadata"] and
               (hi is None or r["longest_gap_metadata"] <= hi) for r in selected)}
              for lo, hi in ((50, 99), (100, None))]
    return selected, strata


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
            raise ValueError("prior original-source ledger changed")
        ids = {r["entry_id"] for r in json.loads(path.read_text())["rows"]}
        if ids & excluded:
            raise ValueError("previous source acquisitions overlap")
        excluded.update(ids)
        bindings[str(path)] = digest(path.read_bytes())
    rows, strata = select(json.loads(data_path.read_text())["rows"], excluded)
    for path in (Path(__file__), Path(acquisition.__file__), Path(acquisition.inventory_mmcif.__code__.co_filename)):
        bindings[str(path)] = digest(path.read_bytes())
    output.mkdir(parents=True, exist_ok=False)
    save(output / "declaration.json", {"grouped_directory": str(grouped), "bindings": bindings,
        "selected_rows": rows, "strata": strata, "count": len(rows), "excluded_prior_source_ids": sorted(excluded),
        "purpose": "replenish authentic long-gap strata; source screening, not sampling/admission",
        "max_source_count": 200, "max_evidence_bytes": 2*1024**3, "serial_requests": True,
        "timeout_seconds": 60, "max_HTTP_attempts": 3,
        "per_source_compressed_max": acquisition.MAX_COMPRESSED,
        "per_source_decompressed_max": acquisition.MAX_DECOMPRESSED,
        "upper_gap_length_cap": None, "new_neural_tasks": 0, "final_cohort_selected": False})
    print(json.dumps({"hard_original_sources_declared": len(rows), "strata": strata}, indent=2))


def run(output):
    plan = json.loads((output / "declaration.json").read_text())
    for name, expected in plan["bindings"].items():
        if digest(Path(name).read_bytes()) != expected:
            raise ValueError(f"hard-source declaration changed: {name}")
    rows = json.loads((Path(plan["grouped_directory"]) / "grouped-census.json").read_text())["rows"]
    selected, strata = select(rows, set(plan["excluded_prior_source_ids"]))
    if selected != plan["selected_rows"] or strata != plan["strata"]:
        raise ValueError("hard-source selection changed")
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
