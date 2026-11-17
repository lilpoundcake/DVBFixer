#!/usr/bin/env python3
"""Declare nonduplicating original-source recruitment from the existing census.

This adds source screening only, not neural tasks or scientific admission.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts import acquire_real_gap_sources as acquisition
from scripts.collect_real_gap_candidates import BINS, digest, save

TARGETS = (60, 100, 60, 0, 0)
SALT = "real-gap-additional-original-sources-v1"


def select(rows, excluded):
    previous = [r for r in rows if r["entry_id"] in excluded]
    used = Counter(r["conservative_entry_component"] for r in previous)
    selected, strata = [], []
    for (lo, hi), target in zip(BINS, TARGETS, strict=True):
        pool = [r for r in rows if r["entry_id"] not in excluded and r.get("natural_gap_count_metadata", 0) > 0
                and not r.get("small_control_not_primary") and lo <= r["longest_gap_metadata"]
                and (hi is None or r["longest_gap_metadata"] <= hi)]
        available = len(pool)
        for _ in range(min(target, available)):
            row = min(pool, key=lambda r: (used[r["conservative_entry_component"]],
                digest(f"{SALT}:{r['entry_id']}".encode())))
            selected.append(row)
            used[row["conservative_entry_component"]] += 1
            pool.remove(row)
        strata.append({"longest_span": [lo, hi], "reserve_target": target, "available": available,
                       "selected": min(target, available)})
    if len({r["entry_id"] for r in selected}) != len(selected):
        raise ValueError("duplicate source recruitment")
    return selected, strata


def declare(grouped, previous, output):
    data_path = grouped / "grouped-census.json"
    complete = json.loads((grouped / "completed.json").read_text())
    if digest(data_path.read_bytes()) != complete["grouped_census_sha256"]:
        raise ValueError("grouped census changed")
    previous_path = previous / "source-results.json"
    previous_complete = json.loads((previous / "completed.json").read_text())
    if digest(previous_path.read_bytes()) != previous_complete["source_results_sha256"]:
        raise ValueError("prior original-source ledger changed")
    excluded = {r["entry_id"] for r in json.loads(previous_path.read_text())["rows"]}
    selected, strata = select(json.loads(data_path.read_text())["rows"], excluded)
    files = [Path(__file__), Path(acquisition.__file__), Path(acquisition.inventory_mmcif.__code__.co_filename),
             data_path, previous_path]
    output.mkdir(parents=True, exist_ok=False)
    save(output / "declaration.json", {"grouped_directory": str(grouped),
        "bindings": {str(p): digest(p.read_bytes()) for p in files}, "selected_rows": selected, "strata": strata,
        "excluded_prior_source_ids": sorted(excluded), "count": len(selected),
        "reason": "101 ready projections cannot satisfy200-entry cohort; recruit unused existing-census originals",
        "new_metadata_entries": 0, "new_ML_tasks": 0, "scientific_admission": False,
        "whole_context_review_and_long_strata_replenishment_still_required": True,
        "no_long_case_removed_or_replaced": True,
        "budgets": {"max_sources": sum(TARGETS), "serial_network": True, "timeout_seconds_per_attempt": 60,
                    "attempts": 3, "max_evidence_bytes": 2*1024**3,
                    "max_source_compressed_bytes": acquisition.MAX_COMPRESSED,
                    "max_source_decompressed_bytes": acquisition.MAX_DECOMPRESSED}})
    print(json.dumps({"additional_original_sources_declared": len(selected), "strata": strata}, indent=2))


def run(output):
    plan = json.loads((output / "declaration.json").read_text())
    for name, expected in plan["bindings"].items():
        if digest(Path(name).read_bytes()) != expected:
            raise ValueError(f"source expansion declaration changed: {name}")
    selected = plan["selected_rows"]
    if any(r["entry_id"] in plan["excluded_prior_source_ids"] for r in selected):
        raise ValueError("would repeat an original source acquisition")
    original = acquisition.select_reserves
    try:
        # Reuse unchanged source download, byte/hash and inventory guards.
        acquisition.select_reserves = lambda rows: (selected, plan["strata"])
        acquisition.run(Path(plan["grouped_directory"]), output / "acquisition")
    finally:
        acquisition.select_reserves = original


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--grouped", type=Path)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.grouped.resolve(), args.previous.resolve(), args.output.resolve())
    else:
        run(args.output.resolve())
