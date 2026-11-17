#!/usr/bin/env python3
"""Recompute union-wide leakage components, never concatenate local group IDs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.collect_real_gap_candidates import BINS, digest, save
from scripts.group_real_gap_candidates import components


def regroup(series):
    rows = [row for source in series for row in source]
    ids = [row["entry_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("grouped metadata sources overlap; explicit reconciliation required")
    rows = sorted(rows, key=lambda r: r["entry_id"])
    groups = components(rows)
    membership = {entry_id: i for i, group in enumerate(groups) for entry_id in group}
    return [{**row, "conservative_entry_component": membership[row["entry_id"]],
             "split_assignment": None} for row in rows], groups


def merge(inputs, output):
    series, bindings = [], {}
    for root in inputs:
        path = root / "grouped-census.json"
        complete_path = root / "completed.json"
        complete = json.loads(complete_path.read_text())
        raw = path.read_bytes()
        if digest(raw) != complete["grouped_census_sha256"]:
            raise ValueError("grouped metadata changed")
        rows = json.loads(raw)["rows"]
        if len(rows) != complete["entries"]:
            raise ValueError("incomplete grouped source")
        series.append(rows)
        bindings[str(path)] = digest(raw)
        bindings[str(complete_path)] = digest(complete_path.read_bytes())
    rows, groups = regroup(series)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "grouped-census.json", {"rows": rows, "components": groups,
        "source_bindings": bindings, "merger_sha256": digest(Path(__file__).read_bytes()),
        "final_split_assignment": None, "scientific_admission_complete": False,
        "limitation": "sequence clusters are grouping controls, not proof of independent folds or training absence"})
    summary = {"entries": len(rows), "components": len(groups), "largest_component": max(map(len, groups)),
        "entries_with_complete_30pct_groups": sum(r["modern_30pct_grouping_complete"] for r in rows),
        "grouped_census_sha256": digest((output / "grouped-census.json").read_bytes()),
        "longest_gap_bins": {f"{lo}-{hi or 'plus'}": sum(lo <= r.get("longest_gap_metadata", 0) and
            (hi is None or r.get("longest_gap_metadata", 0) <= hi) for r in rows) for lo, hi in BINS},
        "final_cohort_selected": False, "scientific_admission_complete": False}
    save(output / "completed.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grouped", required=True, type=Path, action="append")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    merge([p.resolve() for p in args.grouped], args.output.resolve())
