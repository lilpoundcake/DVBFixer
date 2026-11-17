#!/usr/bin/env python3
"""Bounded additional metadata screening for underfilled authentic hard strata."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts import collect_real_gap_candidates as census

STRATA = ((50, 99, 250), (100, 199, 150), (200, None, 100))
SALT = "natural-gap-hard-strata-recruitment-v1"


def rank_new(ids, excluded, stratum):
    return sorted(set(ids)-set(excluded), key=lambda value: census.digest(f"{SALT}:{stratum}:{value}".encode()))


def declare(previous, output):
    raw = previous / "census.json"
    complete = json.loads((previous / "completed.json").read_text())
    if census.digest(raw.read_bytes()) != complete["census_sha256"]:
        raise ValueError("previous metadata census changed")
    excluded = [r["entry_id"] for r in json.loads(raw.read_text())["rows"]]
    if len(excluded) != len(set(excluded)) or len(excluded) != complete["entries"]:
        raise ValueError("incomplete previous census")
    output.mkdir(parents=True, exist_ok=False)
    census.save(output / "declaration.json", {"excluded_prior_metadata_ids": excluded,
        "bindings": {str(p): census.digest(p.read_bytes()) for p in (Path(__file__), Path(census.__file__), raw)},
        "additional_metadata_entries": sum(n for _, _, n in STRATA), "entry_total_missing_proxy_strata": STRATA,
        "reason": "current50-99 source-ready13 is below40 quota; existing longest50+ metadata originals exhausted",
        "not_a_gap_length_ceiling": True, "selection_salt": SALT,
        "new_checkpoint_jobs": 0, "coordinates_downloaded": False,
        "no_final_split": True, "must_regroup_union_with_prior_entries_before_splitting": True,
        "budgets": {"serial_network": True, "timeout_seconds": 60, "max_HTTP_attempts": 3,
                    "max_response_bytes": 32*1024**2, "max_evidence_bytes": 512*1024**2,
                    "min_free_disk_bytes": 2*1024**3}})


def run(output):
    plan = json.loads((output / "declaration.json").read_text())
    for name, expected in plan["bindings"].items():
        if census.digest(Path(name).read_bytes()) != expected:
            raise ValueError(f"additional metadata declaration changed: {name}")
    work = output / "acquisition"
    work.mkdir(parents=True, exist_ok=False)
    (work / "http").mkdir()
    excluded = set(plan["excluded_prior_metadata_ids"])
    selected, populations = [], []
    for index, (lo, hi, count) in enumerate(STRATA):
        response = census.post(work, f"search-{index}", census.SEARCH, census.search_query(lo, hi))
        ids = census.identifiers(response)
        candidates = rank_new(ids, excluded | set(selected), index)
        if len(candidates) < count:
            raise ValueError("proxy stratum insufficient; no substitution or duplicate adoption")
        selected.extend(candidates[:count])
        populations.append({"total_missing_proxy": [lo, hi], "available_new": len(candidates), "selected": count})
    census.save(work / "selection.json", {"entry_ids": selected, "populations": populations,
        "not_final_cohort_or_longest_span_admission": True})
    rows = []
    for batch, start in enumerate(range(0, len(selected), 12)):
        ids = selected[start:start+12]
        response = census.post(work, f"metadata-{batch:03d}", census.GRAPHQL,
                               {"query": census.QUERY, "variables": {"ids": ids}})
        entries = response["data"]["entries"]
        found = {e["rcsb_id"]: e for e in entries if e is not None}
        if len(entries) != len(ids) or set(found) != set(ids):
            raise ValueError("missing/repeated metadata entry")
        for entry_id in ids:
            try:
                row = census.summarize(found[entry_id])
            except (ValueError, KeyError, TypeError) as exc:
                row = {"entry_id": entry_id, "status": "metadata_inventory_refused", "error": str(exc),
                       "scientific_admission_complete": False}
            rows.append(row)
        print(json.dumps({"additional_metadata_read": len(rows), "declared": len(selected)}), flush=True)
    census.save(work / "census.json", {"rows": rows, "declaration_sha256": census.digest((output / "declaration.json").read_bytes()),
        "not_scientific_admission": True, "union_grouping_before_split_required": True})
    summary = {"entries": len(rows), "census_sha256": census.digest((work / "census.json").read_bytes()),
        "metadata_refusals": sum(r.get("status") == "metadata_inventory_refused" for r in rows),
        "longest_gap_bins": {f"{lo}-{hi or 'plus'}": sum(lo <= r.get("longest_gap_metadata", 0) and
            (hi is None or r.get("longest_gap_metadata", 0) <= hi) for r in rows) for lo, hi in census.BINS},
        "prior_entries_repeated": 0, "new_checkpoint_jobs": 0, "final_cohort_selected": False}
    census.save(work / "completed.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.previous.resolve(), args.output.resolve())
    else:
        run(args.output.resolve())
