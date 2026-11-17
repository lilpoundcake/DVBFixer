#!/usr/bin/env python3
"""Immutable cumulative snapshot of settled source attempts and pending tasks.

Does not inspect a running executor's outputs or infer readiness from directories.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save


def account(original_ids, rows, pending_ids):
    settled = [r["entry_id"] for r in rows]
    combined = settled + list(pending_ids)
    if len(combined) != len(set(combined)) or set(combined) != set(original_ids):
        raise ValueError("settled outcomes and pending plan must partition original source tasks exactly once")
    return dict(Counter(r["status"] for r in rows))


def summarize(original, settled_series, pending, output):
    bindings = {}

    def read(path):
        raw = path.read_bytes()
        bindings[str(path)] = digest(path)
        return json.loads(raw)

    original_plan = read(original / "plan.json")
    original_jobs = {j["case_id"]: j for j in original_plan["jobs"]}
    if len(original_jobs) != original_plan["count"]:
        raise ValueError("original source plan contains duplicate or missing jobs")
    rows = []
    for series in settled_series:
        plan = read(series / "plan.json")
        # These series have stopped. Unexecuted cases belong only to pending.
        for job in plan["jobs"]:
            if job["case_id"] not in original_jobs or job["case"] != original_jobs[job["case_id"]]["case"]:
                raise ValueError("continuation changed the source identity")
            work = series / "execution" / f"{job['case_id']}-seed-{job['seed']}"
            receipt_path = work / "execution.json"
            if not receipt_path.exists():
                if work.exists():
                    raise ValueError("uncommitted settled workdir cannot be adopted")
                continue
            receipt = read(receipt_path)
            row = {"entry_id": job["case_id"], "series": str(series), "execution_receipt": str(receipt_path)}
            if not receipt["completed"]:
                if receipt["guard_stop_reason"] != "wall-time-limit":
                    raise ValueError("unclassified technical/resource stop")
                row.update(status="source_time_blocked", elapsed_seconds=receipt["elapsed_seconds"],
                           peak_rss_bytes=receipt["peak_rss_bytes"], not_retried=True)
            else:
                if receipt["returncode"] != 0 or receipt["guard_stop_reason"] is not None:
                    raise ValueError("inconsistent completed execution")
                result_path = work / "result.json"
                result = read(result_path)
                row.update(status=result["status"], result_path=str(result_path))
                if result["status"] == "private_source_projection_prepared":
                    for path_key, hash_key in (("source_input", "source_input_sha256"), ("request_path", "request_sha256")):
                        path = Path(result[path_key])
                        if not path.resolve().is_relative_to(work.resolve()) or digest(path) != result[hash_key]:
                            raise ValueError("ready source artifact changed or escaped its reservation")
                        bindings[str(path)] = result[hash_key]
                    row.update(target_residues=result["target_residues"], longest_gap=max(result["gap_lengths"]),
                               gap_count=result["gap_count"], request_path=result["request_path"])
                elif result["status"] == "source_preparation_refused":
                    row.update(stage=result["stage"], error=result["error"])
                else:
                    raise ValueError("unknown completed source outcome")
            rows.append(row)
    # Read an active plan only, never its running output directories.
    pending_ids = [] if pending is None else [j["case_id"] for j in read(pending / "plan.json")["jobs"]]
    counts = account(original_jobs, rows, pending_ids)
    ready = [r for r in rows if r["status"] == "private_source_projection_prepared"]
    summary = {"original_tasks": len(original_jobs), "settled_outcomes": len(rows),
        "status_counts": counts, "pending_at_snapshot": len(pending_ids),
        "ready_longest_gap_bins": {f"{lo}-{hi or 'plus'}": sum(lo <= r["longest_gap"] and
            (hi is None or r["longest_gap"] <= hi) for r in ready)
            for lo, hi in ((1, 8), (9, 24), (25, 49), (50, 99), (100, None))},
        "not_scientific_backend_admission_or_final_cohort": True, "new_sampling": False}
    output.mkdir(parents=True, exist_ok=False)
    save(output / "ledger.json", {"rows": rows, "pending_entry_ids_at_snapshot": pending_ids,
        "bindings": bindings, "analyzer_sha256": digest(Path(__file__)), "summary": summary})
    save(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", required=True, type=Path)
    parser.add_argument("--settled", required=True, type=Path, action="append")
    parser.add_argument("--pending", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    summarize(args.original.resolve(), [p.resolve() for p in args.settled],
              args.pending.resolve() if args.pending else None, args.output.resolve())
