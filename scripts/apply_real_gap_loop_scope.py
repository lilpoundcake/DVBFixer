#!/usr/bin/env python3
"""Apply current user scope to retained readiness, never admit from length alone."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save


def scoped_rows(rows, excluded):
    result = []
    seen = set()
    for row in rows:
        entry_id = row["entry_id"]
        if entry_id in seen:
            raise ValueError("source readiness contains duplicate entries")
        seen.add(entry_id)
        if entry_id in excluded:
            status, reason = "excluded_by_user", "Explicit exclusion from current test cohort; retain research evidence."
        elif row["small_assembly_control"]:
            status, reason = "small_control_only", "Retained control, not large-structure cohort membership."
        else:
            status = "loop_scope_and_context_review_pending"
            reason = "Input readiness does not establish loop-scale construct/context or backend admission."
        result.append({**row, "current_scope_status": status, "current_scope_reason": reason,
            "current_test_admitted": False, "new_sampling_authorized": False,
            "final_split_assignment": None})
    return result


def apply(capacity, policy_path, output):
    inventory_path = capacity / "capacity-inventory.json"
    data = json.loads(inventory_path.read_text())
    for name, expected in data["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"capacity input changed: {name}")
    policy = json.loads(policy_path.read_text())
    if (policy["new_universal_numeric_gap_cutoff"] is not None
            or policy["historical_outcomes_and_frozen_requests_unchanged"] is not True):
        raise ValueError("scope requires explicit review, without numeric shortcut or mutation")
    rows = scoped_rows(data["rows"], set(policy["explicit_active_cohort_exclusions"]))
    summary = {"retained_ready_projections": len(rows),
        "current_scope_counts": dict(Counter(r["current_scope_status"] for r in rows)),
        "current_test_admitted": 0, "new_sampling_authorized": False,
        "historical_quota_table_is_not_current_execution_requirement": True,
        "replacement_loop_quota_table": None, "new_recruitment_allowed": False,
        "input_requests_modified": False, "historical_preparation_outcomes_unchanged": True}
    output.mkdir(parents=True, exist_ok=False)
    save(output / "loop-scope-inventory.json", {"rows": rows, "summary": summary,
        "source_capacity_sha256": digest(inventory_path), "user_policy_sha256": digest(policy_path),
        "scope_application_sha256": digest(Path(__file__))})
    save(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capacity", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    apply(args.capacity.resolve(), args.policy.resolve(), args.output.resolve())
