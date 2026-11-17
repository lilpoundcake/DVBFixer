#!/usr/bin/env python3
"""Read-only quota/construct/grouping inventory, not final scientific admission."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import DiffusionRequest

BINS = ((1, 8, 30), (9, 24, 50), (25, 49, 50), (50, 99, 40), (100, None, 30))


def task_inventory(request):
    lengths = {t.chain: len(t.sequence) for t in request.target_sequences}
    spans = []
    for gap in request.gaps:
        start, stop = gap.target_interval.start, gap.target_interval.stop
        kind = ("whole-chain" if start == 0 and stop == lengths[gap.chain] else
                "N-terminal" if start == 0 else "C-terminal" if stop == lengths[gap.chain] else "internal")
        spans.append({"chain": gap.chain, "start": start, "stop": stop, "length": stop-start, "kind": kind})
    observed = sum(len(p.observed_residues) for p in request.sequence_placements)
    targets = sum(lengths.values())
    return {"target_residues": targets, "observed_residues": observed, "chains": len(lengths),
        "small_assembly_control": observed < 300 and targets < 400,
        "gaps": spans, "longest_gap": max((g["length"] for g in spans), default=0),
        "terminal_construct_review_required": any(g["kind"] != "internal" for g in spans),
        "internal_gaps_present": any(g["kind"] == "internal" for g in spans),
        "generated_atoms_if_sampled": len(request.generated_atoms), "observed_atoms": len(request.fixed_atoms)}


def capacity(rows):
    result = []
    for lo, hi, quota in BINS:
        subset = [r for r in rows if lo <= r["longest_gap"] and (hi is None or r["longest_gap"] <= hi)]
        large = [r for r in subset if not r["small_assembly_control"]]
        result.append({"longest_gap_bin": [lo, hi], "required_entries": quota,
            "ready_projections": len(subset), "not_small_controls": len(large),
            "minimum_additional_projections_before_other_admission": max(0, quota-len(large)),
            "union_components_represented": len({r["union_component"] for r in large}),
            "complete_30pct_membership_entries": sum(r["group_membership_complete"] for r in large),
            "terminal_review_needed": sum(r["terminal_construct_review_required"] for r in large)})
    return result


def audit(grouped, ledgers, placements, output):
    bindings = {}

    def read(path):
        bindings[str(path)] = digest(path)
        return json.loads(path.read_text())

    complete = read(grouped / "completed.json")
    grouped_path = grouped / "grouped-census.json"
    if digest(grouped_path) != complete["grouped_census_sha256"]:
        raise ValueError("union grouping changed")
    groups = {r["entry_id"]: r for r in read(grouped_path)["rows"]}
    checked = {}
    for directory in placements:
        for row in read(directory / "placement-readback.json")["rows"]:
            if row["status"] != "deposited-axis-checked":
                continue
            if row["entry_id"] in checked or not row["all_placements_match_deposited_axis"]:
                raise ValueError("repeated or inconsistent independent placement evidence")
            checked[row["entry_id"]] = row
    rows, seen = [], set()
    for directory in ledgers:
        ledger = read(directory / "ledger.json")
        for name, expected in ledger["bindings"].items():
            if digest(Path(name)) != expected:
                raise ValueError(f"settled ledger artifact changed: {name}")
        for row in ledger["rows"]:
            entry_id = row["entry_id"]
            if entry_id in seen:
                raise ValueError("settled ledgers repeat an original source task")
            seen.add(entry_id)
            if row["status"] != "private_source_projection_prepared":
                continue
            evidence = checked.get(entry_id)
            if evidence is None:
                raise ValueError("ready projection lacks independent deposited-axis audit")
            result = read(Path(row["result_path"]))
            request_path = Path(result["request_path"])
            source_path = Path(result["source_input"])
            if (digest(request_path) != result["request_sha256"] or digest(source_path) != result["source_input_sha256"]
                    or digest(request_path) != evidence["request_sha256"] or digest(source_path) != evidence["source_sha256"]):
                raise ValueError("ready source/request differs from independent audit")
            bindings[str(request_path)] = digest(request_path)
            bindings[str(source_path)] = digest(source_path)
            request = DiffusionRequest.from_json(request_path.read_text())
            if request.normalized_pdb.sha256 != digest(source_path):
                raise ValueError("request bound to other source bytes")
            inventory = task_inventory(request)
            if [g["length"] for g in inventory["gaps"]] != result["gap_lengths"]:
                raise ValueError("result gap lengths differ from source-bound request")
            if inventory["small_assembly_control"] != evidence["small_assembly_control_not_large_context_evidence"]:
                raise ValueError("actual assembly size differs from independent readback")
            group = groups[entry_id]
            rows.append({"entry_id": entry_id, **inventory, "union_component": group["conservative_entry_component"],
                "group_membership_complete": group["modern_30pct_grouping_complete"],
                "prepared_source_derivative": result["prepared_source_derivative"],
                "original_context_review_complete": False, "scientific_backend_admission_complete": False,
                "request_path": str(request_path), "request_sha256": digest(request_path)})
    large = [r for r in rows if not r["small_assembly_control"]]
    quotas = capacity(rows)
    if sum(q["ready_projections"] for q in quotas) != len(rows):
        raise ValueError("longest-gap bins do not partition ready projections exactly once")
    summary = {"settled_tasks_accounted": len(seen), "independently_checked_ready_projections": len(rows),
        "small_assembly_controls": len(rows)-len(large), "not_small_controls": len(large),
        "union_components_among_not_small": len({r["union_component"] for r in large}),
        "incomplete_group_membership_not_small": sum(not r["group_membership_complete"] for r in large),
        "terminal_construct_review_needed_not_small": sum(r["terminal_construct_review_required"] for r in large),
        "multichain_not_small": sum(r["chains"] > 1 for r in large),
        "multigap_not_small": sum(len(r["gaps"]) > 1 for r in large),
        "at_least1000_observed_residues_not_small": sum(r["observed_residues"] >= 1000 for r in large),
        "quotas": quotas, "final_cohort_selected": False, "split_assignment": None,
        "counts_are_upper_bounds_before_context_construct_and_backend_admission": True}
    output.mkdir(parents=True, exist_ok=False)
    save(output / "capacity-inventory.json", {"rows": rows, "summary": summary, "bindings": bindings,
        "auditor_sha256": digest(Path(__file__)),
        "component_counts_are_descriptive_not_one_entry_per_component_restriction": True})
    save(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grouped", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, action="append", required=True)
    parser.add_argument("--placement", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.grouped.resolve(), [p.resolve() for p in args.ledger],
          [p.resolve() for p in args.placement], args.output.resolve())
