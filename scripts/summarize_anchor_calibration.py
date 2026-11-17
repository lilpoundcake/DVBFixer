#!/usr/bin/env python3
"""Build a paired descriptive map only from completed independent grid readback."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save


def maximum(rows, name):
    if rows is None:
        return None  # Rejected/no final state is not zero observed motion.
    return max((r[name] for r in rows), default=0.)


def make_row(job, original, audited, execution):
    cell = audited["cell"]
    environment = original.get("environment", {})
    moves = audited.get("native_anchor_movements")
    serialized = audited.get("serialized_anchor_movements")
    return {"entry": audited["entry"], "cell": cell["id"], "window": cell["window"],
        "envelope": json.dumps(cell["envelope"]), "positional_k": cell["strength"],
        "final_candidate_available": audited["final_candidate_available"],
        "rejection": ";".join(audited.get("rejection", [])),
        "protein_geometry_passed": audited.get("protein_geometry_passed"),
        "protein_admission_failures": ";".join(audited.get("protein_admission_failures", [])),
        "movement_passed": audited.get("movement_passed"),
        "mask_and_motion_passed": audited.get("restoration_mask_and_motion_passed"),
        "whole_cli_no_new_findings": audited.get("full_cli_no_new_findings"),
        "current_cli_errors": audited.get("full_cli_counts", {}).get("ERROR", 0) if audited["final_candidate_available"] else None,
        "introduced_cli_findings": len(audited["introduced_cli_findings"]) if "introduced_cli_findings" in audited else None,
        "terminal_preparation_findings": len(audited["terminal_preparation_findings"]) if "terminal_preparation_findings" in audited else None,
        "max_anchor_backbone_rmsd": maximum(moves, "backbone_rmsd"),
        "max_anchor_ca_displacement": maximum(moves, "ca"),
        "max_anchor_heavy_displacement": maximum(moves, "heavy_max"),
        "serialized_max_anchor_heavy": maximum(serialized, "heavy_max"),
        "generated_rms_movement": original.get("generated_rms_movement_angstrom"),
        "generated_max_movement": original.get("generated_max_movement_angstrom"),
        "lost_observed_contacts": len(environment["lost_contacts"]) if "lost_contacts" in environment else None,
        "new_observed_contacts": len(environment["new_contacts"]) if "new_contacts" in environment else None,
        "intrasegment_C3_changes": len(environment["secondary_changes"]) if "secondary_changes" in environment else None,
        "execution_seconds": execution["elapsed_seconds"], "shared_preparation_seconds": job["shared_preparation_seconds"],
        "charged_seconds": execution["elapsed_seconds"]+job["shared_preparation_seconds"],
        "peak_rss_bytes": execution["peak_rss_bytes"],
        "main_reporter_last_index": next((c["last_iteration_index"] for c in original.get("minimizer_calls", [])
                                          if c["phase"] == "full-forcefield"), None),
        "internal_generated_issue_count": sum(i["touches_internal_generated"] for i in audited["issues"]) if "issues" in audited else None,
        "terminal_generated_issue_count": sum(i["touches_terminal_generated"] for i in audited["issues"]) if "issues" in audited else None}


def run(grid, readback, output):
    receipt = json.loads((readback / "receipt.json").read_text())
    plan = json.loads((grid / "plan.json").read_text())
    if receipt["grid_plan_sha256"] != digest(grid / "plan.json"):
        raise ValueError("readback does not bind this grid")
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for job in plan["jobs"]:
        work = grid / f"execution/{job['case_id']}-seed-{job['seed']}"
        audited = json.loads((readback / job["case_id"] / "result.json").read_text())
        if audited["original_result_sha256"] != digest(work / "result.json"):
            raise ValueError("original result changed after independent readback")
        original = json.loads((work / "result.json").read_text())
        execution = json.loads((work / "execution.json").read_text())
        rows.append(make_row(job, original, audited, execution))
    for row in rows:
        zero = next(r for r in rows if r["entry"] == row["entry"] and r["cell"] == "fixed")
        row["zero_final_available"] = zero["final_candidate_available"]
        for field in ("lost_observed_contacts", "new_observed_contacts", "intrasegment_C3_changes"):
            row[f"delta_vs_zero_{field}"] = row[field]-zero[field] if row[field] is not None and zero[field] is not None else None
    with (output / "tradeoffs.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save(output / "summary.json", {"cases": receipt["cases"], "declared_cells": len(rows), "blocked_inputs": plan["blocked"],
        "grid_plan_sha256": digest(grid / "plan.json"), "readback_receipt_sha256": digest(readback / "receipt.json"),
        "tradeoffs_csv_sha256": digest(output / "tradeoffs.csv"), "summary_code_sha256": digest(Path(__file__)),
        "geometry_thresholds_changed": False, "policy_selected": False, "optimizer_calls": 0,
        "limitations": ["four development-exposed structures, one retained proposal each, not100independentstructures",
                        "native generator fixation failures remain separate from source-restored input comparisons",
                        "intrasegment C3 excludes intersegment/interchain sheet assignment",
                        "a normal optimizer return does not demonstrate convergence",
                        "mask/motion pass is not whole-entry chemical completeness or MD readiness",
                        "no independent prospective validation cohort yet"]})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--readback", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.grid.resolve(), args.readback.resolve(), args.output.resolve())
