#!/usr/bin/env python3
"""Check shared request placements against deposited CIF label/author evidence."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import DiffusionRequest


def placement_evidence(request, inventory):
    targets = {target.chain: target.sequence for target in request.target_sequences}
    reports = []
    for placement in request.sequence_placements:
        observed_ids = [(r.residue_number, r.insertion_code) for r in placement.observed_residues]
        possible = []
        for chain in inventory["chains"]:
            if chain["sequence"] != targets[placement.chain]:
                continue
            axis = chain["observed_author_axis"]
            if [(r["author_residue_number"], r["insertion_code"]) for r in axis] == observed_ids:
                possible.append((chain["label_chain"], tuple(r["label_seq_id"]-1 for r in axis)))
        expected = {indices for _, indices in possible}
        actual = tuple(placement.observed_target_indices)
        reports.append({"output_chain": placement.chain, "source_label_candidates": [c for c, _ in possible],
            "unambiguous_deposited_position_evidence": len(expected) == 1,
            "request_agrees_with_deposited_axis": len(expected) == 1 and actual in expected,
            "observed_residues": len(observed_ids), "actual_observed_target_indices": actual,
            "expected_target_axes": sorted(expected)})
    return reports


def audit(preparation, sources, output):
    source_path = sources / "source-results.json"
    source_complete = json.loads((sources / "completed.json").read_text())
    if digest(source_path) != source_complete["source_results_sha256"]:
        raise ValueError("original source inventory changed")
    inventory = {r["entry_id"]: r["inventory"] for r in json.loads(source_path.read_text())["rows"]}
    plan = json.loads((preparation / "plan.json").read_text())
    rows = []
    for job in plan["jobs"]:
        work = preparation / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        receipt_path = work / "execution.json"
        if not receipt_path.exists():
            rows.append({"entry_id": job["case_id"], "status": "not_executed", "not_imputed_as_pass": True})
            continue
        receipt = json.loads(receipt_path.read_text())
        if not receipt["completed"]:
            rows.append({"entry_id": job["case_id"], "status": "technical_or_resource_blocked",
                "guard_stop_reason": receipt["guard_stop_reason"], "receipt_sha256": digest(receipt_path)})
            continue
        if receipt["returncode"] != 0 or receipt["guard_stop_reason"] is not None:
            raise ValueError("inconsistent source-execution receipt")
        result_path = work / "result.json"
        result = json.loads(result_path.read_text())
        if result["status"] != "private_source_projection_prepared":
            rows.append({"entry_id": job["case_id"], "status": result["status"], "result_sha256": digest(result_path)})
            continue
        source, request_path = Path(result["source_input"]), Path(result["request_path"])
        if (not source.resolve().is_relative_to(work.resolve()) or not request_path.resolve().is_relative_to(work.resolve())
                or digest(source) != result["source_input_sha256"] or digest(request_path) != result["request_sha256"]):
            raise ValueError("prepared source/request bytes or containment differ")
        request = DiffusionRequest.from_json(request_path.read_text())
        if request.normalized_pdb.sha256 != digest(source):
            raise ValueError("request bound to different source bytes")
        chains = placement_evidence(request, inventory[job["case_id"]])
        observed = sum(c["observed_residues"] for c in chains)
        targets = sum(len(t.sequence) for t in request.target_sequences)
        rows.append({"entry_id": job["case_id"], "status": "deposited-axis-checked",
            "chains": chains, "all_placements_match_deposited_axis": all(c["request_agrees_with_deposited_axis"] for c in chains),
            "observed_residues_in_selected_assembly": observed, "target_residues_in_selected_assembly": targets,
            "small_assembly_control_not_large_context_evidence": observed < 300 and targets < 400,
            "gap_lengths": result["gap_lengths"], "source_sha256": digest(source), "request_sha256": digest(request_path),
            "not_backend_admission_or_generated_gap_success": True})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "placement-readback.json", {"rows": rows, "plan_sha256": digest(preparation / "plan.json"),
        "source_results_sha256": digest(source_path), "auditor_sha256": digest(Path(__file__)), "partial_cohort": True})
    checked = [r for r in rows if r["status"] == "deposited-axis-checked"]
    summary = {"declared_source_tasks": len(rows), "deposited_axes_checked": len(checked),
        "consistent": sum(r["all_placements_match_deposited_axis"] for r in checked),
        "small_assembly_controls": sum(r["small_assembly_control_not_large_context_evidence"] for r in checked),
        "blocked": sum(r["status"] == "technical_or_resource_blocked" for r in rows),
        "unexecuted": sum(r["status"] == "not_executed" for r in rows), "model_sampling": False}
    save(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation", required=True, type=Path)
    parser.add_argument("--sources", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    audit(args.preparation.resolve(), args.sources.resolve(), args.output.resolve())
