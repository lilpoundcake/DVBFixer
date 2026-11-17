#!/usr/bin/env python3
"""Complete readback/export of a settled bounded batch, never replay generation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.audit_retained_full_context_controls import audit_saved_analysis
from scripts.export_retained_context_review import export
from scripts.mac_reconstruction_comparison import digest, save
from scripts.recover_retained_full_context_analysis import run as analyze
from scripts.retained_evidence import verify_bindings


def declare_analysis(previous, output):
    plan_path = previous / "plan.json"
    plan = json.loads(plan_path.read_text())
    resolved = verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    if digest(Path(plan["python"])) != plan["python_sha256"]:
        raise ValueError("declared interpreter changed")
    files = [plan_path, *resolved.values(), Path(__file__), Path(__file__).with_name("recover_retained_full_context_analysis.py"),
             Path(__file__).with_name("retained_modeller_full_context.py")]
    jobs, failures = [], []
    for job in plan["jobs"]:
        work = previous / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        receipt_path = work / "execution.json"
        if not receipt_path.is_file():
            failures.append({"entry_id": job["case_id"], "status": "not_executed"})
            continue
        execution = json.loads(receipt_path.read_text())
        files.append(receipt_path)
        if execution["command"] != job["command"] or execution["case_id"] != job["case_id"]:
            raise ValueError("native execution does not match declaration")
        if not execution["completed"] or execution["returncode"] != 0 or execution["guard_stop_reason"] is not None:
            failures.append({"entry_id": job["case_id"], "status": "retained_failure", "execution": str(receipt_path)})
            continue
        stages_path = work / "native/stages.json"
        stages = json.loads(stages_path.read_text())
        files += [stages_path, work / "native/native-interface-calls-final.json"]
        for row in stages["stages"]:
            label = row["stage"]
            filename = {"raw": "native.pdb", "generated-only-CG": "generated-only-CG-native.pdb"}[label]
            native = work / "native" / filename
            files.append(native)
            jobs.append({"entry_id": job["case_id"], "seed": job["seed"], "stage": label, "case": job["case"],
                         "native": str(native), "original_execution": str(receipt_path)})
    if not jobs:
        raise ValueError("no completed native stages; do not manufacture analysis success")
    output.mkdir(parents=True, exist_ok=False)
    save(output / "plan.json", {"jobs": jobs, "bindings": {str(path): digest(path) for path in files},
         "settled_execution_failures": failures, "model_calls": 0, "optimizer_calls": 0,
         "original_declaration": str(plan_path), "no_failed_case_replayed": True})


def finish(previous, output):
    output.mkdir(parents=True, exist_ok=False)
    analysis, readback = output / "analysis", output / "independent-readback"
    declare_analysis(previous, analysis)
    analyze(analysis)
    audit_saved_analysis(analysis, readback)
    rows = json.loads((readback / "receipt.json").read_text())["rows"]
    # One immutable final stage per case, never retrospectively select a donor.
    accepted = [row for row in rows if row["accepted"]]
    if accepted:
        export(analysis, readback, output / "review")
    save(output / "completion.json", {"original_run": str(previous), "independent_rows": [
        {"entry_id": row["entry_id"], "stage": row["stage"], "accepted": row["accepted"],
         "native_observations_changed": row["changed_native_observations"], "new_geometry_findings": len(row["new_geometry_findings"])}
        for row in rows], "settled_execution_failures": json.loads((analysis / "plan.json").read_text())["settled_execution_failures"],
        "exports_created": bool(accepted), "model_calls": 0, "optimizer_calls": 0,
        "public_capability_changed": False, "readback_sha256": digest(readback / "receipt.json")})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    finish(args.run.resolve(), args.output.resolve())
