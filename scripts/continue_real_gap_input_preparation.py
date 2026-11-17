#!/usr/bin/env python3
"""Declare unexecuted source-preparation jobs; never retry a timed-out case."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save


def declare(previous, output):
    original_plan = previous / "plan.json"
    plan = json.loads(original_plan.read_text())
    for name, expected in plan["source_code_and_input_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"initial source-preparation code/input changed: {name}")
    if digest(Path(plan["python"])) != plan["python_sha256"]:
        raise ValueError("initial source-preparation interpreter changed")
    reused, blocked, jobs, bindings = [], [], [], dict(plan["source_code_and_input_sha256"])
    bindings[str(original_plan)] = digest(original_plan)
    bindings[str(Path(__file__))] = digest(Path(__file__))
    for job in plan["jobs"]:
        work = previous / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        receipt_path = work / "execution.json"
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            bindings[str(receipt_path)] = digest(receipt_path)
            if receipt["completed"] and receipt["returncode"] == 0 and receipt["guard_stop_reason"] is None:
                result_path = work / "result.json"
                row = json.loads(result_path.read_text())
                bindings[str(result_path)] = digest(result_path)
                if row["status"] == "private_source_projection_prepared":
                    for path_key, hash_key in (("source_input", "source_input_sha256"), ("request_path", "request_sha256")):
                        path = Path(row[path_key])
                        if not path.resolve().is_relative_to(work.resolve()) or digest(path) != row[hash_key]:
                            raise ValueError("completed source-preparation artifact differs or escapes its workdir")
                        bindings[str(path)] = row[hash_key]
                reused.append({"entry_id": job["case_id"], "result_path": str(result_path),
                               "result_sha256": digest(result_path), "status": row["status"]})
            elif receipt["guard_stop_reason"] == "wall-time-limit" and not receipt["completed"]:
                blocked.append({"entry_id": job["case_id"], "receipt": str(receipt_path),
                    "receipt_sha256": digest(receipt_path), "not_retried": True,
                    "classification": "source-preparation time limit, not a neural gap-length failure"})
            else:
                raise ValueError("unexpected prior technical/resource outcome needs a separate amendment")
        else:
            if work.exists():
                raise ValueError("uncommitted prior workdir must not be adopted")
            command = list(job["command"])
            command[command.index("--output")+1] = str(output / "execution" / work.name)
            jobs.append({**job, "command": command})
    if len(reused)+len(blocked)+len(jobs) != plan["count"]:
        raise ValueError("continuation does not account for every original source task")
    output.mkdir(parents=True, exist_ok=False)
    save(output / "amendment.json", {"previous_plan": str(original_plan), "previous_plan_sha256": digest(original_plan),
        "completed_outcomes_reused": reused, "resource_blocked_cases_retained": blocked,
        "remaining_first_executions": len(jobs), "no_completed_case_repeated": True,
        "no_timeout_budget_increase": True, "no_neural_sampling": True,
        "reason": "continue unchanged unexecuted source tasks after an isolated retained1800s preparation timeout"})
    save(output / "plan.json", {**plan, "jobs": jobs, "count": len(jobs), "source_code_and_input_sha256": bindings,
        "continuation_of": str(previous), "reused_count": len(reused), "resource_blocked_count": len(blocked)})
    print(json.dumps({"reused": len(reused), "blocked_not_retried": len(blocked), "new_first_execution_jobs": len(jobs)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    declare(args.previous.resolve(), args.output.resolve())
