#!/usr/bin/env python3
"""Serial retained-data continuation; blocking lock waits, no status polling."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import signal
import subprocess
import sys
from pathlib import Path

from scripts.continue_real_gap_input_preparation import declare as declare_continuation
from scripts.mac_reconstruction_comparison import digest, save
from scripts.prepare_real_gap_inputs import declare as declare_preparation
from scripts.run_mac_modeller_comparison import execute

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"
SHARED_LOCK = PROJECT / ".artifacts/mac-reconstruction-2026-11-11/ml.lock"


def wait_for_executor(output, stage):
    # This is scheduler coordination, not repeated reads of a running job's files.
    print(json.dumps({"stage": stage, "waiting_on_shared_executor": True}), flush=True)
    with SHARED_LOCK.open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
    # The unchanged executor acquires its own lock. A competing reservation wins
    # fail-closed; never bypass its lock or adopt a resulting partial workdir.
    print(json.dumps({"stage": stage, "shared_executor_released": True}), flush=True)


def settled_and_pending(plan_root):
    plan = json.loads((plan_root / "plan.json").read_text())
    settled, pending = [], []
    for job in plan["jobs"]:
        work = plan_root / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        receipt = work / "execution.json"
        if not receipt.exists():
            if work.exists():
                raise ValueError("uncommitted workdir cannot be adopted by queue")
            pending.append(job["case_id"])
            continue
        row = json.loads(receipt.read_text())
        if row["case_id"] != job["case_id"] or row["command"] != job["command"]:
            raise ValueError("execution receipt does not match declared source task")
        if row["guard_stop_reason"] == "wall-time-limit" and not row["completed"]:
            status = "retained_time_block"
        elif row["completed"] and row["returncode"] == 0 and row["guard_stop_reason"] is None:
            if not (work / "result.json").is_file():
                raise ValueError("completion receipt lacks result")
            if row.get("result") != json.loads((work / "result.json").read_text()):
                raise ValueError("retained result differs from execution receipt")
            status = "completed_outcome"
        else:
            raise RuntimeError(f"technical/non-time resource failure requires investigation: {receipt}")
        settled.append({"entry_id": job["case_id"], "status": status, "receipt_sha256": digest(receipt)})
    return settled, pending


def execute_source_series(plan_root, output, label):
    current = plan_root
    history = []
    initial_count = json.loads((current / "plan.json").read_text())["count"]
    for index in range(initial_count+1):
        wait_for_executor(output, label)
        if not (current / "execution").exists():
            try:
                execute(current)
            except RuntimeError:
                # Only a retained per-case wall-time stop permits unchanged
                # first executions to continue. All other errors remain blockers.
                settled_and_pending(current)
        settled, pending = settled_and_pending(current)
        history.append({"plan": str(current), "settled": settled, "unexecuted": pending})
        save(output / f"{label}-accounting-{index:02d}.json", history[-1])
        if not pending:
            return history
        next_root = output / f"{label}-continuation-{index+1:02d}"
        declare_continuation(current, next_root)
        next_count = json.loads((next_root / "plan.json").read_text())["count"]
        if next_count != len(pending) or next_count >= initial_count:
            raise RuntimeError("continuation made no strict progress; no automatic retry")
        initial_count = next_count
        current = next_root
    raise RuntimeError("finite source continuation bound exhausted")


def execute_declared_stage(plan_root, output, label):
    wait_for_executor(output, label)
    if not (plan_root / "execution").exists():
        execute(plan_root)
        return
    completed_path = plan_root / "completed.json"
    if not completed_path.is_file():
        raise ValueError("partial non-source stage requires an explicit recovery amendment")
    plan = json.loads((plan_root / "plan.json").read_text())
    completed = json.loads(completed_path.read_text())
    if completed["completed"] != plan["count"] or completed["executed"] != plan["count"]:
        raise ValueError("existing non-source stage is incomplete")
    for name, expected in plan["source_code_and_input_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError("existing completed stage binding changed")
    settled, pending = settled_and_pending(plan_root)
    if pending or any(row["status"] != "completed_outcome" for row in settled):
        raise ValueError("non-source completion contains resource or missing outcomes")
    save(output / f"{label}-verified-reuse.json", {"plan": str(plan_root),
        "completion_sha256": digest(completed_path), "no_stage_repeated": True})


def declare(output):
    files = [Path(__file__), PROJECT / "scripts/continue_real_gap_input_preparation.py",
             PROJECT / "scripts/prepare_real_gap_inputs.py", PROJECT / "scripts/run_mac_modeller_comparison.py"]
    plans = {name: ROOT / name for name in ("source-preparation-hard-attempt01",
        "assembly-preflight-medium-attempt01", "retained-loop-protenix-native-preflight-attempt01")}
    files += [path / "plan.json" for path in plans.values()]
    output.mkdir(parents=True, exist_ok=False)
    save(output / "queue-plan.json", {"plans": {name: str(path) for name, path in plans.items()},
        "bindings": {str(path): digest(path) for path in files}, "python": sys.executable,
        "python_sha256": digest(Path(sys.executable)), "maximum_new_medium_source_tasks": 141,
        "source_continuation": "unchanged first executions only; retain wall-time blocks, never retry them",
        "per_child_guards": "existing serial1800s/RSS8GiB/pressure/headroom/10GiB executor",
        "full_fast_suite_timeout_seconds": 1800, "new_sampling_calls": 0, "larger_context_trials": 0,
        "new_source_acquisitions": 0, "new_model_assets": 0,
        "scientific_admission_or_final_cohort_not_implied": True})


def run(output):
    plan = json.loads((output / "queue-plan.json").read_text())
    for name, expected in plan["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"queue declaration changed: {name}")
    if digest(Path(plan["python"])) != plan["python_sha256"]:
        raise ValueError("queue interpreter changed")
    os.environ["PATH"] = str(Path(plan["python"]).parent)+os.pathsep+os.environ.get("PATH", "")
    # Exclusive scheduler reservation. A crash-left queue is not resumed silently.
    (output / "running-reservation").mkdir(exist_ok=False)
    try:
        paths = {name: Path(value) for name, value in plan["plans"].items()}
        history = execute_source_series(paths["source-preparation-hard-attempt01"], output, "hard-source")
        save(output / "hard-source-finished.json", {"series": history, "no_timeout_retried": True})
        medium = paths["assembly-preflight-medium-attempt01"]
        execute_declared_stage(medium, output, "medium-assembly")
        # Correctness preflight uses its own pinned interpreter and existing data.
        native = paths["retained-loop-protenix-native-preflight-attempt01"]
        execute_declared_stage(native, output, "native-data-preflight")
        source_root = output / "medium-source-preparation"
        declare_preparation(medium, source_root)
        if json.loads((source_root / "plan.json").read_text())["count"] > plan["maximum_new_medium_source_tasks"]:
            raise ValueError("source preparation exceeds declared queue task bound")
        history = execute_source_series(source_root, output, "medium-source")
        save(output / "medium-source-finished.json", {"series": history, "no_timeout_retried": True})
        # Test execution holds the shared lock, avoiding competition with studies.
        with SHARED_LOCK.open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            with (output / "full-fast-suite.log").open("xb") as log:
                process = subprocess.Popen([plan["python"], "-m", "pytest", "-m", "not slow", "-q"], cwd=PROJECT,
                    stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=plan["full_fast_suite_timeout_seconds"])
                except BaseException:
                    if process.poll() is None:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
                    raise
        save(output / "full-fast-suite-result.json", {"returncode": code,
            "log_sha256": digest(output / "full-fast-suite.log"), "optional_skips_are_reported_in_log": True})
        if code:
            raise RuntimeError("full fast suite failed; retained log requires investigation")
        save(output / "completed.json", {"planned_stages_completed": True, "new_sampling_calls": 0,
            "new_recruitment": False, "new_limits": False, "scientific_cohort_admission": False})
    except Exception as exc:
        save(output / "queue-blocked.json", {"type": type(exc).__name__, "message": str(exc),
            "remaining_stages_not_launched": True, "automatic_retry": False})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    (declare if args.mode == "declare" else run)(args.output.resolve())
