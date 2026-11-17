import json
import sys
from pathlib import Path

import pytest
from scripts import run_retained_loop_correctness_queue as queue
from scripts.mac_reconstruction_comparison import digest


def plan(root):
    root.mkdir()
    jobs = [{"case_id": name, "case": "retained.pdb", "seed": 0,
        "command": [sys.executable, "worker", "--output", str(root / "execution" / f"{name}-seed-0")]}
        for name in ("finished", "pending")]
    (root / "plan.json").write_text(json.dumps({"python": sys.executable,
        "python_sha256": digest(Path(sys.executable)), "source_code_and_input_sha256": {},
        "jobs": jobs, "count": 2}))
    receipt(root, jobs[0])
    return jobs


def receipt(root, job, *, blocked=False, technical=False):
    work = root / "execution" / f"{job['case_id']}-seed-{job['seed']}"
    work.mkdir(parents=True)
    if not blocked and not technical:
        (work / "result.json").write_text(json.dumps({"status": "source_preparation_refused", "entry_id": job["case_id"]}))
    result = {"result": json.loads((work / "result.json").read_text())} if not blocked and not technical else {}
    (work / "execution.json").write_text(json.dumps({**job, **result, "completed": not blocked and not technical,
        "returncode": -15 if blocked else 1 if technical else 0,
        "guard_stop_reason": "wall-time-limit" if blocked else None}))


def test_queue_continues_only_unexecuted_tasks_without_repeating_finished(tmp_path, monkeypatch):
    original, output = tmp_path / "original", tmp_path / "queue"
    plan(original)
    output.mkdir()
    calls = []
    monkeypatch.setattr(queue, "wait_for_executor", lambda *args: None)

    def execute(root):
        calls.append(root)
        jobs = json.loads((root / "plan.json").read_text())["jobs"]
        assert [j["case_id"] for j in jobs] == ["pending"]
        for job in jobs:
            receipt(root, job)

    monkeypatch.setattr(queue, "execute", execute)
    history = queue.execute_source_series(original, output, "test")
    assert len(calls) == 1
    assert not history[-1]["unexecuted"]
    assert history[0]["settled"][0]["entry_id"] == "finished"


def test_time_block_is_retained_not_retried(tmp_path):
    root = tmp_path / "series"
    jobs = plan(root)
    receipt(root, jobs[1], blocked=True)
    settled, pending = queue.settled_and_pending(root)
    assert not pending
    assert settled[-1]["status"] == "retained_time_block"


def test_technical_failure_or_uncommitted_directory_is_not_adopted(tmp_path):
    root = tmp_path / "series"
    jobs = plan(root)
    receipt(root, jobs[1], technical=True)
    with pytest.raises(RuntimeError, match="technical"):
        queue.settled_and_pending(root)
    other = tmp_path / "other"
    jobs = plan(other)
    (other / "execution/pending-seed-0").mkdir()
    with pytest.raises(ValueError, match="uncommitted"):
        queue.settled_and_pending(other)
