import json
import sys
from pathlib import Path

import pytest
from scripts.continue_real_gap_input_preparation import declare
from scripts.mac_reconstruction_comparison import digest


def previous_plan(tmp_path):
    previous = tmp_path / "previous"
    previous.mkdir()
    jobs = [{"case_id": name, "case": "input.pdb", "seed": 0,
             "command": [sys.executable, "worker", "--output", str(previous / "execution" / f"{name}-seed-0")]}
            for name in ("finished", "timeout", "unexecuted")]
    (previous / "plan.json").write_text(json.dumps({"python": sys.executable,
        "python_sha256": digest(Path(sys.executable)), "jobs": jobs, "count": 3,
        "source_code_and_input_sha256": {}}))
    completed = previous / "execution/finished-seed-0"
    completed.mkdir(parents=True)
    source, request = completed / "source.pdb", completed / "request.json"
    source.write_text("retained source")
    request.write_text("retained request")
    (completed / "result.json").write_text(json.dumps({"status": "private_source_projection_prepared",
        "source_input": str(source), "source_input_sha256": digest(source),
        "request_path": str(request), "request_sha256": digest(request)}))
    (completed / "execution.json").write_text(json.dumps({"completed": True, "returncode": 0, "guard_stop_reason": None}))
    blocked = previous / "execution/timeout-seed-0"
    blocked.mkdir()
    (blocked / "execution.json").write_text(json.dumps({"completed": False, "returncode": -15,
                                                      "guard_stop_reason": "wall-time-limit"}))
    return previous, source


def test_continuation_neither_repeats_completed_nor_retries_timeout(tmp_path):
    previous, source = previous_plan(tmp_path)
    output = tmp_path / "continuation"
    declare(previous, output)
    plan = json.loads((output / "plan.json").read_text())
    assert [j["case_id"] for j in plan["jobs"]] == ["unexecuted"]
    assert plan["reused_count"] == 1 and plan["resource_blocked_count"] == 1
    assert plan["jobs"][0]["command"][-1] == str(output / "execution/unexecuted-seed-0")
    assert source.read_text() == "retained source"
    assert json.loads((output / "amendment.json").read_text())["no_timeout_budget_increase"]


def test_changed_complete_artifact_and_uncommitted_workdir_fail_closed(tmp_path):
    previous, source = previous_plan(tmp_path)
    source.write_text("changed")
    with pytest.raises(ValueError, match="artifact differs"):
        declare(previous, tmp_path / "must-not-adopt")
    source.write_text("retained source")
    (previous / "execution/unexecuted-seed-0").mkdir()
    with pytest.raises(ValueError, match="uncommitted"):
        declare(previous, tmp_path / "must-not-overwrite")
