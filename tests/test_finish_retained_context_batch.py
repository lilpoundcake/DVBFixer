import json
import sys
from pathlib import Path

import pytest
from scripts import finish_retained_context_batch as batch
from scripts.mac_reconstruction_comparison import digest


def plan(tmp_path):
    previous = tmp_path / "previous"
    previous.mkdir()
    job = {"case_id": "3C7D", "seed": 104, "case": str(previous / "cases/3C7D"), "command": ["declared-worker"]}
    (previous / "plan.json").write_text(json.dumps({"python": sys.executable,
        "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {}, "jobs": [job]}))
    return previous, job


def test_unexecuted_case_does_not_become_an_analysis_success(tmp_path):
    previous, _ = plan(tmp_path)
    output = tmp_path / "analysis"
    with pytest.raises(ValueError, match="no completed native stages"):
        batch.declare_analysis(previous, output)
    assert not output.exists()


def test_execution_receipt_must_match_the_declared_command(tmp_path):
    previous, job = plan(tmp_path)
    work = previous / "execution/3C7D-seed-104"
    work.mkdir(parents=True)
    (work / "execution.json").write_text(json.dumps({**job, "command": ["other-worker"]}))
    with pytest.raises(ValueError, match="does not match declaration"):
        batch.declare_analysis(previous, tmp_path / "analysis")


@pytest.mark.parametrize("accepted", [False, True])
def test_finish_preserves_order_and_exports_only_independent_acceptance(tmp_path, monkeypatch, accepted):
    calls = []

    def declare(previous, output):
        calls.append("declare")
        output.mkdir()
        (output / "plan.json").write_text(json.dumps({"settled_execution_failures": []}))

    def analyze(output):
        calls.append("analyze")

    def audit(analysis, output):
        calls.append("audit")
        output.mkdir()
        (output / "receipt.json").write_text(json.dumps({"rows": [{"entry_id": "3C7D", "stage": "raw",
            "accepted": accepted, "changed_native_observations": 0, "new_geometry_findings": [] if accepted else [{}]}]}))

    def export(analysis, readback, output):
        calls.append("export")
        assert (readback / "receipt.json").is_file()
        output.mkdir()

    monkeypatch.setattr(batch, "declare_analysis", declare)
    monkeypatch.setattr(batch, "analyze", analyze)
    monkeypatch.setattr(batch, "audit_saved_analysis", audit)
    monkeypatch.setattr(batch, "export", export)
    output = tmp_path / "finished"
    batch.finish(tmp_path / "previous", output)
    assert calls == ["declare", "analyze", "audit"] + (["export"] if accepted else [])
    result = json.loads((output / "completion.json").read_text())
    assert result["exports_created"] is accepted
    assert result["model_calls"] == result["optimizer_calls"] == 0
    assert (output / "review").exists() is accepted
