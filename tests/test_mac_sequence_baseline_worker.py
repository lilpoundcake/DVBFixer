import json

import pytest
from scripts import match_mac_sequence_baseline as baseline


def test_worker_uses_executor_reservation_without_adopting_existing_payload(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    (source / "source-restored.pdb").write_text("retained raw bytes\n")
    output = tmp_path / "worker"
    output.mkdir()
    for name in ("command.json", "output.log"):
        (output / name).write_text("")
    calls = []

    def clean_no_repair(audit, destination, *, platform):
        calls.append((audit, platform))
        destination.mkdir()

    monkeypatch.setattr(baseline, "refine", clean_no_repair)
    baseline.worker(source, output)
    result = json.loads((output / "result.json").read_text())
    assert calls == [(source, "CPU")]
    assert result["final_path"] == str(source / "source-restored.pdb")
    assert result["resampled"] is False
    assert result["repair_failure"] is None
    with pytest.raises(ValueError, match="fresh executor-owned"):
        baseline.worker(source, output)


def test_worker_refuses_unreserved_or_payload_directory(tmp_path):
    with pytest.raises(ValueError, match="fresh executor-owned"):
        baseline.worker(tmp_path, tmp_path / "absent")
    (tmp_path / "refined.pdb").write_text("must not adopt\n")
    with pytest.raises(ValueError, match="fresh executor-owned"):
        baseline.worker(tmp_path, tmp_path)
