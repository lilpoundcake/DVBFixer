"""Monitoring must stay independent of scientific output and success."""

import json
import os

import pytest

from dvbfixer.job_progress import emit_job_progress, job_progress_descriptor


def test_optional_progress_does_not_write_logs(monkeypatch, capsys):
    monkeypatch.delenv("DVBFIXER_JOB_PROGRESS_FD", raising=False)
    emit_job_progress("preparation")
    assert capsys.readouterr() == ("", "")


def test_versioned_events_use_only_the_configured_pipe(monkeypatch, capsys):
    reader, writer = os.pipe()
    try:
        monkeypatch.setenv("DVBFIXER_JOB_PROGRESS_FD", str(writer))
        emit_job_progress("runner")
        emit_job_progress("validation", completed=2, total=3, seed=19)
        os.close(writer)
        writer = -1
        events = [json.loads(line) for line in os.read(reader, 4096).splitlines()]
        assert events == [
            {"version": 1, "stage": "runner"},
            {"version": 1, "stage": "validation", "completed": 2, "total": 3, "seed": 19},
        ]
        assert capsys.readouterr() == ("", "")
    finally:
        os.close(reader)
        if writer >= 0:
            os.close(writer)


@pytest.mark.parametrize("descriptor", ["-1", "1", "2", "invalid", "99999999", str(2**100)])
def test_unavailable_or_log_descriptors_do_not_interrupt_execution(monkeypatch, capsys, descriptor):
    monkeypatch.setenv("DVBFIXER_JOB_PROGRESS_FD", descriptor)
    emit_job_progress("publication")
    assert capsys.readouterr() == ("", "")


@pytest.mark.parametrize("completed,total", [(1, None), (None, 3), (-1, 3), (4, 3), (0, 0), (True, 2), (1, 2**53)])
def test_invalid_counts_are_rejected(completed, total):
    with pytest.raises(ValueError):
        emit_job_progress("sampling", completed=completed, total=total)


def test_closed_monitor_pipe_does_not_fail_calculation(monkeypatch):
    reader, writer = os.pipe()
    os.close(reader)
    try:
        monkeypatch.setenv("DVBFIXER_JOB_PROGRESS_FD", str(writer))
        emit_job_progress("validation")
    finally:
        os.close(writer)


def test_regular_files_are_not_accepted_as_monitor_pipes(tmp_path, monkeypatch):
    target = tmp_path / "must-not-change"
    with target.open("wb") as stream:
        monkeypatch.setenv("DVBFIXER_JOB_PROGRESS_FD", str(stream.fileno()))
        assert job_progress_descriptor() is None
        emit_job_progress("sampling")
    assert target.read_bytes() == b""


@pytest.mark.parametrize("seed", [-1, True, 0.5, 2**53])
def test_invalid_seed_identifiers_are_rejected(seed):
    with pytest.raises(ValueError):
        emit_job_progress("sampling", seed=seed)
