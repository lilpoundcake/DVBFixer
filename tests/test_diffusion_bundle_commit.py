"""Manifest commit must survive unsupported rename flags, races and interruption."""

import json
import os
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from dvbfixer.model.diffusion import bundle_commit
from dvbfixer.model.diffusion.bundle_commit import commit_manifest_last


def _staging(root: Path, name: str = "stage") -> Path:
    staging = root / name
    staging.mkdir()
    (staging / "structure.pdb").write_text(name)
    (staging / "evidence.json").write_text('{"complete": true}')
    (staging / "bundle.json").write_text(json.dumps({"publisher": name}))
    return staging


def test_commit_index_is_last_and_every_payload_has_one_link(tmp_path, monkeypatch):
    staging = _staging(tmp_path)
    destination = tmp_path / "published"
    transfer = bundle_commit._transfer
    observed = []

    def inspect(source, name, descriptor, linked):
        assert not (destination / "bundle.json").exists()
        if name == "bundle.json":
            assert (destination / "structure.pdb").read_text() == "stage"
            assert json.loads((destination / "evidence.json").read_text())["complete"]
            assert all(p.stat().st_nlink == 1 for p in destination.iterdir())
        transfer(source, name, descriptor, linked)
        observed.append(name)

    monkeypatch.setattr(bundle_commit, "_transfer", inspect)
    commit_manifest_last(staging, destination)
    assert observed[-1] == "bundle.json"
    assert all(p.stat().st_nlink == 1 for p in destination.iterdir())
    assert not staging.exists()


@pytest.mark.parametrize("existing", ["empty", "nonempty", "file", "symlink"])
def test_existing_destination_is_never_adopted_or_replaced(tmp_path, existing):
    staging = _staging(tmp_path)
    destination = tmp_path / "published"
    if existing in {"empty", "nonempty"}:
        destination.mkdir()
        if existing == "nonempty":
            (destination / "keep").write_text("original")
    elif existing == "file":
        destination.write_text("original")
    else:
        destination.symlink_to(tmp_path / "missing")
    identity = destination.lstat()
    with pytest.raises(FileExistsError):
        commit_manifest_last(staging, destination)
    assert destination.lstat() == identity
    assert (staging / "structure.pdb").read_text() == "stage"


@pytest.mark.parametrize("fail_at", ["structure.pdb", "evidence.json", "bundle.json"])
def test_transfer_failure_rolls_back_only_our_destination(tmp_path, monkeypatch, fail_at):
    staging = _staging(tmp_path)
    destination = tmp_path / "published"
    transfer = bundle_commit._transfer

    def fail(source, name, descriptor, linked):
        if name == fail_at:
            raise OSError("injected transfer failure")
        transfer(source, name, descriptor, linked)

    monkeypatch.setattr(bundle_commit, "_transfer", fail)
    with pytest.raises(OSError, match="injected"):
        commit_manifest_last(staging, destination)
    assert not destination.exists()


@pytest.mark.parametrize("fail_call", [1, 2])
def test_fsync_failure_cannot_leave_an_index_without_payloads(tmp_path, monkeypatch, fail_call):
    staging = _staging(tmp_path)
    destination = tmp_path / "published"
    fsync = os.fsync
    calls = 0

    def fail(descriptor):
        nonlocal calls
        calls += 1
        if calls == fail_call:
            raise OSError("injected fsync failure")
        fsync(descriptor)

    monkeypatch.setattr(bundle_commit.os, "fsync", fail)
    with pytest.raises(OSError, match="injected"):
        commit_manifest_last(staging, destination)
    assert not destination.exists()


def test_competing_publishers_have_exactly_one_winner(tmp_path):
    first, second = _staging(tmp_path, "first"), _staging(tmp_path, "second")
    destination = tmp_path / "published"
    barrier = threading.Barrier(2)

    def commit(staging):
        barrier.wait(timeout=10)
        try:
            commit_manifest_last(staging, destination)
            return staging.name
        except FileExistsError:
            return None

    with ThreadPoolExecutor(max_workers=2) as executor:
        winners = list(executor.map(commit, (first, second)))
    winner, = [name for name in winners if name]
    assert (destination / "structure.pdb").read_text() == winner
    assert json.loads((destination / "bundle.json").read_text())["publisher"] == winner
    loser = second if winner == "first" else first
    assert (loser / "structure.pdb").read_text() == loser.name


def test_unrelated_injected_entry_is_preserved_on_failure(tmp_path, monkeypatch):
    staging = _staging(tmp_path)
    destination = tmp_path / "published"
    transfer = bundle_commit._transfer

    def inject(source, name, descriptor, linked):
        if name == "bundle.json":
            (destination / "keep").write_text("not owned by the publisher")
            raise KeyboardInterrupt
        transfer(source, name, descriptor, linked)

    monkeypatch.setattr(bundle_commit, "_transfer", inject)
    with pytest.raises(KeyboardInterrupt):
        commit_manifest_last(staging, destination)
    assert (destination / "keep").read_text() == "not owned by the publisher"
    assert not (destination / "bundle.json").exists()


def test_process_crash_leaves_unready_reservation_and_retry_preserves_it(tmp_path):
    staging = _staging(tmp_path)
    destination = tmp_path / "published"
    script = """
import os,sys
from pathlib import Path
from dvbfixer.model.diffusion import bundle_commit as b
transfer = b._transfer
def crash(source, name, descriptor, linked):
    if name == 'bundle.json':
        os._exit(23)
    transfer(source, name, descriptor, linked)
b._transfer = crash
b.commit_manifest_last(Path(sys.argv[1]), Path(sys.argv[2]))
"""
    env = {**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src")}
    result = subprocess.run([sys.executable, "-c", script, str(staging), str(destination)],
                            env=env, timeout=20, check=False)
    assert result.returncode == 23
    assert (destination / "structure.pdb").read_text() == "stage"
    assert not (destination / "bundle.json").exists()
    fresh = _staging(tmp_path, "retry")
    with pytest.raises(FileExistsError):
        commit_manifest_last(fresh, destination)
    assert (destination / "structure.pdb").read_text() == "stage"
