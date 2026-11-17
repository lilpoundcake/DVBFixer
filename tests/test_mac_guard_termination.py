import subprocess
from types import SimpleNamespace

import pytest
from scripts import run_mac_reconstruction_guarded as guard


def test_dead_child_permission_race_is_reaped_and_reported(monkeypatch):
    def denied(*args):
        raise PermissionError("dead group permission race")

    waits = []
    monkeypatch.setattr(guard.os, "killpg", denied)
    process = SimpleNamespace(pid=123, wait=lambda **kw: waits.append(kw) or -15)
    assert "permission race" in guard.stop_owned_process(process)
    assert waits == [{"timeout": 10}]


def test_permission_denial_for_live_child_is_not_treated_as_success(monkeypatch):
    def denied(*args):
        raise PermissionError("denied")

    def still_running(**kwargs):
        raise subprocess.TimeoutExpired("owned-child", 10)

    monkeypatch.setattr(guard.os, "killpg", denied)
    with pytest.raises(RuntimeError, match="remains unconfirmed"):
        guard.stop_owned_process(SimpleNamespace(pid=123, wait=still_running))


def test_unresponsive_owned_child_uses_bounded_sigkill(monkeypatch):
    signals, waits = [], []
    monkeypatch.setattr(guard.os, "killpg", lambda pid, sig: signals.append((pid, sig)))

    def wait(**kwargs):
        waits.append(kwargs)
        if len(waits) == 1:
            raise subprocess.TimeoutExpired("owned-child", 10)
        return -9

    assert guard.stop_owned_process(SimpleNamespace(pid=123, wait=wait)) is None
    assert signals == [(123, guard.signal.SIGTERM), (123, guard.signal.SIGKILL)]
    assert waits == [{"timeout": 10}, {"timeout": 10}]
