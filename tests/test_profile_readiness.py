"""Selected-profile failures and deferred license state must remain distinct."""

import subprocess
from pathlib import Path

import pytest

from dvbfixer import doctor, profile_readiness


def simulated_profile(monkeypatch, *, missing=None, unlicensed=False):
    monkeypatch.setattr(Path, "glob", lambda *args: iter(()))
    monkeypatch.setattr(profile_readiness.importlib.util, "find_spec",
                        lambda name: None if name in {"torch", "protenix"} else object())
    monkeypatch.setattr(profile_readiness.shutil, "which", lambda name: None if name == missing else f"/env/bin/{name}")

    def run(argv, **kwargs):
        failed = unlicensed and argv[-1] == "import modeller"
        return subprocess.CompletedProcess(argv, int(failed), b"", b"private diagnostic")

    monkeypatch.setattr(profile_readiness.subprocess, "run", run)
    return profile_readiness.collect_profile_readiness("standard")


def test_unselected_optional_tools_are_not_profile_errors(monkeypatch):
    report = simulated_profile(monkeypatch, missing="probe")
    assert report["state"] == "ready"
    assert "probe" not in {entry["name"] for entry in report["checks"]}


def test_modeller_license_is_action_required_without_raw_diagnostics(monkeypatch):
    report = simulated_profile(monkeypatch, unlicensed=True)
    assert report["state"] == "action-required"
    assert "private diagnostic" not in str(report)


def test_required_executable_failure_dominates_deferred_license(monkeypatch):
    report = simulated_profile(monkeypatch, missing="reduce", unlicensed=True)
    assert report["state"] == "error"
    assert next(entry for entry in report["checks"] if entry["name"] == "reduce")["state"] == "error"


@pytest.mark.parametrize("state,code", [("error", 1), ("action-required", 2)])
def test_profile_doctor_exit_code_and_json_state(monkeypatch, capsys, state, code):
    monkeypatch.setattr(doctor, "collect_capabilities", lambda **kwargs: {})
    monkeypatch.setattr(profile_readiness, "collect_profile_readiness",
                        lambda *args: {"state": state, "checks": []})
    with pytest.raises(SystemExit) as exc:
        doctor.main(["--profile", "standard", "--format", "json"])
    assert exc.value.code == code
    assert f'"state": "{state}"' in capsys.readouterr().out
