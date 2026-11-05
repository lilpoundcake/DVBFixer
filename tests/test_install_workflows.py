"""Installer UX, explicit updates and observable long-running commands."""

import json
import subprocess
import sys
import threading

import pytest

from dvbfixer import install, install_discovery, runtime_install
from dvbfixer.install_progress import logged_run, stage


def interactive(monkeypatch, answers):
    prompts = []
    iterator = iter(answers)
    monkeypatch.setattr(sys.stdin, "isatty", lambda: True)

    def answer(prompt):
        prompts.append(prompt)
        return next(iterator)

    monkeypatch.setattr("builtins.input", answer)
    return prompts


def manager(monkeypatch, tmp_path):
    monkeypatch.setattr(install.shutil, "which", lambda _: "/tools/micromamba")
    monkeypatch.setattr(install.Manager, "validate", lambda _: "2.9.0")
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([tmp_path / "envs"], []))
    discover = install_discovery.discover_installations
    monkeypatch.setattr(install_discovery, "discover_installations",
                        lambda dirs, envs, **kwargs: discover(dirs, envs, repository=tmp_path))
    monkeypatch.chdir(tmp_path)


def existing(tmp_path, profile="standard", addons=()):
    prefix = tmp_path / "installed"
    install.reserve_prefix(prefix, profile, addons, False)
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta/history").write_text("existing history")
    return prefix


def test_initial_menu_separates_profile_from_backend(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    prompts = interactive(monkeypatch, ["install", "standard", "micromamba", "new", ""])
    assert install.main(["--dry-run"]) == 0
    assert prompts[0] == "Choice: "
    assert prompts[1] == "Choice: "
    output = capsys.readouterr().out
    assert "1. New Standard/Desktop installation" in output
    assert "1. Standard" in output
    assert all("runtime [" not in prompt and "checkpoint" not in prompt for prompt in prompts)
    assert all("directory" not in prompt.lower() and "path" not in prompt.lower() for prompt in prompts)
    assert "Directory:" in output
    assert "PREVIEW COMPLETED —" in output
    assert "SUCCESS —" not in output
    assert not (tmp_path / "new.dvbfixer-install.json").exists()


def test_new_name_rejects_paths_and_uses_manager_location(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    prompts = interactive(monkeypatch, ["install", "standard", "micromamba", "not/a/name", "dvbfixer", ""])
    assert install.main(["--dry-run"]) == 0
    assert str(tmp_path / "envs/dvbfixer") in capsys.readouterr().out
    assert not (tmp_path / "envs").exists()
    assert all("directory" not in prompt.lower() for prompt in prompts)


def test_install_profile_and_manager_are_selectable_using_only_numbers(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    interactive(monkeypatch, ["1", "2", "1", "desktop-numeric", "2", "0"])
    assert install.main(["--dry-run"]) == 0
    output = capsys.readouterr().out
    assert "Install DVBFixer desktop" in output
    assert "1. micromamba" in output
    assert "1. Yes" in output and "2. No" in output
    assert not (tmp_path / "envs").exists()


def test_backend_selection_is_numbered_and_separate_from_profiles(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    interactive(monkeypatch, ["4", "4", "1", "mps-numeric", ""])
    calls = []

    def provision(**kwargs):
        calls.append(kwargs)
        return {"runtime": kwargs["name"], "prefix": str(kwargs["prefix"]), "state": "planned"}

    monkeypatch.setattr(runtime_install, "provision_runtime", provision)
    assert install.main(["--dry-run"]) == 0
    assert calls[0]["name"] == "protpardelle-mps"
    assert calls[0]["dry_run"] is True
    assert "4. protpardelle-mps" in capsys.readouterr().out


def test_resume_infers_backend_without_new_directory_or_profile_prompt(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    prefix = existing(tmp_path, "runtime:protpardelle-mps")
    prompts = interactive(monkeypatch, ["backend", "micromamba", "1", ""])
    calls = []

    def provision(**kwargs):
        calls.append(kwargs)
        return {"runtime": kwargs["name"], "prefix": str(prefix), "state": "action-required",
                "actions": ["Provide checkpoint"], "scientific_acceptance": "not-tested"}

    monkeypatch.setattr(runtime_install, "provision_runtime", provision)
    assert install.main(["--resume"]) == 2
    assert calls[0]["name"] == "protpardelle-mps"
    assert calls[0]["resume"] is True
    assert not any("NEW" in prompt or "Profile" in prompt or "profile [" in prompt for prompt in prompts)
    assert all("directory" not in prompt.lower() and "path" not in prompt.lower() for prompt in prompts if "checkpoint" not in prompt)
    output = capsys.readouterr().out
    assert "ACTION REQUIRED" in output and "not ready" in output
    assert "Provide checkpoint" in output
    summary = output.split("=== Installation result ===")[1]
    assert "ACTION REQUIRED —" in summary and "Provide checkpoint" in summary
    assert "managed download" in summary
    assert "SUCCESS —" not in output
    assert '"schema_version"' not in output


def test_resume_infers_profile_and_addons_noninteractively(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    prefix = existing(tmp_path, addons=("msa-alternatives",))
    assert install.main(["--resume", "--prefix", str(prefix), "--manager", "micromamba",
                         "--non-interactive", "--dry-run"]) == 0
    assert "Resume DVBFixer standard" in capsys.readouterr().out


def test_update_infers_composition_preserves_license_and_user_data(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    prefix = existing(tmp_path, addons=("msa-alternatives",))
    user_data = prefix / "user-project.json"
    user_data.write_text("untouched")
    before = prefix.with_name(prefix.name + ".dvbfixer-install.json").read_bytes()
    calls, licenses = [], []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        output = "old-version" if "from importlib.metadata import version; print(version('dvbfixer'))" in argv else ""
        if install.MODELLER_READ in argv:
            output = "protected-license"
        return subprocess.CompletedProcess(argv, 0, output, "")

    monkeypatch.setattr(install.subprocess, "run", run)
    monkeypatch.setattr(install, "run_checked", lambda argv, **kwargs: calls.append((argv, kwargs)))
    monkeypatch.setattr(install, "configure_modeller", lambda *args, **kwargs: licenses.append((args, kwargs)))
    assert install.main(["--update", "--manager", "micromamba", "--prefix", str(prefix), "--non-interactive"]) == 0
    assert licenses[0][0][2] == "protected-license"
    assert licenses[0][1]["verify_import"] is True
    assert any(argv[1] == "install" and "--file" in argv for argv, _ in calls)
    assert any("pip" in argv and "--upgrade" in argv for argv, _ in calls)
    assert any("--addon" in argv and "msa-alternatives" in argv for argv, _ in calls)
    assert user_data.read_text() == "untouched"
    assert prefix.with_name(prefix.name + ".dvbfixer-install.json").read_bytes() == before
    captured = capsys.readouterr()
    output = captured.out + captured.err
    assert "protected-license" not in output
    assert "old-version ->" in output and "READY:" in output


@pytest.mark.parametrize("kind", ["unowned", "runtime", "profile-mismatch", "addon-mismatch", "incomplete"])
def test_update_refuses_wrong_targets_before_any_command(tmp_path, monkeypatch, kind):
    if kind == "unowned":
        prefix = tmp_path / "manual"
        prefix.mkdir()
    else:
        prefix = existing(tmp_path, "runtime:protpardelle-mps" if kind == "runtime" else "standard")
    extra = {"profile-mismatch": ["--profile", "desktop"], "addon-mismatch": ["--addon", "msa-alternatives"]}.get(kind, [])
    if kind == "incomplete":
        (prefix / "conda-meta/history").unlink()
    monkeypatch.setattr(install.shutil, "which", lambda _: pytest.fail("must refuse before commands"))
    assert install.main(["--update", "--manager", "micromamba", "--prefix", str(prefix), "--non-interactive", *extra]) == 1
    assert not prefix.with_name(prefix.name + ".installation-logs").exists()


def test_update_confirmation_cancel_does_not_write(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    prefix = existing(tmp_path)
    prompts = interactive(monkeypatch, ["micromamba", "1", "no"])
    monkeypatch.setattr(install.subprocess, "run", lambda argv, **kwargs: subprocess.CompletedProcess(argv, 0, "0.1", ""))
    assert install.main(["--update"]) == 0
    assert not any("NEW" in prompt or "Add-ons" in prompt or "Profile" in prompt for prompt in prompts)
    assert "Cancelled" in capsys.readouterr().out
    assert not prefix.with_name(prefix.name + ".installation-logs").exists()


def test_progress_emits_heartbeat_without_fake_percent(capsys):
    with stage("silent solver", interval=0.01):
        threading.Event().wait(0.05)
    output = capsys.readouterr().out
    assert "Starting: silent solver" in output
    assert "Still running: silent solver" in output
    assert "Completed: silent solver" in output
    assert "%" not in output


def test_real_child_output_and_failure_are_retained(tmp_path, capsys):
    log = tmp_path / "package.log"
    with pytest.raises(subprocess.CalledProcessError):
        logged_run([sys.executable, "-c", "import sys; print('details'); print('error', file=sys.stderr); sys.exit(3)"],
                   log=log, env=install.child_environment())
    assert "details" in log.read_text() and "error" in log.read_text()
    output = capsys.readouterr().out
    assert str(log) in output and "Failed/interrupted:" in output
    assert "Completed:" not in output


def test_bad_record_is_not_adopted(tmp_path):
    prefix = existing(tmp_path)
    marker = prefix.with_name(prefix.name + ".dvbfixer-install.json")
    marker.write_text(json.dumps([]))
    with pytest.raises(ValueError, match="Invalid"):
        install.recorded_installation(prefix)


def test_desktop_update_rebuilds_only_selected_prefix_and_keeps_workspace(tmp_path, monkeypatch):
    manager(monkeypatch, tmp_path)
    prefix = existing(tmp_path, "desktop")
    gui = prefix / "share/dvbfixer/gui"
    workspace = gui / "structures/projects/example/manifest.json"
    workspace.parent.mkdir(parents=True)
    workspace.write_text("owned workspace")
    (gui / ".env").write_text("operator settings")
    calls = []

    def run(argv, **kwargs):
        text = "license" if install.MODELLER_READ in argv else "0.10.2"
        return subprocess.CompletedProcess(argv, 0, text, "")

    monkeypatch.setattr(install.subprocess, "run", run)
    monkeypatch.setattr(install, "configure_modeller", lambda *args, **kwargs: None)
    monkeypatch.setattr(install, "run_checked", lambda argv, **kwargs: calls.append((argv, kwargs)))
    # Exercise real copytree and its structures/.env exclusions without copying
    # an entire GUI tree; the helper's root is controlled through __file__.
    source = tmp_path / "source"
    module = source / "src/dvbfixer/install.py"
    module.parent.mkdir(parents=True)
    (source / "pyproject.toml").write_text("[project]\nname='test'\n")
    (source / "gui/structures").mkdir(parents=True)
    (source / "gui/structures/unrelated").write_text("do not import")
    (source / "gui/.env").write_text("do not replace operator settings")
    (source / "gui/package.json").write_text("{}")
    monkeypatch.setattr(install, "__file__", str(module))
    assert install.main(["--update", "--prefix", str(prefix), "--manager", "micromamba", "--non-interactive"]) == 0
    assert workspace.read_text() == "owned workspace"
    assert (gui / ".env").read_text() == "operator settings"
    assert not (gui / "structures/unrelated").exists()
    gui_calls = [(argv, kwargs) for argv, kwargs in calls if "npm" in argv]
    assert len(gui_calls) == 2
    assert all(kwargs["cwd"] == gui for _, kwargs in gui_calls)
    assert any(argv[-3:] == ["npm", "run", "build"] for argv, _ in gui_calls)


def test_update_cannot_read_license_refuses_before_package_mutation(tmp_path, monkeypatch, capsys):
    manager(monkeypatch, tmp_path)
    prefix = existing(tmp_path)

    def run(argv, **kwargs):
        return subprocess.CompletedProcess(argv, 1 if install.MODELLER_READ in argv else 0, "", "sensitive output")

    monkeypatch.setattr(install.subprocess, "run", run)
    monkeypatch.setattr(install, "run_checked", lambda *args, **kwargs: pytest.fail("must not update packages"))
    assert install.main(["--update", "--manager", "micromamba", "--prefix", str(prefix), "--non-interactive"]) == 1
    captured = capsys.readouterr()
    assert "Cannot safely preserve MODELLER" in captured.err
    assert "sensitive output" not in captured.out + captured.err


def test_interrupted_prompt_has_no_traceback_or_installation(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: (_ for _ in ()).throw(KeyboardInterrupt()))
    assert install.main([]) == 130
    assert "Cancelled" in capsys.readouterr().err


def test_invalid_interactive_backend_reprompts_and_can_cancel(tmp_path, monkeypatch):
    prompts = interactive(monkeypatch, ["backend", "bogus", "0"])
    assert install.main(["--dry-run"]) == 130
    assert prompts == ["Choice: "] * 3


@pytest.mark.parametrize("interrupted", [False, True])
def test_failed_update_attempts_private_license_recovery(tmp_path, monkeypatch, capsys, interrupted):
    manager(monkeypatch, tmp_path)
    prefix = existing(tmp_path)
    recovered = []

    def run(argv, **kwargs):
        output = "protected-key" if install.MODELLER_READ in argv else "0.1"
        return subprocess.CompletedProcess(argv, 0, output, "")

    def fail(*args, **kwargs):
        if interrupted:
            raise KeyboardInterrupt()
        raise subprocess.CalledProcessError(1, ["native-package-manager"])

    monkeypatch.setattr(install.subprocess, "run", run)
    monkeypatch.setattr(install, "run_checked", fail)
    monkeypatch.setattr(install, "configure_modeller", lambda *args, **kwargs: recovered.append((args, kwargs)))
    assert install.main(["--update", "--prefix", str(prefix), "--manager", "micromamba", "--non-interactive"]) == (
        130 if interrupted else 1
    )
    assert recovered[0][0][2] == "protected-key"
    assert recovered[0][1] == {"verify_import": False}
    captured = capsys.readouterr()
    assert "protected-key" not in captured.out + captured.err
    assert "READY:" not in captured.out


def test_version_snapshot_log_can_be_retained_without_dumping_json(tmp_path, capsys):
    log = tmp_path / "versions.log"
    with stage("Package inventory", log, follow=False):
        log.write_text('{"packages": "large snapshot"}')
    assert "large snapshot" not in capsys.readouterr().out
    assert "large snapshot" in log.read_text()
