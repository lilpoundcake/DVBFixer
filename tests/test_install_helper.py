"""Protect existing environments and secrets; portable manager command contracts."""

import os
from pathlib import Path

import pytest

from dvbfixer import install
from dvbfixer.install import (
    MODELLER_SETUP,
    Manager,
    child_environment,
    configure_modeller,
    reserve_prefix,
    validate_prefix_path,
)


@pytest.mark.parametrize("name,subcommand", [("micromamba", ["create"]),
                                           ("mamba", ["env", "create"]),
                                           ("conda", ["env", "create"])])
def test_manager_argv_keeps_paths_with_spaces(name, subcommand):
    manager = Manager(name, f"/tools with spaces/{name}")
    prefix = Path("/env with spaces/standard")
    spec = Path("/tmp with spaces/spec.yml")
    channel_flags = ["--override-channels"] if name in ("micromamba", "mamba") else []
    assert manager.create(prefix, spec) == [manager.executable, *subcommand, *channel_flags, "--yes", "--prefix",
                                           str(prefix), "--file", str(spec)]
    flags = ["--no-capture-output"] if name == "conda" else []
    assert manager.run(prefix, ["python", "-m", "pip"]) == [manager.executable, "run", *flags,
                                                           "--prefix", str(prefix), "python", "-m", "pip"]


def test_existing_prefix_is_never_adopted(tmp_path):
    prefix = tmp_path / "existing"
    prefix.mkdir()
    with pytest.raises(ValueError, match="already exists"):
        reserve_prefix(prefix, "standard", (), False)
    with pytest.raises(ValueError, match="no installation record"):
        reserve_prefix(prefix, "standard", (), True)


def test_resume_binds_profile_and_prefix(tmp_path):
    prefix = tmp_path / "new standard"
    marker = reserve_prefix(prefix, "standard", (), False)
    assert reserve_prefix(prefix, "standard", (), True) == marker
    with pytest.raises(ValueError, match="differ"):
        reserve_prefix(prefix, "desktop", (), True)
    with pytest.raises(FileExistsError):
        reserve_prefix(prefix, "standard", (), False)


@pytest.mark.parametrize("name,subcommand", [("micromamba", ["install", "--yes"]),
                                           ("mamba", ["env", "update"]),
                                           ("conda", ["env", "update"])])
def test_resume_reconciles_interrupted_native_transaction(name, subcommand):
    manager = Manager(name, name)
    channel_flags = ["--override-channels"] if name in ("micromamba", "mamba") else []
    assert manager.update(Path("/target"), Path("/spec.yml")) == [
        name, *subcommand, *channel_flags, "--prefix", "/target", "--file", "/spec.yml"
    ]


def test_target_symlink_rejected(tmp_path):
    prefix = tmp_path / "link"
    prefix.symlink_to(tmp_path / "other")
    with pytest.raises(ValueError, match="symlink"):
        reserve_prefix(prefix, "standard", (), False)


def test_prefix_whitespace_fails_before_creating_an_environment():
    with pytest.raises(ValueError, match="MAFFT"):
        validate_prefix_path(Path("/target with spaces/profile"))
    validate_prefix_path(Path("/target/profile"))


def test_child_does_not_import_another_worktree_or_inherit_license(monkeypatch):
    for name in ("PYTHONPATH", "PYTHONHOME", "KEY_MODELLER"):
        monkeypatch.setenv(name, "must-not-propagate")
    env = child_environment()
    assert not {"PYTHONPATH", "PYTHONHOME", "KEY_MODELLER"} & env.keys()
    assert env["PYTHONNOUSERSITE"] == "1"


def test_child_starts_without_caller_conda_activation(monkeypatch):
    keys = ("CONDA_PREFIX", "CONDA_PREFIX_1", "CONDA_SHLVL", "CONDA_DEFAULT_ENV", "CONDA_PROMPT_MODIFIER")
    for key in keys:
        monkeypatch.setenv(key, "caller-environment")
    monkeypatch.setenv("CONDA_PKGS_DIRS", "/chosen/package-cache")
    env = child_environment()
    assert not set(keys) & env.keys()
    assert env["CONDA_PKGS_DIRS"] == "/chosen/package-cache"
    # Clearing the child activation must not deactivate the user's shell.
    assert all(os.environ[key] == "caller-environment" for key in keys)


def test_child_enforces_native_transaction_safety_without_changing_parent(monkeypatch):
    monkeypatch.setenv("CONDA_SAFETY_CHECKS", "disabled")
    monkeypatch.setenv("CONDA_PATH_CONFLICT", "clobber")
    env = child_environment()
    assert env["CONDA_SAFETY_CHECKS"] == "enabled"
    assert env["CONDA_PATH_CONFLICT"] == "prevent"
    assert os.environ["CONDA_SAFETY_CHECKS"] == "disabled"
    assert os.environ["CONDA_PATH_CONFLICT"] == "clobber"


def test_native_probe_managed_addon_fails_before_prefix_reservation(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(install, "platform_subdir", lambda: "osx-arm64")
    monkeypatch.setattr(install.shutil, "which", lambda _: pytest.fail("Manager must not be invoked"))
    prefix = tmp_path / "new-standard"
    result = install.main(["--profile", "standard", "--manager", "micromamba", "--prefix", str(prefix),
                           "--addon", "probe", "--non-interactive", "--skip-modeller"])
    assert result == 1
    message = capsys.readouterr().err
    assert "Managed add-on 'probe' is unavailable on osx-arm64" in message
    assert "docs/installation.md#native-apple-probe-alternative" in message
    assert not prefix.exists()
    assert not prefix.with_name(prefix.name + ".dvbfixer-install.json").exists()


@pytest.mark.parametrize("subdir,addon", [("linux-64", "probe"), ("osx-arm64", "msa-alternatives")])
def test_supported_managed_addon_dry_run_is_unchanged(tmp_path, monkeypatch, subdir, addon):
    monkeypatch.setattr(install, "platform_subdir", lambda: subdir)
    monkeypatch.setattr(install.shutil, "which", lambda _: "/tools/micromamba")
    monkeypatch.setattr(Manager, "validate", lambda _: "2.9.0")
    assert install.main(["--profile", "standard", "--manager", "micromamba", "--prefix", str(tmp_path / "new"),
                         "--addon", addon, "--non-interactive", "--skip-modeller", "--dry-run"]) == 0


def test_license_failure_output_is_not_disclosed(monkeypatch, capsys):
    import subprocess

    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 1, "secret", "secret")

    monkeypatch.setattr(subprocess, "run", run)
    with pytest.raises(ValueError, match="details suppressed") as exc:
        configure_modeller(Manager("micromamba", "/tools/micromamba"), Path("/target"), "secret")
    assert "secret" not in str(exc.value)
    assert "secret" not in " ".join(calls[0][0])
    assert calls[0][1]["input"] == "secret"
    assert "secret" not in capsys.readouterr().out


def test_license_configuration_replaces_only_target_and_breaks_hardlinks(tmp_path):
    import os
    import stat
    import subprocess
    import sys

    config = tmp_path / "lib/modeller-test/modlib/modeller/config.py"
    config.parent.mkdir(parents=True)
    config.write_text("# target config\nlicense = 'XXXX'\ninstall_dir = 'unchanged'\n")
    original = tmp_path / "package-cache-config.py"
    os.link(config, original)
    subprocess.run([sys.executable, "-c",
                    f"import sys; sys.prefix = {str(tmp_path)!r}; exec({MODELLER_SETUP!r})"],
                   input="test-value", text=True, check=True)
    assert "license = 'test-value'" in config.read_text()
    assert "install_dir = 'unchanged'" in config.read_text()
    assert "license = 'XXXX'" in original.read_text()
    assert stat.S_IMODE(config.stat().st_mode) == 0o600
