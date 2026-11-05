"""Name-based provisioning and read-only discovery of helper-owned installations."""

import json
import subprocess
from pathlib import Path

import pytest

from dvbfixer import install, install_discovery
from dvbfixer.install import Manager
from dvbfixer.install_discovery import (
    Installation,
    creation_directory,
    discover_installations,
    installed_version,
    manager_inventory,
    new_named_prefix,
    protected_environment,
    select_existing,
    validate_name,
)


def record(prefix, profile="standard", version="0.10.3"):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    install.reserve_prefix(prefix, profile, (), False)
    if not profile.startswith("runtime:"):
        metadata = prefix / f"lib/python3.12/site-packages/dvbfixer-{version}.dist-info/METADATA"
        metadata.parent.mkdir(parents=True)
        metadata.write_text(f"Name: dvbfixer\nVersion: {version}\n")
    return prefix


@pytest.mark.parametrize("name,info", [
    ("micromamba", {"envs directories": ["/manager/envs"]}),
    ("mamba", {"envs_dirs": ["/manager/envs", "/another/envs"]}),
    ("conda", {"envs_dirs": ["/manager/envs"]}),
    ("micromamba", {"base environment": "/manager"}),
    ("conda", {"root_prefix": "/manager"}),
])
def test_inventory_supports_native_manager_json_without_mutation(monkeypatch, name, info):
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        data = info if argv[1] == "info" else {"envs": ["/custom/desktop", "/custom/backend/sampling"]}
        return subprocess.CompletedProcess(argv, 0, json.dumps(data), "")

    monkeypatch.setattr(subprocess, "run", run)
    directories, envs = manager_inventory(Manager(name, f"/bin/{name}"))
    assert directories[0] == Path("/manager/envs")
    assert envs == [Path("/custom/desktop"), Path("/custom/backend/sampling")]
    assert [argv[1:] for argv, _ in calls] == [["info", "--json"], ["env", "list", "--json"]]
    assert all(kwargs["capture_output"] and kwargs["timeout"] == 30 for _, kwargs in calls)


@pytest.mark.parametrize("payload", [{}, {"envs_dirs": ["relative"]}, {"envs_dirs": "wrong"}])
def test_invalid_manager_layout_fails_without_guessing_a_directory(monkeypatch, payload):
    monkeypatch.setattr(subprocess, "run", lambda argv, **kwargs: subprocess.CompletedProcess(
        argv, 0, json.dumps(payload if argv[1] == "info" else {"envs": []}), ""))
    with pytest.raises(ValueError, match="directories"):
        manager_inventory(Manager("conda", "conda"))


def test_discovery_finds_versions_custom_prefixes_composite_and_interrupted_roots(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    managed = tmp_path / "manager/envs"
    standard = record(managed / "standard-v1", version="0.10.2")
    desktop = record(tmp_path / "custom/desktop", profile="desktop")
    backend = record(tmp_path / "custom/backend", profile="runtime:protpardelle-mps")
    interrupted = record(managed / "interrupted", profile="runtime:psi4")
    unowned = tmp_path / "manual-dvbfixer"
    unowned.mkdir()
    before = sorted(str(path) for path in tmp_path.rglob("*"))
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: pytest.fail("Do not execute discovered environments"))
    found = discover_installations([managed], [standard, desktop, backend / "sampling", backend / "refinement", unowned],
                                   repository=tmp_path)
    assert {entry.prefix for entry in found} == {standard, desktop, backend, interrupted}
    assert {entry.version for entry in found if not entry.backend} == {"0.10.2", "0.10.3"}
    assert len([entry for entry in found if entry.prefix == backend]) == 1
    assert sorted(str(path) for path in tmp_path.rglob("*")) == before


def test_discovery_rejects_corrupt_and_redirected_records(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    valid = record(tmp_path / "envs/valid")
    corrupt = tmp_path / "envs/corrupt.dvbfixer-install.json"
    corrupt.write_text('{"prefix": "not this installation"}')
    alias = tmp_path / "envs/alias"
    alias.symlink_to(valid, target_is_directory=True)
    alias.with_name(alias.name + ".dvbfixer-install.json").write_text(
        valid.with_name(valid.name + ".dvbfixer-install.json").read_text())
    found = discover_installations([valid.parent], [], repository=tmp_path)
    assert [entry.prefix for entry in found] == [valid]
    assert "Ignoring unusable installation record" in capsys.readouterr().out


@pytest.mark.parametrize("name", ["../existing", "/tmp/target", "nested/name", "name with space", "", "base", "root", "x\\y", ".."])
def test_environment_names_are_not_paths(name):
    with pytest.raises(ValueError):
        validate_name(name)


def test_new_name_suggestion_avoids_existing_names_and_crash_left_reservations(tmp_path, monkeypatch):
    directory = tmp_path / "manager/envs"
    directory.mkdir(parents=True)
    (directory / "dvbfixer-standard-0.10.3").mkdir()
    (directory / "dvbfixer-standard-0.10.3-2.dvbfixer-install.json").write_text("reservation")
    prompts = []
    monkeypatch.setattr("builtins.input", lambda prompt: prompts.append(prompt) or "")
    prefix = new_named_prefix(directory, [Path("/other/dvbfixer-standard-0.10.3-3")],
                              default="dvbfixer-standard-0.10.3", name=None, interactive=True)
    assert prefix == directory / "dvbfixer-standard-0.10.3-4"
    assert not prefix.exists()
    assert all("Environment name" in prompt for prompt in prompts)


def test_explicit_occupied_name_fails_instead_of_reusing_prefix(tmp_path):
    (tmp_path / "existing").mkdir()
    with pytest.raises(ValueError, match="already in use"):
        new_named_prefix(tmp_path, [], default="default", name="existing", interactive=False)


def test_number_selection_distinguishes_duplicate_names_without_path_input(tmp_path, monkeypatch, capsys):
    entries = [Installation(tmp_path / "one/desktop", "desktop", "0.10.2"),
               Installation(tmp_path / "two/desktop", "desktop", "0.10.3")]
    answers = iter(["99", "/not/a/selection", "2"])
    prompts = []

    def answer(prompt):
        prompts.append(prompt)
        return next(answers)

    monkeypatch.setattr("builtins.input", answer)
    assert select_existing(entries) == entries[1]
    assert prompts == ["Installation number: "] * 3
    output = capsys.readouterr().out
    assert "0.10.2" in output and "0.10.3" in output
    assert str(entries[0].prefix) in output and str(entries[1].prefix) in output


def test_unknown_version_is_not_fabricated(tmp_path):
    assert installed_version(tmp_path) == "(version unknown)"


def test_noninteractive_named_install_never_prompts(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(install.shutil, "which", lambda _: "/tools/micromamba")
    monkeypatch.setattr(Manager, "validate", lambda _: "2.9.0")
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([tmp_path / "envs"], []))
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("Non-interactive installation must not prompt"))
    assert install.main(["--profile", "desktop", "--manager", "micromamba", "--name", "desktop-v1",
                         "--non-interactive", "--dry-run"]) == 0
    assert str(tmp_path / "envs/desktop-v1") in capsys.readouterr().out
    assert not (tmp_path / "envs").exists()


def test_explicit_name_and_prefix_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        install.main(["--name", "a", "--prefix", "/b", "--non-interactive"])


def test_new_location_follows_first_writable_manager_directory(tmp_path, monkeypatch):
    first, second = tmp_path / "readonly/envs", tmp_path / "writable/envs"
    first.parent.mkdir()
    second.parent.mkdir()
    monkeypatch.setattr(install_discovery.os, "access", lambda path, mode: path == second.parent)
    assert creation_directory([first, second]) == second
    assert not second.exists()


def test_no_writable_manager_location_does_not_invent_a_fallback(tmp_path, monkeypatch):
    monkeypatch.setattr(install_discovery.os, "access", lambda *args: False)
    with pytest.raises(ValueError, match="no writable"):
        creation_directory([tmp_path / "envs"])


def test_update_chooses_one_numbered_version_and_never_lists_backends(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    directory = tmp_path / "envs"
    first = record(directory / "alpha", version="0.10.1")
    second = record(directory / "zeta", version="0.10.2")
    record(directory / "model-backend", profile="runtime:protpardelle-mps")
    for prefix in (first, second):
        (prefix / "conda-meta").mkdir()
        (prefix / "conda-meta/history").write_text("existing")
    monkeypatch.setattr(install.shutil, "which", lambda _: "/tools/micromamba")
    monkeypatch.setattr(Manager, "validate", lambda _: "2.9.0")
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([directory], []))
    prompts, calls = [], []
    answers = iter(["micromamba", "2", "yes"])

    def answer(prompt):
        prompts.append(prompt)
        return next(answers)

    def run(argv, **kwargs):
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 0, "0.10.2", "")

    monkeypatch.setattr("builtins.input", answer)
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr(install.subprocess, "run", run)
    before = sorted(str(path) for path in tmp_path.rglob("*"))
    assert install.main(["--update", "--dry-run"]) == 0
    output = capsys.readouterr().out
    assert "1. DVBFixer Standard 0.10.1 — alpha" in output
    assert "2. DVBFixer Standard 0.10.2 — zeta" in output
    assert "Additional backend" not in output
    assert all("directory" not in prompt.lower() and "path" not in prompt.lower() for prompt in prompts)
    assert all(str(second) in argv for argv in calls)
    assert sorted(str(path) for path in tmp_path.rglob("*")) == before


def test_ambiguous_name_requires_number_selection_instead_of_guessing(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    first = record(tmp_path / "one/desktop", profile="desktop")
    second = record(tmp_path / "two/desktop", profile="desktop")
    monkeypatch.setattr(install.shutil, "which", lambda _: "/tools/micromamba")
    monkeypatch.setattr(Manager, "validate", lambda _: "2.9.0")
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([first.parent, second.parent], []))
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: pytest.fail("Ambiguous selection must not execute a target"))
    assert install.main(["--update", "--name", "desktop", "--manager", "micromamba", "--non-interactive"]) == 1
    assert "ambiguous" in capsys.readouterr().err


def test_no_installations_does_not_fall_back_to_path_prompt(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(install.shutil, "which", lambda _: "/tools/micromamba")
    monkeypatch.setattr(Manager, "validate", lambda _: "2.9.0")
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([tmp_path / "envs"], []))
    monkeypatch.setattr(install_discovery, "discover_installations", lambda *args, **kwargs: [])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    prompts = []

    def answer(prompt):
        prompts.append(prompt)
        assert prompt == "Choice: "
        return "micromamba"

    monkeypatch.setattr("builtins.input", answer)
    assert install.main(["--update"]) == 1
    assert "No matching DVBFixer installations" in capsys.readouterr().err
    assert len(prompts) == 1
    assert not (tmp_path / "envs").exists()


def test_extension_protects_foreign_manager_bases_and_application_bundles(tmp_path):
    base = tmp_path / "micromamba"
    foreign = tmp_path / "miniforge3"
    envs = [foreign, foreign / "envs/science", base / "envs/dvbfixer"]
    assert protected_environment(base, base, envs)
    assert protected_environment(foreign, base, envs)
    assert protected_environment(tmp_path / "PyMOL.app/Contents", base, envs)
    assert protected_environment(tmp_path / "PyMOL.app/Contents/envs/python", base, envs)
    assert not protected_environment(base / "envs/dvbfixer", base, envs)
    assert not protected_environment(base / "envs/science.app", base, envs)


def test_other_environments_are_collapsed_by_default(tmp_path, monkeypatch, capsys):
    primary = Installation(tmp_path / "dvbfixer", "manual", "0.10.2", managed=False)
    other = Installation(tmp_path / "unrelated-science", "manual", "(version unknown)", managed=False, package_present=False)
    monkeypatch.setattr("builtins.input", lambda _: "1")
    assert select_existing([primary], other_entries=[other]) == primary
    output = capsys.readouterr().out
    assert "2. Other existing environments" in output
    assert "unrelated-science" not in output


def test_other_environments_can_be_selected_explicitly_by_number(tmp_path, monkeypatch, capsys):
    primary = Installation(tmp_path / "dvbfixer", "manual", "0.10.2", managed=False)
    other = Installation(tmp_path / "science", "manual", "(version unknown)", managed=False, package_present=False)
    answers = iter(["2", "1"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    assert select_existing([primary], other_entries=[other]) == other
    assert "DVBFixer not installed" in capsys.readouterr().out
