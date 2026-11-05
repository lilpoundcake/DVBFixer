"""Pre-profile installs remain discoverable and explicitly updatable, without adoption."""

import json
import subprocess

import pytest

from dvbfixer import install, install_discovery, install_update
from dvbfixer.install_discovery import discover_installations, package_details


def legacy(tmp_path, *, name="existing", layout="wheel"):
    prefix = tmp_path / name
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta/history").write_text("existing native packages")
    site = prefix / "lib/python3.12/site-packages"
    site.mkdir(parents=True)
    if layout in {"wheel", "editable"}:
        dist = site / "DVBFixer-0.10.2.dist-info"
        dist.mkdir()
        (dist / "METADATA").write_text("Name: DVBFixer\nVersion: 0.10.2\n")
        if layout == "editable":
            (dist / "direct_url.json").write_text(json.dumps({"dir_info": {"editable": True}, "url": "file:///source"}))
    elif layout == "egg-info":
        dist = site / "dvbfixer.egg-info"
        dist.mkdir()
        (dist / "PKG-INFO").write_text("Name: dvbfixer\nVersion: 0.8.6\n")
    elif layout == "egg-link":
        source = tmp_path / "editable-source"
        dist = source / "dvbfixer.egg-info"
        dist.mkdir(parents=True)
        (dist / "PKG-INFO").write_text("Name: dvbfixer\nVersion: 0.8.6\n")
        (site / "dvbfixer.egg-link").write_text(str(source) + "\n.\n")
    elif layout in {"console-only", "editable-marker"}:
        (prefix / "bin").mkdir()
        (prefix / "bin/dvbfixer").write_text("# preserved CLI entry point; metadata missing\n")
        if layout == "editable-marker":
            (site / "__editable__.dvbfixer-0.10.2.pth").write_text("raise RuntimeError('must never execute while discovering')\n")
    return prefix


def fake_manager(monkeypatch, tmp_path, envs):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(install.shutil, "which", lambda _: "/tools/micromamba")
    monkeypatch.setattr(install.Manager, "validate", lambda _: "2.9.0")
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([tmp_path], envs))


@pytest.mark.parametrize("layout,version,editable", [
    ("wheel", "0.10.2", False), ("editable", "0.10.2", True),
    ("egg-info", "0.8.6", False), ("egg-link", "0.8.6", True),
    ("console-only", "(version unknown)", False), ("editable-marker", "(version unknown)", True),
])
def test_old_and_metadata_incomplete_installs_are_detected_without_execution(tmp_path, monkeypatch, layout, version, editable):
    prefix = legacy(tmp_path, layout=layout)
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: pytest.fail("Do not execute discovered interpreters"))
    assert package_details(prefix) == (version, editable, True)


def test_manual_installations_are_selectable_for_update_but_not_resume(tmp_path, monkeypatch):
    prefix = legacy(tmp_path)
    monkeypatch.chdir(tmp_path)
    found = discover_installations([tmp_path], [prefix], repository=tmp_path, include_manual=True)
    assert len(found) == 1 and found[0].prefix == prefix
    assert not found[0].managed
    assert "existing/manual" in found[0].label
    assert discover_installations([tmp_path], [prefix], repository=tmp_path) == []


def test_console_only_dvbfixer_is_shown_in_numbered_update_menu(tmp_path, monkeypatch, capsys):
    prefix = legacy(tmp_path, name="dvbfixer", layout="console-only")
    fake_manager(monkeypatch, tmp_path, [prefix])
    prompts = []
    answers = iter(["micromamba", "1", "1"])

    def answer(prompt):
        prompts.append(prompt)
        return next(answers)

    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", answer)
    before = sorted(str(path) for path in tmp_path.rglob("*"))
    assert install.main(["--update", "--dry-run"]) == 0
    output = capsys.readouterr().out
    assert "1. DVBFixer (version unknown) — dvbfixer" in output
    assert "Package-only update: dvbfixer" in output
    assert "No matching" not in output
    assert all("directory" not in prompt.lower() and "path" not in prompt.lower() for prompt in prompts)
    assert sorted(str(path) for path in tmp_path.rglob("*")) == before


@pytest.mark.parametrize("editable", [False, True])
def test_confirmed_manual_update_changes_only_package_not_dependencies_or_license(tmp_path, monkeypatch, capsys, editable):
    prefix = legacy(tmp_path, layout="editable" if editable else "wheel")
    fake_manager(monkeypatch, tmp_path, [prefix])
    config = prefix / "lib/modeller-10.8/modlib/modeller/config.py"
    config.parent.mkdir(parents=True)
    config.write_text("license = 'protected-key'\n")
    workspace = prefix / "share/dvbfixer/gui/structures/projects/p/manifest.json"
    workspace.parent.mkdir(parents=True)
    workspace.write_text("user workspace")
    answers = iter(["micromamba", "1", "1", "yes"])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    calls = []
    monkeypatch.setattr(install_update, "run_checked", lambda argv, **kwargs: calls.append(argv))
    monkeypatch.setattr(install_update, "logged_run", lambda argv, **kwargs: calls.append(argv) or subprocess.CompletedProcess(argv, 0))
    assert install.main(["--update"]) == 0
    assert len(calls) == 2
    package = calls[0]
    assert "pip" in package and "--no-deps" in package
    assert ("-e" in package) is editable
    assert str(prefix) in package
    assert all(argv[1] == "run" for argv in calls)
    assert all("doctor" in argv or "pip" in argv for argv in calls)
    assert config.read_text() == "license = 'protected-key'\n"
    assert workspace.read_text() == "user workspace"
    assert not prefix.with_name(prefix.name + ".dvbfixer-install.json").exists()
    output = capsys.readouterr().out
    assert "protected-key" not in output
    assert "no Standard/Desktop profile readiness is claimed" in output


def test_manual_update_cancel_does_not_create_record_or_logs(tmp_path, monkeypatch):
    prefix = legacy(tmp_path)
    fake_manager(monkeypatch, tmp_path, [prefix])
    answers = iter(["micromamba", "1", "1", "no"])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    monkeypatch.setattr(install_update, "run_checked", lambda *args, **kwargs: pytest.fail("cancel must not install"))
    before = sorted(str(path) for path in tmp_path.rglob("*"))
    assert install.main(["--update"]) == 0
    assert sorted(str(path) for path in tmp_path.rglob("*")) == before


def test_noninteractive_manual_update_requires_explicit_package_only_consent(tmp_path, monkeypatch, capsys):
    prefix = legacy(tmp_path)
    fake_manager(monkeypatch, tmp_path, [prefix])
    assert install.main(["--update", "--name", prefix.name, "--manager", "micromamba", "--non-interactive"]) == 1
    assert "Add --package-only" in capsys.readouterr().err
    assert install.main(["--update", "--package-only", "--name", prefix.name, "--manager", "micromamba",
                         "--non-interactive", "--dry-run"]) == 0


def test_runtime_child_is_not_misidentified_as_manual_application(tmp_path, monkeypatch):
    root = tmp_path / "frozen"
    install.reserve_prefix(root, "runtime:protpardelle-mps", (), False)
    sampling = legacy(root, name="sampling")
    monkeypatch.chdir(tmp_path)
    entries = discover_installations([tmp_path], [sampling], repository=tmp_path, include_manual=True)
    assert [entry.prefix for entry in entries] == [root]
    assert entries[0].backend


def test_package_only_cannot_be_used_to_mutate_an_unrelated_environment(tmp_path, monkeypatch):
    prefix = legacy(tmp_path)
    for metadata in (prefix / "lib/python3.12/site-packages").glob("*.dist-info/METADATA"):
        metadata.write_text("Name: unrelated\nVersion: 1\n")
    fake_manager(monkeypatch, tmp_path, [prefix])
    monkeypatch.setattr(install_update, "run_checked", lambda *args, **kwargs: pytest.fail("no existing dvbfixer"))
    assert install.main(["--update", "--package-only", "--prefix", str(prefix), "--manager", "micromamba", "--non-interactive"]) == 1


def test_corrupt_editable_metadata_does_not_make_existing_package_disappear(tmp_path, monkeypatch):
    prefix = legacy(tmp_path)
    direct = next(prefix.glob("lib/python*/site-packages/*.dist-info")) / "direct_url.json"
    direct.write_text("not json")
    monkeypatch.chdir(tmp_path)
    entries = discover_installations([tmp_path], [prefix], repository=tmp_path, include_manual=True)
    assert len(entries) == 1 and entries[0].version == "0.10.2"


def test_bad_ownership_record_cannot_be_bypassed_as_a_manual_install(tmp_path, monkeypatch):
    prefix = legacy(tmp_path)
    prefix.with_name(prefix.name + ".dvbfixer-install.json").write_text("invalid record")
    monkeypatch.chdir(tmp_path)
    assert discover_installations([tmp_path], [prefix], repository=tmp_path, include_manual=True) == []


def test_custom_manual_location_also_offers_package_only_without_adopting(tmp_path, monkeypatch):
    prefix = legacy(tmp_path)
    fake_manager(monkeypatch, tmp_path, [prefix])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda prompt: pytest.fail("Explicit manager and dry run need no more input"))
    assert install.main(["--update", "--package-only", "--prefix", str(prefix), "--manager", "micromamba", "--dry-run"]) == 0
    assert not prefix.with_name(prefix.name + ".dvbfixer-install.json").exists()
