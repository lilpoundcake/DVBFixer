"""Explicit profile attachment: planning, refusal, preservation and commit ordering."""

import argparse
import json
import subprocess

import pytest

from dvbfixer import desktop_shortcuts, install, install_discovery, install_extend
from dvbfixer.install import Manager
from dvbfixer.install_extend import extension_composition, prepare_extension
from dvbfixer.install_ui import addons, choose, confirm


def environment(tmp_path, *, package=True):
    prefix = tmp_path / "envs/existing"
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta/history").write_text("original transaction history")
    for name, version in (("python", "3.12.8"), ("numpy", "1.26.4")):
        (prefix / f"conda-meta/{name}.json").write_text(json.dumps({"name": name, "version": version, "build": "b"}))
    if package:
        metadata = prefix / "lib/python3.12/site-packages/dvbfixer-0.10.2.dist-info"
        metadata.mkdir(parents=True)
        (metadata / "METADATA").write_text("Name: dvbfixer\nVersion: 0.10.2\n")
    return prefix


def arguments(prefix, *, gui=False, profile="standard"):
    return argparse.Namespace(prefix=prefix, profile="desktop" if gui else profile, gui_only=gui,
                              addon=[], modeller_key_fd=None)


def harness(tmp_path, monkeypatch, *, native_plan=None, pip_failure=False, doctor=0):
    source = tmp_path / "source"
    module = source / "src/dvbfixer/install_extend.py"
    module.parent.mkdir(parents=True)
    (source / "pyproject.toml").write_text("[project]\nname='dvbfixer'\nversion='0.10.3'\n")
    (source / "gui").mkdir()
    (source / "gui/package.json").write_text("{}")
    (source / "scripts").mkdir()
    (source / "scripts/install_diffusion_backend.py").write_text("# Fixture marker; runtime provisioning is mocked.\n")
    monkeypatch.setattr(install_extend, "__file__", str(module))
    monkeypatch.setattr(install, "__file__", str(module.with_name("install.py")))
    monkeypatch.setattr(install.shutil, "which", lambda _: "/tools/micromamba")
    monkeypatch.setattr(Manager, "validate", lambda _: "2.9.0")
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        payload = ""
        code = 0
        if argv[1:3] == ["info", "--json"]:
            payload = json.dumps({"root_prefix": str(tmp_path / "base")})
        elif "--dry-run" in argv and "--json" in argv:
            payload = json.dumps(native_plan or {"success": True, "actions": {"LINK": [], "UNLINK": []}})
        elif "--report" in argv:
            if pip_failure:
                return subprocess.CompletedProcess(argv, 1, "", "")
            report = install_extend.Path(argv[argv.index("--report") + 1])
            report.write_text(json.dumps({"install": [{"metadata": {"name": "dvbfixer", "version": "0.10.3"}}]}))
        elif install.MODELLER_READ in argv:
            payload = "protected-key"
        elif "doctor" in argv:
            code = doctor
        if kwargs.get("stdout") is not None:
            kwargs["stdout"].write("mock command output\n")
        return subprocess.CompletedProcess(argv, code, payload, "")

    monkeypatch.setattr(install.subprocess, "run", run)
    return calls, source


def test_numbered_choices_confirmations_and_addons(monkeypatch, capsys):
    answers = iter(["9", "2", "", "3,1", "0"])
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    assert choose("Profile", [("standard", "Standard"), ("desktop", "Desktop")]) == "desktop"
    assert confirm("Apply?") is False
    assert addons([("a", "A"), ("b", "B"), ("c", "C")]) == ["c", "a"]
    with pytest.raises(EOFError):
        choose("Cancel", [("a", "A")])
    assert "1. Standard" in capsys.readouterr().out


def test_preview_and_cancel_leave_environment_and_ownership_unchanged(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    calls, source = harness(tmp_path, monkeypatch)
    before = {str(path.relative_to(prefix)): path.read_bytes() for path in prefix.rglob("*") if path.is_file()}
    plan = prepare_extension(arguments(prefix), Manager("micromamba", "/tools/micromamba"))
    try:
        assert any("--override-channels" in argv and "--dry-run" in argv for argv, _ in calls)
        assert any("--report" in argv and "pyscf" in argv for argv, _ in calls)
        assert "conda-forge::mdanalysis" in plan.native_command
        assert "conda-forge::python=3.12.*" in plan.native_command
        assert not plan.marker.exists()
        assert before == {str(path.relative_to(prefix)): path.read_bytes() for path in prefix.rglob("*") if path.is_file()}
    finally:
        plan.close()


@pytest.mark.parametrize("name", ["micromamba", "mamba", "conda"])
def test_extension_native_plan_uses_channel_isolated_install_for_each_manager(tmp_path, monkeypatch, name):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    plan = prepare_extension(arguments(prefix), Manager(name, f"/tools/{name}"))
    try:
        assert plan.native_command[:2] == [f"/tools/{name}", "install"]
        assert "--override-channels" in plan.native_command
        assert "--clobber" not in plan.native_command
    finally:
        plan.close()


@pytest.mark.parametrize("failure", ["native", "pip"])
def test_conflicts_fail_before_any_package_or_ownership_change(tmp_path, monkeypatch, failure):
    prefix = environment(tmp_path)
    calls, _ = harness(tmp_path, monkeypatch,
        native_plan={"success": False, "message": "incompatible Python"} if failure == "native" else None,
        pip_failure=failure == "pip")
    with pytest.raises(ValueError, match="plan failed|preview failed"):
        prepare_extension(arguments(prefix), Manager("micromamba", "/tools/micromamba"))
    assert not prefix.with_name(prefix.name + ".dvbfixer-install.json").exists()
    assert all("--dry-run" in argv or argv[1:3] == ["info", "--json"] for argv, _ in calls)


def test_environment_drift_is_rejected_before_apply(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    plan = prepare_extension(arguments(prefix), Manager("micromamba", "/tools/micromamba"))
    try:
        (prefix / "conda-meta/history").write_text("another process changed this")
        with pytest.raises(ValueError, match="changed after preview"):
            plan.verify_before_apply()
        assert not plan.marker.exists()
    finally:
        plan.close()


def test_gui_only_does_not_request_scientific_dependencies_or_change_ownership(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    plan = prepare_extension(arguments(prefix, gui=True), Manager("micromamba", "/tools/micromamba"))
    try:
        assert any("nodejs" in arg for arg in plan.native_command)
        assert "python=3.12.8=b" in plan.native_command
        assert "numpy=1.26.4=b" in plan.native_command
        assert all(not value.startswith("conda-forge::python") for value in plan.native_command)
        assert not any("openmm" in arg or "pyscf" in arg or "ambertools" in arg for arg in plan.native_command)
        assert "--no-deps" in plan.pip_command
        plan.commit()
        assert not plan.marker.exists()
    finally:
        plan.close()


def test_gui_only_refuses_solver_changes_to_existing_science(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch, native_plan={"success": True, "actions": {
        "UNLINK": [{"name": "numpy", "version": "1.26.4", "build": "b"}],
        "LINK": [{"name": "numpy", "version": "2.4.0", "build": "b"}],
    }})
    with pytest.raises(ValueError, match="GUI-only plan would change"):
        prepare_extension(arguments(prefix, gui=True), Manager("micromamba", "/tools/micromamba"))


def test_profile_commit_retains_existing_desktop_and_addons(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    marker = prefix.with_name(prefix.name + ".dvbfixer-install.json")
    marker.write_text(json.dumps({"schema_version": 1, "prefix": str(prefix), "profile": "desktop", "addons": ["msa-alternatives"]}))
    args = arguments(prefix)
    profile, selected = extension_composition(args)
    assert profile == "desktop" and selected == ["msa-alternatives"]
    args.profile, args.addon = profile, selected
    harness(tmp_path, monkeypatch)
    plan = prepare_extension(args, Manager("micromamba", "/tools/micromamba"))
    try:
        assert json.loads(marker.read_text())["profile"] == "desktop"
        plan.commit()
        assert json.loads(marker.read_text())["addons"] == ["msa-alternatives"]
    finally:
        plan.close()


@pytest.mark.parametrize("kind", ["base", "runtime", "openbabel", "torch"])
def test_unsafe_profile_targets_are_refused(tmp_path, monkeypatch, kind):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    if kind == "runtime":
        (prefix / "runtime-plan.json").write_text("{}")
        with pytest.raises(ValueError, match="Frozen backend"):
            extension_composition(arguments(prefix))
        return
    if kind == "openbabel":
        (prefix / "conda-meta/openbabel.json").write_text(json.dumps({"name": "openbabel", "version": "3.1.1", "build": "b"}))
    elif kind == "torch":
        (prefix / "conda-meta/torch.json").write_text(json.dumps({"name": "torch", "version": "2.6.0", "build": "b"}))
    else:
        original = install.subprocess.run

        def run(argv, **kwargs):
            if argv[1:3] == ["info", "--json"]:
                return subprocess.CompletedProcess(argv, 0, json.dumps({"root_prefix": str(prefix)}), "")
            return original(argv, **kwargs)

        monkeypatch.setattr(install.subprocess, "run", run)
    with pytest.raises(ValueError, match="manager base|Open Babel|Keep PyTorch"):
        prepare_extension(arguments(prefix), Manager("micromamba", "/tools/micromamba"))


@pytest.mark.parametrize("gui,doctor", [(False, 0), (False, 2), (False, 1), (True, 0)])
def test_complete_extension_preserves_files_and_commits_only_after_checks(tmp_path, monkeypatch, gui, doctor, capsys):
    prefix = environment(tmp_path)
    calls, source = harness(tmp_path, monkeypatch, doctor=doctor)
    config = prefix / "lib/modeller-10.8/modlib/modeller/config.py"
    config.parent.mkdir(parents=True)
    config.write_text("license = 'protected-key'\n")
    workspace = prefix / "share/dvbfixer/gui/structures/projects/p/manifest.json"
    workspace.parent.mkdir(parents=True)
    workspace.write_text("existing project")
    (workspace.parents[3] / ".env").write_text("local settings")
    applied = []
    monkeypatch.setattr(install, "run_checked", lambda argv, **kwargs: applied.append(argv))
    flags = ["--gui-only"] if gui else ["--profile", "desktop"]
    assert install.main(["--extend", *flags, "--prefix", str(prefix), "--manager", "micromamba", "--non-interactive"]) == doctor
    marker = prefix.with_name(prefix.name + ".dvbfixer-install.json")
    assert marker.exists() is (not gui and doctor in (0, 2))
    assert workspace.read_text() == "existing project"
    assert config.read_text() == "license = 'protected-key'\n"
    assert any("npm" in argv and "build" in argv for argv in applied)
    assert all("--clobber" not in argv for argv in applied)
    if gui:
        assert all("smoke_installation_profile.py" not in " ".join(argv) for argv in applied)
        assert "--no-deps" in next(argv for argv in applied if "pip" in argv)
        assert all(install.MODELLER_READ not in argv for argv, _ in calls)
        assert all("import modeller" not in argv for argv, _ in calls)
    output = capsys.readouterr().out
    assert "protected-key" not in output
    if gui:
        assert "Extend DVBFixer GUI only" in output
    assert (workspace.parents[3] / ".env").read_text() == "local settings"


def test_number_only_interactive_extension_can_select_empty_environment_and_cancel(tmp_path, monkeypatch, capsys):
    prefix = environment(tmp_path, package=False)
    harness(tmp_path, monkeypatch)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([prefix.parent], [prefix]))
    answers = iter(["3", "1", "1", "2", "2", "0", "2"])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    monkeypatch.setattr(install, "run_checked", lambda *args, **kwargs: pytest.fail("cancel must not install"))
    assert install.main([]) == 0
    output = capsys.readouterr().out
    assert "DVBFixer not installed" in output
    assert "2. Desktop" in output
    assert "Cancelled" in output
    assert not prefix.with_name(prefix.name + ".dvbfixer-install.json").exists()


def test_build_failure_does_not_replace_prior_profile_record(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    marker = prefix.with_name(prefix.name + ".dvbfixer-install.json")
    original = json.dumps({"schema_version": 1, "prefix": str(prefix), "profile": "standard", "addons": []})
    marker.write_text(original)
    harness(tmp_path, monkeypatch)

    def run(argv, **kwargs):
        if "npm" in argv and "build" in argv:
            raise subprocess.CalledProcessError(1, argv)

    monkeypatch.setattr(install, "run_checked", run)
    assert install.main(["--extend", "--profile", "desktop", "--prefix", str(prefix),
                         "--manager", "micromamba", "--non-interactive", "--skip-modeller"]) == 1
    assert marker.read_text() == original
    assert not list(prefix.parent.glob("*.installation-logs/*/extension-completed.json"))


def test_new_unconfirmed_pip_plan_stops_after_native_before_package_publication(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    plan = prepare_extension(arguments(prefix), Manager("micromamba", "/tools/micromamba"))
    try:
        monkeypatch.setattr(plan, "preview_pip", lambda _: {("unexpected-package", "100")})
        with pytest.raises(ValueError, match="Pip plan changed"):
            plan.verify_pip_after_native()
        assert not plan.marker.exists()
    finally:
        plan.close()


def test_missing_or_unexpected_native_output_is_not_accepted(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    plan = prepare_extension(arguments(prefix), Manager("micromamba", "/tools/micromamba"))
    try:
        (prefix / "conda-meta/numpy.json").unlink()
        with pytest.raises(ValueError, match="Applied native packages differ"):
            plan.verify_pip_after_native()
        assert not plan.marker.exists()
    finally:
        plan.close()


def test_profile_record_is_not_overwritten_if_changed_during_install(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    plan = prepare_extension(arguments(prefix), Manager("micromamba", "/tools/micromamba"))
    try:
        plan.marker.write_text("unrelated user record")
        with pytest.raises(ValueError, match="Ownership record changed"):
            plan.commit()
        assert plan.marker.read_text() == "unrelated user record"
    finally:
        plan.close()


def test_new_native_transitive_python_package_is_pinned_before_pip_apply(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch, native_plan={"success": True, "actions": {
        "UNLINK": [], "LINK": [{"name": "h5py", "version": "3.15.0", "build": "b"}],
    }})
    plan = prepare_extension(arguments(prefix), Manager("micromamba", "/tools/micromamba"))
    try:
        (prefix / "conda-meta/h5py.json").write_text(json.dumps({"name": "h5py", "version": "3.15.0", "build": "b"}))
        metadata = prefix / "lib/python3.12/site-packages/h5py-3.15.0.dist-info"
        metadata.mkdir()
        (metadata / "METADATA").write_text("Name: h5py\nVersion: 3.15.0\n")

        def preview(label):
            assert "h5py==3.15.0" in plan.constraints.read_text()
            return plan.pip_plan

        monkeypatch.setattr(plan, "preview_pip", preview)
        plan.verify_pip_after_native()
    finally:
        plan.close()


@pytest.mark.parametrize("apply", [False, True])
def test_openbabel_profile_conflict_offers_explicit_gui_only_recovery(tmp_path, monkeypatch, capsys, apply):
    prefix = environment(tmp_path)
    openbabel = prefix / "conda-meta/openbabel.json"
    original = json.dumps({"name": "openbabel", "version": "3.1.1", "build": "b"})
    openbabel.write_text(original)
    calls, _ = harness(tmp_path, monkeypatch)
    answers = iter(["1" if apply else "2", "0", "1", "1" if apply else "2"])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    applied = []
    monkeypatch.setattr(install, "run_checked", lambda argv, **kwargs: applied.append(argv))
    shortcuts = []
    monkeypatch.setattr(desktop_shortcuts, "install_desktop_shortcut", lambda selected: shortcuts.append(selected) or (tmp_path / "DVBFixer.app"))
    before = {str(path.relative_to(prefix)): path.read_bytes() for path in prefix.rglob("*") if path.is_file()}
    assert install.main(["--extend", "--profile", "desktop", "--prefix", str(prefix), "--manager", "micromamba"]) == 0
    output = capsys.readouterr().out
    assert "No packages or profile ownership were changed" in output
    assert "1. Add GUI only" in output
    assert "Scope changed explicitly to GUI only" in output
    assert openbabel.read_text() == original
    assert not prefix.with_name(prefix.name + ".dvbfixer-install.json").exists()
    if apply:
        pip = next(argv for argv in applied if "pip" in argv)
        assert "--no-deps" in pip
        assert all("openbabel-wheel" not in arg and "pyscf" not in arg for argv in applied for arg in argv)
        assert all(install.MODELLER_READ not in argv and "import modeller" not in argv for argv, _ in calls)
        assert "GUI only (scientific stack unchanged)" in output
        assert shortcuts == [prefix]
    else:
        assert not applied
        assert not shortcuts
        assert before == {str(path.relative_to(prefix)): path.read_bytes() for path in prefix.rglob("*") if path.is_file()}


def test_openbabel_conflict_can_stop_without_changes(tmp_path, monkeypatch, capsys):
    prefix = environment(tmp_path)
    (prefix / "conda-meta/openbabel.json").write_text(json.dumps({"name": "openbabel", "version": "3.1.1", "build": "b"}))
    harness(tmp_path, monkeypatch)
    answers = iter(["2", "0", "2"])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    monkeypatch.setattr(install, "run_checked", lambda *args, **kwargs: pytest.fail("cancel must not apply"))
    assert install.main(["--extend", "--profile", "desktop", "--prefix", str(prefix), "--manager", "micromamba"]) == 0
    assert "Cancelled" in capsys.readouterr().out
    assert not prefix.with_name(prefix.name + ".dvbfixer-install.json").exists()


def test_noninteractive_conflict_never_changes_scope_automatically(tmp_path, monkeypatch, capsys):
    prefix = environment(tmp_path)
    (prefix / "conda-meta/openbabel.json").write_text(json.dumps({"name": "openbabel", "version": "3.1.1", "build": "b"}))
    calls, _ = harness(tmp_path, monkeypatch)
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("non-interactive must not prompt"))
    assert install.main(["--extend", "--profile", "desktop", "--prefix", str(prefix), "--manager", "micromamba", "--non-interactive"]) == 1
    assert "--gui-only" in capsys.readouterr().err
    assert all(argv[1:3] == ["info", "--json"] for argv, _ in calls)


def test_extension_menu_hides_bases_bundles_and_unrelated_environments(tmp_path, monkeypatch, capsys):
    prefix = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    base, foreign, bundle = tmp_path / "base", tmp_path / "miniforge3", tmp_path / "PyMOL.app/Contents"
    other, foreign_child = prefix.parent / "unrelated-science", foreign / "envs/other"
    for candidate in (base, foreign, bundle, other, foreign_child):
        (candidate / "conda-meta").mkdir(parents=True)
        (candidate / "conda-meta/history").write_text("fixture history")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([prefix.parent], [prefix, base, foreign, bundle, other, foreign_child]))
    answers = iter(["1", "1", "3", "2"])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    assert install.main(["--extend", "--dry-run"]) == 0
    output = capsys.readouterr().out
    assert "2. Other existing environments" in output
    assert "PyMOL.app" not in output and "miniforge3" not in output
    assert "unrelated-science" not in output


def test_gui_only_preserves_known_editable_package_mode(tmp_path, monkeypatch):
    prefix = environment(tmp_path)
    dist = next(prefix.glob("lib/python*/site-packages/dvbfixer-*.dist-info"))
    (dist / "direct_url.json").write_text(json.dumps({"dir_info": {"editable": True}, "url": "file:///old-source"}))
    harness(tmp_path, monkeypatch)
    plan = prepare_extension(arguments(prefix, gui=True), Manager("micromamba", "/tools/micromamba"))
    try:
        assert "-e" in plan.pip_command
        assert "--no-deps" in plan.pip_command
    finally:
        plan.close()


def test_second_update_scope_option_previews_complete_profile_not_cancellation(tmp_path, monkeypatch, capsys):
    prefix = environment(tmp_path)
    calls, _ = harness(tmp_path, monkeypatch)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([prefix.parent], [prefix]))
    answers = iter(["1", "1", "2", "2", "2"])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    monkeypatch.setattr(install, "run_checked", lambda *args, **kwargs: pytest.fail("dry run must not install"))
    assert install.main(["--update", "--dry-run"]) == 0
    output = capsys.readouterr().out
    assert "2. Full installation" in output
    assert "Update profile DVBFixer desktop" in output
    assert "Package-only update:" not in output
    assert any("pyscf" in argv and "--report" in argv for argv, _ in calls)
    assert not prefix.with_name(prefix.name + ".dvbfixer-install.json").exists()


def test_diffusion_is_available_from_selected_application_as_separate_runtime(tmp_path, monkeypatch, capsys):
    from dvbfixer import runtime_install

    application = environment(tmp_path)
    harness(tmp_path, monkeypatch)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(install, "platform_subdir", lambda: "osx-arm64")
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([application.parent], [application]))
    # Manager, scientific installation, diffusion component, MPS model, new
    # runtime, manager default, runtime name, deferred checkpoint.
    answers = iter(["1", "1", "4", "1", "1", "", "new-mps-runtime", ""])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    calls = []

    def provision(**kwargs):
        calls.append(kwargs)
        return {"runtime": kwargs["name"], "prefix": str(kwargs["prefix"]), "state": "planned"}

    monkeypatch.setattr(runtime_install, "provision_runtime", provision)
    before = {str(path.relative_to(application)): path.read_bytes() for path in application.rglob("*") if path.is_file()}
    assert install.main(["--extend", "--dry-run"]) == 0
    assert len(calls) == 1
    assert calls[0]["name"] == "protpardelle-mps"
    assert calls[0]["prefix"] == application.parent / "new-mps-runtime"
    assert calls[0]["prefix"] != application
    assert calls[0]["dry_run"] is True
    assert calls[0]["resume"] is False
    assert before == {str(path.relative_to(application)): path.read_bytes() for path in application.rglob("*") if path.is_file()}
    output = capsys.readouterr().out
    assert "4. Diffusion backend" in output
    assert "per-user discovery" in output
    assert "not PyTorch inside Standard/Desktop" in output


def test_diffusion_existing_runtime_can_be_resumed_without_reusing_app_prefix(tmp_path, monkeypatch):
    from dvbfixer import runtime_install

    application = environment(tmp_path)
    runtime = application.parent / "existing-mps"
    install.reserve_prefix(runtime, "runtime:protpardelle-mps", (), False)
    (runtime / "sampling").mkdir(parents=True)
    harness(tmp_path, monkeypatch)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(install, "platform_subdir", lambda: "osx-arm64")
    monkeypatch.setattr(install_discovery, "manager_inventory", lambda _: ([application.parent], [application, runtime / "sampling"]))
    answers = iter(["1", "1", "4", "1", "2", "", "1", ""])
    monkeypatch.setattr(install.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    calls = []

    def provision(**kwargs):
        calls.append(kwargs)
        return {"runtime": kwargs["name"], "prefix": str(kwargs["prefix"]), "state": "planned"}

    monkeypatch.setattr(runtime_install, "provision_runtime", provision)
    assert install.main(["--extend", "--dry-run"]) == 0
    assert calls[0]["prefix"] == runtime
    assert calls[0]["resume"] is True
