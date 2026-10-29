"""Runtime isolation, operator assets, source integrity and no-clobber publication."""

from __future__ import annotations

import hashlib
import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from dvbfixer import runtime_install
from dvbfixer.install import Manager
from dvbfixer.runtime_install import (
    copy_checked_asset,
    prepare_source,
    provision_runtime,
    publish_launcher,
    runtime_environment,
    verify_assets,
)

ROOT = Path(__file__).resolve().parents[1]


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


@pytest.fixture
def engine_inputs(tmp_path):
    source = tmp_path / "operator engine"
    source.mkdir()
    env = {**os.environ, "GIT_AUTHOR_NAME": "Installer test", "GIT_COMMITTER_NAME": "Installer test",
           "GIT_AUTHOR_EMAIL": "installer@example.invalid", "GIT_COMMITTER_EMAIL": "installer@example.invalid"}
    subprocess.run(["git", "init", str(source)], check=True, capture_output=True, env=env)
    (source / "models.py").write_text("original\n")
    (source / "other.py").write_text("untouched\n")
    subprocess.run(["git", "-C", str(source), "add", "models.py", "other.py"], check=True, env=env)
    subprocess.run(["git", "-C", str(source), "commit", "-m", "Fixture engine"],
                   check=True, capture_output=True, env=env)
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    repository = tmp_path / "application"
    repository.mkdir()
    (repository / "callback.patch").write_text(
        "diff --git a/models.py b/models.py\n--- a/models.py\n+++ b/models.py\n"
        "@@ -1 +1 @@\n-original\n+patched\n")
    entry = {"repository": str(source), "revision": revision, "patch": "callback.patch"}
    requirements = {"sources": {"models.py": digest(b"patched\n")}}
    return source, repository, entry, requirements


def test_engine_checkout_and_resume_preserve_operator_sources(tmp_path, engine_inputs):
    source, repository, entry, requirements = engine_inputs
    engine = tmp_path / "new engine"
    log = tmp_path / "source.log"
    prepare_source(entry, repository, engine, source, log, os.environ.copy(), requirements)
    assert (source / "models.py").read_text() == "original\n"
    assert (engine / "models.py").read_text() == "patched\n"
    prepare_source(entry, repository, engine, source, log, os.environ.copy(), requirements)
    assert (engine / "models.py").read_text() == "patched\n"


@pytest.mark.parametrize("tamper", ["other.py", "models.py", "untracked.py"])
def test_resume_refuses_foreign_source_changes_without_cleaning_them(tmp_path, engine_inputs, tamper):
    source, repository, entry, requirements = engine_inputs
    engine = tmp_path / "engine"
    log = tmp_path / "source.log"
    prepare_source(entry, repository, engine, source, log, os.environ.copy(), requirements)
    target = engine / tamper
    target.write_text(target.read_text() + "foreign\n" if target.exists() else "foreign\n")
    with pytest.raises((ValueError, subprocess.CalledProcessError)):
        prepare_source(entry, repository, engine, source, log, os.environ.copy(), requirements)
    assert "foreign" in target.read_text()
    assert (source / "models.py").read_text() == "original\n"


def test_missing_checkpoint_is_deferred_but_bad_or_symlinked_checkpoint_is_error(tmp_path):
    requirements = {"checkpoint": {"name": "weights.pth", "sha256": digest(b"weights")}}
    assert verify_assets(requirements, checkpoint=None, kalign=None, auxiliary=None)
    checkpoint = tmp_path / "weights.pth"
    checkpoint.write_bytes(b"wrong")
    with pytest.raises(ValueError, match="SHA256"):
        verify_assets(requirements, checkpoint=checkpoint, kalign=None, auxiliary=None)
    checkpoint.write_bytes(b"weights")
    assert verify_assets(requirements, checkpoint=checkpoint, kalign=None, auxiliary=None) == []
    link = tmp_path / "linked.pth"
    link.symlink_to(checkpoint)
    with pytest.raises(ValueError, match="regular file"):
        verify_assets(requirements, checkpoint=link, kalign=None, auxiliary=None)


def test_launcher_publication_is_no_clobber_and_explicit_replace_is_atomic(tmp_path):
    candidate = tmp_path / "launcher"
    candidate.write_text("#!/bin/sh\nnew\n")
    directory = tmp_path / "runners"
    directory.mkdir()
    existing = directory / candidate.name
    existing.write_text("existing\n")
    with pytest.raises(FileExistsError):
        publish_launcher(candidate, directory)
    assert existing.read_text() == "existing\n"
    assert publish_launcher(candidate, directory, replace=True) == existing
    assert existing.read_text() == candidate.read_text()
    assert stat.S_IMODE(existing.stat().st_mode) == 0o700
    assert sorted(path.name for path in directory.iterdir()) == [candidate.name]


def test_launcher_symlink_is_not_replaced_even_when_replace_is_selected(tmp_path):
    candidate = tmp_path / "launcher"
    candidate.write_text("new")
    directory = tmp_path / "runners"
    directory.mkdir()
    (directory / candidate.name).symlink_to(candidate)
    with pytest.raises(ValueError, match="symlink"):
        publish_launcher(candidate, directory, replace=True)
    assert candidate.read_text() == "new"


def test_identical_private_launcher_is_idempotent_on_resume(tmp_path):
    candidate = tmp_path / "launcher"
    candidate.write_text("#!/bin/sh\nnew\n")
    target = publish_launcher(candidate, tmp_path / "runners")
    before = target.stat()
    assert publish_launcher(candidate, tmp_path / "runners") == target
    assert target.stat().st_ino == before.st_ino


def test_interrupted_asset_copy_leaves_no_destination_or_partial_resume_input(monkeypatch, tmp_path):
    source = tmp_path / "source"
    source.write_bytes(b"complete")
    destination = tmp_path / "data/common/file"

    def interrupted(source, target):
        Path(target).write_bytes(b"partial")
        raise OSError("interrupted copy")

    monkeypatch.setattr(runtime_install.shutil, "copyfile", interrupted)
    with pytest.raises(OSError, match="interrupted copy"):
        copy_checked_asset(source, destination, digest(b"complete"))
    assert not destination.exists()
    assert list(destination.parent.iterdir()) == []


def test_asset_copy_verifies_digest_and_never_overwrites_existing_data(tmp_path):
    source = tmp_path / "source"
    source.write_bytes(b"complete")
    destination = tmp_path / "data/file"
    with pytest.raises(ValueError, match="SHA256"):
        copy_checked_asset(source, destination, digest(b"wrong"))
    assert not destination.exists()
    copy_checked_asset(source, destination, digest(b"complete"))
    assert destination.read_bytes() == b"complete"
    destination.write_bytes(b"existing")
    with pytest.raises(FileExistsError):
        copy_checked_asset(source, destination, digest(b"complete"))
    assert destination.read_bytes() == b"existing"


def test_runtime_assets_can_be_inspected_without_loading_scientific_namespace():
    program = (
        "import sys; from pathlib import Path; from dvbfixer.installation_profiles import load_catalog; "
        "from dvbfixer.runtime_install import asset_requirements; "
        f"r=asset_requirements(load_catalog()['external']['protpardelle-cuda'],Path({str(ROOT)!r})); "
        "assert r['checkpoint']['name']=='cc94_epoch3100.pth'; "
        "assert 'numpy' not in sys.modules and 'torch' not in sys.modules"
    )
    subprocess.run([sys.executable, "-c", program], check=True,
                   env={**os.environ, "PYTHONPATH": str(ROOT / "src")})


def test_runtime_environment_cleans_mps_tuning_and_never_inherits_license(monkeypatch, tmp_path):
    for key in ("KEY_MODELLER", "PYTORCH_MPS_FAST_MATH", "PYTORCH_MPS_HIGH_WATERMARK_RATIO"):
        monkeypatch.setenv(key, "must-not-propagate")
    env = runtime_environment({"device": "mps", "source_subdir": "src"}, ROOT, tmp_path)
    assert env["PYTORCH_ENABLE_MPS_FALLBACK"] == "0"
    assert env["PYTHONDONTWRITEBYTECODE"] == "1"
    assert not {"KEY_MODELLER", "PYTORCH_MPS_FAST_MATH", "PYTORCH_MPS_HIGH_WATERMARK_RATIO"} & env.keys()
    assert env["PYTHONPATH"].split(os.pathsep) == [str(tmp_path / "src"), str(ROOT / "src")]


def test_dry_run_and_unsupported_platform_do_not_create_prefix(monkeypatch, tmp_path):
    monkeypatch.setattr(Manager, "validate", lambda self: "2.9.0")
    monkeypatch.setattr("dvbfixer.runtime_install.platform_subdir", lambda: "linux-64")
    prefix = tmp_path / "new"
    result = provision_runtime(name="protpardelle-cuda", prefix=prefix,
                               manager=Manager("micromamba", "micromamba"), repository=ROOT, dry_run=True)
    assert result["state"] == "planned" and result["actions"]
    assert not prefix.exists() and list(tmp_path.iterdir()) == []
    with pytest.raises(ValueError, match="osx-arm64"):
        provision_runtime(name="protpardelle-mps", prefix=prefix,
                          manager=Manager("micromamba", "micromamba"), repository=ROOT)
    assert list(tmp_path.iterdir()) == []


def test_qm_runtime_plan_is_isolated_and_rejects_diffusion_inputs(monkeypatch, tmp_path):
    monkeypatch.setattr(Manager, "validate", lambda self: "2.9.0")
    monkeypatch.setattr(runtime_install, "platform_subdir", lambda: "linux-64")
    prefix = tmp_path / "qm"
    options = dict(name="psi4", prefix=prefix, manager=Manager("micromamba", "micromamba"), repository=ROOT)
    report = provision_runtime(**options, dry_run=True)
    assert report["state"] == "planned"
    assert report["environment_prefix"] == str(prefix / "environment")
    assert {"pydantic<2", "libint=2.9.0"} <= set(report["dependencies"])
    assert list(tmp_path.iterdir()) == []
    with pytest.raises(ValueError, match="no engine/checkpoint/launcher inputs"):
        provision_runtime(**options, checkpoint=tmp_path / "missing.pth")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("assets,passed,expected", [
    (False, False, "action-required"), (True, False, "error"), (True, True, "ready"),
])
def test_only_complete_assets_and_successful_canonical_preflight_publish(
    monkeypatch, tmp_path, assets, passed, expected,
):
    checkpoint = tmp_path / "weights.pth"
    checkpoint.write_bytes(b"weights")
    requirements = {"checkpoint": {"name": checkpoint.name, "sha256": digest(b"weights")}}
    monkeypatch.setattr(runtime_install, "asset_requirements", lambda *args: requirements)
    monkeypatch.setattr(runtime_install, "platform_subdir", lambda: "linux-64")
    monkeypatch.setattr(Manager, "validate", lambda self: "2.9.0")
    monkeypatch.setattr(runtime_install, "prepare_source", lambda *args: None)
    published = tmp_path / "runners"
    calls = []
    launcher_name = "dvbfixer-diffusion-protpardelle-1c-cc94-cuda-v1"

    def run(argv, **kwargs):
        calls.append(argv)
        if "--runner-directory" in argv:
            staging = Path(argv[argv.index("--runner-directory") + 1])
            staging.mkdir(parents=True)
            (staging / launcher_name).write_text("#!/bin/sh\nexit 0\n")
        payload = {"preflight": {"passed": passed}, "runner_directory": str(published)}
        return subprocess.CompletedProcess(argv, 0, json.dumps(payload), "")

    monkeypatch.setattr(runtime_install.subprocess, "run", run)
    options = dict(name="protpardelle-cuda", prefix=tmp_path / "runtime",
                   manager=Manager("micromamba", "micromamba"), repository=ROOT,
                   checkpoint=checkpoint if assets else None, runner_directory=published)
    if expected == "error":
        with pytest.raises(ValueError, match="no launcher was published"):
            provision_runtime(**options)
        report = json.loads((tmp_path / "runtime/evidence/installation.json").read_text())
    else:
        report = provision_runtime(**options)
    assert report["state"] == expected
    assert (published / launcher_name).exists() == (expected == "ready")
    assert any("--runner-directory" in call for call in calls) == assets
    assert report["scientific_acceptance"] == "not-tested"
