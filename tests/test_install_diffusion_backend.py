"""Tests for private diffusion-launcher registration."""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parent.parent
INSTALLER = ROOT / "scripts/install_diffusion_backend.py"


def _load_installer() -> ModuleType:
    spec = importlib.util.spec_from_file_location("install_diffusion_backend", INSTALLER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_protenix_launcher_binds_current_source_engine_and_kalign(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    installer = _load_installer()
    repository = tmp_path / "repository"
    scripts = repository / "scripts"
    runner = repository / "deploy/protenix-v1/production_runner.py"
    source = repository / "src"
    engine = tmp_path / "Protenix"
    destination = tmp_path / "runners"
    scripts.mkdir(parents=True)
    runner.parent.mkdir(parents=True)
    source.mkdir()
    engine.mkdir()
    runner.write_text("# runner\n")
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"checkpoint")
    kalign = tmp_path / "kalign-bin/kalign"
    kalign.parent.mkdir()
    kalign.write_bytes(b"kalign")
    kalign.chmod(0o700)

    monkeypatch.setattr(installer, "__file__", str(scripts / INSTALLER.name))
    monkeypatch.setattr(installer.sys, "executable", "/backend/bin/python")
    monkeypatch.setattr(installer, "diffusion_runner_directory", lambda: destination)
    monkeypatch.setattr(
        installer,
        "BACKENDS",
        {
            "protenix": installer.Backend(
                "deploy/protenix-v1/production_runner.py",
                "dvbfixer-diffusion-protenix-v1-cuda",
                _sha256(checkpoint),
                ".",
                kalign_sha256=_sha256(kalign),
            )
        },
    )
    monkeypatch.setattr(
        installer.sys,
        "argv",
        [
            str(INSTALLER),
            "protenix",
            str(checkpoint),
            "--engine-root",
            str(engine),
            "--kalign",
            str(kalign),
        ],
    )

    installer.main()

    launcher = destination / "dvbfixer-diffusion-protenix-v1-cuda"
    text = launcher.read_text()
    assert launcher.stat().st_mode & 0o700 == 0o700
    assert "export PYTHONNOUSERSITE=1" in text
    assert "export PYTHONDONTWRITEBYTECODE=1" in text
    assert f"export PYTHONPATH={engine}:{source}" in text
    assert f"export PATH={kalign.parent}:$PATH" in text
    assert str(runner) in text
    assert str(checkpoint) in text
    assert "/backend/bin/python" in text


def test_protenix_registration_requires_kalign(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    installer = _load_installer()
    repository = tmp_path / "repository"
    scripts = repository / "scripts"
    runner = repository / "deploy/protenix-v1/production_runner.py"
    engine = tmp_path / "Protenix"
    scripts.mkdir(parents=True)
    runner.parent.mkdir(parents=True)
    engine.mkdir()
    runner.write_text("# runner\n")
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"checkpoint")

    monkeypatch.setattr(installer, "__file__", str(scripts / INSTALLER.name))
    monkeypatch.setattr(
        installer,
        "BACKENDS",
        {
            "protenix": installer.Backend(
                "deploy/protenix-v1/production_runner.py",
                "dvbfixer-diffusion-protenix-v1-cuda",
                _sha256(checkpoint),
                ".",
                kalign_sha256="0" * 64,
            )
        },
    )
    monkeypatch.setattr(
        installer.sys,
        "argv",
        [
            str(INSTALLER),
            "protenix",
            str(checkpoint),
            "--engine-root",
            str(engine),
        ],
    )

    with pytest.raises(SystemExit, match="2"):
        installer.main()
