"""Tests for the Apple Silicon diffusion inventory collector."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_inventory() -> ModuleType:
    path = ROOT / "scripts/collect_apple_diffusion_inventory.py"
    spec = importlib.util.spec_from_file_location("collect_apple_diffusion_inventory", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_operator_comparison() -> ModuleType:
    path = ROOT / "deploy/protpardelle-1c/compare_operator_smokes.py"
    spec = importlib.util.spec_from_file_location("compare_operator_smokes", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_smoke(path: Path, device: str, coordinates: np.ndarray) -> None:
    coordinate_path = path.with_name(f"{path.stem}.coordinates.npy")
    np.save(coordinate_path, coordinates, allow_pickle=False)
    path.write_text(
        json.dumps(
            {
                "status": "passed",
                "device": device,
                "coordinate_artifact": coordinate_path.name,
                "coordinate_artifact_sha256": hashlib.sha256(
                    coordinate_path.read_bytes()
                ).hexdigest(),
                "shape": list(coordinates.shape),
                "config_sha256": "1" * 64,
                "checkpoint_sha256": "2" * 64,
                "length": coordinates.shape[1],
                "steps": 1,
                "seed": 7,
                "sample_seconds": 1.0,
            }
        ),
        encoding="utf-8",
    )


def test_inventory_rejects_non_apple_host_before_importing_torch(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    inventory = _load_inventory()
    monkeypatch.setattr(inventory.platform, "system", lambda: "Linux")
    monkeypatch.setattr(inventory.platform, "machine", lambda: "x86_64")

    with pytest.raises(RuntimeError, match="native macOS arm64"):
        inventory.collect(tmp_path / "missing-lock", ROOT)


def test_inventory_requires_fallback_disabled_before_importing_torch(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    inventory = _load_inventory()
    monkeypatch.setattr(inventory.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(inventory.platform, "machine", lambda: "arm64")
    monkeypatch.delenv("PYTORCH_ENABLE_MPS_FALLBACK", raising=False)

    with pytest.raises(RuntimeError, match="PYTORCH_ENABLE_MPS_FALLBACK=0"):
        inventory.collect(tmp_path / "missing-lock", ROOT)


def test_operator_comparison_verifies_pair_and_reports_coordinate_divergence(
    tmp_path: Path,
) -> None:
    comparison = _load_operator_comparison()
    cpu = np.zeros((1, 2, 1, 3), dtype=np.float32)
    mps = cpu.copy()
    mps[0, 1, 0, 0] = 2.0
    cpu_path = tmp_path / "cpu.json"
    mps_path = tmp_path / "mps.json"
    _write_smoke(cpu_path, "cpu", cpu)
    _write_smoke(mps_path, "mps:0", mps)

    result = comparison.compare(cpu_path, mps_path)

    assert result["coordinate_byte_identical"] is False
    assert result["coordinate_maximum_displacement"] == 2.0
    assert result["coordinate_mean_displacement"] == 1.0
    assert result["coordinate_rmsd"] == pytest.approx(2**0.5)


def test_operator_smoke_requires_explicit_mps_fallback_policy() -> None:
    source = (ROOT / "deploy/protpardelle-1c/operator_smoke.py").read_text(
        encoding="utf-8"
    )

    assert 'os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK") != "0"' in source
    assert 'coordinates.device.type != device.type' in source
    assert '"record_trajectory": False' in source
