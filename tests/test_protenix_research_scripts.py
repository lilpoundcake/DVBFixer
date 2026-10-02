"""Focused tests for the private Protenix research wrappers."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    AtomIdentity,
    BackendProvenance,
    DiffusionStatus,
    RunnerDiagnostics,
    RunnerResult,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
CHECKPOINT_SMOKE = REPO_ROOT / "deploy/protenix-v1/checkpoint_gap_smoke.py"


def _load_refinement_script() -> ModuleType:
    path = REPO_ROOT / "deploy/protenix-v1/refine_candidate.py"
    spec = importlib.util.spec_from_file_location("protenix_refine_candidate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_production_runner() -> ModuleType:
    path = REPO_ROOT / "deploy/protenix-v1/production_runner.py"
    spec = importlib.util.spec_from_file_location("protenix_production_runner", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _atom_line(serial: int, residue_number: int, atom_name: str, x: float) -> str:
    return (
        f"ATOM  {serial:5d}  {atom_name:<3} ALA A{residue_number:4d}    "
        f"{x:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00  0.00           C\n"
    )


def test_refinement_rewrites_only_requested_coordinate_columns() -> None:
    module = _load_refinement_script()
    fixed = _atom_line(1, 1, "CA", 1.0)
    generated = _atom_line(2, 2, "CA", 4.0)
    text = "HEADER preserved\n" + fixed + generated + "END\n"
    identity = AtomIdentity("A", "2", "", "CA")

    rewritten = module._rewrite_generated_coordinates(
        text,
        {identity: np.asarray([7.0, 8.0, 9.0])},
    )

    original_lines = text.splitlines(keepends=True)
    rewritten_lines = rewritten.splitlines(keepends=True)
    assert rewritten_lines[0] == original_lines[0]
    assert rewritten_lines[1] == original_lines[1]
    assert rewritten_lines[2][:30] == original_lines[2][:30]
    assert rewritten_lines[2][30:54] == "   7.000   8.000   9.000"
    assert rewritten_lines[2][54:] == original_lines[2][54:]
    assert rewritten_lines[3] == original_lines[3]


def test_refinement_rejects_missing_or_duplicate_generated_atoms() -> None:
    module = _load_refinement_script()
    generated = _atom_line(2, 2, "CA", 4.0)
    identity = AtomIdentity("A", "2", "", "CA")
    coordinates = {identity: np.asarray([7.0, 8.0, 9.0])}

    with pytest.raises(ValueError, match="omits 1 generated atoms"):
        module._rewrite_generated_coordinates("END\n", coordinates)
    with pytest.raises(ValueError, match="duplicate generated atom"):
        module._rewrite_generated_coordinates(generated + generated, coordinates)


def test_checkpoint_smoke_reports_peak_vram() -> None:
    smoke = CHECKPOINT_SMOKE.read_text(encoding="utf-8")

    assert "torch.cuda.reset_peak_memory_stats()" in smoke
    assert "peak_vram_bytes=torch.cuda.max_memory_allocated()" in smoke
    assert '"peak_vram_bytes": runner_result.resource_metrics.peak_vram_bytes' in smoke


def test_production_runner_writes_protocol_result(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    production = _load_production_runner()
    checkpoint = tmp_path / "protenix_base_default_v1.0.0.pt"
    checkpoint.write_bytes(b"checkpoint")
    (tmp_path / "request.json").write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(production.shutil, "which", lambda _name: "/usr/bin/true")

    expected = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.FAILED,
        candidates=(),
        runner_diagnostics=RunnerDiagnostics(exit_code=1, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="fake-protenix",
            runner_protocol_version=3,
            engine_repository="builtin://test",
            engine_revision="test",
            device="cuda:0",
            framework_version="2.13.0+cu129",
            cuda_version="12.9",
        ),
        message="test result",
    )

    class FakeBuilder:
        @staticmethod
        def build_input(_request: Path, output: Path) -> Path:
            output.mkdir(parents=True)
            path = output / "input.json"
            path.write_text("{}")
            return path

    class FakeAdapter:
        @staticmethod
        def run(*args: object, **kwargs: object) -> tuple[RunnerResult, Path]:
            assert kwargs["cycles"] == 1
            assert kwargs["steps"] == 200
            assert (
                kwargs["ablation_mode"].value
                == "template-conditioning-plus-reinjection"
            )
            return expected, Path("summary.json")

    observed = production.run(
        profile="protenix-v1-cuda",
        checkpoint=checkpoint,
        scripts=(FakeBuilder(), FakeAdapter()),
        refiner=lambda _request, result, _output: result,
        preflight=lambda: None,
    )

    assert observed.status is DiffusionStatus.FAILED
    assert observed.message == expected.message
    assert observed.backend_provenance.source_license == "Apache-2.0"
    assert observed.backend_provenance.environment_identity.startswith("linux-amd64")
    assert RunnerResult.from_json((tmp_path / "result.json").read_text()) == observed


def test_production_runner_rejects_wrong_profile(tmp_path: Path) -> None:
    production = _load_production_runner()
    with pytest.raises(ValueError, match="unsupported Protenix profile"):
        production.run(
            profile="protpardelle-1c-mps",
            checkpoint=tmp_path / "checkpoint.pt",
        )
