"""Focused tests for the private Protenix research wrappers."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity

REPO_ROOT = Path(__file__).resolve().parent.parent
CHECKPOINT_SMOKE = REPO_ROOT / "deploy/protenix-v1/checkpoint_gap_smoke.py"


def _load_refinement_script() -> ModuleType:
    path = REPO_ROOT / "deploy/protenix-v1/refine_candidate.py"
    spec = importlib.util.spec_from_file_location("protenix_refine_candidate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
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
