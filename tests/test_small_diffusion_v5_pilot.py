"""Tests for the resumable compact-model v5 pilot driver."""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_driver() -> ModuleType:
    path = ROOT / "scripts/run_small_diffusion_v5_pilot.py"
    spec = importlib.util.spec_from_file_location("run_small_diffusion_v5_pilot", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_completion_marker_requires_current_input_and_output_digests(tmp_path: Path) -> None:
    driver = _load_driver()
    source = tmp_path / "source"
    output = tmp_path / "output"
    marker = tmp_path / "stage.json"
    source.write_text("input", encoding="utf-8")
    output.write_text("result", encoding="utf-8")
    marker.write_text(
        json.dumps(
            {
                "status": "complete",
                "inputs": {"source": driver._sha256(source)},
                "outputs": {"output": driver._sha256(output)},
            }
        ),
        encoding="utf-8",
    )

    assert driver._completion_valid(marker, {"source": source}, {"output": output})
    source.write_text("changed", encoding="utf-8")
    assert not driver._completion_valid(marker, {"source": source}, {"output": output})


def test_preregistered_manifest_is_balanced_unique_and_ordered() -> None:
    value = json.loads(
        (ROOT / "docs/research/small-diffusion-v5-expanded-pilot.json").read_text(
            encoding="utf-8"
        )
    )
    cases = value["cases"]

    assert len(cases) == 24
    assert len({case["case_id"] for case in cases}) == 24
    assert sum(case["gap_length"] == 5 for case in cases) == 12
    assert sum(case["gap_length"] == 10 for case in cases) == 12
    for gap_length in (5, 10):
        indices = [
            case["screening_index"] for case in cases if case["gap_length"] == gap_length
        ]
        assert indices == sorted(indices)
    assert not set(value["pilot_case_ids_excluded"]) & {
        case["case_id"] for case in cases
    }


def test_selection_retains_valid_raw_when_refinement_regresses() -> None:
    driver = _load_driver()

    assert driver._selection_source(
        {"validation": {"passed": True}},
        {"validation": {"passed": False}},
    ) == "raw"
    assert driver._selection_source(
        {"validation": {"passed": False}},
        {"validation": {"passed": True}},
    ) == "refined"
    assert (
        driver._selection_source(
            {"validation": {"passed": False}},
            {"validation": {"passed": False}},
        )
        is None
    )


def test_output_namespace_is_one_contained_component(tmp_path: Path) -> None:
    driver = _load_driver()

    assert driver._case_output_root(tmp_path, "reinjection-followup") == (
        tmp_path / "reinjection-followup"
    )
    for invalid in ("", "../escape", "nested/output", str(tmp_path / "absolute")):
        with pytest.raises(driver.PilotError, match="one relative path component"):
            driver._case_output_root(tmp_path, invalid)


def test_driver_exposes_isolated_per_step_mode() -> None:
    source = (ROOT / "scripts/run_small_diffusion_v5_pilot.py").read_text(encoding="utf-8")

    assert 'run_parser.add_argument("--per-step-reinjection", action="store_true")' in source
    assert 'run_parser.add_argument("--output-subdir"' in source
    assert '"git",' in source and '"--reverse",' in source and '"--check",' in source


def test_failed_stage_writes_only_failure_marker(tmp_path: Path) -> None:
    driver = _load_driver()
    source = tmp_path / "source"
    output = tmp_path / "output"
    marker = tmp_path / "stage.json"
    source.write_text("input", encoding="utf-8")

    completed = driver._run_stage(
        marker,
        (sys.executable, "-c", "raise SystemExit(7)"),
        {"source": source},
        {"output": output},
        env=os.environ.copy(),
    )

    assert completed is False
    assert not marker.exists()
    failure = json.loads((tmp_path / "stage.failure.json").read_text(encoding="utf-8"))
    assert failure["status"] == "failed"
    assert failure["exit_status"] == 7
