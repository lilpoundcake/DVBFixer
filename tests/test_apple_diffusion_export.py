"""Tests for the Apple diffusion input/output structure exporter."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parents[1]


def _load_exporter() -> ModuleType:
    path = ROOT / "scripts/export_apple_diffusion_structures.py"
    spec = importlib.util.spec_from_file_location("export_apple_diffusion", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_result(directory: Path, content: str, *, passed: bool) -> None:
    directory.mkdir(parents=True)
    candidate = directory / "candidate.pdb"
    candidate.write_text(content, encoding="utf-8")
    summary = {
        "candidate_sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
        "validation": {
            "passed": passed,
            "hard_gate_failures": [] if passed else ["example-failure"],
        },
    }
    (directory / "summary.json").write_text(json.dumps(summary), encoding="utf-8")


def test_export_uses_validation_first_selection_and_case_id_names(tmp_path: Path) -> None:
    exporter = _load_exporter()
    cohort = tmp_path / "cohort"
    for case_id in ("case-a", "case-b"):
        input_dir = cohort / "cases" / case_id / "workspace" / "input"
        input_dir.mkdir(parents=True)
        (input_dir / "normalized.pdb").write_text(
            f"REMARK input {case_id}\n", encoding="utf-8"
        )
    first = cohort / "cases/case-a/workspace"
    _write_result(first / "prefix-raw", "REMARK raw a\n", passed=True)
    _write_result(first / "prefix-refined", "REMARK refined a\n", passed=True)
    second = cohort / "cases/case-b/workspace"
    _write_result(second / "remaining-raw", "REMARK raw b\n", passed=True)
    _write_result(second / "remaining-refined", "REMARK refined b\n", passed=False)
    manifest = {
        "reused_prefix": {
            "case_count": 1,
            "raw_output_subdir": "prefix-raw",
            "refined_output_subdir": "prefix-refined",
        },
        "sampling": {"remaining_output_subdir": "remaining-raw"},
        "refinement": {"remaining_output_subdir": "remaining-refined"},
        "cases": [
            {"case_id": "case-a", "screening_index": 1},
            {"case_id": "case-b", "screening_index": 2},
        ],
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    output = tmp_path / "export"

    result = exporter.export_structures(
        manifest_path, cohort, output, require_complete=True
    )

    assert result["exported_count"] == 2
    assert (output / "case-a/case-a_input.pdb").is_file()
    assert (output / "case-a/case-a_output.pdb").read_text() == "REMARK refined a\n"
    assert (output / "case-b/case-b_output.pdb").read_text() == "REMARK raw b\n"
    metadata = json.loads((output / "case-b/case-b_metadata.json").read_text())
    assert metadata["output"]["selection"] == "raw"
    assert metadata["output"]["validation_passed"] is True
