"""Tests for the resumable public-CLI Protenix cohort replay."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_runner() -> ModuleType:
    path = ROOT / "scripts/run_protenix_public_cohort.py"
    spec = importlib.util.spec_from_file_location("protenix_public_cohort", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifest(path: Path, cases: list[dict[str, object]]) -> Path:
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "case_count": len(cases),
                "source_cohort_sha256": "c" * 64,
                "cases": cases,
            }
        )
    )
    return path


def _materialized_case(root: Path, case: dict[str, object]) -> Path:
    case_root = root / "cases" / str(case["case_id"])
    workspace = case_root / "workspace"
    (workspace / "input").mkdir(parents=True)
    artifacts: dict[str, dict[str, object]] = {}
    for name, data in (
        ("input/normalized.pdb", b"END\n"),
        ("target.fasta", b">A\nAAAAA\n"),
        ("request.json", b"{}"),
    ):
        artifact = workspace / name
        artifact.write_bytes(data)
        artifacts[name] = {"sha256": _sha256(artifact), "bytes": len(data)}
    accepted = {
        "case_id": case["case_id"],
        "screening_index": case["screening_index"],
        "selected_for_inference": True,
        "artifacts": artifacts,
    }
    accepted_path = case_root / "accepted.json"
    accepted_path.write_text(json.dumps(accepted))
    (case_root / ".complete.json").write_text(
        json.dumps(
            {
                "status": "accepted",
                "case_id": case["case_id"],
                "cohort_sha256": "c" * 64,
                "record": "accepted.json",
                "record_sha256": _sha256(accepted_path),
            }
        )
    )
    return workspace


def _bundle(path: Path, *, seed: int = 7, passed: bool = True) -> None:
    path.mkdir(parents=True)
    references = {}
    for label, suffix in (
        ("pdb", "pdb"),
        ("dat", "dat"),
        ("sampler_trace", "sampler-trace.json"),
        ("provenance", "diffusion.json"),
    ):
        artifact = path / f"candidate.{suffix}"
        artifact.write_text(label)
        references[label] = {"path": artifact.name, "sha256": _sha256(artifact)}
    (path / "bundle.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "requested_profile": "protenix-v1-cuda",
                "status": "success" if passed else "validation_failed",
                "candidates": [
                    {
                        "candidate_id": "candidate",
                        "seed": seed,
                        "validation_passed": passed,
                        "hard_gate_failures": [] if passed else ["geometry"],
                        **references,
                    }
                ],
            }
        )
    )


def test_load_manifest_rejects_membership_drift(tmp_path: Path) -> None:
    runner = _load_runner()
    manifest = _manifest(
        tmp_path / "manifest.json",
        [{"case_id": "one", "screening_index": 1, "gap_length": 5, "target_length": 10}],
    )
    runner.EXPECTED_CASES = 1
    runner.EXPECTED_STRATA = {5: 1}
    _document, cases = runner.load_manifest(manifest)
    assert cases[0]["case_id"] == "one"

    raw = json.loads(manifest.read_text())
    raw["cases"][0]["target_length"] = 513
    manifest.write_text(json.dumps(raw))
    with pytest.raises(runner.PublicCohortError, match="invalid frozen cohort case"):
        runner.load_manifest(manifest)


def test_replay_resumes_valid_bundle_and_runs_missing_case(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    cases = [
        {"case_id": "one", "screening_index": 1, "gap_length": 5, "target_length": 10},
        {"case_id": "two", "screening_index": 2, "gap_length": 10, "target_length": 20},
    ]
    manifest = _manifest(tmp_path / "manifest.json", cases)
    cohort = tmp_path / "cohort"
    for case in cases:
        _materialized_case(cohort, case)
    output = tmp_path / "output"
    output.mkdir()
    _bundle(output / "cases/one/bundle")
    executable = tmp_path / "dvbfixer"
    executable.write_text("#!/bin/sh\n")
    work_parent = tmp_path / "work"
    work_parent.mkdir()
    observed: list[str] = []

    def fake_run_case(**kwargs: object) -> int:
        bundle = kwargs["bundle"]
        assert isinstance(bundle, Path)
        observed.append(bundle.parent.name)
        log = kwargs["log"]
        assert isinstance(log, Path)
        log.write_text("complete")
        _bundle(bundle, passed=False)
        return 0

    monkeypatch.setattr(runner, "EXPECTED_CASES", 2)
    monkeypatch.setattr(runner, "EXPECTED_STRATA", {5: 1, 10: 1})
    monkeypatch.setattr(runner, "_run_case", fake_run_case)
    report = runner.replay(
        manifest_path=manifest,
        cohort_root=cohort,
        output_root=output,
        work_parent=work_parent,
        executable=executable,
        seed=7,
        timeout=900,
    )

    assert observed == ["two"]
    assert report["complete"] is True
    assert report["success_count"] == 1
    assert report["validation_failed_count"] == 1
    assert [case["resumed"] for case in report["cases"]] == [True, False]
    assert json.loads((output / "report.json").read_text()) == report


def test_replay_stops_and_records_operational_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    case = {"case_id": "one", "screening_index": 1, "gap_length": 5, "target_length": 10}
    manifest = _manifest(tmp_path / "manifest.json", [case])
    cohort = tmp_path / "cohort"
    _materialized_case(cohort, case)
    output = tmp_path / "output"
    output.mkdir()
    executable = tmp_path / "dvbfixer"
    executable.write_text("#!/bin/sh\n")
    work_parent = tmp_path / "work"
    work_parent.mkdir()

    def fake_run_case(**kwargs: object) -> int:
        log = kwargs["log"]
        assert isinstance(log, Path)
        log.write_text("failed")
        return 9

    monkeypatch.setattr(runner, "EXPECTED_CASES", 1)
    monkeypatch.setattr(runner, "EXPECTED_STRATA", {5: 1})
    monkeypatch.setattr(runner, "_run_case", fake_run_case)
    with pytest.raises(runner.PublicCohortError, match="failed operationally"):
        runner.replay(
            manifest_path=manifest,
            cohort_root=cohort,
            output_root=output,
            work_parent=work_parent,
            executable=executable,
            seed=7,
            timeout=900,
        )

    failure = json.loads((output / "cases/one/failure.json").read_text())
    assert failure["exit_status"] == 9
    assert not (output / "cases/one/bundle").exists()


def test_dry_run_does_not_claim_completion(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    case = {"case_id": "one", "screening_index": 1, "gap_length": 5, "target_length": 10}
    manifest = _manifest(tmp_path / "manifest.json", [case])
    cohort = tmp_path / "cohort"
    _materialized_case(cohort, case)
    output = tmp_path / "output"
    output.mkdir()
    executable = tmp_path / "dvbfixer"
    executable.write_text("#!/bin/sh\n")
    work_parent = tmp_path / "work"
    work_parent.mkdir()
    monkeypatch.setattr(runner, "EXPECTED_CASES", 1)
    monkeypatch.setattr(runner, "EXPECTED_STRATA", {5: 1})

    report = runner.replay(
        manifest_path=manifest,
        cohort_root=cohort,
        output_root=output,
        work_parent=work_parent,
        executable=executable,
        seed=7,
        timeout=900,
        dry_run=True,
    )

    assert report["completed_case_count"] == 0
    assert report["complete"] is False
    assert not (output / "report.json").exists()
