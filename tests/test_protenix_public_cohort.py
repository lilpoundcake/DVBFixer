"""Tests for the resumable public-CLI Protenix cohort replay."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

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
        (
            "request.json",
            json.dumps(
                {
                    "candidate_count": 1,
                    "seeds": [7],
                    "backend_options": [
                        {"name": "profile", "value": "protenix-v1-cuda"}
                    ],
                }
            ).encode(),
        ),
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


def test_public_command_passes_request_without_construction_options(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    observed = []

    def fake_run(command, **kwargs):
        observed.extend(command)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    assert runner._run_case(
        executable=tmp_path / "dvbfixer",
        workspace=tmp_path,
        bundle=tmp_path / "bundle",
        work_parent=tmp_path,
        seed=7,
        timeout=3600,
        log=tmp_path / "command.log",
    ) == 0
    assert "--diffusion-request" in observed
    assert "--fasta" not in observed
    assert "--diffusion-seed" not in observed


def test_migration_preserves_frozen_masks_and_original_request(tmp_path: Path) -> None:
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

    runner = _load_runner()
    source = tmp_path / "source"
    (source / "input").mkdir(parents=True)
    fixture = ROOT / "tests/fixtures/8cz8/8cz8_a_u.pdb"
    retained = {61, 62, 63, 64, 70, 71, 72, 73}
    text = "".join(
        line for line in fixture.read_text().splitlines(keepends=True)
        if line.startswith("ATOM  ") and line[21] == "C"
        and int(line[22:26]) in retained
    ) + "TER\nEND\n"
    pdb = source / "input/normalized.pdb"
    pdb.write_text(text)
    request = build_cli_diffusion_request(
        pdb, {"C": "SNRFSGSKSGNTA"}, seeds=(7,), profile="protenix-v1-cuda",
    )
    raw = request.to_dict()
    raw["normalized_pdb"]["path"] = "input/normalized.pdb"
    raw["schema_version"] = 3
    (source / "request.json").write_text(json.dumps(raw))
    before = (source / "request.json").read_bytes()
    destination = runner.prepare_replay_request(source, tmp_path / "replay-input")
    migrated = DiffusionRequest.from_json((destination / "request.json").read_text())
    assert migrated.schema_version == 4
    assert migrated.gaps == request.gaps
    assert migrated.sequence_placements == request.sequence_placements
    assert migrated.seeds == request.seeds
    assert (source / "request.json").read_bytes() == before
    assert (destination / "input/normalized.pdb").read_bytes() == pdb.read_bytes()


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
    monkeypatch.setattr(runner, "prepare_replay_request", lambda source, _destination: source)
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
    assert report["timeout_seconds"] == 900
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
    monkeypatch.setattr(runner, "prepare_replay_request", lambda source, _destination: source)
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


def test_archive_survives_loss_of_temporary_publication_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    import shutil

    runner = _load_runner()
    case = {"case_id": "one", "screening_index": 1, "gap_length": 5, "target_length": 10}
    manifest = _manifest(tmp_path / "manifest.json", [case])
    cohort = tmp_path / "cohort"
    _materialized_case(cohort, case)
    output = tmp_path / "temporary"
    case_root = output / "cases/one"
    _bundle(case_root / "bundle")
    (case_root / "dvbfixer.log").write_text("accepted")
    archive = tmp_path / "durable"
    runner.archive_case(case_root, archive, case_id="one", seed=7)
    assert (archive / "cases/one/dvbfixer.log").read_text() == "accepted"
    shutil.rmtree(output)
    executable = tmp_path / "dvbfixer"
    executable.write_text("#!/bin/sh\n")
    monkeypatch.setattr(runner, "EXPECTED_CASES", 1)
    monkeypatch.setattr(runner, "EXPECTED_STRATA", {5: 1})
    monkeypatch.setattr(runner, "_run_case", lambda **_kwargs: pytest.fail("archive must resume"))
    report = runner.replay(
        manifest_path=manifest, cohort_root=cohort, output_root=output,
        work_parent=tmp_path, executable=executable, seed=7, timeout=3600,
        archive_root=archive,
    )
    assert report["complete"] is True
    assert report["cases"][0]["archived"] is True
    assert json.loads((archive / "report.json").read_text()) == report


def test_corrupt_bundle_is_never_archived(tmp_path: Path) -> None:
    runner = _load_runner()
    case_root = tmp_path / "case"
    _bundle(case_root / "bundle")
    (case_root / "bundle/candidate.pdb").write_text("corrupt")
    archive = tmp_path / "archive"
    with pytest.raises(runner.PublicCohortError, match="digest mismatch"):
        runner.archive_case(case_root, archive, case_id="one", seed=7)
    assert not (archive / "cases/one").exists()
