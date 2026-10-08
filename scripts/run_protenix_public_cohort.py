#!/usr/bin/env python3
"""Run the frozen 231-case Protenix cohort through the public DVBFixer CLI."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import stat
import subprocess
import time
from pathlib import Path
from typing import Any

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    DiffusionContractError,
    DiffusionRequest,
)
from dvbfixer.model.diffusion.scope import assess_diffusion_scope

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "docs/research/small-diffusion-full-followup-v2.json"
FROZEN_MANIFEST_SHA256 = "fa2bd455ef44c0a413cc2a174a4dfd017af51a99617fb69eb7bdc364eba0850f"
PROFILE = "protenix-v1-cuda"
EXPECTED_CASES = 231
EXPECTED_STRATA = {5: 128, 10: 103}
REPORT_SCHEMA_VERSION = 1
BUNDLE_SCHEMA_VERSION = 2
DEFAULT_TIMEOUT_SECONDS = 3600


class PublicCohortError(RuntimeError):
    """Raised when replay inputs, execution, or published evidence are invalid."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise PublicCohortError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PublicCohortError(f"JSON document must be an object: {path}")
    return value


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    data = (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode()
    try:
        with temporary.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    finally:
        temporary.unlink(missing_ok=True)


def load_manifest(path: Path) -> tuple[dict[str, Any], tuple[dict[str, Any], ...]]:
    manifest = _read_json(path)
    cases = manifest.get("cases")
    if manifest.get("schema_version") != 1 or not isinstance(cases, list):
        raise PublicCohortError("unsupported frozen cohort manifest")
    if manifest.get("case_count") != EXPECTED_CASES or len(cases) != EXPECTED_CASES:
        raise PublicCohortError(f"frozen cohort must contain exactly {EXPECTED_CASES} cases")

    normalized: list[dict[str, Any]] = []
    for case in cases:
        if not isinstance(case, dict):
            raise PublicCohortError("frozen cohort case must be an object")
        case_id = case.get("case_id")
        screening_index = case.get("screening_index")
        gap_length = case.get("gap_length")
        target_length = case.get("target_length")
        if (
            not isinstance(case_id, str)
            or not case_id
            or not isinstance(screening_index, int)
            or screening_index < 0
            or gap_length not in EXPECTED_STRATA
            or not isinstance(target_length, int)
            or not 0 < target_length <= 512
        ):
            raise PublicCohortError(f"invalid frozen cohort case: {case!r}")
        normalized.append(case)

    case_ids = [case["case_id"] for case in normalized]
    indices = [case["screening_index"] for case in normalized]
    if len(set(case_ids)) != len(case_ids) or len(set(indices)) != len(indices):
        raise PublicCohortError("frozen cohort repeats a case ID or screening index")
    if indices != sorted(indices):
        raise PublicCohortError("frozen cohort is not in screening-index order")
    observed_strata = {
        gap_length: sum(case["gap_length"] == gap_length for case in normalized)
        for gap_length in EXPECTED_STRATA
    }
    if observed_strata != EXPECTED_STRATA:
        raise PublicCohortError(f"frozen cohort strata changed: {observed_strata!r}")
    return manifest, tuple(normalized)


def _regular_file(path: Path) -> bool:
    try:
        return stat.S_ISREG(path.lstat().st_mode)
    except OSError:
        return False


def validate_materialized_case(
    cohort_root: Path,
    case: dict[str, Any],
    *,
    source_cohort_sha256: str,
    seed: int,
) -> Path:
    case_root = cohort_root / "cases" / case["case_id"]
    completion = _read_json(case_root / ".complete.json")
    if (
        completion.get("status") != "accepted"
        or completion.get("case_id") != case["case_id"]
        or completion.get("cohort_sha256") != source_cohort_sha256
        or completion.get("record") != "accepted.json"
    ):
        raise PublicCohortError(f"invalid materialization marker for {case['case_id']}")
    accepted_path = case_root / "accepted.json"
    if not _regular_file(accepted_path) or _sha256(accepted_path) != completion.get("record_sha256"):
        raise PublicCohortError(f"materialization record digest mismatch for {case['case_id']}")
    accepted = _read_json(accepted_path)
    if (
        accepted.get("case_id") != case["case_id"]
        or accepted.get("screening_index") != case["screening_index"]
        or accepted.get("selected_for_inference") is not True
    ):
        raise PublicCohortError(f"materialization metadata mismatch for {case['case_id']}")
    artifacts = accepted.get("artifacts")
    if not isinstance(artifacts, dict):
        raise PublicCohortError(f"missing materialization artifacts for {case['case_id']}")
    for name in ("input/normalized.pdb", "target.fasta", "request.json"):
        metadata = artifacts.get(name)
        artifact = case_root / "workspace" / name
        if (
            not isinstance(metadata, dict)
            or not _regular_file(artifact)
            or metadata.get("sha256") != _sha256(artifact)
        ):
            raise PublicCohortError(
                f"materialized artifact digest mismatch for {case['case_id']}: {name}"
            )
    request = _read_json(case_root / "workspace/request.json")
    profile = next(
        (
            option.get("value")
            for option in request.get("backend_options", [])
            if isinstance(option, dict) and option.get("name") == "profile"
        ),
        None,
    )
    if (
        request.get("candidate_count") != 1
        or request.get("seeds") != [seed]
        or profile not in {None, PROFILE}
    ):
        raise PublicCohortError(f"materialized request policy mismatch for {case['case_id']}")
    return case_root / "workspace"


def _contained_artifact(bundle: Path, value: object, label: str) -> tuple[str, str]:
    if not isinstance(value, dict):
        raise PublicCohortError(f"bundle {label} reference is missing")
    relative = value.get("path")
    digest = value.get("sha256")
    if not isinstance(relative, str) or not isinstance(digest, str):
        raise PublicCohortError(f"bundle {label} reference is invalid")
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise PublicCohortError(f"bundle {label} path escapes its directory")
    artifact = bundle / relative_path
    if not _regular_file(artifact) or _sha256(artifact) != digest:
        raise PublicCohortError(f"bundle {label} digest mismatch")
    return relative, digest


def prepare_replay_request(source: Path, destination: Path) -> Path:
    """Explicitly migrate frozen v3 internal-gap requests into a private v4 copy.

    V4 added optional terminal/context/ownership fields. V3's persisted internal
    placements, masks and seeds are retained verbatim and validated by the current
    contract and scope admission; original materialization is never overwritten.
    """
    raw = _read_json(source / "request.json")
    version = raw.get("schema_version")
    if version not in {3, DIFFUSION_SCHEMA_VERSION}:
        raise PublicCohortError(f"unsupported frozen request schema {version}")
    if version == 3:
        raw["schema_version"] = DIFFUSION_SCHEMA_VERSION
    try:
        request = DiffusionRequest.from_dict(raw)
    except DiffusionContractError as exc:
        raise PublicCohortError(f"invalid migrated request: {exc}") from exc
    source_bytes = (source / "input/normalized.pdb").read_bytes()
    admission = assess_diffusion_scope(request, source_bytes)
    if not admission.supported:
        raise PublicCohortError(f"migrated request is outside scope: {admission.reasons}")
    if request.normalized_pdb.path != "input/normalized.pdb":
        raise PublicCohortError("frozen request must reference input/normalized.pdb")
    (destination / "input").mkdir(parents=True, exist_ok=True)
    (destination / "input/normalized.pdb").write_bytes(source_bytes)
    _atomic_json(destination / "request.json", request.to_dict())
    _atomic_json(
        destination / "migration.json",
        {
            "source_schema_version": version,
            "schema_version": DIFFUSION_SCHEMA_VERSION,
            "source_request_sha256": _sha256(source / "request.json"),
            "replay_request_sha256": _sha256(destination / "request.json"),
        },
    )
    return destination


def validate_bundle(bundle: Path, *, case_id: str, seed: int) -> dict[str, Any]:
    index_path = bundle / "bundle.json"
    if not _regular_file(index_path):
        raise PublicCohortError(f"case {case_id} has no regular bundle.json")
    index = _read_json(index_path)
    if (
        index.get("schema_version") != BUNDLE_SCHEMA_VERSION
        or index.get("requested_profile") != PROFILE
        or index.get("status") not in {"success", "validation_failed"}
    ):
        raise PublicCohortError(f"case {case_id} has an invalid public bundle")
    candidates = index.get("candidates")
    if not isinstance(candidates, list) or len(candidates) != 1:
        raise PublicCohortError(f"case {case_id} must publish exactly one candidate")
    candidate = candidates[0]
    if not isinstance(candidate, dict) or candidate.get("seed") != seed:
        raise PublicCohortError(f"case {case_id} published the wrong seed")
    artifacts = {
        label: _contained_artifact(bundle, candidate.get(label), label)
        for label in ("pdb", "dat", "sampler_trace", "provenance")
    }
    passed = candidate.get("validation_passed")
    failures = candidate.get("hard_gate_failures")
    if not isinstance(passed, bool) or not isinstance(failures, list) or not all(
        isinstance(item, str) for item in failures
    ):
        raise PublicCohortError(f"case {case_id} has invalid validation evidence")
    if (index["status"] == "success") != passed or passed == bool(failures):
        raise PublicCohortError(f"case {case_id} has inconsistent validation status")
    return {
        "case_id": case_id,
        "status": index["status"],
        "candidate_id": candidate.get("candidate_id"),
        "seed": seed,
        "validation_passed": passed,
        "hard_gate_failures": failures,
        "bundle_json_sha256": _sha256(index_path),
        "artifacts": {
            label: {"path": path, "sha256": digest}
            for label, (path, digest) in artifacts.items()
        },
    }


def _run_case(
    *,
    executable: Path,
    workspace: Path,
    bundle: Path,
    work_parent: Path,
    seed: int,
    timeout: int,
    log: Path,
) -> int:
    command = (
        str(executable),
        "model",
        str(workspace / "input/normalized.pdb"),
        "--backend",
        "diffusion",
        "--diffusion-model",
        "protenix",
        "--diffusion-request",
        str(workspace / "request.json"),
        "--diffusion-timeout",
        str(timeout),
        "--diffusion-work-parent",
        str(work_parent),
        "-o",
        str(bundle),
    )
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("xb") as stream:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return completed.returncode


def _report(
    *,
    manifest_path: Path,
    manifest_sha256: str,
    cohort_root: Path,
    seed: int,
    timeout: int,
    cases: tuple[dict[str, Any], ...],
) -> dict[str, Any]:
    successful = sum(case["status"] == "success" for case in cases)
    rejected = sum(case["status"] == "validation_failed" for case in cases)
    completed = successful + rejected
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "profile": PROFILE,
        "seed": seed,
        "timeout_seconds": timeout,
        "manifest": str(manifest_path),
        "manifest_sha256": manifest_sha256,
        "cohort_root": str(cohort_root),
        "expected_case_count": EXPECTED_CASES,
        "completed_case_count": completed,
        "success_count": successful,
        "validation_failed_count": rejected,
        "complete": completed == EXPECTED_CASES,
        "cases": list(cases),
    }


def replay(
    *,
    manifest_path: Path,
    cohort_root: Path,
    output_root: Path,
    work_parent: Path,
    executable: Path,
    seed: int,
    timeout: int,
    selected_case_ids: frozenset[str] = frozenset(),
    max_cases: int | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    manifest, manifest_cases = load_manifest(manifest_path)
    source_cohort_sha256 = manifest.get("source_cohort_sha256")
    if not isinstance(source_cohort_sha256, str):
        raise PublicCohortError("manifest has no source cohort digest")
    source_aggregate_sha256 = manifest.get("source_aggregate_sha256")
    if isinstance(source_aggregate_sha256, str):
        source_aggregate = cohort_root / "aggregate.json"
        if not _regular_file(source_aggregate) or _sha256(source_aggregate) != source_aggregate_sha256:
            raise PublicCohortError("source aggregate digest does not match the frozen manifest")
    if not work_parent.is_dir() or not output_root.parent.is_dir():
        raise PublicCohortError("work parent and output parent must already exist")
    if seed < 0 or timeout <= 0:
        raise PublicCohortError("seed must be non-negative and timeout must be positive")
    if output_root.is_symlink():
        raise PublicCohortError("output root must not be a symlink")
    output_root.mkdir(exist_ok=True)
    if not _regular_file(executable):
        raise PublicCohortError(f"DVBFixer executable is unavailable: {executable}")

    cases = tuple(
        case
        for case in manifest_cases
        if not selected_case_ids or case["case_id"] in selected_case_ids
    )
    if selected_case_ids - {case["case_id"] for case in cases}:
        raise PublicCohortError("requested case is not in the frozen cohort")
    if max_cases is not None:
        if max_cases <= 0:
            raise PublicCohortError("max_cases must be positive")
        cases = cases[:max_cases]

    results: list[dict[str, Any]] = []
    for case in cases:
        case_id = case["case_id"]
        workspace = validate_materialized_case(
            cohort_root,
            case,
            source_cohort_sha256=source_cohort_sha256,
            seed=seed,
        )
        case_root = output_root / "cases" / case_id
        bundle = case_root / "bundle"
        if bundle.exists():
            result = validate_bundle(bundle, case_id=case_id, seed=seed)
            results.append({**case, **result, "resumed": True})
            continue
        if dry_run:
            results.append({**case, "status": "pending", "resumed": False})
            continue

        case_root.mkdir(parents=True, exist_ok=True)
        workspace = prepare_replay_request(workspace, case_root / "replay-input")
        log = case_root / "dvbfixer.log"
        failure = case_root / "failure.json"
        log.unlink(missing_ok=True)
        failure.unlink(missing_ok=True)
        started = time.time()
        exit_status = _run_case(
            executable=executable,
            workspace=workspace,
            bundle=bundle,
            work_parent=work_parent,
            seed=seed,
            timeout=timeout,
            log=log,
        )
        if exit_status != 0:
            record = {
                "schema_version": REPORT_SCHEMA_VERSION,
                "case_id": case_id,
                "exit_status": exit_status,
                "started_unix_seconds": started,
                "wall_time_seconds": time.time() - started,
                "log_sha256": _sha256(log),
            }
            _atomic_json(failure, record)
            raise PublicCohortError(
                f"case {case_id} failed operationally with exit status {exit_status}"
            )
        result = validate_bundle(bundle, case_id=case_id, seed=seed)
        results.append(
            {**case, **result, "resumed": False, "log_sha256": _sha256(log)}
        )
        _atomic_json(
            output_root / "report.json",
            _report(
                manifest_path=manifest_path,
                manifest_sha256=_sha256(manifest_path),
                cohort_root=cohort_root,
                seed=seed,
                timeout=timeout,
                cases=tuple(results),
            ),
        )

    report = _report(
        manifest_path=manifest_path,
        manifest_sha256=_sha256(manifest_path),
        cohort_root=cohort_root,
        seed=seed,
        timeout=timeout,
        cases=tuple(results),
    )
    if not dry_run:
        _atomic_json(output_root / "report.json", report)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohort_root", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--work-parent", required=True, type=Path)
    parser.add_argument("--dvbfixer", default="dvbfixer", type=Path)
    parser.add_argument("--seed", default=7, type=int)
    parser.add_argument("--timeout", default=DEFAULT_TIMEOUT_SECONDS, type=int)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--max-cases", type=int)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        manifest_path = args.manifest.absolute()
        if manifest_path == DEFAULT_MANIFEST and _sha256(manifest_path) != FROZEN_MANIFEST_SHA256:
            raise PublicCohortError("default frozen cohort manifest digest changed")
        executable_name = str(args.dvbfixer)
        resolved_executable = shutil.which(executable_name)
        executable = (
            Path(resolved_executable)
            if resolved_executable is not None
            else args.dvbfixer.absolute()
        )
        report = replay(
            manifest_path=manifest_path,
            cohort_root=args.cohort_root.absolute(),
            output_root=args.output_root.absolute(),
            work_parent=args.work_parent.absolute(),
            executable=executable,
            seed=args.seed,
            timeout=args.timeout,
            selected_case_ids=frozenset(args.case),
            max_cases=args.max_cases,
            dry_run=args.dry_run,
        )
    except PublicCohortError as exc:
        print(f"ERROR: {exc}")
        return 1
    print(json.dumps({key: report[key] for key in report if key != "cases"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
