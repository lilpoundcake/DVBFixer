#!/usr/bin/env python3
"""Run and aggregate the preregistered Protpardelle v5 pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import subprocess
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dvbfixer.model.diffusion.contract import DiffusionRequest

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "docs/research/small-diffusion-v5-expanded-pilot.json"
ADAPTER = ROOT / "deploy/protpardelle-1c/checkpoint_gap_smoke.py"
CALLBACK_PATCH = ROOT / "deploy/protpardelle-1c/per-step-callback.patch"
REFINER = ROOT / "deploy/protenix-v1/refine_candidate.py"
BOUNDARY_REFINEMENT = ROOT / "src/dvbfixer/model/diffusion/boundary_refinement.py"
PROTPARDELLE_REVISION = "ee378400f25b801fa481028000f9060183d7fb4c"
ENGINE_REPOSITORY = "https://github.com/ProteinDesignLab/protpardelle-1c"


class PilotError(RuntimeError):
    """Raised when pilot provenance or a stage is invalid."""


@dataclass(frozen=True, slots=True)
class PilotCase:
    case_id: str
    screening_index: int
    gap_length: int
    target_length: int


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_manifest(path: Path, cohort_root: Path) -> tuple[dict[str, Any], tuple[PilotCase, ...]]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema_version") != 1:
        raise PilotError("unsupported pilot manifest schema")
    expected_digest = value.get("source_aggregate_sha256")
    aggregate = cohort_root / "aggregate.json"
    if not aggregate.is_file() or _sha256(aggregate) != expected_digest:
        raise PilotError("source aggregate digest mismatch")
    cases = tuple(PilotCase(**item) for item in value.get("cases", []))
    expected_count = int(value.get("case_count", 24))
    if len(cases) != expected_count or len({case.case_id for case in cases}) != len(cases):
        raise PilotError(f"pilot manifest must contain {expected_count} unique cases")
    expected_strata = value.get("gap_strata", {"5": 12, "10": 12})
    actual_strata = Counter(str(case.gap_length) for case in cases)
    if actual_strata != Counter({str(key): int(count) for key, count in expected_strata.items()}):
        raise PilotError("pilot manifest gap strata differ from the frozen counts")
    if any(case.target_length > 512 for case in cases):
        raise PilotError("pilot manifest contains a target above the cc89 limit")
    return value, cases


def _verify_request(case: PilotCase, workspace: Path) -> Path:
    request_path = workspace / "request.json"
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if len(request.target_sequences) != 1 or len(request.gaps) != 1:
        raise PilotError(f"{case.case_id}: request is outside the single-gap scope")
    if len(request.target_sequences[0].sequence) != case.target_length:
        raise PilotError(f"{case.case_id}: target length differs from manifest")
    gap = request.gaps[0]
    if gap.target_interval.stop - gap.target_interval.start != case.gap_length:
        raise PilotError(f"{case.case_id}: gap length differs from manifest")
    return request_path


def _case_output_root(workspace: Path, output_subdir: str) -> Path:
    relative = Path(output_subdir)
    if not output_subdir or relative.is_absolute() or len(relative.parts) != 1:
        raise PilotError("output subdirectory must be one relative path component")
    return workspace / relative


def _completion_valid(
    marker: Path,
    inputs: dict[str, Path],
    outputs: dict[str, Path],
) -> bool:
    if not marker.is_file() or any(not path.is_file() for path in (*inputs.values(), *outputs.values())):
        return False
    try:
        value = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return value.get("status") == "complete" and value.get("inputs") == {
        name: _sha256(path) for name, path in inputs.items()
    } and value.get("outputs") == {name: _sha256(path) for name, path in outputs.items()}


def _run_stage(
    marker: Path,
    command: tuple[str, ...],
    inputs: dict[str, Path],
    outputs: dict[str, Path],
    *,
    env: dict[str, str],
) -> bool:
    if _completion_valid(marker, inputs, outputs):
        return True
    if any(path.exists() for path in outputs.values()):
        raise PilotError(f"invalid unmarked stage output blocks resume: {marker.stem}")
    start = time.perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True, env=env, check=False)
    record: dict[str, Any] = {
        "command": list(command),
        "exit_status": completed.returncode,
        "inputs": {name: _sha256(path) for name, path in inputs.items()},
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "wall_time_seconds": time.perf_counter() - start,
    }
    if completed.returncode or any(not path.is_file() for path in outputs.values()):
        record["status"] = "failed"
        _write_json(marker.with_name(f"{marker.stem}.failure.json"), record)
        return False
    record["status"] = "complete"
    record["outputs"] = {name: _sha256(path) for name, path in outputs.items()}
    _write_json(marker, record)
    return True


def run_pilot(
    manifest_path: Path,
    cohort_root: Path,
    protpardelle_python: Path,
    protpardelle_checkout: Path,
    config: Path,
    checkpoint: Path,
    dvbfixer_python: Path,
    *,
    selected_case_ids: frozenset[str] = frozenset(),
    output_subdir: str = "small-diffusion-expanded",
    per_step_reinjection: bool = False,
) -> None:
    manifest, cases = load_manifest(manifest_path, cohort_root)
    sampling = manifest["sampling"]
    if _sha256(config) != sampling["config_sha256"]:
        raise PilotError("Protpardelle config digest mismatch")
    if _sha256(checkpoint) != sampling["checkpoint_sha256"]:
        raise PilotError("Protpardelle checkpoint digest mismatch")
    revision = subprocess.run(
        ("git", "-C", str(protpardelle_checkout), "rev-parse", "HEAD"),
        capture_output=True,
        text=True,
        check=False,
    )
    if revision.returncode or revision.stdout.strip() != PROTPARDELLE_REVISION:
        raise PilotError("Protpardelle source revision mismatch")
    if per_step_reinjection:
        patch_check = subprocess.run(
            (
                "git",
                "-C",
                str(protpardelle_checkout),
                "apply",
                "--reverse",
                "--check",
                str(CALLBACK_PATCH),
            ),
            capture_output=True,
            text=True,
            check=False,
        )
        if patch_check.returncode:
            raise PilotError("Protpardelle per-step callback patch is not applied")
    unknown = selected_case_ids - {case.case_id for case in cases}
    if unknown:
        raise PilotError(f"requested cases are absent from manifest: {sorted(unknown)}")

    for case in cases:
        if selected_case_ids and case.case_id not in selected_case_ids:
            continue
        workspace = cohort_root / "cases" / case.case_id / "workspace"
        request = _verify_request(case, workspace)
        output_root = _case_output_root(workspace, output_subdir)
        raw = output_root / "protpardelle-raw"
        refinement = manifest["refinement"]
        refined = output_root / str(
            refinement.get("output_subdir", "protpardelle-refined-cpu")
        )
        markers = output_root / ".stages"
        raw_command = (
            str(protpardelle_python),
            str(ADAPTER),
            str(request),
            str(raw),
            "--config",
            str(config),
            "--checkpoint",
            str(checkpoint),
            "--steps",
            str(sampling["steps"]),
            "--step-scale",
            str(sampling["step_scale"]),
            "--s-churn",
            str(sampling["s_churn"]),
        ) + (("--per-step-reinjection",) if per_step_reinjection else ())
        raw_env = os.environ.copy()
        raw_env["PYTHONPATH"] = os.pathsep.join(
            (str(ROOT / "src"), str(protpardelle_checkout / "src"))
        )
        raw_inputs = {
            "manifest": manifest_path,
            "request": request,
            "config": config,
            "checkpoint": checkpoint,
            "adapter": ADAPTER,
        }
        if per_step_reinjection:
            raw_inputs["callback_patch"] = CALLBACK_PATCH
        if not _run_stage(
            markers / "raw.json",
            raw_command,
            raw_inputs,
            {"candidate": raw / "candidate.pdb", "summary": raw / "summary.json"},
            env=raw_env,
        ):
            raise PilotError(f"{case.case_id}: raw Protpardelle stage failed")

        raw_summary = json.loads((raw / "summary.json").read_text(encoding="utf-8"))
        raw_digest = _sha256(raw / "candidate.pdb")
        if raw_summary.get("candidate_sha256") != raw_digest:
            raise PilotError(f"{case.case_id}: raw summary candidate digest mismatch")
        if bool(raw_summary.get("per_step_reinjection", False)) != per_step_reinjection:
            raise PilotError(f"{case.case_id}: raw summary reinjection mode mismatch")
        if per_step_reinjection:
            callback = raw_summary.get("per_step_callback", {})
            if callback.get("callback_count") != sampling["steps"]:
                raise PilotError(f"{case.case_id}: raw callback count mismatch")
            if callback.get("maximum_post_projection_error_angstrom") != 0.0:
                raise PilotError(f"{case.case_id}: raw callback projection was not exact")
        backend = "protpardelle-1c-cc89"
        if per_step_reinjection:
            backend += "-per-step-reinjection"
        refine_command = (
            str(dvbfixer_python),
            str(REFINER),
            str(request),
            str(raw / "candidate.pdb"),
            str(refined),
            "--platform",
            str(refinement["platform"]),
            "--restart-count",
            str(refinement["restart_count"]),
            "--expected-candidate-sha256",
            raw_digest,
            "--backend",
            backend,
            "--engine-repository",
            ENGINE_REPOSITORY,
            "--engine-revision",
            PROTPARDELLE_REVISION,
            "--checkpoint-sha256",
            str(sampling["checkpoint_sha256"]),
        )
        refine_env = os.environ.copy()
        if "openmm_cpu_threads" in refinement:
            refine_env["OPENMM_CPU_THREADS"] = str(refinement["openmm_cpu_threads"])
        refinement_inputs = {
            "manifest": manifest_path,
            "request": request,
            "raw_candidate": raw / "candidate.pdb",
            "refiner": REFINER,
            "boundary_refinement": BOUNDARY_REFINEMENT,
        }
        if not _run_stage(
            markers / "refinement.json",
            refine_command,
            refinement_inputs,
            {
                "candidate": refined / "candidate.pdb",
                "summary": refined / "summary.json",
            },
            env=refine_env,
        ):
            raise PilotError(f"{case.case_id}: boundary refinement stage failed")
        print(case.case_id)


def _median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def _selection_source(
    raw: dict[str, Any] | None,
    refined: dict[str, Any] | None,
) -> str | None:
    if refined is not None and refined["validation"]["passed"]:
        return "refined"
    if raw is not None and raw["validation"]["passed"]:
        return "raw"
    return None


def aggregate(
    manifest_path: Path,
    cohort_root: Path,
    output: Path,
    *,
    output_subdir: str = "small-diffusion-expanded",
) -> dict[str, Any]:
    manifest, cases = load_manifest(manifest_path, cohort_root)
    existing = json.loads((cohort_root / "aggregate.json").read_text(encoding="utf-8"))
    comparator_maps = {
        backend: {item["case_id"]: item for item in items}
        for backend, items in existing["case_results"].items()
    }
    rows: list[dict[str, Any]] = []
    for case in cases:
        workspace = cohort_root / "cases" / case.case_id / "workspace"
        root = _case_output_root(workspace, output_subdir)
        raw_path = root / "protpardelle-raw/summary.json"
        refined_path = root / str(
            manifest["refinement"].get("output_subdir", "protpardelle-refined-cpu")
        ) / "summary.json"
        raw_stage_path = root / ".stages/raw.json"
        refinement_stage_path = root / ".stages/refinement.json"
        raw = json.loads(raw_path.read_text(encoding="utf-8")) if raw_path.is_file() else None
        refined = (
            json.loads(refined_path.read_text(encoding="utf-8"))
            if refined_path.is_file()
            else None
        )
        selection_source = _selection_source(raw, refined)
        if raw is not None and raw.get("per_step_reinjection"):
            callback = raw.get("per_step_callback", {})
            if callback.get("callback_count") != raw.get("diffusion_steps"):
                raise PilotError(f"{case.case_id}: aggregate callback count mismatch")
            if callback.get("maximum_post_projection_error_angstrom") != 0.0:
                raise PilotError(f"{case.case_id}: aggregate callback projection was not exact")
        raw_stage = (
            json.loads(raw_stage_path.read_text(encoding="utf-8"))
            if raw_stage_path.is_file()
            else None
        )
        refinement_stage = (
            json.loads(refinement_stage_path.read_text(encoding="utf-8"))
            if refinement_stage_path.is_file()
            else None
        )
        rows.append(
            {
                "case_id": case.case_id,
                "screening_index": case.screening_index,
                "gap_length": case.gap_length,
                "target_length": case.target_length,
                "raw": raw,
                "refined": refined,
                "selection_source": selection_source,
                "operational": {
                    "raw_stage_wall_time_seconds": (
                        raw_stage["wall_time_seconds"] if raw_stage is not None else None
                    ),
                    "refinement_stage_wall_time_seconds": (
                        refinement_stage["wall_time_seconds"]
                        if refinement_stage is not None
                        else None
                    ),
                    "total_wall_time_seconds": (
                        raw_stage["wall_time_seconds"] + refinement_stage["wall_time_seconds"]
                        if raw_stage is not None and refinement_stage is not None
                        else None
                    ),
                },
                "protenix": comparator_maps["protenix-v1"][case.case_id],
                "boltz": comparator_maps["boltz-2-proxy-only"][case.case_id],
                "modeller": comparator_maps["modeller-10.8"][case.case_id],
            }
        )

    def candidate_item(row: dict[str, Any], key: str) -> dict[str, Any] | None:
        if key != "selected":
            return row[key]
        source = row["selection_source"]
        return row[source] if source is not None else None

    def summarize(key: str, gap_length: int | None = None) -> dict[str, Any]:
        selected = [row for row in rows if gap_length is None or row["gap_length"] == gap_length]
        available = [item for row in selected if (item := candidate_item(row, key)) is not None]
        result: dict[str, Any] = {
            "case_count": len(selected),
            "completed_count": len(available),
            "passed_count": sum(bool(item["validation"]["passed"]) for item in available),
            "hard_gate_failure_counts": dict(
                sorted(
                    Counter(
                        failure
                        for item in available
                        for failure in item["validation"]["hard_gate_failures"]
                    ).items()
                )
            ),
            "median_gap_backbone_rmsd_angstrom": _median(
                [float(item["quality"]["gap_backbone_rmsd_angstrom"]) for item in available]
            ),
            "median_gap_all_heavy_rmsd_angstrom": _median(
                [float(item["quality"]["gap_all_heavy_rmsd_angstrom"]) for item in available]
            ),
        }
        if key != "selected":
            result["median_wall_time_seconds"] = _median(
                [float(item["wall_time_seconds"]) for item in available]
            )
        if key == "raw":
            result.update(
                median_native_conditioning_rmsd_angstrom=_median(
                    [
                        float(item["native_conditioning_fit"]["rmsd_angstrom"])
                        for item in available
                    ]
                ),
                median_native_conditioning_max_displacement_angstrom=_median(
                    [
                        float(item["native_conditioning_fit"]["max_displacement_angstrom"])
                        for item in available
                    ]
                ),
                peak_vram_bytes=max(int(item["peak_vram_bytes"]) for item in available),
            )
        return result

    def summarize_comparator(key: str, gap_length: int | None = None) -> dict[str, Any]:
        selected = [row for row in rows if gap_length is None or row["gap_length"] == gap_length]
        return {
            "case_count": len(selected),
            "passed_count": sum(bool(row[key]["validation_passed"]) for row in selected),
            "median_gap_backbone_rmsd_angstrom": _median(
                [float(row[key]["gap_backbone_rmsd_angstrom"]) for row in selected]
            ),
            "median_wall_time_seconds": _median(
                [float(row[key]["wall_time_seconds"]) for row in selected]
            ),
        }

    def paired_candidate_comparison(candidate_key: str, key: str) -> dict[str, Any]:
        def candidate_passed(row: dict[str, Any]) -> bool:
            candidate = candidate_item(row, candidate_key)
            return candidate is not None and bool(candidate["validation"]["passed"])

        def candidate_backbone_rmsd(row: dict[str, Any]) -> float:
            candidate = candidate_item(row, candidate_key)
            if candidate is None:
                raise PilotError("passing candidate is unexpectedly absent")
            return float(candidate["quality"]["gap_backbone_rmsd_angstrom"])

        both_passing = [
            row
            for row in rows
            if candidate_passed(row) and row[key]["validation_passed"]
        ]
        return {
            "both_passed_count": len(both_passing),
            "candidate_only_passed_count": sum(
                candidate_passed(row) and not bool(row[key]["validation_passed"])
                for row in rows
            ),
            "comparator_only_passed_count": sum(
                not candidate_passed(row) and bool(row[key]["validation_passed"])
                for row in rows
            ),
            "neither_passed_count": sum(
                not candidate_passed(row) and not bool(row[key]["validation_passed"])
                for row in rows
            ),
            "median_paired_backbone_rmsd_difference_angstrom": _median(
                [
                    candidate_backbone_rmsd(row)
                    - float(row[key]["gap_backbone_rmsd_angstrom"])
                    for row in both_passing
                ]
            ),
        }

    result = {
        "schema_version": 1,
        "output_subdir": output_subdir,
        "manifest_sha256": _sha256(manifest_path),
        "source_aggregate_sha256": manifest["source_aggregate_sha256"],
        "summary": {
            key: {
                "all": summarize(key),
                "gap_5": summarize(key, 5),
                "gap_10": summarize(key, 10),
            }
            for key in ("raw", "refined", "selected")
        },
        "comparators": {
            key: {
                "all": summarize_comparator(key),
                "gap_5": summarize_comparator(key, 5),
                "gap_10": summarize_comparator(key, 10),
            }
            for key in ("protenix", "boltz", "modeller")
        },
        "paired_refined_comparisons": {
            key: paired_candidate_comparison("refined", key)
            for key in ("protenix", "boltz", "modeller")
        },
        "paired_selected_comparisons": {
            key: paired_candidate_comparison("selected", key)
            for key in ("protenix", "boltz", "modeller")
        },
        "selection": {
            "refined_count": sum(row["selection_source"] == "refined" for row in rows),
            "raw_count": sum(row["selection_source"] == "raw" for row in rows),
            "failed_count": sum(row["selection_source"] is None for row in rows),
            "failed_case_ids": [
                row["case_id"] for row in rows if row["selection_source"] is None
            ],
        },
        "operations": {
            "completed_case_count": sum(
                row["operational"]["total_wall_time_seconds"] is not None for row in rows
            ),
            "median_total_wall_time_seconds": _median(
                [
                    float(row["operational"]["total_wall_time_seconds"])
                    for row in rows
                    if row["operational"]["total_wall_time_seconds"] is not None
                ]
            ),
            "maximum_total_wall_time_seconds": max(
                (
                    float(row["operational"]["total_wall_time_seconds"])
                    for row in rows
                    if row["operational"]["total_wall_time_seconds"] is not None
                ),
                default=None,
            ),
        },
        "cases": rows,
    }
    _write_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    run_parser.add_argument("--cohort-root", required=True, type=Path)
    run_parser.add_argument("--protpardelle-python", required=True, type=Path)
    run_parser.add_argument("--protpardelle-checkout", required=True, type=Path)
    run_parser.add_argument("--config", required=True, type=Path)
    run_parser.add_argument("--checkpoint", required=True, type=Path)
    run_parser.add_argument("--dvbfixer-python", required=True, type=Path)
    run_parser.add_argument("--case", action="append", default=[])
    run_parser.add_argument("--output-subdir", default="small-diffusion-expanded")
    run_parser.add_argument("--per-step-reinjection", action="store_true")
    aggregate_parser = subparsers.add_parser("aggregate")
    aggregate_parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    aggregate_parser.add_argument("--cohort-root", required=True, type=Path)
    aggregate_parser.add_argument("--output", required=True, type=Path)
    aggregate_parser.add_argument("--output-subdir", default="small-diffusion-expanded")
    args = parser.parse_args()
    if args.command == "run":
        run_pilot(
            args.manifest,
            args.cohort_root,
            args.protpardelle_python,
            args.protpardelle_checkout,
            args.config,
            args.checkpoint,
            args.dvbfixer_python,
            selected_case_ids=frozenset(args.case),
            output_subdir=args.output_subdir,
            per_step_reinjection=args.per_step_reinjection,
        )
    else:
        aggregate(
            args.manifest,
            args.cohort_root,
            args.output,
            output_subdir=args.output_subdir,
        )


if __name__ == "__main__":
    main()
