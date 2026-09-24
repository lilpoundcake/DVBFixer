#!/usr/bin/env python3
"""Compare MODELLER outputs against a withheld-coordinate benchmark reference."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from pathlib import Path

from dvbfixer.model.diffusion.benchmark import candidate_quality, load_pdb_coordinates
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResult,
)
from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION
from dvbfixer.model.diffusion.validate import build_validated_result


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(workspace: Path, models: tuple[Path, ...]) -> dict[str, object]:
    """Return per-model and median MODELLER comparison metrics."""
    workspace = workspace.expanduser().resolve()
    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    reference = load_pdb_coordinates(
        workspace,
        ArtifactReference("reference.pdb", _sha256(workspace / "reference.pdb")),
    )
    rows: list[dict[str, object]] = []
    for index, raw_path in enumerate(models, start=1):
        path = raw_path.expanduser().resolve()
        relative = path.relative_to(workspace).as_posix()
        runner_candidate = RunnerCandidate(
            candidate_id=f"modeller-{index:04d}",
            seed=index,
            coordinate_artifact=ArtifactReference(relative, _sha256(path)),
            generated_atoms=request.generated_atoms,
            generated_residues=request.gaps[0].generated_residues,
            raw_backend_score=None,
            score_provenance="modeller-10.8:molpdf-not-cross-backend-comparable",
        )
        runner = RunnerResult(
            schema_version=DIFFUSION_SCHEMA_VERSION,
            status=DiffusionStatus.SUCCESS,
            candidates=(runner_candidate,),
            runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
            backend_provenance=BackendProvenance(
                backend="modeller-comparator",
                runner_protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
                engine_repository="https://salilab.org/modeller/",
                engine_revision="10.8",
                source_license="MODELLER academic license",
            ),
        )
        validated = build_validated_result(request, runner, workspace=workspace)
        summary = validated.validation_summaries[0]
        coordinates = load_pdb_coordinates(workspace, runner_candidate.coordinate_artifact)
        quality = candidate_quality(
            runner_candidate.candidate_id,
            reference,
            coordinates,
            generated_residues=set(request.gaps[0].generated_residues),
            fixed_atoms=set(request.fixed_atoms),
            anchor_residues={request.gaps[0].left_anchor, request.gaps[0].right_anchor},
            closure_passed=(
                "junction-peptide-connectivity" not in summary.hard_gate_failures
            ),
            validation_passed=summary.passed,
        )
        rows.append(
            {
                "path": relative,
                "validation_passed": summary.passed,
                "hard_gate_failures": summary.hard_gate_failures,
                "junction_passed": quality.closure_passed,
                "gap_backbone_rmsd_angstrom": quality.gap_backbone_rmsd_angstrom,
                "gap_all_heavy_rmsd_angstrom": quality.gap_all_heavy_rmsd_angstrom,
                "fixed_heavy_rmsd_angstrom": quality.fixed_heavy_rmsd_angstrom,
            }
        )
    return {
        "candidate_count": len(rows),
        "median_gap_backbone_rmsd_angstrom": statistics.median(
            float(row["gap_backbone_rmsd_angstrom"]) for row in rows
        ),
        "median_gap_all_heavy_rmsd_angstrom": statistics.median(
            float(row["gap_all_heavy_rmsd_angstrom"]) for row in rows
        ),
        "median_fixed_heavy_rmsd_angstrom": statistics.median(
            float(row["fixed_heavy_rmsd_angstrom"]) for row in rows
        ),
        "junction_pass_rate": sum(bool(row["junction_passed"]) for row in rows)
        / len(rows),
        "validation_pass_count": sum(bool(row["validation_passed"]) for row in rows),
        "candidates": rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("workspace", type=Path)
    parser.add_argument("models", type=Path, nargs="+")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(json.dumps(analyze(args.workspace, tuple(args.models)), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
