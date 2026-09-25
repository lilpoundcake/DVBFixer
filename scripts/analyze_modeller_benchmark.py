#!/usr/bin/env python3
"""Compare MODELLER outputs against a withheld-coordinate benchmark reference."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import tempfile
from pathlib import Path

from dvbfixer.model.cli import AA3TO1
from dvbfixer.model.diffusion.benchmark import candidate_quality, load_pdb_coordinates
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    ResidueIdentity,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResult,
)
from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION
from dvbfixer.model.diffusion.validate import build_validated_result


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _normalize_target_identities(request: DiffusionRequest, pdb_text: str) -> str:
    """Map a complete MODELLER target chain back to request residue identities."""
    if len(request.target_sequences) != 1 or len(request.sequence_placements) != 1:
        raise ValueError("MODELLER benchmark normalization requires one target chain")
    if len(request.gaps) != 1:
        raise ValueError("MODELLER benchmark normalization requires one gap")

    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    gap = request.gaps[0]
    identities_by_index = dict(
        zip(placement.observed_target_indices, placement.observed_residues)
    )
    identities_by_index.update(
        zip(
            range(gap.target_interval.start, gap.target_interval.stop),
            gap.generated_residues,
        )
    )
    if set(identities_by_index) != set(range(len(target.sequence))):
        raise ValueError("diffusion request does not map every target sequence ordinal")

    ordered_keys: list[tuple[str, str]] = []
    residue_names: dict[tuple[str, str], str] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 27 or line[21] != target.chain:
            continue
        key = (line[22:26].strip(), line[26].strip())
        if key not in residue_names:
            ordered_keys.append(key)
            residue_names[key] = line[17:20].strip()
    if len(ordered_keys) != len(target.sequence):
        raise ValueError(
            "MODELLER target residue count does not match the diffusion request"
        )

    remapping: dict[tuple[str, str], ResidueIdentity] = {}
    for index, key in enumerate(ordered_keys):
        residue_name = residue_names[key]
        if AA3TO1.get(residue_name) != target.sequence[index]:
            raise ValueError(
                f"MODELLER target sequence mismatch at ordinal {index + 1}"
            )
        remapping[key] = identities_by_index[index]

    output: list[str] = []
    for line in pdb_text.splitlines(keepends=True):
        if line.startswith("ATOM  ") and len(line) >= 27 and line[21] == target.chain:
            key = (line[22:26].strip(), line[26].strip())
            identity = remapping[key]
            line = (
                line[:21]
                + identity.chain
                + f"{int(identity.residue_number):4d}"
                + (identity.insertion_code or " ")
                + line[27:]
            )
        output.append(line)
    return "".join(output)


def analyze(workspace: Path, models: tuple[Path, ...]) -> dict[str, object]:
    """Return per-model and median MODELLER comparison metrics."""
    workspace = workspace.expanduser().resolve()
    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    reference = load_pdb_coordinates(
        workspace,
        ArtifactReference("reference.pdb", _sha256(workspace / "reference.pdb")),
    )
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix=".modeller-benchmark-", dir=workspace) as temp:
        temporary_root = Path(temp)
        for index, raw_path in enumerate(models, start=1):
            path = raw_path.expanduser().resolve()
            source_relative = path.relative_to(workspace).as_posix()
            source_text = path.read_text(encoding="utf-8")
            normalized_text = _normalize_target_identities(request, source_text)
            normalized_path = temporary_root / f"candidate-{index:04d}.pdb"
            normalized_path.write_text(normalized_text, encoding="utf-8")
            normalized_relative = normalized_path.relative_to(workspace).as_posix()
            runner_candidate = RunnerCandidate(
                candidate_id=f"modeller-{index:04d}",
                seed=index,
                coordinate_artifact=ArtifactReference(
                    normalized_relative,
                    _sha256(normalized_path),
                ),
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
            coordinates = load_pdb_coordinates(
                workspace,
                runner_candidate.coordinate_artifact,
            )
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
                    "path": source_relative,
                    "identity_normalized": normalized_text != source_text,
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
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rendered = json.dumps(analyze(args.workspace, tuple(args.models)), indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
