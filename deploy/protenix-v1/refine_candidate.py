"""Apply DVBFixer boundary refinement to a materialized Protenix candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from dvbfixer.model.diffusion.benchmark import candidate_quality, load_pdb_coordinates
from dvbfixer.model.diffusion.boundary_refinement import (
    BOUNDARY_REFINEMENT_REVISION,
    refine_generated_region,
)
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResourceMetrics,
    RunnerResult,
)
from dvbfixer.model.diffusion.validate import validate_runner_result


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _atom_identity(line: str) -> AtomIdentity:
    return AtomIdentity(
        line[21:22],
        line[22:26].strip(),
        line[26:27].strip(),
        line[12:16].strip(),
    )


def _rewrite_generated_coordinates(
    pdb_text: str,
    coordinates: dict[AtomIdentity, np.ndarray],
) -> str:
    expected = set(coordinates)
    seen: set[AtomIdentity] = set()
    output: list[str] = []
    for line in pdb_text.splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            output.append(line)
            continue
        identity = _atom_identity(line)
        if identity not in expected:
            output.append(line)
            continue
        if identity in seen:
            raise ValueError(f"candidate contains duplicate generated atom: {identity}")
        xyz = np.asarray(coordinates[identity], dtype=np.float64)
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise ValueError(f"refinement returned invalid coordinates for {identity}")
        fields = tuple(f"{float(value):8.3f}" for value in xyz)
        if any(len(field) != 8 for field in fields):
            raise ValueError(f"refined coordinate cannot be represented in PDB: {identity}")
        output.append(line[:30] + "".join(fields) + line[54:])
        seen.add(identity)
    missing = expected - seen
    if missing:
        raise ValueError(f"candidate omits {len(missing)} generated atoms")
    return "".join(output)


def _relative_artifact(workspace: Path, path: Path) -> ArtifactReference:
    try:
        relative = path.resolve().relative_to(workspace.resolve())
    except ValueError as exc:
        raise ValueError("refined output must be inside the request workspace") from exc
    return ArtifactReference(relative.as_posix(), _sha256(path))


def run(
    request_path: Path,
    candidate_path: Path,
    output_dir: Path,
    *,
    platform: str,
    restart_count: int,
    expected_candidate_sha256: str,
) -> None:
    request_path = request_path.resolve()
    workspace = request_path.parent
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if len(request.gaps) != 1 or request.candidate_count != 1:
        raise ValueError("Protenix refinement smoke requires one gap and candidate")
    if restart_count < 0:
        raise ValueError("restart count must be non-negative")

    candidate_path = candidate_path.resolve()
    candidate_sha256 = _sha256(candidate_path)
    if expected_candidate_sha256 and candidate_sha256 != expected_candidate_sha256:
        raise ValueError("raw candidate digest mismatch")
    candidate_text = candidate_path.read_text(encoding="ascii")
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)

    gap = request.gaps[0]
    refinement_kwargs: dict[str, object] = {}
    if len(gap.generated_residues) <= 5:
        refinement_kwargs.update(max_iterations=500, perturbation_angstrom=0.25)
    start = time.perf_counter()
    refinement = refine_generated_region(
        candidate_text,
        generated_residues=gap.generated_residues,
        generated_atoms=request.generated_atoms,
        restart_count=restart_count,
        random_seed=request.seeds[0],
        platform_name=platform,
        **refinement_kwargs,
    )
    refined_text = _rewrite_generated_coordinates(
        candidate_text,
        refinement.coordinates_angstrom,
    )
    output_path = output_dir / "candidate.pdb"
    output_path.write_text(refined_text, encoding="ascii")

    candidate = RunnerCandidate(
        candidate_id=f"protenix-v1-refined-seed-{request.seeds[0]}",
        seed=request.seeds[0],
        coordinate_artifact=_relative_artifact(workspace, output_path),
        generated_atoms=request.generated_atoms,
        generated_residues=gap.generated_residues,
        raw_backend_score=None,
        score_provenance=f"{BOUNDARY_REFINEMENT_REVISION}:unranked",
    )
    runner_result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=(candidate,),
        runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="protenix-v1-hook-spike+boundary-refinement",
            runner_protocol_version=3,
            engine_repository="https://github.com/bytedance/Protenix",
            engine_revision="85767b811c40ed46e73a9b39519cf6bfca8701ba",
            checkpoint_sha256=(
                "2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04"
            ),
            device="cpu",
            precision="float64",
            framework="OpenMM",
            deterministic_algorithms=True,
            deterministic_flags=(
                f"platform={refinement.platform}",
                f"random_seed={request.seeds[0]}",
                f"restart_count={restart_count}",
            ),
        ),
        resource_metrics=RunnerResourceMetrics(wall_time_seconds=time.perf_counter() - start),
    )
    validation = validate_runner_result(request, runner_result, workspace=workspace)[0]
    reference_path = workspace / "reference.pdb"
    if not reference_path.is_file():
        raise ValueError("refinement smoke workspace requires reference.pdb")
    reference = load_pdb_coordinates(
        workspace,
        ArtifactReference("reference.pdb", _sha256(reference_path)),
    )
    candidate_coordinates = load_pdb_coordinates(workspace, candidate.coordinate_artifact)
    quality = candidate_quality(
        candidate.candidate_id,
        reference,
        candidate_coordinates,
        generated_residues=set(gap.generated_residues),
        fixed_atoms=set(request.fixed_atoms),
        anchor_residues={gap.left_anchor, gap.right_anchor},
        closure_passed=(
            "junction-peptide-connectivity" not in validation.summary.hard_gate_failures
        ),
        validation_passed=validation.summary.passed,
    )
    summary = {
        "boundary_refinement_revision": BOUNDARY_REFINEMENT_REVISION,
        "raw_candidate_sha256": candidate_sha256,
        "candidate_sha256": candidate.coordinate_artifact.sha256,
        "pre_coordinate_sha256": refinement.pre_coordinate_sha256,
        "post_coordinate_sha256": refinement.post_coordinate_sha256,
        "initial_energy_kj_mol": refinement.initial_energy_kj_mol,
        "final_energy_kj_mol": refinement.final_energy_kj_mol,
        "chirality_repairs": [asdict(item) for item in refinement.chirality_repairs],
        "platform": refinement.platform,
        "restart_count": restart_count,
        "quality": {
            "gap_backbone_rmsd_angstrom": quality.gap_backbone_rmsd_angstrom,
            "gap_all_heavy_rmsd_angstrom": quality.gap_all_heavy_rmsd_angstrom,
            "fixed_heavy_rmsd_angstrom": quality.fixed_heavy_rmsd_angstrom,
        },
        "validation": asdict(validation.summary),
        "ranking_metrics": [asdict(metric) for metric in validation.ranking_metrics],
        "wall_time_seconds": runner_result.resource_metrics.wall_time_seconds,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--platform", default="Reference")
    parser.add_argument("--restart-count", default=8, type=int)
    parser.add_argument("--expected-candidate-sha256", default="")
    args = parser.parse_args()
    if args.expected_candidate_sha256 and (
        len(args.expected_candidate_sha256) != 64
        or any(character not in "0123456789abcdef" for character in args.expected_candidate_sha256)
    ):
        parser.error("--expected-candidate-sha256 must be 64 lowercase hexadecimal characters")
    run(
        args.request,
        args.candidate,
        args.output_dir,
        platform=args.platform,
        restart_count=args.restart_count,
        expected_candidate_sha256=args.expected_candidate_sha256,
    )


if __name__ == "__main__":
    main()
