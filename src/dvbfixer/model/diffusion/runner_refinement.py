"""Shared generated-region refinement for external diffusion runners."""

from __future__ import annotations

import hashlib
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from dvbfixer.model.diffusion.boundary_refinement import (
    BOUNDARY_REFINEMENT_MAX_ITERATIONS,
    BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM,
    BOUNDARY_REFINEMENT_RESTART_COUNT,
    BOUNDARY_REFINEMENT_REVISION,
    refine_generated_region,
)
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    RunnerCandidate,
    RunnerResourceMetrics,
    RunnerResult,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _atom_identity(line: str) -> AtomIdentity:
    return AtomIdentity(
        line[21:22],
        line[22:26].strip(),
        line[26:27].strip(),
        line[12:16].strip(),
    )


def rewrite_generated_coordinates(
    pdb_text: str,
    coordinates: dict[AtomIdentity, np.ndarray],
) -> str:
    """Replace only generated coordinate columns in a candidate PDB."""
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
        fields = tuple(f"{float(value):8.3f}" for value in xyz)
        if xyz.shape != (3,) or not np.isfinite(xyz).all() or any(
            len(field) != 8 for field in fields
        ):
            raise ValueError(f"invalid refined coordinate for {identity}")
        output.append(line[:30] + "".join(fields) + line[54:])
        seen.add(identity)
    if seen != expected:
        raise ValueError(f"candidate omits {len(expected - seen)} generated atoms")
    return "".join(output)


def refine_runner_result(
    request_path: Path,
    raw_result: RunnerResult,
    output_dir: Path,
    *,
    platform_name: str,
) -> RunnerResult:
    """Refine one runner candidate while preserving its sampling provenance."""
    request_path = request_path.resolve()
    workspace = request_path.parent
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if len(raw_result.candidates) != 1 or len(request.gaps) != 1:
        raise ValueError("production refinement requires one candidate and one gap")
    raw_candidate = raw_result.candidates[0]
    raw_path = workspace.joinpath(*Path(raw_candidate.coordinate_artifact.path).parts)
    raw_text = raw_path.read_text(encoding="ascii")
    gap = request.gaps[0]
    short_gap = len(gap.generated_residues) <= 5

    start = time.perf_counter()
    refinement = refine_generated_region(
        raw_text,
        generated_residues=gap.generated_residues,
        generated_atoms=request.generated_atoms,
        max_iterations=500 if short_gap else BOUNDARY_REFINEMENT_MAX_ITERATIONS,
        restart_count=BOUNDARY_REFINEMENT_RESTART_COUNT,
        perturbation_angstrom=(
            0.25 if short_gap else BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM
        ),
        platform_name=platform_name,
        random_seed=request.seeds[0],
    )
    output_dir.mkdir(parents=True, exist_ok=False)
    refined_path = output_dir / "candidate.pdb"
    refined_path.write_text(
        rewrite_generated_coordinates(raw_text, refinement.coordinates_angstrom),
        encoding="ascii",
    )
    relative_path = refined_path.resolve().relative_to(workspace).as_posix()
    candidate = RunnerCandidate(
        candidate_id=f"{raw_candidate.candidate_id}-refined",
        seed=raw_candidate.seed,
        coordinate_artifact=ArtifactReference(relative_path, _sha256(refined_path)),
        generated_atoms=raw_candidate.generated_atoms,
        generated_residues=raw_candidate.generated_residues,
        raw_backend_score=raw_candidate.raw_backend_score,
        score_provenance=f"{BOUNDARY_REFINEMENT_REVISION}:unranked",
        warnings=raw_candidate.warnings,
    )
    raw_metrics = raw_result.resource_metrics
    elapsed = time.perf_counter() - start
    platform_label = refinement.platform
    return replace(
        raw_result,
        candidates=(candidate,),
        backend_provenance=replace(
            raw_result.backend_provenance,
            backend=(
                f"{raw_result.backend_provenance.backend}"
                f"+{platform_label.lower()}-boundary-refinement"
            ),
            deterministic_flags=(
                *raw_result.backend_provenance.deterministic_flags,
                f"refinement-platform={platform_label}",
                f"refinement={BOUNDARY_REFINEMENT_REVISION}",
                f"refinement-restart-count={BOUNDARY_REFINEMENT_RESTART_COUNT}",
            ),
            known_nondeterministic_operations=(
                *raw_result.backend_provenance.known_nondeterministic_operations,
                f"OpenMM {platform_label} minimization",
            ),
        ),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=(raw_metrics.wall_time_seconds or 0.0) + elapsed,
            model_load_seconds=raw_metrics.model_load_seconds,
            peak_ram_bytes=raw_metrics.peak_ram_bytes,
            peak_vram_bytes=raw_metrics.peak_vram_bytes,
        ),
    )
