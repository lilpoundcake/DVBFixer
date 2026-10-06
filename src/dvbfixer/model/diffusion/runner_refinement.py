"""Shared generated-region refinement for external diffusion runners."""

from __future__ import annotations

import hashlib
import os
import stat
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from dvbfixer.model.diffusion.boundary_refinement import (
    BOUNDARY_REFINEMENT_MAX_ITERATIONS,
    BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM,
    BOUNDARY_REFINEMENT_RESTART_COUNT,
    BOUNDARY_REFINEMENT_REVISION,
    BoundaryRefinementError,
    refine_generated_region,
)
from dvbfixer.model.diffusion.contract import (
    MAX_TRACE_BYTES,
    ArtifactReference,
    AtomIdentity,
    BackendOption,
    DiffusionRequest,
    ResidueIdentity,
    RunnerCandidate,
    RunnerResourceMetrics,
    RunnerResult,
    SamplerTrace,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_private_artifact(
    workspace: Path,
    artifact: ArtifactReference,
    *,
    label: str,
    max_bytes: int,
    encoding: str,
) -> str:
    path = workspace.joinpath(*Path(artifact.path).parts)
    current = workspace
    for part in Path(artifact.path).parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f"raw {label} path contains a symlink")
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise ValueError(f"raw {label} does not exist") from exc
    if not resolved.is_relative_to(workspace):
        raise ValueError(f"raw {label} escapes the request workspace")
    file_stat = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(file_stat.st_mode) or file_stat.st_nlink != 1:
        raise ValueError(f"raw {label} is not a private regular file")
    if file_stat.st_size > max_bytes:
        raise ValueError(f"raw {label} exceeds the size limit")
    data = path.read_bytes()
    if len(data) > max_bytes:
        raise ValueError(f"raw {label} exceeds the size limit")
    if hashlib.sha256(data).hexdigest() != artifact.sha256:
        raise ValueError(f"raw {label} digest mismatch")
    try:
        return data.decode(encoding)
    except UnicodeDecodeError as exc:
        raise ValueError(f"raw {label} is not valid {encoding}") from exc


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
    if len(raw_result.candidates) != 1:
        raise ValueError("production refinement requires one candidate")
    raw_candidate = raw_result.candidates[0]
    raw_text = _read_private_artifact(
        workspace,
        raw_candidate.coordinate_artifact,
        label="coordinate artifact",
        max_bytes=1_000_000_000,
        encoding="ascii",
    )
    raw_trace = SamplerTrace.from_json(
        _read_private_artifact(
            workspace,
            raw_candidate.sampler_trace_artifact,
            label="sampler trace",
            max_bytes=MAX_TRACE_BYTES,
            encoding="utf-8",
        )
    )
    movable_residues = tuple(dict.fromkeys(
        residue for gap in request.gaps for residue in gap.movable_junction_residues
    ))
    short_gaps = all(len(gap.generated_residues) <= 5 for gap in request.gaps)

    start = time.perf_counter()
    refined_text = raw_text
    refinement_platform = platform_name
    try:
        for gap_index, gap in enumerate(request.gaps):
            gap_residues = set(gap.generated_residues)
            gap_atoms = tuple(
                atom
                for atom in request.generated_atoms
                if ResidueIdentity(
                    atom.chain,
                    atom.residue_number,
                    atom.insertion_code,
                ) in gap_residues
            )
            short_gap = len(gap.generated_residues) <= 5
            refinement_input = (
                _protein_only_refinement_text(refined_text)
                if request.heterogen_contexts
                else refined_text
            )
            refinement = refine_generated_region(
                refinement_input,
                generated_residues=gap.generated_residues,
                generated_atoms=gap_atoms,
                max_iterations=(
                    500 if short_gap else BOUNDARY_REFINEMENT_MAX_ITERATIONS
                ),
                restart_count=BOUNDARY_REFINEMENT_RESTART_COUNT,
                perturbation_angstrom=(
                    0.25 if short_gap else BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM
                ),
                platform_name=platform_name,
                random_seed=request.seeds[0] + gap_index,
            )
            refinement_platform = refinement.platform
            refined_text = rewrite_generated_coordinates(
                refined_text,
                refinement.coordinates_angstrom,
            )
    except BoundaryRefinementError:
        elapsed = time.perf_counter() - start
        rejected = replace(
            raw_candidate,
            warnings=(
                *raw_candidate.warnings,
                "Localized OpenMM boundary refinement failed; the raw sampler "
                "candidate was retained for inspection.",
            ),
            postprocessing_failures=(
                *raw_candidate.postprocessing_failures,
                "localized-openmm-boundary-refinement-failed",
            ),
        )
        raw_metrics = raw_result.resource_metrics
        return replace(
            raw_result,
            candidates=(rejected,),
            resource_metrics=RunnerResourceMetrics(
                wall_time_seconds=(raw_metrics.wall_time_seconds or 0.0) + elapsed,
                model_load_seconds=raw_metrics.model_load_seconds,
                peak_ram_bytes=raw_metrics.peak_ram_bytes,
                peak_vram_bytes=raw_metrics.peak_vram_bytes,
            ),
            message=(
                "localized OpenMM boundary refinement failed; raw sampler "
                "candidate retained"
            ),
        )
    output_dir.mkdir(parents=True, exist_ok=False)
    refined_path = output_dir / "candidate.pdb"
    refined_path.write_text(refined_text, encoding="ascii")
    trace_path = output_dir / "sampler-trace.json"
    trace = replace(
        raw_trace,
        refinement_mode="localized-openmm-boundary-refinement",
        refinement_parameters=(
            BackendOption("platform", refinement_platform),
            BackendOption("gap_count", str(len(request.gaps))),
            BackendOption(
                "max_iterations",
                str(500 if short_gaps else BOUNDARY_REFINEMENT_MAX_ITERATIONS),
            ),
            BackendOption("restart_count", str(BOUNDARY_REFINEMENT_RESTART_COUNT)),
            BackendOption(
                "perturbation_angstrom",
                str(0.25 if short_gaps else BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM),
            ),
        ),
        localized_refinement_residues=movable_residues,
        final_heavy_coordinate_operations=(
            *raw_trace.final_heavy_coordinate_operations,
            "localized-openmm-boundary-refinement",
        ),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=(raw_trace.resource_metrics.wall_time_seconds or 0.0)
            + (time.perf_counter() - start),
            model_load_seconds=raw_trace.resource_metrics.model_load_seconds,
            peak_ram_bytes=raw_trace.resource_metrics.peak_ram_bytes,
            peak_vram_bytes=raw_trace.resource_metrics.peak_vram_bytes,
        ),
    )
    trace_path.write_text(trace.to_json(), encoding="utf-8")
    relative_path = refined_path.resolve().relative_to(workspace).as_posix()
    candidate = RunnerCandidate(
        candidate_id=f"{raw_candidate.candidate_id}-refined",
        seed=raw_candidate.seed,
        coordinate_artifact=ArtifactReference(relative_path, _sha256(refined_path)),
        sampler_trace_artifact=ArtifactReference(
            trace_path.resolve().relative_to(workspace).as_posix(),
            _sha256(trace_path),
        ),
        generated_atoms=raw_candidate.generated_atoms,
        generated_residues=raw_candidate.generated_residues,
        raw_backend_score=raw_candidate.raw_backend_score,
        score_provenance=f"{BOUNDARY_REFINEMENT_REVISION}:unranked",
        warnings=raw_candidate.warnings,
        postprocessing_failures=raw_candidate.postprocessing_failures,
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
                raw_result.backend_provenance.known_nondeterministic_operations
                if platform_label == "Reference"
                else (
                    *raw_result.backend_provenance.known_nondeterministic_operations,
                    f"OpenMM {platform_label} minimization",
                )
            ),
        ),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=(raw_metrics.wall_time_seconds or 0.0) + elapsed,
            model_load_seconds=raw_metrics.model_load_seconds,
            peak_ram_bytes=raw_metrics.peak_ram_bytes,
            peak_vram_bytes=raw_metrics.peak_vram_bytes,
        ),
    )


def _protein_only_refinement_text(pdb_text: str) -> str:
    """Remove retained chemistry from the temporary protein FF refinement input."""
    return "".join(
        line
        for line in pdb_text.splitlines(keepends=True)
        if not line.startswith(("HETATM", "CONECT", "LINK  "))
    )
