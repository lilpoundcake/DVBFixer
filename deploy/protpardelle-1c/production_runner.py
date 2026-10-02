#!/usr/bin/env python3
"""Protocol wrapper for the frozen Apple Protpardelle-1c profile."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import platform
import subprocess
import sys
import time
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from types import ModuleType

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

PROFILE = "protpardelle-1c-mps"
ENGINE_REVISION = "ee378400f25b801fa481028000f9060183d7fb4c"
APPLE_PATCH_SHA256 = "a87fe2e9f0c143102441d6a39ff181bd0c2b1c411cdfdf65236d7baf38a8d858"
PATCHED_MODELS_SHA256 = "3ad9efdc4e1086e14dbba941d88ca62521e956f13a1df88bba5fc6edec81c19d"
_UNSET_MPS_VARIABLES = (
    "PYTORCH_MPS_FAST_MATH",
    "PYTORCH_MPS_PREFER_METAL",
    "PYTORCH_MPS_HIGH_WATERMARK_RATIO",
    "PYTORCH_MPS_LOW_WATERMARK_RATIO",
)


def _load_adapter() -> ModuleType:
    path = Path(__file__).with_name("checkpoint_gap_smoke.py")
    spec = importlib.util.spec_from_file_location(
        "dvbfixer_protpardelle_checkpoint_gap_smoke",
        path,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load the Protpardelle adapter")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _config_for_checkpoint(checkpoint: Path) -> Path:
    config = checkpoint.parent.parent / "configs" / "cc89.yaml"
    if not config.is_file():
        raise FileNotFoundError(
            "cc89 config was not found beside the checkpoint; expected "
            f"{config}"
        )
    return config


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _verify_profile_environment() -> None:
    if sys.platform != "darwin" or platform.machine() != "arm64":
        raise RuntimeError("protpardelle-1c-mps requires native macOS arm64")
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError("protpardelle-1c-mps requires Python 3.12")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        raise RuntimeError("PYTHONNOUSERSITE=1 is required")
    configured = [name for name in _UNSET_MPS_VARIABLES if os.environ.get(name) is not None]
    if configured:
        raise RuntimeError("frozen MPS profile requires unset variables: " + ", ".join(configured))

    patch_path = Path(__file__).with_name("apple-portability.patch")
    if _sha256(patch_path) != APPLE_PATCH_SHA256:
        raise RuntimeError("Apple portability patch digest mismatch")
    spec = importlib.util.find_spec("protpardelle")
    if spec is None or spec.origin is None:
        raise RuntimeError("protpardelle is not installed in the runner environment")
    source_root = Path(spec.origin).resolve().parents[2]
    revision = subprocess.run(
        ["git", "-C", str(source_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if revision != ENGINE_REVISION:
        raise RuntimeError(f"Protpardelle revision mismatch: {revision}")
    models_path = source_root / "src/protpardelle/core/models.py"
    if _sha256(models_path) != PATCHED_MODELS_SHA256:
        raise RuntimeError("Protpardelle source does not match the frozen Apple patch state")


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


def _refine_result(
    request_path: Path,
    raw_result: RunnerResult,
    output_dir: Path,
) -> RunnerResult:
    request_path = request_path.resolve()
    workspace = request_path.parent
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if len(raw_result.candidates) != 1 or len(request.gaps) != 1:
        raise ValueError("Apple production refinement requires one candidate and one gap")
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
        platform_name="CPU",
        random_seed=request.seeds[0],
    )
    output_dir.mkdir(parents=True, exist_ok=False)
    refined_path = output_dir / "candidate.pdb"
    refined_path.write_text(
        _rewrite_generated_coordinates(raw_text, refinement.coordinates_angstrom),
        encoding="ascii",
    )
    relative_path = refined_path.resolve().relative_to(workspace).as_posix()
    candidate = RunnerCandidate(
        candidate_id=f"protpardelle-1c-cc89-refined-seed-{raw_candidate.seed}",
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
    return replace(
        raw_result,
        candidates=(candidate,),
        backend_provenance=replace(
            raw_result.backend_provenance,
            backend=f"{raw_result.backend_provenance.backend}+cpu-boundary-refinement",
            deterministic_flags=(
                *raw_result.backend_provenance.deterministic_flags,
                "refinement-platform=CPU",
                f"refinement={BOUNDARY_REFINEMENT_REVISION}",
                f"refinement-restart-count={BOUNDARY_REFINEMENT_RESTART_COUNT}",
            ),
            known_nondeterministic_operations=(
                *raw_result.backend_provenance.known_nondeterministic_operations,
                "OpenMM CPU minimization",
            ),
        ),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=(raw_metrics.wall_time_seconds or 0.0) + elapsed,
            model_load_seconds=raw_metrics.model_load_seconds,
            peak_ram_bytes=raw_metrics.peak_ram_bytes,
            peak_vram_bytes=raw_metrics.peak_vram_bytes,
        ),
    )


def run(
    *,
    profile: str,
    checkpoint: Path,
    adapter: ModuleType | None = None,
    refiner: Callable[[Path, RunnerResult, Path], RunnerResult] = _refine_result,
    preflight: Callable[[], None] = _verify_profile_environment,
) -> RunnerResult:
    if profile != PROFILE:
        raise ValueError(f"unsupported Protpardelle profile: {profile}")
    if os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK") != "0":
        raise RuntimeError("PYTORCH_ENABLE_MPS_FALLBACK=0 is required")
    preflight()

    request_path = Path(os.environ.get("DVBFIXER_DIFFUSION_REQUEST", "request.json"))
    result_path = Path(os.environ.get("DVBFIXER_DIFFUSION_RESULT", "result.json"))
    checkpoint = checkpoint.expanduser().resolve()
    config = _config_for_checkpoint(checkpoint)
    output_root = Path("candidates") / PROFILE
    raw_output_dir = output_root / "raw"
    raw_output_dir.mkdir(parents=True, exist_ok=False)

    active_adapter = adapter or _load_adapter()
    raw_result, _summary_path = active_adapter._run_staged(
        request_path,
        raw_output_dir,
        config,
        checkpoint,
        steps=500,
        step_scale=1.2,
        s_churn=200.0,
        per_step_reinjection=False,
        device_name="mps",
        mps_profile=False,
    )
    if (
        raw_result.backend_provenance.device != "mps"
        or raw_result.backend_provenance.framework_version != "2.6.0"
    ):
        raise RuntimeError("runner did not use the frozen PyTorch 2.6.0 MPS profile")
    result = refiner(request_path, raw_result, output_root / "refined")
    result_path.write_text(result.to_json(), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", required=True)
    parser.add_argument("--checkpoint", required=True, type=Path)
    args = parser.parse_args()
    run(profile=args.profile, checkpoint=args.checkpoint)


if __name__ == "__main__":
    main()
