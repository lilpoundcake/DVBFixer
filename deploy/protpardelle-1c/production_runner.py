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
from collections.abc import Callable
from dataclasses import replace
from functools import partial
from pathlib import Path
from types import ModuleType

from dvbfixer.model.diffusion.contract import RunnerResult
from dvbfixer.model.diffusion.runner_refinement import refine_runner_result

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


def run(
    *,
    profile: str,
    checkpoint: Path,
    adapter: ModuleType | None = None,
    refiner: Callable[[Path, RunnerResult, Path], RunnerResult] = partial(
        refine_runner_result,
        platform_name="CPU",
    ),
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
    raw_result = replace(
        raw_result,
        backend_provenance=replace(
            raw_result.backend_provenance,
            source_license="MIT",
            environment_identity="macos-arm64;python=3.12;torch=2.6.0",
            known_nondeterministic_operations=("PyTorch MPS sampling",),
        ),
    )
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
