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

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    DiffusionStatus,
    RunnerResult,
)
from dvbfixer.model.diffusion.preflight import (
    AdapterPreflightCode,
    AdapterPreflightIssue,
    RunnerPreflightReport,
    runtime_preflight_facts,
)
from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION
from dvbfixer.model.diffusion.runner_refinement import refine_runner_result

PROFILE = "protpardelle-1c-mps"
ENGINE_REVISION = "ee378400f25b801fa481028000f9060183d7fb4c"
APPLE_PATCH_SHA256 = "627891e28d5055cb0d903f542af695569d7133b0480cab8f25d154dc8ccc78c9"
PATCHED_MODELS_SHA256 = "8513fac0d18d3080b14673bb9b6ed3b231e1d31d30449d3b8c6cc09f0445800b"
CHECKPOINT_SHA256 = "dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483"
CONFIG_SHA256 = "e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d"
PATCH_IDENTITY = f"sha256:{APPLE_PATCH_SHA256}"
REFINEMENT_PLATFORM = "OpenMM CPU"
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
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _source_tree_is_frozen(root: Path) -> bool:
    changed = subprocess.run(
        ["git", "-C", str(root), "diff", "--name-only", "HEAD", "--"],
        check=True,
        capture_output=True,
        text=True,
        timeout=10,
    ).stdout.splitlines()
    untracked_modules = subprocess.run(
        [
            "git", "-C", str(root), "ls-files", "--others", "--",
            "*.py", "*.pyc", "*.so", "*.pyd", "*.dylib",
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=10,
    ).stdout.splitlines()
    return set(changed) == {"src/protpardelle/core/models.py"} and not untracked_modules


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
    if not _source_tree_is_frozen(source_root):
        raise RuntimeError("Protpardelle source tree contains changes outside the frozen patch")
    models_path = source_root / "src/protpardelle/core/models.py"
    if _sha256(models_path) != PATCHED_MODELS_SHA256:
        raise RuntimeError("Protpardelle source does not match the frozen Apple patch state")


def preflight_report(profile: str, checkpoint: Path) -> RunnerPreflightReport:
    """Verify the frozen environment without importing adapters or loading weights."""
    issues: list[AdapterPreflightIssue] = []

    def issue(code: AdapterPreflightCode, message: str) -> None:
        issues.append(AdapterPreflightIssue(code, message))

    if profile != PROFILE:
        issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "runner does not support the requested profile")
    if sys.platform != "darwin":
        issue(AdapterPreflightCode.INCOMPATIBLE_PLATFORM, "profile requires macOS")
    if platform.machine() != "arm64":
        issue(AdapterPreflightCode.INCOMPATIBLE_ARCHITECTURE, "profile requires arm64")
    if sys.version_info[:2] != (3, 12):
        issue(AdapterPreflightCode.INCOMPATIBLE_PYTHON, "profile requires Python 3.12")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "PYTHONNOUSERSITE must be enabled")
    configured = [name for name in _UNSET_MPS_VARIABLES if os.environ.get(name) is not None]
    if configured:
        issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "frozen MPS tuning variables must be unset")
    if importlib.util.find_spec("openmm") is None:
        issue(AdapterPreflightCode.MISSING_RESOURCE, "required OpenMM refinement runtime is unavailable")
    fallback_disabled = os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK") == "0"
    if not fallback_disabled:
        issue(AdapterPreflightCode.FALLBACK_ENABLED, "MPS fallback must be explicitly disabled")

    checkpoint_sha256 = ""
    if checkpoint.is_symlink() or not checkpoint.is_file():
        issue(AdapterPreflightCode.MISSING_CHECKPOINT, "required checkpoint is unavailable")
    else:
        try:
            checkpoint_sha256 = _sha256(checkpoint)
        except OSError:
            issue(AdapterPreflightCode.MISSING_CHECKPOINT, "required checkpoint could not be read")
        if checkpoint_sha256 and checkpoint_sha256 != CHECKPOINT_SHA256:
            issue(AdapterPreflightCode.CHECKPOINT_DIGEST_MISMATCH, "checkpoint digest does not match the frozen profile")
        try:
            config = _config_for_checkpoint(checkpoint)
            if _sha256(config) != CONFIG_SHA256:
                issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "Protpardelle config does not match the frozen profile")
        except (FileNotFoundError, OSError):
            issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "frozen Protpardelle config is unavailable")

    revision = ""
    try:
        patch_path = Path(__file__).with_name("apple-portability.patch")
        if _sha256(patch_path) != APPLE_PATCH_SHA256:
            issue(AdapterPreflightCode.INCOMPATIBLE_PATCH, "Apple portability patch digest is invalid")
        spec = importlib.util.find_spec("protpardelle")
        if spec is None or spec.origin is None:
            raise RuntimeError("source unavailable")
        source_root = Path(spec.origin).resolve().parents[2]
        revision = subprocess.run(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.strip()
        if revision != ENGINE_REVISION:
            issue(AdapterPreflightCode.INCOMPATIBLE_SOURCE, "Protpardelle revision does not match the frozen profile")
        if not _source_tree_is_frozen(source_root):
            issue(AdapterPreflightCode.INCOMPATIBLE_SOURCE, "Protpardelle source tree contains unapproved changes")
        if _sha256(source_root / "src/protpardelle/core/models.py") != PATCHED_MODELS_SHA256:
            issue(AdapterPreflightCode.INCOMPATIBLE_PATCH, "Protpardelle patch state does not match the frozen profile")
    except (OSError, RuntimeError, subprocess.SubprocessError):
        issue(AdapterPreflightCode.INCOMPATIBLE_SOURCE, "Protpardelle source identity could not be verified")

    framework_version = ""
    accelerator_available = False
    effective_device = ""
    try:
        import torch

        framework_version = torch.__version__.split("+")[0]
        if framework_version != "2.6.0":
            issue(AdapterPreflightCode.INCOMPATIBLE_FRAMEWORK, "PyTorch version does not match the frozen profile")
        accelerator_available = bool(torch.backends.mps.is_available())
        if not accelerator_available:
            issue(AdapterPreflightCode.MISSING_MPS, "required MPS device is unavailable")
        else:
            effective_device = str(torch.device("mps"))
    except (ImportError, RuntimeError):
        issue(AdapterPreflightCode.INCOMPATIBLE_FRAMEWORK, "frozen PyTorch MPS runtime is unavailable")

    return RunnerPreflightReport(
        profile=profile,
        runner_protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
        contract_schema_version=DIFFUSION_SCHEMA_VERSION,
        facts=runtime_preflight_facts(
            engine_revision=revision,
            patch_identity=PATCH_IDENTITY,
            environment_identity="macos-arm64;python=3.12;torch=2.6.0",
            checkpoint_sha256=checkpoint_sha256,
            framework_version=framework_version,
            accelerator_available=accelerator_available,
            effective_device=effective_device,
            fallback_disabled=fallback_disabled,
            sampling_platform="native macOS arm64 / MPS",
            refinement_platform=REFINEMENT_PLATFORM,
        ),
        issues=tuple(issues),
    )


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
    if os.environ.get("DVBFIXER_DIFFUSION_PROTOCOL_VERSION") != str(
        DIFFUSION_RUNNER_PROTOCOL_VERSION
    ):
        raise RuntimeError("diffusion execution protocol version is incompatible")
    preflight()
    sys.dont_write_bytecode = True

    request_path = Path(os.environ.get("DVBFIXER_DIFFUSION_REQUEST", "request.json"))
    result_path = Path(os.environ.get("DVBFIXER_DIFFUSION_RESULT", "result.json"))
    checkpoint = checkpoint.expanduser().resolve()
    config = _config_for_checkpoint(checkpoint)
    os.environ["PROTPARDELLE_MODEL_PARAMS"] = str(checkpoint.parent.parent)
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
    raw_result = replace(
        raw_result,
        backend_provenance=replace(
            raw_result.backend_provenance,
            source_license="MIT",
            environment_identity="macos-arm64;python=3.12;torch=2.6.0",
            known_nondeterministic_operations=("PyTorch MPS sampling",),
        ),
    )
    if raw_result.status is not DiffusionStatus.SUCCESS:
        result_path.write_text(raw_result.to_json(), encoding="utf-8")
        return raw_result
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
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--profile", required=True)
    parser.add_argument("--checkpoint", required=True, type=Path)
    args = parser.parse_args()
    if args.preflight:
        print(preflight_report(args.profile, args.checkpoint).to_json(), end="")
        return
    run(profile=args.profile, checkpoint=args.checkpoint)


if __name__ == "__main__":
    main()
