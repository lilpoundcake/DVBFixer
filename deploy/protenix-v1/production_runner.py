#!/usr/bin/env python3
"""Protocol wrapper for the frozen Linux/NVIDIA Protenix v1 profile."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import platform
import shutil
import subprocess
import sys
from collections.abc import Callable
from dataclasses import replace
from functools import partial
from pathlib import Path
from types import ModuleType

from dvbfixer.model.diffusion.contract import DIFFUSION_SCHEMA_VERSION, RunnerResult
from dvbfixer.model.diffusion.preflight import (
    AdapterPreflightCode,
    AdapterPreflightIssue,
    RunnerPreflightReport,
    runtime_preflight_facts,
)
from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION
from dvbfixer.model.diffusion.runner_refinement import refine_runner_result
from dvbfixer.model.diffusion.sampler import SamplingAblationMode

PROFILE = "protenix-v1-cuda"
ENGINE_REVISION = "85767b811c40ed46e73a9b39519cf6bfca8701ba"
PATCH_SHA256 = "244cfe4fc876fd71df8737b805fb9b290069211fecc5178ef02029424b2895d0"
TORCH_VERSION = "2.13.0"
CUDA_VERSION = "12.9"
CHECKPOINT_SHA256 = "2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04"
PATCH_IDENTITY = f"sha256:{PATCH_SHA256}"
REFINEMENT_PLATFORM = "OpenMM Reference"


def _sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _load_script(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path.name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _source_root() -> Path:
    spec = importlib.util.find_spec("protenix")
    if spec is None or spec.origin is None:
        raise RuntimeError("Protenix is not installed in the runner environment")
    return Path(spec.origin).resolve().parent.parent


def _normalized_patch(text: str) -> str:
    return "".join(line for line in text.splitlines(keepends=True) if not line.startswith("index "))


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
    return set(changed) == {
        "protenix/model/generator.py",
        "protenix/model/protenix.py",
    } and not untracked_modules


def _verify_profile_environment() -> None:
    if sys.platform != "linux" or platform.machine() not in {"x86_64", "amd64"}:
        raise RuntimeError("protenix-v1-cuda requires Linux amd64")
    if sys.version_info[:2] != (3, 13):
        raise RuntimeError("protenix-v1-cuda requires Python 3.13")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        raise RuntimeError("PYTHONNOUSERSITE=1 is required")

    root = _source_root()
    revision = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if revision != ENGINE_REVISION:
        raise RuntimeError(f"Protenix revision mismatch: {revision}")
    if not _source_tree_is_frozen(root):
        raise RuntimeError("Protenix source tree contains changes outside the maintained patch")
    patch_path = Path(__file__).with_name("per-step-callback.patch")
    if _sha256(patch_path) != PATCH_SHA256:
        raise RuntimeError("Protenix callback patch digest mismatch")
    applied = subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "diff",
            "HEAD",
            "--",
            "protenix/model/generator.py",
            "protenix/model/protenix.py",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if _normalized_patch(applied) != _normalized_patch(patch_path.read_text()):
        raise RuntimeError("Protenix source does not match the maintained callback patch")

    import torch

    if torch.__version__.split("+")[0] != TORCH_VERSION:
        raise RuntimeError(f"protenix-v1-cuda requires PyTorch {TORCH_VERSION}")
    if str(torch.version.cuda or "") != CUDA_VERSION or not torch.cuda.is_available():
        raise RuntimeError(f"protenix-v1-cuda requires available CUDA {CUDA_VERSION}")
    if not torch.cuda.is_bf16_supported():
        raise RuntimeError("protenix-v1-cuda requires GPU bfloat16 support")


def preflight_report(profile: str, checkpoint: Path) -> RunnerPreflightReport:
    """Verify the frozen environment without importing adapters or loading weights."""
    issues: list[AdapterPreflightIssue] = []

    def issue(code: AdapterPreflightCode, message: str) -> None:
        issues.append(AdapterPreflightIssue(code, message))

    if profile != PROFILE:
        issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "runner does not support the requested profile")
    if sys.platform != "linux":
        issue(AdapterPreflightCode.INCOMPATIBLE_PLATFORM, "profile requires Linux")
    if platform.machine() not in {"x86_64", "amd64"}:
        issue(AdapterPreflightCode.INCOMPATIBLE_ARCHITECTURE, "profile requires amd64")
    if sys.version_info[:2] != (3, 13):
        issue(AdapterPreflightCode.INCOMPATIBLE_PYTHON, "profile requires Python 3.13")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "PYTHONNOUSERSITE must be enabled")
    if shutil.which("kalign") is None:
        issue(AdapterPreflightCode.MISSING_RESOURCE, "required kalign executable is unavailable")
    if importlib.util.find_spec("openmm") is None:
        issue(AdapterPreflightCode.MISSING_RESOURCE, "required OpenMM refinement runtime is unavailable")

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

    revision = ""
    try:
        root = _source_root()
        revision = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.strip()
        if revision != ENGINE_REVISION:
            issue(AdapterPreflightCode.INCOMPATIBLE_SOURCE, "Protenix revision does not match the frozen profile")
        patch_path = Path(__file__).with_name("per-step-callback.patch")
        if _sha256(patch_path) != PATCH_SHA256:
            issue(AdapterPreflightCode.INCOMPATIBLE_PATCH, "maintained callback patch digest is invalid")
        else:
            if not _source_tree_is_frozen(root):
                issue(AdapterPreflightCode.INCOMPATIBLE_SOURCE, "Protenix source tree contains unapproved changes")
            applied = subprocess.run(
                ["git", "-C", str(root), "diff", "HEAD", "--", "protenix/model/generator.py", "protenix/model/protenix.py"],
                check=True,
                capture_output=True,
                text=True,
                timeout=10,
            ).stdout
            if _normalized_patch(applied) != _normalized_patch(patch_path.read_text()):
                issue(AdapterPreflightCode.INCOMPATIBLE_PATCH, "Protenix patch state does not match the frozen profile")
    except (OSError, RuntimeError, subprocess.SubprocessError):
        issue(AdapterPreflightCode.INCOMPATIBLE_SOURCE, "Protenix source identity could not be verified")

    framework_version = ""
    accelerator_available = False
    effective_device = ""
    try:
        import torch

        framework_version = torch.__version__.split("+")[0]
        if framework_version != TORCH_VERSION or str(torch.version.cuda or "") != CUDA_VERSION:
            issue(AdapterPreflightCode.INCOMPATIBLE_FRAMEWORK, "PyTorch or CUDA version does not match the frozen profile")
        accelerator_available = bool(torch.cuda.is_available())
        if not accelerator_available:
            issue(AdapterPreflightCode.MISSING_CUDA, "required CUDA device is unavailable")
        else:
            effective_device = f"cuda:{torch.cuda.current_device()}"
            if not torch.cuda.is_bf16_supported():
                issue(AdapterPreflightCode.INCOMPATIBLE_DEVICE, "CUDA device lacks required bfloat16 support")
    except (ImportError, RuntimeError):
        issue(AdapterPreflightCode.INCOMPATIBLE_FRAMEWORK, "frozen PyTorch CUDA runtime is unavailable")

    return RunnerPreflightReport(
        profile=profile,
        runner_protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
        contract_schema_version=DIFFUSION_SCHEMA_VERSION,
        facts=runtime_preflight_facts(
            engine_revision=revision,
            patch_identity=PATCH_IDENTITY,
            environment_identity="linux-amd64;python=3.13;torch=2.13.0;cuda=12.9",
            checkpoint_sha256=checkpoint_sha256,
            framework_version=framework_version,
            accelerator_available=accelerator_available,
            effective_device=effective_device,
            fallback_disabled=True,
            sampling_platform="Linux amd64 / NVIDIA CUDA 12.9",
            refinement_platform=REFINEMENT_PLATFORM,
        ),
        issues=tuple(issues),
    )


def _load_production_scripts() -> tuple[ModuleType, ModuleType]:
    directory = Path(__file__).parent
    reinjection = _load_script("reinjection", directory / "reinjection.py")
    sys.modules["reinjection"] = reinjection
    builder = _load_script("dvbfixer_protenix_input_builder", directory / "build-template-input.py")
    adapter = _load_script("dvbfixer_protenix_checkpoint_gap_smoke", directory / "checkpoint_gap_smoke.py")
    return builder, adapter


def run(
    *,
    profile: str,
    checkpoint: Path,
    scripts: tuple[ModuleType, ModuleType] | None = None,
    refiner: Callable[[Path, RunnerResult, Path], RunnerResult] = partial(
        refine_runner_result,
        platform_name="Reference",
    ),
    preflight: Callable[[], None] = _verify_profile_environment,
) -> RunnerResult:
    if profile != PROFILE:
        raise ValueError(f"unsupported Protenix profile: {profile}")
    preflight()
    sys.dont_write_bytecode = True

    request_path = Path(os.environ.get("DVBFIXER_DIFFUSION_REQUEST", "request.json"))
    result_path = Path(os.environ.get("DVBFIXER_DIFFUSION_RESULT", "result.json"))
    checkpoint = checkpoint.expanduser().resolve()
    kalign = shutil.which("kalign")
    if kalign is None:
        raise RuntimeError("protenix-v1-cuda requires the kalign executable")
    root = _source_root() if scripts is None else Path.cwd()
    os.environ["PROTENIX_ROOT_DIR"] = str(root)
    os.environ["LAYERNORM_TYPE"] = "torch"

    builder, adapter = scripts or _load_production_scripts()
    output_root = Path("candidates") / PROFILE
    input_json = builder.build_input(request_path, output_root / "input")
    raw_result, _summary_path = adapter.run(
        request_path,
        input_json,
        output_root / "raw",
        Path(kalign),
        checkpoint,
        cycles=1,
        steps=200,
        ablation_mode=SamplingAblationMode.REINJECTION,
    )
    provenance = raw_result.backend_provenance
    if (
        not provenance.device.startswith("cuda")
        or provenance.framework_version.split("+")[0] != TORCH_VERSION
        or provenance.cuda_version != CUDA_VERSION
    ):
        raise RuntimeError("runner did not use the frozen PyTorch/CUDA profile")
    raw_result = replace(
        raw_result,
        backend_provenance=replace(
            provenance,
            source_license="Apache-2.0",
            checkpoint_license="Apache-2.0 (upstream claim)",
            environment_identity="linux-amd64;python=3.13;torch=2.13.0;cuda=12.9",
        ),
    )
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
