#!/usr/bin/env python3
"""Protocol wrapper for the frozen Linux/NVIDIA Protenix v1 profile."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
from collections.abc import Callable
from dataclasses import replace
from functools import partial
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from types import ModuleType

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    DiffusionRequest,
    DiffusionStatus,
    RunnerResourceMetrics,
    RunnerResult,
    SamplerTrace,
)
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
PATCH_SHA256 = "cc4153be3dfd241124ea183d592884799300046b6ea7d3eccae8409b9fe21aa0"
TORCH_VERSION = "2.13.0"
CUDA_VERSION = "12.9"
CHECKPOINT_SHA256 = "2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04"
CHECKPOINT_NAME = "protenix_base_default_v1.0.0.pt"
PYTHON_VERSION = (3, 13, 15)
GPU_NAME = "NVIDIA A100-SXM4-40GB"
GPU_COMPUTE_CAPABILITY = (8, 0)
MAX_RAM_BYTES = 150_323_855_360
PROFILE_LOCK_SHA256 = "3fb3661be2b6748650b90bf9b0aaa7a09cbba9accc13253a6c29f644001b5d2e"
PATCH_IDENTITY = f"sha256:{PATCH_SHA256}"
REFINEMENT_PLATFORM = "OpenMM Reference"


def _sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _profile_lock() -> dict[str, object]:
    path = Path(__file__).with_name("profile-lock.json")
    if _sha256(path) != PROFILE_LOCK_SHA256:
        raise RuntimeError("Protenix profile lock digest mismatch")
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise RuntimeError("Protenix profile lock must be a JSON object")
    return raw


def _verify_locked_packages(lock: dict[str, object]) -> None:
    packages = lock.get("packages")
    if not isinstance(packages, dict):
        raise RuntimeError("Protenix profile lock omits package pins")
    for distribution, expected in packages.items():
        if not isinstance(distribution, str) or not isinstance(expected, str):
            raise RuntimeError("Protenix profile lock contains an invalid package pin")
        try:
            observed = version(distribution)
        except PackageNotFoundError as exc:
            raise RuntimeError(f"required package is unavailable: {distribution}") from exc
        if observed != expected:
            raise RuntimeError(
                f"package version mismatch for {distribution}: {observed} != {expected}"
            )


def _verify_kalign(path: Path, lock: dict[str, object]) -> Path:
    identity = lock.get("kalign")
    if not isinstance(identity, dict):
        raise RuntimeError("Protenix profile lock omits the kalign identity")
    expected_digest = identity.get("sha256")
    expected_version = identity.get("version")
    if not isinstance(expected_digest, str) or not isinstance(expected_version, str):
        raise RuntimeError("Protenix profile lock contains an invalid kalign identity")
    resolved = path.resolve()
    if not resolved.is_file() or _sha256(resolved) != expected_digest:
        raise RuntimeError("kalign executable digest does not match the frozen profile")
    observed_version = subprocess.run(
        [str(resolved), "--version"],
        check=True,
        capture_output=True,
        text=True,
        timeout=10,
    ).stdout.strip()
    if observed_version != expected_version:
        raise RuntimeError(
            f"kalign version mismatch: {observed_version!r} != {expected_version!r}"
        )
    return resolved


def _verify_checkpoint(checkpoint: Path) -> Path:
    expanded = checkpoint.expanduser()
    if expanded.is_symlink() or not expanded.is_file():
        raise RuntimeError("Protenix checkpoint must be a regular, non-symlink file")
    resolved = expanded.resolve()
    if resolved.name != CHECKPOINT_NAME:
        raise RuntimeError(f"Protenix checkpoint must be named {CHECKPOINT_NAME}")
    if _sha256(resolved) != CHECKPOINT_SHA256:
        raise RuntimeError("Protenix checkpoint digest mismatch")
    return resolved


def _cgroup_ram_limit() -> int | None:
    for path in (
        Path("/sys/fs/cgroup/memory.max"),
        Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    ):
        try:
            value = path.read_text(encoding="ascii").strip()
        except OSError:
            continue
        if value != "max":
            try:
                return int(value)
            except ValueError:
                return None
    return None


def _verify_ram_bound() -> None:
    limit = _cgroup_ram_limit()
    if limit is None or limit <= 0 or limit > MAX_RAM_BYTES:
        raise RuntimeError("protenix-v1-cuda requires a bounded cgroup RAM limit")


def _workspace_path(environment_name: str, default: str) -> Path:
    value = Path(os.environ.get(environment_name, default))
    if value.is_absolute() or ".." in value.parts:
        raise RuntimeError(f"{environment_name} must be a workspace-relative path")
    return value


def _write_result_atomic(path: Path, result: RunnerResult) -> None:
    if path.exists():
        raise RuntimeError("diffusion result path already exists")
    temporary = path.with_name(f".{path.name}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write(result.to_json())
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _container_digest() -> str:
    digest = os.environ.get("DVBFIXER_PROTENIX_IMAGE_DIGEST", "")
    payload = digest.removeprefix("sha256:")
    if digest and (
        not digest.startswith("sha256:")
        or len(payload) != 64
        or any(character not in "0123456789abcdef" for character in payload)
    ):
        raise RuntimeError("DVBFIXER_PROTENIX_IMAGE_DIGEST is not a SHA-256 digest")
    return digest


def _verify_artifact(workspace: Path, relative_path: str, expected_sha256: str) -> Path:
    path = workspace.joinpath(*Path(relative_path).parts)
    current = workspace
    for part in Path(relative_path).parts:
        current /= part
        if current.is_symlink():
            raise RuntimeError("Protenix result artifact path contains a symlink")
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise RuntimeError("Protenix result artifact is missing") from exc
    if not resolved.is_relative_to(workspace) or path.is_symlink():
        raise RuntimeError("Protenix result artifact escapes the private workspace")
    file_stat = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(file_stat.st_mode) or file_stat.st_nlink != 1:
        raise RuntimeError("Protenix result artifact is not a private regular file")
    if _sha256(path) != expected_sha256:
        raise RuntimeError("Protenix result artifact digest mismatch")
    return path


def _read_trace(workspace: Path, result: RunnerResult) -> SamplerTrace:
    if len(result.candidates) != 1:
        raise RuntimeError("successful Protenix result must contain exactly one candidate")
    artifact = result.candidates[0].sampler_trace_artifact
    path = _verify_artifact(workspace, artifact.path, artifact.sha256)
    return SamplerTrace.from_json(path.read_text(encoding="utf-8"))


def _verify_success_result(result: RunnerResult, workspace: Path, *, refined: bool) -> None:
    if result.status is not DiffusionStatus.SUCCESS:
        return
    provenance = result.backend_provenance
    if (
        provenance.engine_revision != ENGINE_REVISION
        or provenance.checkpoint_sha256 != CHECKPOINT_SHA256
        or not provenance.device.startswith("cuda:")
        or provenance.precision != "bf16"
        or provenance.framework != "torch"
        or provenance.framework_version.split("+")[0] != TORCH_VERSION
        or provenance.cuda_version != CUDA_VERSION
        or provenance.deterministic_algorithms is not True
        or provenance.known_nondeterministic_operations
    ):
        raise RuntimeError("runner did not prove the frozen Protenix execution profile")
    metrics = result.resource_metrics
    if not metrics.wall_time_seconds or not metrics.peak_ram_bytes or not metrics.peak_vram_bytes:
        raise RuntimeError("runner did not provide complete nonzero resource evidence")
    if len(result.candidates) != 1 or result.candidates[0].postprocessing_failures:
        raise RuntimeError("successful Protenix result contains postprocessing failures")
    coordinate = result.candidates[0].coordinate_artifact
    _verify_artifact(workspace, coordinate.path, coordinate.sha256)
    trace = _read_trace(workspace, result)
    if (
        trace.profile != PROFILE
        or trace.device != provenance.device
        or not trace.fallback_disabled
    ):
        raise RuntimeError("runner trace does not prove CUDA execution with fallback disabled")
    expected_refinement = "localized-openmm-boundary-refinement" if refined else "none"
    if trace.refinement_mode != expected_refinement:
        raise RuntimeError("runner trace does not prove the expected refinement state")


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


def _verify_patched_sources(root: Path, lock: dict[str, object]) -> None:
    identities = lock.get("patched_source_sha256")
    if not isinstance(identities, dict) or not identities:
        raise RuntimeError("Protenix profile lock omits patched source identities")
    for relative_path, expected_digest in identities.items():
        if not isinstance(relative_path, str) or not isinstance(expected_digest, str):
            raise RuntimeError("Protenix profile lock contains an invalid source identity")
        path = root.joinpath(*Path(relative_path).parts)
        if path.is_symlink() or not path.is_file() or _sha256(path) != expected_digest:
            raise RuntimeError(f"patched Protenix source mismatch: {relative_path}")


def _verify_auxiliary_artifacts(root: Path, lock: dict[str, object]) -> None:
    identities = lock.get("auxiliary_artifact_sha256")
    if not isinstance(identities, dict) or not identities:
        raise RuntimeError("Protenix profile lock omits auxiliary artifact identities")
    for relative_path, expected_digest in identities.items():
        if not isinstance(relative_path, str) or not isinstance(expected_digest, str):
            raise RuntimeError("Protenix profile lock contains an invalid artifact identity")
        path = root.joinpath(*Path(relative_path).parts)
        if path.is_symlink() or not path.is_file() or _sha256(path) != expected_digest:
            raise RuntimeError(f"Protenix auxiliary artifact mismatch: {relative_path}")


def _verify_profile_environment() -> None:
    if sys.platform != "linux" or platform.machine() not in {"x86_64", "amd64"}:
        raise RuntimeError("protenix-v1-cuda requires Linux amd64")
    if sys.version_info[:3] != PYTHON_VERSION:
        raise RuntimeError("protenix-v1-cuda requires Python 3.13.15")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        raise RuntimeError("PYTHONNOUSERSITE=1 is required")

    lock = _profile_lock()
    _verify_locked_packages(lock)
    kalign = shutil.which("kalign")
    if kalign is None:
        raise RuntimeError("protenix-v1-cuda requires the kalign executable")
    _verify_kalign(Path(kalign), lock)
    _verify_ram_bound()
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
    _verify_patched_sources(root, lock)
    _verify_auxiliary_artifacts(root, lock)
    patch_path = Path(__file__).with_name("per-step-callback.patch")
    if _sha256(patch_path) != PATCH_SHA256:
        raise RuntimeError("Protenix callback patch digest mismatch")

    import torch

    if torch.__version__.split("+")[0] != TORCH_VERSION:
        raise RuntimeError(f"protenix-v1-cuda requires PyTorch {TORCH_VERSION}")
    if str(torch.version.cuda or "") != CUDA_VERSION or not torch.cuda.is_available():
        raise RuntimeError(f"protenix-v1-cuda requires available CUDA {CUDA_VERSION}")
    if not torch.cuda.is_bf16_supported():
        raise RuntimeError("protenix-v1-cuda requires GPU bfloat16 support")
    device = torch.cuda.current_device()
    if (
        torch.cuda.get_device_name(device) != GPU_NAME
        or torch.cuda.get_device_capability(device) != GPU_COMPUTE_CAPABILITY
    ):
        raise RuntimeError("protenix-v1-cuda requires the frozen NVIDIA A100 device class")


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
    if sys.version_info[:3] != PYTHON_VERSION:
        issue(AdapterPreflightCode.INCOMPATIBLE_PYTHON, "profile requires Python 3.13.15")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "PYTHONNOUSERSITE must be enabled")
    kalign = shutil.which("kalign")
    if kalign is None:
        issue(AdapterPreflightCode.MISSING_RESOURCE, "required kalign executable is unavailable")
    if importlib.util.find_spec("openmm") is None:
        issue(AdapterPreflightCode.MISSING_RESOURCE, "required OpenMM refinement runtime is unavailable")
    try:
        lock = _profile_lock()
        _verify_locked_packages(lock)
        if kalign is not None:
            _verify_kalign(Path(kalign), lock)
    except (
        OSError,
        RuntimeError,
        ValueError,
        json.JSONDecodeError,
        subprocess.SubprocessError,
    ):
        issue(
            AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT,
            "profile lock or package versions do not match",
        )
    try:
        _verify_ram_bound()
    except RuntimeError:
        issue(
            AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT,
            "runner does not have the frozen cgroup RAM bound",
        )

    checkpoint_sha256 = ""
    if checkpoint.is_symlink() or not checkpoint.is_file():
        issue(AdapterPreflightCode.MISSING_CHECKPOINT, "required checkpoint is unavailable")
    else:
        if checkpoint.name != CHECKPOINT_NAME:
            issue(
                AdapterPreflightCode.MISSING_CHECKPOINT,
                f"checkpoint must be named {CHECKPOINT_NAME}",
            )
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
            try:
                _verify_patched_sources(root, _profile_lock())
            except RuntimeError:
                issue(
                    AdapterPreflightCode.INCOMPATIBLE_PATCH,
                    "Protenix patch state does not match the frozen profile",
                )
            try:
                _verify_auxiliary_artifacts(root, _profile_lock())
            except RuntimeError:
                issue(
                    AdapterPreflightCode.MISSING_RESOURCE,
                    "Protenix auxiliary artifacts do not match the frozen profile",
                )
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
            device = torch.cuda.current_device()
            if (
                torch.cuda.get_device_name(device) != GPU_NAME
                or torch.cuda.get_device_capability(device) != GPU_COMPUTE_CAPABILITY
            ):
                issue(
                    AdapterPreflightCode.INCOMPATIBLE_DEVICE,
                    "CUDA device is not the frozen NVIDIA A100 class",
                )
    except (ImportError, RuntimeError):
        issue(AdapterPreflightCode.INCOMPATIBLE_FRAMEWORK, "frozen PyTorch CUDA runtime is unavailable")

    return RunnerPreflightReport(
        profile=profile,
        runner_protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
        contract_schema_version=DIFFUSION_SCHEMA_VERSION,
        facts=runtime_preflight_facts(
            engine_revision=revision,
            patch_identity=PATCH_IDENTITY,
            environment_identity=(
                f"linux-amd64;profile-lock-sha256={PROFILE_LOCK_SHA256}"
            ),
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


def _profile_result(raw_result: RunnerResult) -> RunnerResult:
    return replace(
        raw_result,
        backend_provenance=replace(
            raw_result.backend_provenance,
            source_license="Apache-2.0",
            checkpoint_license="Apache-2.0 (upstream claim)",
            environment_identity=(f"linux-amd64;profile-lock-sha256={PROFILE_LOCK_SHA256}"),
            container_digest=_container_digest(),
        ),
    )


def _merge_candidate_results(results: tuple[RunnerResult, ...]) -> RunnerResult:
    if not results:
        raise RuntimeError("Protenix multi-seed execution produced no results")
    first = results[0]
    if any(result.status is not DiffusionStatus.SUCCESS for result in results):
        raise RuntimeError("Protenix multi-seed execution contains a failed result")
    if any(result.backend_provenance != first.backend_provenance for result in results[1:]):
        raise RuntimeError("Protenix multi-seed provenance changed between candidates")

    metrics = tuple(result.resource_metrics for result in results)
    wall_times = tuple(
        metric.wall_time_seconds for metric in metrics if metric.wall_time_seconds is not None
    )
    model_load_times = tuple(
        metric.model_load_seconds for metric in metrics if metric.model_load_seconds is not None
    )
    peak_ram = tuple(
        metric.peak_ram_bytes for metric in metrics if metric.peak_ram_bytes is not None
    )
    peak_vram = tuple(
        metric.peak_vram_bytes for metric in metrics if metric.peak_vram_bytes is not None
    )
    return replace(
        first,
        candidates=tuple(candidate for result in results for candidate in result.candidates),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=sum(wall_times) if wall_times else None,
            model_load_seconds=sum(model_load_times) if model_load_times else None,
            peak_ram_bytes=max(peak_ram) if peak_ram else None,
            peak_vram_bytes=max(peak_vram) if peak_vram else None,
        ),
        message=f"generated {len(results)} deterministic Protenix candidates",
    )


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
    if os.environ.get("DVBFIXER_DIFFUSION_PROTOCOL_VERSION") != str(
        DIFFUSION_RUNNER_PROTOCOL_VERSION
    ):
        raise RuntimeError("diffusion execution protocol version is incompatible")
    checkpoint = _verify_checkpoint(checkpoint)
    preflight()
    sys.dont_write_bytecode = True

    request_path = _workspace_path("DVBFIXER_DIFFUSION_REQUEST", "request.json")
    result_path = _workspace_path("DVBFIXER_DIFFUSION_RESULT", "result.json")
    kalign = shutil.which("kalign")
    if kalign is None:
        raise RuntimeError("protenix-v1-cuda requires the kalign executable")
    kalign_path = _verify_kalign(Path(kalign), _profile_lock())
    root = _source_root() if scripts is None else Path.cwd()
    os.environ["PROTENIX_ROOT_DIR"] = str(root)
    os.environ["LAYERNORM_TYPE"] = "torch"

    builder, adapter = scripts or _load_production_scripts()
    output_root = Path("candidates") / PROFILE
    request_data = json.loads(request_path.read_text(encoding="utf-8"))
    request_seeds = request_data.get("seeds") if isinstance(request_data, dict) else None
    if not isinstance(request_seeds, list) or len(request_seeds) <= 1:
        candidate_requests = ((request_path, output_root),)
    else:
        request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
        candidate_requests_list: list[tuple[Path, Path]] = []
        for seed in request.seeds:
            candidate_request_path = request_path.with_name(f".protenix-seed-{seed}.request.json")
            with candidate_request_path.open("x", encoding="utf-8") as stream:
                stream.write(replace(request, candidate_count=1, seeds=(seed,)).to_json())
            candidate_requests_list.append((candidate_request_path, output_root / f"seed-{seed}"))
        candidate_requests = tuple(candidate_requests_list)

    candidate_results: list[RunnerResult] = []
    for candidate_request_path, candidate_root in candidate_requests:
        input_json = builder.build_input(candidate_request_path, candidate_root / "input")
        raw_result, _summary_path = adapter.run(
            candidate_request_path,
            input_json,
            candidate_root / "raw",
            kalign_path,
            checkpoint,
            cycles=1,
            steps=200,
            ablation_mode=SamplingAblationMode.REINJECTION,
        )
        raw_result = _profile_result(raw_result)
        if raw_result.status is not DiffusionStatus.SUCCESS:
            _write_result_atomic(result_path, raw_result)
            return raw_result
        _verify_success_result(raw_result, Path.cwd().resolve(), refined=False)
        refined_result = refiner(
            candidate_request_path,
            raw_result,
            candidate_root / "refined",
        )
        _verify_success_result(refined_result, Path.cwd().resolve(), refined=True)
        candidate_results.append(refined_result)

    result = _merge_candidate_results(tuple(candidate_results))
    _write_result_atomic(result_path, result)
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
