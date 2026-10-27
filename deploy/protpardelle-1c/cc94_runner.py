#!/usr/bin/env python3
"""Production cc94 ODE runner: native relative chains, fixed sequence, no fallback."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import platform
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

from dvbfixer.model.diffusion import cc94
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    DiffusionRequest,
    RunnerResourceMetrics,
    RunnerResult,
)
from dvbfixer.model.diffusion.preflight import (
    AdapterPreflightCode,
    AdapterPreflightIssue,
    RunnerPreflightReport,
    runtime_preflight_facts,
)
from dvbfixer.model.diffusion.protein_context import validate_cc94_request
from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION
from dvbfixer.model.diffusion.scope import assess_diffusion_scope
from dvbfixer.model.diffusion.validate import validate_runner_result


def _load(name):
    spec = importlib.util.spec_from_file_location(f"cc94_{name}", Path(__file__).with_name(name + ".py"))
    if spec is None or spec.loader is None:
        raise RuntimeError("Protpardelle adapter is unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def refinement_environment(python):
    return {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1",
        "OPENMM_CPU_THREADS": "1", "PATH": str(python.parent) + os.pathsep + os.defpath}


def refinement_preflight(profile):
    if profile == cc94.MPS_PROFILE:
        _load("mps_runner").refinement_preflight()
        return
    import numpy as np
    import openmm
    if openmm.__version__ != "8.6.1" or np.__version__ != "2.4.6" or sys.version_info[:2] != (3, 13):
        raise RuntimeError("cc94 Reference refinement requires Python 3.13, OpenMM 8.6.1, NumPy 2.4.6")
    system = openmm.System()
    system.addParticle(12)
    force = openmm.CustomExternalForce("x*x+y*y+z*z")
    force.addParticle(0, [])
    system.addForce(force)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    context.setPositions([[1, 1, 1]])
    assert np.isclose(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole), 3)


def preflight_report(profile, checkpoint, refinement_python):
    issues = []
    revision, device, framework, checkpoint_digest = "", "", "", ""
    mps = profile == cc94.MPS_PROFILE

    def issue(code, message):
        issues.append(AdapterPreflightIssue(code, message))

    if profile not in cc94.PROFILES or (sys.platform != "darwin" if mps else sys.platform != "linux"):
        issue(AdapterPreflightCode.INCOMPATIBLE_PLATFORM, "cc94 profile does not match the host")
    if platform.machine().lower() not in ({"arm64"} if mps else {"x86_64", "amd64"}):
        issue(AdapterPreflightCode.INCOMPATIBLE_ARCHITECTURE, "cc94 requires native arm64 MPS or amd64 CUDA")
    if sys.version_info[:2] != (3, 12) or os.environ.get("PYTHONNOUSERSITE") != "1":
        issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "cc94 sampling requires isolated Python 3.12")
    try:
        if checkpoint.is_symlink() or not checkpoint.is_file() or checkpoint.name != cc94.CHECKPOINT_NAME:
            raise ValueError("cc94 checkpoint must be a regular file with the pinned name")
        checkpoint_digest = _sha(checkpoint)
        if checkpoint_digest != cc94.CHECKPOINT_SHA256:
            issue(AdapterPreflightCode.CHECKPOINT_DIGEST_MISMATCH, "cc94 checkpoint digest mismatch; cc89 is not a fallback")
        if _sha(checkpoint.parent.parent / "configs/cc94.yaml") != cc94.CONFIG_SHA256:
            issue(AdapterPreflightCode.INCOMPATIBLE_ENVIRONMENT, "cc94 config digest mismatch")
    except (OSError, ValueError) as exc:
        issue(AdapterPreflightCode.MISSING_CHECKPOINT, str(exc))
    try:
        if mps:
            _load("production_runner")._verify_profile_environment()
            revision = cc94.ENGINE_REVISION
        else:
            spec = importlib.util.find_spec("protpardelle")
            if spec is None or spec.origin is None:
                raise ValueError("Protpardelle source unavailable")
            root = Path(spec.origin).resolve().parents[2]
            revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True, timeout=10).strip()
            changed = subprocess.check_output(["git", "-C", str(root), "diff", "--name-only", "HEAD", "--"], text=True, timeout=10).splitlines()
            untracked = subprocess.check_output(["git", "-C", str(root), "ls-files", "--others", "--", "*.py", "*.pyc", "*.so"], text=True, timeout=10).splitlines()
            hashes = _load("cuda_runner").SOURCE_HASHES
            if revision != cc94.ENGINE_REVISION or set(changed) != {"src/protpardelle/core/models.py"} or untracked or any(_sha(root/name) != expected for name, expected in hashes.items()):
                raise ValueError("engine source differs from pinned callback patch")
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        issue(AdapterPreflightCode.INCOMPATIBLE_SOURCE, str(exc))
    fallback_disabled = not mps or os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK") == "0"
    if not fallback_disabled:
        issue(AdapterPreflightCode.FALLBACK_ENABLED, "MPS fallback must be disabled")
    try:
        import torch
        framework = torch.__version__.split("+")[0]
        if torch.__version__ != ("2.6.0" if mps else "2.6.0+cu124"):
            issue(AdapterPreflightCode.INCOMPATIBLE_FRAMEWORK, "cc94 requires pinned PyTorch 2.6.0")
        available = torch.backends.mps.is_available() if mps else torch.cuda.is_available()
        if available:
            device = "mps" if mps else f"cuda:{torch.cuda.current_device()}"
            # Exercise real kernels before starting the sample.
            tensor = torch.ones((2, 2), device=device)
            if float((tensor @ tensor).sum().item()) != 8:
                raise RuntimeError("accelerator tensor probe failed")
        else:
            issue(AdapterPreflightCode.MISSING_MPS if mps else AdapterPreflightCode.MISSING_CUDA, "cc94 accelerator unavailable")
    except (ImportError, RuntimeError) as exc:
        issue(AdapterPreflightCode.INCOMPATIBLE_FRAMEWORK, str(exc))
    try:
        subprocess.run([str(refinement_python), str(Path(__file__).resolve()), "--profile", profile, "--refinement-preflight"],
            env=refinement_environment(refinement_python), check=True, capture_output=True, timeout=90)
    except (OSError, subprocess.SubprocessError):
        issue(AdapterPreflightCode.INCOMPATIBLE_REFINEMENT_PLATFORM, "cc94 isolated refinement kernel probe failed")
    return RunnerPreflightReport(profile, DIFFUSION_RUNNER_PROTOCOL_VERSION, DIFFUSION_SCHEMA_VERSION,
        runtime_preflight_facts(engine_revision=revision,
            patch_identity="cc94-native-relative-v1;source-physical-v1;" + ("apple-patch-v1" if mps else "callback-patch-v1"),
            environment_identity=cc94.MPS_ENVIRONMENT if mps else cc94.CUDA_ENVIRONMENT,
            checkpoint_sha256=checkpoint_digest, framework_version=framework,
            accelerator_available=bool(device), effective_device=device, fallback_disabled=fallback_disabled,
            sampling_platform="native macOS arm64 / MPS" if mps else "Linux amd64 / NVIDIA CUDA 12.4",
            refinement_platform="OpenMM CPU" if mps else "OpenMM Reference"), tuple(issues))


def run(args):
    if os.environ.get("DVBFIXER_DIFFUSION_PROTOCOL_VERSION") != str(DIFFUSION_RUNNER_PROTOCOL_VERSION):
        raise RuntimeError("incompatible execution protocol")
    report = preflight_report(args.profile, args.checkpoint, args.refinement_python)
    if not report.passed:
        raise RuntimeError("; ".join(i.message for i in report.issues))
    path = Path(os.environ.get("DVBFIXER_DIFFUSION_REQUEST", "request.json")).resolve()
    request = DiffusionRequest.from_json(path.read_text())
    validate_cc94_request(request)
    if next(o.value for o in request.backend_options if o.name == "profile") != args.profile:
        raise ValueError("cc94 request and launcher profiles differ")
    admission = assess_diffusion_scope(request, (path.parent/request.normalized_pdb.path).read_bytes())
    if not admission.supported:
        raise ValueError("unsupported cc94 input: " + ", ".join(admission.reasons))
    adapter = _load("checkpoint_gap_smoke")
    os.environ["PROTPARDELLE_MODEL_PARAMS"] = str(args.checkpoint.parent.parent)
    mps = args.profile == cc94.MPS_PROFILE
    if not mps:
        import torch
        torch.cuda.init()
    results = []
    for seed in request.seeds:
        single = replace(request, seeds=(seed,), candidate_count=1)
        single_path = path.parent/f"sampling-seed-{seed}.json"
        single_path.write_text(single.to_json())
        directory = path.parent/"candidates"/args.profile/f"seed-{seed}"

        def sample(request_path, output):
            raw, _ = adapter._run_staged(request_path, output,
                args.checkpoint.parent.parent/"configs/cc94.yaml", args.checkpoint,
                steps=500, step_scale=1.0, s_churn=0.0, per_step_reinjection=False,
                device_name="mps" if mps else "cuda:0", mps_profile=False,
                chain_phase_policy="native", position_feature_policy="source-physical-v1", checkpoint_model="cc94")
            return replace(raw, backend_provenance=replace(raw.backend_provenance,
                source_license="MIT", environment_identity=cc94.MPS_ENVIRONMENT if mps else cc94.CUDA_ENVIRONMENT))

        if single.gap_contexts is not None:
            from dvbfixer.model.diffusion.independent_gaps import sample_independent_gaps
            raw = sample_independent_gaps(single_path, directory/"raw", sample)
        else:
            (directory/"raw").mkdir(parents=True, exist_ok=False)
            raw = sample(single_path, directory/"raw")
        raw_path = directory/"raw-result.json"
        raw_path.write_text(raw.to_json())
        refined_path = directory/"refined-result.json"
        subprocess.run([str(args.refinement_python), str(Path(__file__).resolve()), "--profile", args.profile,
            "--refine-worker", "--request", str(single_path), "--raw-result", str(raw_path),
            "--output", str(directory/"refined"), "--result", str(refined_path)],
            env=refinement_environment(args.refinement_python), check=True)
        refined = RunnerResult.from_json(refined_path.read_text())
        validate_runner_result(single, refined, workspace=path.parent)
        results.append(refined)
    final = replace(results[0], candidates=tuple(c for r in results for c in r.candidates),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=sum(r.resource_metrics.wall_time_seconds or 0 for r in results),
            peak_ram_bytes=max(r.resource_metrics.peak_ram_bytes or 0 for r in results),
            peak_vram_bytes=max(r.resource_metrics.peak_vram_bytes or 0 for r in results)))
    Path(os.environ.get("DVBFIXER_DIFFUSION_RESULT", "result.json")).write_text(final.to_json())


def main():
    sys.dont_write_bytecode = True
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=sorted(cc94.PROFILES), default=cc94.MPS_PROFILE if sys.platform == "darwin" else cc94.CUDA_PROFILE)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--refinement-python", type=Path)
    for name in ("preflight", "refine-worker", "refinement-preflight"):
        parser.add_argument("--" + name, action="store_true")
    for name in ("request", "raw-result", "output", "result"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    if args.refinement_preflight:
        refinement_preflight(args.profile)
    elif args.refine_worker:
        from dvbfixer.model.diffusion.runner_refinement import refine_runner_result
        refinement_preflight(args.profile)
        raw = RunnerResult.from_json(args.raw_result.read_text())
        refined = refine_runner_result(args.request, raw, args.output,
            platform_name="CPU" if args.profile == cc94.MPS_PROFILE else "Reference")
        args.result.write_text(refined.to_json())
    elif args.preflight:
        report = preflight_report(args.profile, args.checkpoint, args.refinement_python)
        print(report.to_json())
        raise SystemExit(0 if report.passed else 1)
    else:
        run(args)


if __name__ == "__main__":
    main()
