#!/usr/bin/env python3
"""Two declared private MPS context controls; no frozen/public runner changes."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.protein_context import (
    context_evidence,
    validate_context_capability,
    validate_context_source,
)

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"
RUNTIME = PROJECT / ".artifacts/mac-reconstruction-2026-11-11"
ENCODED = ROOT / "retained-loop-protenix-encoder-attempt01"
NATIVE = ROOT / "retained-loop-protenix-native-preflight-attempt01"
CALLBACK = ROOT / "retained-context-callback-controls-attempt01"
CLONE = RUNTIME / "engines/protenix-v1"


@contextmanager
def context_wiring(module, request):
    """Narrow private adapter: selected observations overwrite; others preserve."""
    plan = request.protein_context
    if plan is None:
        raise ValueError("private context wiring requires an explicit plan")
    selected = tuple(item.atom for item in plan.fixed_conditioning)
    if set(selected) & set(plan.preserved_only) or set(selected) | set(plan.preserved_only) != set(request.fixed_atoms):
        raise ValueError("context roles do not partition immutable source observations")
    original_read, original_class = module._read_fixed_coordinates, module.MPSReinjector

    def read_selected(path, fixed):
        if tuple(fixed) != request.fixed_atoms:
            raise ValueError("legacy caller changed the original fixed mask")
        result = original_read(path, selected)
        if set(result) != set(selected) or any(not np.array_equal(result[item.atom], item.coordinate_angstrom)
                                               for item in plan.fixed_conditioning):
            raise ValueError("selected source references differ from the frozen context plan")
        return result

    class ContextReinjector(original_class):
        def __init__(self, axis, coordinates):
            if set(coordinates) != set(selected):
                raise ValueError("context callback received an incorrect observation mask")
            evidence = context_evidence(request, axis)
            super().__init__(axis, coordinates, alignment_atoms=evidence.alignment_atoms)

    module._read_fixed_coordinates, module.MPSReinjector = read_selected, ContextReinjector
    try:
        yield
    finally:
        module._read_fixed_coordinates, module.MPSReinjector = original_read, original_class


def worker(case, output):
    sys.path.insert(0, str(CLONE))
    sys.path.insert(0, str(PROJECT / "deploy/protenix-v1"))
    import mps_research_smoke as legacy

    request_path = case / "context-request.json"
    request = DiffusionRequest.from_json(request_path.read_text())
    source = case / request.normalized_pdb.path
    if digest(source) != request.normalized_pdb.sha256 or request.seeds != (7,) or request.candidate_count != 1:
        raise ValueError("private control source/candidate/seed differs from declaration")
    validate_context_source(request, source.read_bytes())
    validate_context_capability(request, "protenix")
    if case.name not in {"1NON", "2BE1"}:
        raise ValueError("undeclared development control")
    if output.exists():
        raise ValueError("native prediction directory must not preexist")
    save(output.parent / "context-wiring-policy.json", {"private_adapter": "retained-loop-mps-context-v1",
        "request_sha256": digest(request_path), "source_sha256": digest(source),
        "represented_fixed_atoms": len(request.protein_context.fixed_conditioning),
        "preserved_only_atoms": len(request.protein_context.preserved_only),
        "all_generated_atoms_joint": len(request.generated_atoms), "legacy_runner_sha256": digest(Path(legacy.__file__)),
        "fixed_source_and_generated_ownership_unchanged": True, "public_runner_changed": False})
    args = SimpleNamespace(request=request_path, input=case / "encoded/input.json", assets=RUNTIME / "assets/protenix",
                           kalign=RUNTIME / "tools/bin/kalign", output=output, seed=7)
    with context_wiring(legacy, request):
        legacy.run(args, allowed_seeds=(7,))


def declare(output, previous=None, unexecuted_after=None):
    if previous is not None and unexecuted_after is not None:
        raise ValueError("startup recovery and unexecuted continuation are distinct")
    files = [Path(__file__), PROJECT / "scripts/run_mac_reconstruction_guarded.py",
        PROJECT / "deploy/protenix-v1/mps_research_smoke.py", PROJECT / "deploy/protenix-v1/reinjection.py",
        PROJECT / "deploy/protenix-v1/checkpoint_gap_smoke.py", PROJECT / "deploy/protenix-v1/profile-lock.json"]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    files += list((CLONE / "configs").glob("*.py"))
    for folder in ("protenix/model", "protenix/data", "runner"):
        files += list((CLONE / folder).rglob("*.py"))
    jobs = []
    for entry in ("1NON", "2BE1"):
        case = ENCODED / entry
        native = NATIVE / "execution" / f"{entry}-native-data-seed-0" / "result.json"
        callback = CALLBACK / "execution" / f"{entry}-callback-seed-0" / "result.json"
        if json.loads(native.read_text())["checkpoint_loaded"] or json.loads(callback.read_text())["diffusion_steps"]:
            raise ValueError("prerequisites are not the declared data-only/synthetic controls")
        files += [native, callback] + [p for p in case.rglob("*") if p.is_file()]
        work = output / entry
        command = [sys.executable, str(PROJECT / "scripts/run_mac_reconstruction_guarded.py"),
            "--python", str(RUNTIME / "envs/protenix/bin/python"), "--module", "scripts.run_retained_loop_context_sampling",
            "--output", str(work), "--resource-policy", "pressure-v2", "--timeout", "1800", "--",
            "worker", "--case", str(case), "--output", str(work / "predictions")]
        jobs.append({"case_id": entry, "command": command, "work": str(work)})
    recovery = None
    continuation = None
    if unexecuted_after is not None:
        retained = unexecuted_after / "resource-failure-readback.json"
        failure = json.loads(retained.read_text())
        if (failure["entry_id"] != "1NON" or not failure["model_call_consumed"]
                or failure["resource_outcome"] != "mps-driver-limit"
                or failure["candidate_written"] or failure["second_case_launched"]
                or (unexecuted_after / "2BE1").exists()):
            raise ValueError("continuation requires retained1NON resource failure and unexecuted2BE1")
        for name, expected in failure["bindings"].items():
            if digest(Path(name)) != expected:
                raise ValueError("retained resource failure evidence changed")
        import psutil

        marker = str(unexecuted_after / "1NON/child-entry.py")
        for process in psutil.process_iter(["cmdline"]):
            try:
                if marker in (process.info["cmdline"] or []):
                    raise RuntimeError("original scientific child remains alive; continuation refused")
            except psutil.NoSuchProcess:
                continue
        files += [retained, unexecuted_after / "run_mac_reconstruction_guarded-executed-snapshot.py",
                  unexecuted_after / "run_retained_loop_context_sampling-executed-snapshot.py"]
        files += [Path(name) for name in failure["bindings"]]
        jobs = [job for job in jobs if job["case_id"] == "2BE1"]
        continuation = {"previous_attempt": str(unexecuted_after), "retained1NON": "consumed seed7; MPS resource refusal; never repeated",
            "new_first_execution_only": "2BE1", "original_two_actual_call_budget_unchanged": True,
            "guard_fix": "reap confirmed-dead child after group-signal permission race; thresholds unchanged"}
    if previous is not None:
        blocked_path = previous / "blocked.json"
        blocked = json.loads(blocked_path.read_text())
        work = previous / blocked["case_id"]
        device_path = work / "device-final.json"
        log_path = work / "output.log"
        device = json.loads(device_path.read_text())
        log = log_path.read_text()
        module_failure = "No module named scripts.run_retained_loop_context_sampling" in log
        layernorm_failure = ("No module named 'fast_layer_norm_cuda_v2'" in log
            and "Ninja is required to load C++ extensions" in log
            and "from runner.inference import InferenceRunner, update_inference_configs" in log)
        if (device["model_mps_calls"] != 0 or (work / "predictions").exists()
                or not (module_failure or layernorm_failure) or (previous / "2BE1").exists()):
            raise ValueError("recovery requires a proven pre-runner import failure and no second launch")
        files += [previous / "plan.json", blocked_path, device_path, log_path, work / "execution.json",
                  previous / "entry-script-before-launcher-fix.py"]
        recovery = {"previous_attempt": str(previous), "reason": (
            "native checkout's concrete scripts package shadowed project namespace" if module_failure else
            "missing established MPS LAYERNORM_TYPE=torch caused CUDA extension import before runner construction"),
            "previous_model_calls": 0, "previous_worker_started": not module_failure,
            "previous_runner_or_sampling_started": False,
            "sampler_seed7_first_consumption_still_pending": True, "unchanged_native_call_budget": 2}
    output.mkdir(parents=True, exist_ok=False)
    save(output / "plan.json", {"jobs": jobs, "bindings": {str(p): digest(p) for p in files},
        "new_native_calls": len(jobs), "seed": 7, "candidate_count_per_case": 1, "diffusion_updates": 200,
        "cycles": 1, "MSA": "none", "layernorm_type": "torch", "profile": "private-retained-loop-mps-context-v1",
        "checkpoint": "existing pinned protenix_base_default_v1.0.0.pt; legacy verifies exact asset hashes",
        "context_residues": [354, 365], "gap_lengths": [[18, 18, 18, 18], [10, 20, 8, 10, 20, 8]],
        "scope": "two prospective development technical controls, not final cohort, ten-entry tranche or independent evaluation",
        "input_admission": "existing canonical, geometry-checked and axis-audited source projections; internal gaps only; omitted heterogens water only",
        "no_new_sources_or_limits": True, "no_donor_approval_implied": True,
        "repair_calls_declared": 0, "automatic_sampling_retries": 0, "technical_startup_recovery": recovery,
        "unexecuted_only_continuation": continuation,
        "resource_guard": "serial pressure-v2 RSS8GiB/MPS6GiB,1800s/call, host headroom, no generic fallback"})


def run(output):
    plan = json.loads((output / "plan.json").read_text())
    for name, expected in plan["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"declared development input/code changed: {name}")
    env = {**os.environ, "PYTHONNOUSERSITE": "1", "PYTORCH_ENABLE_MPS_FALLBACK": "0",
        "PROTENIX_ROOT_DIR": str(RUNTIME / "assets/protenix"), "LAYERNORM_TYPE": "torch",
        # The project scripts namespace must load before adding the upstream
        # checkout: its concrete scripts package otherwise shadows this module.
        "PYTHONPATH": os.pathsep.join((str(PROJECT), str(PROJECT / "src")))}
    (output / "run-reservation").mkdir(exist_ok=False)
    for job in plan["jobs"]:
        result = subprocess.run(job["command"], cwd=PROJECT, env=env, check=False)
        if result.returncode:
            save(output / "blocked.json", {"case_id": job["case_id"], "returncode": result.returncode,
                "remaining_calls_not_launched": True, "automatic_retry": False})
            raise RuntimeError("native/resource failure retained; no retries or silent replacements")
    save(output / "completed.json", {"native_calls_completed": len(plan["jobs"]), "scientific_acceptance": False,
        "independent_raw_geometry_and_custody_audit_required": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--case", type=Path)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--unexecuted-after", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve(), args.previous.resolve() if args.previous else None,
                args.unexecuted_after.resolve() if args.unexecuted_after else None)
    elif args.mode == "run":
        run(args.output.resolve())
    else:
        worker(args.case.resolve(), args.output.resolve())
