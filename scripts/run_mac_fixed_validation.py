#!/usr/bin/env python3
"""Frozen private fixed-source validation: six masks, five reserved seeds.

90 serial ML calls plus 30 one-fast-loop MODELLER comparators. No automatic
backend/candidate selection, anchor motion, remote MSA or public publication.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import traceback
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save
from scripts.mac_validation_inputs import prepare, transport
from scripts.mac_validation_policy import SEEDS, require_policy
from scripts.refine_mac_reconstruction_smoke import run as refine
from scripts.run_mac_modeller_comparison import execute as execute_cpu
from scripts.run_mac_modeller_proposal import run as modeller_proposal

from dvbfixer.model.diffusion.contract import DiffusionRequest

PROJECT = Path(__file__).resolve().parents[1]
OLD = PROJECT / "tmp/mac-reconstruction-2026-11-11"
RUNTIME = PROJECT / ".artifacts/mac-reconstruction-2026-11-11"


def declare(cohort: Path, output: Path):
    candidates = json.loads((cohort / "results.json").read_text())
    if len(candidates) != 6 or any(not r.get("mask_admitted") for r in candidates):
        raise ValueError("all six prospectively declared masks must retain their admission outcomes")
    output.mkdir(parents=True, exist_ok=False)
    assets = {
        RUNTIME / "assets/rfd3.ckpt": "9b3f85923e0d51e9453e15cdd2f8c666e7ce096a60577f57d11bbc54ae6d67c1",
        RUNTIME / "assets/boltz/boltz2_conf.ckpt": "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1",
        RUNTIME / "assets/protenix/checkpoint/protenix_base_default_v1.0.0.pt": "2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04",
    }
    for path, expected in assets.items():
        if digest(path) != expected:
            raise ValueError(f"previously reviewed checkpoint changed: {path}")
    revisions = {"rfd3": "7a9c246fbf140a743d14a35248f5039e7723a23a",
        "boltz": "1431d573c3a02a2cb6b5c3b775e7128c7a98d539", "protenix-v1": "85767b811c40ed46e73a9b39519cf6bfca8701ba"}
    (output / "engine-source").mkdir()
    for name, expected in revisions.items():
        engine = RUNTIME / "engines" / name
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=engine, text=True).strip()
        if revision != expected:
            raise ValueError(f"reviewed engine revision changed: {name}")
        (output / "engine-source" / f"{name}.patch").write_bytes(subprocess.check_output(["git", "diff"], cwd=engine))
    cases = []
    for candidate in candidates:
        case = output / "cases" / candidate["case"]
        prepare(cohort / candidate["case"], case)
        cases.append(case)
    python = Path(sys.executable).absolute()
    jobs, cpu_jobs = [], []
    for seed in SEEDS:
        for case in cases:
            request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
            for model in ("rfd3", "boltz", "protenix"):
                old = OLD / ("rfd3-explicit-seed7" if model == "rfd3" else f"{model}-smoke-attempt02")
                native = json.loads((old / "command.json").read_text())
                arguments = list(native[4:])
                work = output / "runs" / f"{case.name}-{model}-seed{seed}"
                work.mkdir(parents=True)
                def set_arg(key, value, args=arguments):
                    args[args.index(key)+1] = str(value)
                if model == "rfd3":
                    arguments = [f"inputs={case / 'rfd3/rfd3-input.json'}" if a.startswith("inputs=") else
                        f"out_dir={work / 'run/predictions'}" if a.startswith("out_dir=") else
                        f"seed={seed}" if a.startswith("seed=") else a for a in arguments]
                    module = "rfd3_mps_validation_research"
                elif model == "boltz":
                    arguments[1] = str(case / "boltz-2/input.yaml")
                    set_arg("--out_dir", work / "run/predictions")
                    set_arg("--seed", seed)
                    module = native[3]
                else:
                    current = replace(request, seeds=(seed,), candidate_count=1)
                    (work / "request.json").write_text(current.to_json())
                    (work / "source.pdb").write_bytes((case / "source.pdb").read_bytes())
                    inputs = json.loads((case / "protenix-v1/input.json").read_text())
                    inputs[0]["modelSeeds"] = [seed]
                    save(work / "input.json", inputs)
                    set_arg("--request", work / "request.json")
                    set_arg("--input", work / "input.json")
                    set_arg("--output", work / "run/predictions")
                    if "--seed" in arguments:
                        set_arg("--seed", seed)
                    else:
                        arguments += ["--seed", str(seed)]
                    module = "mps_validation_smoke"
                command = [str(python), "scripts/run_mac_reconstruction_guarded.py", "--python", native[0], "--module", module,
                    "--output", str(work / "run"), "--resource-policy", "pressure-v2", "--timeout", "1800", "--", *arguments]
                jobs.append({"case": str(case), "model": model, "seed": seed, "work": str(work), "command": command})
            work = output / "cpu/execution" / f"{case.name}-seed-{seed}"
            cpu_jobs.append({"case": str(case), "case_id": case.name, "seed": seed,
                "command": [str(python), "-m", "scripts.run_mac_fixed_validation", "modeller-worker", "--output", str(work),
                    "--case", str(case), "--seed", str(seed), "--policy", str(output / "policy.json")]})
    files = [PROJECT / "scripts" / name for name in ("run_mac_fixed_validation.py", "mac_validation_inputs.py", "mac_validation_policy.py",
        "rfd3_mps_validation_research.py", "trace_rfd3_seed_transport.py", "run_mac_reconstruction_guarded.py", "preflight_mac_conditioning.py",
        "run_mac_modeller_comparison.py", "run_mac_modeller_proposal.py", "mac_conditioning_contracts.py", "refine_mac_reconstruction_smoke.py",
        "preflight_mac_protenix_validation.py", "preflight_mac_boltz_validation.py")]
    files += [PROJECT / "deploy/protenix-v1" / name for name in ("mps_research_smoke.py", "mps_validation_smoke.py", "checkpoint_gap_smoke.py", "reinjection.py", "profile-lock.json")]
    files.append(PROJECT / "deploy/boltz-2/atom_mapping.py")
    for name in revisions:
        engine = RUNTIME / "engines" / name
        prefixes = [engine / "models/rfd3/src"] if name == "rfd3" else [engine / "src/boltz"] if name == "boltz" else [engine / "protenix", engine / "runner", engine / "configs"]
        files.extend(p for prefix in prefixes for p in prefix.rglob("*.py"))
    files += list((output / "engine-source").glob("*.patch"))
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    files += [p for p in (output / "cases").rglob("*") if p.is_file()]
    files += [p for p in (output / "runs").rglob("*") if p.is_file()]
    bindings = {str(path): digest(path) for path in files}
    bindings.update({str(path): expected for path, expected in assets.items()})
    save(output / "jobs.json", jobs)
    policy = {"protocol": "private-fixed-source-validation-v1", "seeds": SEEDS, "models": ["rfd3", "boltz", "protenix", "modeller"],
        "observed_motion": "none", "automatic_candidate_selection": False, "public_capability_enabled": False,
        "validation_seed_sampling_authorized": True, "context": "full unique protein chain; canonical heavy atoms; declared observed OXT pass-through where absent from native axis",
        "native_fidelity": "RFD3/Protenix represented observations fixed; Boltz/MODELLER native drift measured before explicit restoration",
        "source_scope": "canonical protein-only source projections, exact original occupancies and observations; no source repair",
        "repair": "necessary-only all-generated joint CPU1thread, 1000 iterations, optimizerseed7, zero restarts/preconditioning/anchors",
        "acceptance": "unchanged existing versioned geometry plus exact source custody, finite complete inventory and global L",
        "accuracy": "all-candidate heavy/backbone RMSD after generation/refinement only; no reference-based best-of-five",
        "statistical_units": "two protein families/entries with three masks each, not 30 independent structures",
        "training_membership": "not established; workflow-held-out does not mean absent from model training",
        "resource_policy": "serial pressure-v2, 1800s/ML call, RSS8GiB/MPS6GiB; no generic fallback or remote MSA",
        "storage": "estimated under 2GiB; 10GiB study-output cap checked between jobs; require 2GiB free headroom before each launch",
        "calls": {"ML": 90, "MODELLER": 30}, "production_qualification": False,
        "engine_revisions": revisions, "same_prior_approved_checkpoints": {str(p): value for p, value in assets.items()},
        "cohort_admission_sha256": digest(cohort / "results.json"), "jobs_sha256": digest(output / "jobs.json"), "file_sha256": bindings}
    save(output / "policy.json", policy)
    (output / "cpu").mkdir()
    save(output / "cpu/plan.json", {"python": str(python), "python_sha256": digest(python), "jobs": cpu_jobs, "count": 30,
        "source_code_and_input_sha256": bindings, "policy_sha256": digest(output / "policy.json"), "qualification": False})
    save(output / "plan.json", {"policy_sha256": digest(output / "policy.json"), "cpu_plan_sha256": digest(output / "cpu/plan.json"),
        "calls": 90, "modeller_calls": 30, "preflight_required": "all six RFD3 masks before model weights; native transport/custody checked for every model output",
        "code_and_input_files": bindings, "qualification": False})
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps({"ML_calls": 90, "MODELLER_calls": 30, "policy_sha256": digest(output / "policy.json"), "preflight_before_sampling": True}, indent=2))


def finalize(case, work, *, model, seed):
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    source = (case / "source.pdb").read_text()
    assess = assessment_adapter(request, source)
    raw = (work / "audit/source-restored.pdb").read_text()
    before = assess(raw)
    repair_error = None
    try:
        refine(work / "audit", work / "common-repair", platform="CPU")
    except Exception:
        repair_error = traceback.format_exc()
        save(work / "repair-failure.json", {"traceback": repair_error, "scientific_failure_not_dropped": True})
    final = work / "common-repair/refined.pdb"
    repaired = final.is_file()
    if not repaired:
        final = work / "audit/source-restored.pdb"
    text = final.read_text()
    assessed = assess(text)
    originals = [line for line in source.splitlines(keepends=True) if line.startswith("ATOM  ")]
    source_exact = all(line in text.splitlines(keepends=True) for line in originals)
    if not source_exact:
        raise ValueError("validation changed observed ATOM records, including terminal OXT")
    all_l = True
    try:
        structure = _structure(text, "fixed-validation-final.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
    except Exception:
        all_l = False
    reference, positions = inventory(case / "reference-projection-NOT-CONDITIONING.pdb"), inventory(final)
    errors = {a: float(np.sum((reference[a][2]-positions[a][2])**2)) for a in request.generated_atoms}
    row = {"case": case.name, "backend": model, "seed": seed, "repaired": repaired, "repair_error": repair_error,
        "raw_restored_geometry_passed": before.geometry_passed, "raw_admission_failures": list(before.admission_failures),
        "final_geometry_passed": assessed.geometry_passed, "final_admission_failures": list(assessed.admission_failures),
        "source_atom_records_byte_identical": source_exact, "global_all_l": all_l,
        "final_accepted": repair_error is None and assessed.geometry_passed and not assessed.admission_failures and all_l,
        "heldout_heavy_rmsd": float(np.sqrt(np.mean(list(errors.values())))),
        "heldout_backbone_rmsd": float(np.sqrt(np.mean([e for a, e in errors.items() if a.atom_name in {"N", "CA", "C", "O"}]))),
        "native_transport": json.loads((work / "audit/transport.json").read_text()),
        "final_path": str(final), "final_sha256": digest(final), "assessment": asdict(assessed), "qualification": False}
    save(work / "result.json", row)
    return row


def modeller_worker(case, output, seed, policy):
    require_policy("modeller", seed, policy)
    modeller_proposal(case / "evaluation-mask.json", case / "source.pdb", output / "proposal", seed, validation_policy=policy)
    transport(case, output / "proposal/mapped-native.pdb", output / "audit", model="modeller")
    finalize(case, output, model="modeller", seed=seed)


def run(output):
    plan = json.loads((output / "plan.json").read_text())
    if digest(output / "policy.json") != plan["policy_sha256"] or digest(output / "cpu/plan.json") != plan["cpu_plan_sha256"]:
        raise ValueError("frozen validation policy/CPU plan changed")
    policy = json.loads((output / "policy.json").read_text())
    if digest(output / "jobs.json") != policy["jobs_sha256"]:
        raise ValueError("frozen validation sampling matrix changed")
    for name, expected in plan["code_and_input_files"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"frozen validation code/input changed: {name}")
    for case in sorted((output / "cases").iterdir()):
        preflight = json.loads((case / "rfd3/preflight/preflight.json").read_text())
        if (preflight["manifest_sha256"] != digest(case / "rfd3/manifest.json")
                or preflight["spec_sha256"] != digest(case / "rfd3/rfd3-input.json")
                or not preflight["exact_source_identity_masks_checked"] or not preflight["all_sequence_fixed"]
                or not preflight["unfixed_coordinates_zero"]):
            raise ValueError("validation RFD3 builder preflight changed")
        protenix = json.loads((case / "protenix-v1/preflight/preflight.json").read_text())
        boltz = json.loads((case / "boltz-2/preflight/preflight.json").read_text())
        if (not protenix["all_requested_observations_including_OXT_represented"] or not boltz["exact_token_atom_feature_axis_verified"]
                or protenix["checkpoint_loaded"] or boltz["checkpoint_loaded"]
                or protenix["request_sha256"] != digest(case / "evaluation-mask.json")
                or boltz["request_sha256"] != digest(case / "evaluation-mask.json")):
            raise ValueError("pinned native feature-axis preflight did not pass")
    os.environ.update(DVBFIXER_VALIDATION_POLICY=str(output / "policy.json"), DVBFIXER_VALIDATION_POLICY_SHA256=plan["policy_sha256"])
    save(output / "execution-started.json", {"policy_sha256": plan["policy_sha256"], "qualification": False})
    rows = []
    for job in json.loads((output / "jobs.json").read_text()):
        if not output.is_dir() or shutil.disk_usage(output).free < 2*1024**3:
            raise RuntimeError("validation output root/headroom unavailable; no more jobs launched")
        size = 0
        for path in output.rglob("*"):
            try:
                if path.is_file():
                    size += path.stat().st_size
            except FileNotFoundError:
                if not output.is_dir():
                    raise
        if size > 10*1024**3:
            raise RuntimeError("validation study-output cap exceeded; no more jobs launched")
        case, work, model = Path(job["case"]), Path(job["work"]), job["model"]
        permission = require_policy(model, job["seed"])
        save(work / "validation-permission.json", permission)
        paths = [PROJECT, PROJECT / "src", PROJECT / "scripts"]
        if model == "protenix":
            paths += [PROJECT / "deploy/protenix-v1", PROJECT / ".artifacts/mac-reconstruction-2026-11-11/engines/protenix-v1"]
        env = {**os.environ, "PYTHONPATH": os.pathsep.join(map(str, paths)),
            "PROTENIX_ROOT_DIR": str(PROJECT / ".artifacts/mac-reconstruction-2026-11-11/assets/protenix"), "LAYERNORM_TYPE": "torch"}
        if model == "rfd3":
            manifest = json.loads((case / "rfd3/manifest.json").read_text())
            env.update(DVBFIXER_RFD3_MODEL_SEED=str(job["seed"]), DVBFIXER_RFD3_TRANSPORT_OUTPUT=str(work / "run/transport"),
                DVBFIXER_RFD3_EXPECTED_FIXED=str(manifest["sampler_fixed_atoms"]), DVBFIXER_RFD3_EXPECTED_GENERATED=str(manifest["sampler_unfixed_atoms"]))
        with (work / "orchestration.log").open("xb") as log:
            code = subprocess.run(job["command"], env=env, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT).returncode
        execution_path = work / "run/execution.json"
        execution = json.loads(execution_path.read_text()) if execution_path.exists() else {}
        save(work / "execution.json", execution)
        if code or not execution or execution.get("guard_stop_reason"):
            save(work / "failure.json", {"returncode": code, "execution": execution, "not_retried": True})
            raise RuntimeError("validation technical/resource failure retained; remaining jobs not launched")
        pattern = {"rfd3": "*.cif.gz", "boltz": "*_model_0.pdb", "protenix": "candidate.pdb"}[model]
        native = list((work / "run/predictions").rglob(pattern))
        if len(native) != 1:
            raise ValueError("native output missing or ambiguous")
        try:
            transport(case, native[0], work / "audit", model=model)
            row = finalize(case, work, model=model, seed=job["seed"])
        except Exception:
            save(work / "failure.json", {"traceback": traceback.format_exc(), "unchanged_native_output_retained": str(native[0])})
            raise
        rows.append(row)
        (output / "ML-progress.json").write_text(json.dumps(rows, indent=2)+"\n")
        print(json.dumps({"case": case.name, "model": model, "seed": job["seed"], "accepted": row["final_accepted"]}), flush=True)
    save(output / "ML-completed.json", {"executed": len(rows), "qualified": False})
    execute_cpu(output / "cpu")
    save(output / "completed.json", {"ML_executed": len(rows), "MODELLER_executed": 30, "qualification": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "modeller-worker"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cohort", type=Path)
    parser.add_argument("--case", type=Path)
    parser.add_argument("--seed", type=int, choices=SEEDS)
    parser.add_argument("--policy", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.cohort.resolve(), args.output.resolve())
    elif args.mode == "run":
        run(args.output.resolve())
    else:
        modeller_worker(args.case.resolve(), args.output.resolve(), args.seed, args.policy.resolve())
