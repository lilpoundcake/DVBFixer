#!/usr/bin/env python3
"""Prospective single-case development replication; no held-out seed use.

The three seed-7 pilots are historical controls, not prospectively declared here.
Runs additional approved development seeds serially with the same model settings.
This is not the later multi-case conditioning or anchor-calibration study.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import traceback
from dataclasses import replace
from pathlib import Path

from scripts.audit_mac_reconstruction_smoke import audit

from dvbfixer.model.diffusion.contract import DiffusionRequest

ROOT = Path("tmp/mac-reconstruction-2026-11-11").absolute()
RUNTIME = Path(".artifacts/mac-reconstruction-2026-11-11").absolute()
SEEDS = (42, 101, 102, 103, 104)
ARMS = ("rfd3", "boltz", "protenix")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def declare(output):
    output.mkdir(parents=True, exist_ok=False)
    request = DiffusionRequest.from_json((ROOT / "rfd3-audit/evaluation-mask.json").read_text())
    jobs = []
    for seed in SEEDS:
        for arm in ARMS:
            base = ROOT / f"{arm}-smoke-attempt02"
            original = json.loads((base / "command.json").read_text())
            work = output / f"{arm}-seed{seed}"
            work.mkdir()
            args = [part.replace(str(base), str(work / "run")).replace(
                str(base.relative_to(Path.cwd())), str(work / "run")) for part in original[4:]]
            if arm == "rfd3":
                args = [f"seed={seed}" if part == "seed=7" else part for part in args]
            elif arm == "boltz":
                args[args.index("--seed") + 1] = str(seed)
            else:
                current = replace(request, seeds=(seed,))
                (work / "request.json").write_text(current.to_json())
                (work / current.normalized_pdb.path).write_bytes((ROOT / "rfd3-audit/source.pdb").read_bytes())
                inputs = json.loads((ROOT / "protenix-input/input.json").read_text())
                inputs[0]["modelSeeds"] = [seed]
                (work / "input.json").write_text(json.dumps(inputs) + "\n")
                args[args.index("--request") + 1] = str(work / "request.json")
                args[args.index("--input") + 1] = str(work / "input.json")
                args += ["--seed", str(seed)]
            command = [sys.executable, "scripts/run_mac_reconstruction_guarded.py", "--python", original[0],
                       "--module", original[3], "--output", str(work / "run"), "--resource-policy", "pressure-v2",
                       "--timeout", "1800", "--", *args]
            jobs.append({"arm": arm, "seed": seed, "work": str(work), "command": command})
    files = [Path(p) for p in ("scripts/run_mac_reconstruction_guarded.py", "scripts/audit_mac_reconstruction_smoke.py",
             "scripts/run_mac_reconstruction_development_panel.py", "deploy/protenix-v1/mps_research_smoke.py",
             "deploy/protenix-v1/mps-confidence-dtype.patch")]
    files += list((ROOT / "smoke-input-03").glob("*"))
    files += list((ROOT / "protenix-input").glob("*"))
    for job in jobs:
        files += [p for p in Path(job["work"]).iterdir() if p.is_file()]
    plan = {"case": "explicit 1UBQ A1-72 gap10-13", "seeds": [7, *SEEDS], "new_seeds": list(SEEDS),
            "reserved_seeds_not_used": [211, 223, 227, 229, 233],
            "seed7": "reuse retained original pilots, including technical failures and materialization-only recovery",
            "arm_semantics": {"rfd3": "fixed-sequence reconstruction", "boltz": "template-conditioned proposal",
                              "protenix": "template plus per-update fixed-coordinate reinjection; explicit CPU 3x3 fit"},
            "same_sampler_settings_as_seed7": True, "sampling_steps_requested": 200,
            "serial_ml": True, "automatic_retry": False, "new_refinement_in_this_panel": False,
            "ranking_optimization": False, "qualification": False,
            "evaluation": "same existing inventory, drift, junction and geometry gates; held-out RMSD scored only after proposal",
            "file_sha256": {str(p.absolute()): digest(p) for p in files}, "jobs": jobs}
    (output / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    print(f"Declared {len(jobs)} new serial jobs; seed 7 reused, reserved seeds untouched")


def run(output):
    plan = json.loads((output / "plan.json").read_text())
    for path, expected in plan["file_sha256"].items():
        if digest(Path(path)) != expected:
            raise ValueError(f"declared input/code changed: {path}")
    env = os.environ.copy()
    env.update(PYTHONPATH=os.pathsep.join(str(Path(p).absolute()) for p in
        (".", "src", "deploy/protenix-v1", str(RUNTIME / "engines/protenix-v1"))),
        PROTENIX_ROOT_DIR=str(RUNTIME / "assets/protenix"), LAYERNORM_TYPE="torch")
    results = []
    with (output / "orchestration.log").open("xb") as log:
        for job in plan["jobs"]:
            work = Path(job["work"])
            result = subprocess.run(job["command"], env=env, stdout=log, stderr=subprocess.STDOUT)
            row = {"arm": job["arm"], "seed": job["seed"], "returncode": result.returncode}
            execution_path = work / "run/execution.json"
            execution = json.loads(execution_path.read_text()) if execution_path.exists() else {}
            row["execution"] = execution
            if result.returncode == 0:
                try:
                    pattern = {"rfd3": "*.cif.gz", "boltz": "*_model_0.pdb", "protenix": "candidate.pdb"}[job["arm"]]
                    candidates = list((work / "run/predictions").rglob(pattern))
                    if len(candidates) != 1:
                        raise ValueError("missing/ambiguous prediction")
                    row["audit"] = audit(ROOT / "smoke-input-03", candidates[0], work / "audit", seed=job["seed"])
                except Exception:
                    row["audit_failure"] = traceback.format_exc()
            results.append(row)
            (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
            print(json.dumps({k: row[k] for k in ("arm", "seed", "returncode")}), flush=True)
            if execution.get("guard_stop_reason") or not execution:
                raise RuntimeError("panel stopped on resource/orchestration failure; no automatic retry")
    (output / "completed.json").write_text(json.dumps({"jobs_executed": len(results), "qualified": False}) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    (declare if args.mode == "declare" else run)(args.output.absolute())
