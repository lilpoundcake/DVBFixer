#!/usr/bin/env python3
"""Declare/run one seed per admitted reference case/backend, strictly serially.

This precedes a calibrated mobility envelope or broad qualification. All original
raw failures, native drift, source-restored damage and private artifacts remain.
No validation seeds, reference-conditioned proposals or implicit retries.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import traceback
from pathlib import Path

from scripts.audit_mac_reconstruction_smoke import audit
from scripts.remap_reference_model_output import run as remap

PROJECT = Path(__file__).parents[1]
EVIDENCE = PROJECT / "tmp/mac-reconstruction-2026-11-11"
RUNTIME = PROJECT / ".artifacts/mac-reconstruction-2026-11-11"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def declare(preparation, output):
    rows = json.loads((preparation / "preparation-results.json").read_text())
    jobs = []
    files = []
    output.mkdir(exist_ok=False)
    for row in rows:
        if not row["ready"]:
            continue
        case = preparation / row["case_id"]
        manifest = json.loads((case / "manifest.json").read_text())
        if not manifest["source_admissible"] or not manifest["reference_baseline_accepted"]:
            raise ValueError("case source/reference readiness not established")
        files.extend(p for p in case.rglob("*") if p.is_file())
        for arm in ("rfd3","boltz","protenix"):
            work = output / f"{row['case_id']}-{arm}"
            old = EVIDENCE / ("rfd3-explicit-seed7" if arm == "rfd3" else f"{arm}-smoke-attempt02")
            native = json.loads((old / "command.json").read_text())
            args = list(native[4:])
            def set_arg(key, value, arguments=args):
                arguments[arguments.index(key)+1] = str(value)
            if arm == "rfd3":
                args = [f"inputs={case / 'rfd3/rfd3-input.json'}" if a.startswith("inputs=") else
                        f"out_dir={work / 'run/predictions'}" if a.startswith("out_dir=") else a for a in args]
            elif arm == "boltz":
                args[1] = str(case / "boltz-2/input.yaml")
                set_arg("--out_dir", work / "run/predictions")
            else:
                set_arg("--request", case / "evaluation-mask.json")
                set_arg("--input", case / "protenix-v1/input.json")
                set_arg("--output", work / "run/predictions")
            command = [sys.executable, "scripts/run_mac_reconstruction_guarded.py", "--python", native[0],
                "--module", native[3], "--output", str(work / "run"), "--resource-policy", "pressure-v2", "--timeout", "1800", "--", *args]
            jobs.append({"case":str(case), "arm":arm, "seed":7, "work":str(work), "command":command,
                         "expected_fixed":manifest["observed_atoms"], "expected_generated":manifest["generated_atoms"]})
    files += [PROJECT / p for p in ("scripts/run_reference_reconstruction_pilot.py", "scripts/audit_mac_reconstruction_smoke.py",
        "scripts/remap_reference_model_output.py", "scripts/run_mac_reconstruction_guarded.py", "scripts/trace_rfd3_seed_transport.py",
        "scripts/rfd3_mps_seeded_research.py", "deploy/protenix-v1/mps_research_smoke.py", "deploy/boltz-2/atom_mapping.py",
        "deploy/protenix-v1/mps-confidence-dtype.patch", "deploy/rfd3/source-component-width.patch")]
    plan = {"purpose":"two new reference-case seed-7 transport/science pilots, before envelope calibration",
        "seed":7, "serial_ml":True, "reserved_seeds_used":False, "automatic_retries":False,
        "sampler_settings":"same per-backend 200-step settings as prior tiny-case arms; RFD3 explicit MPS seed contract",
        "boltz_semantics":"template-conditioned proposal, native drift retained before source restoration",
        "protenix_semantics":"per-step exact reinjection; native extra OXT retained but excluded from requested projection",
        "heldout_reference_to_models":False, "qualification":False,
        "excluded_cases":[r for r in rows if not r["ready"]], "file_sha256":{str(p):digest(p) for p in files}, "jobs":jobs}
    (output / "plan.json").write_text(json.dumps(plan,indent=2)+"\n")
    print(f"Declared {len(jobs)} serial reference-case pilots")


def run(output):
    plan = json.loads((output / "plan.json").read_text())
    for name, expected in plan["file_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"declared pilot input/code changed: {name}")
    rows = []
    with (output / "orchestration.log").open("xb") as log:
        for job in plan["jobs"]:
            env = os.environ.copy()
            paths = [".","src","scripts"]
            if job["arm"] == "protenix":
                paths += ["deploy/protenix-v1", ".artifacts/mac-reconstruction-2026-11-11/engines/protenix-v1"]
            env.update(PYTHONPATH=os.pathsep.join(str(PROJECT / p) for p in paths),
                PROTENIX_ROOT_DIR=str(RUNTIME / "assets/protenix"), LAYERNORM_TYPE="torch")
            work = Path(job["work"])
            if job["arm"] == "rfd3":
                env.update(DVBFIXER_RFD3_MODEL_SEED="7",DVBFIXER_RFD3_TRANSPORT_OUTPUT=str(work / "run/transport"),
                    DVBFIXER_RFD3_EXPECTED_FIXED=str(job["expected_fixed"]),DVBFIXER_RFD3_EXPECTED_GENERATED=str(job["expected_generated"]))
            code = subprocess.run(job["command"],env=env,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT).returncode
            execution_path = work / "run/execution.json"
            execution = json.loads(execution_path.read_text()) if execution_path.exists() else {}
            row = {"case":job["case"], "arm":job["arm"], "sampling_seed":7, "returncode":code, "execution":execution}
            if code == 0:
                try:
                    patterns = {"rfd3":"*.cif.gz", "boltz":"*_model_0.pdb", "protenix":"candidate.pdb"}
                    candidates = list((work / "run/predictions").rglob(patterns[job["arm"]]))
                    if len(candidates) != 1:
                        raise ValueError("native output is missing/ambiguous")
                    prediction = candidates[0]
                    if job["arm"] in ("rfd3","boltz"):
                        prediction = remap(Path(job["case"]),prediction,work / "native-mapping",model=job["arm"])
                    row["audit"] = audit(Path(job["case"]),prediction,work / "audit",seed=7)
                except Exception:
                    row["audit_failure"] = traceback.format_exc()
            rows.append(row)
            (output / "results.json").write_text(json.dumps(rows,indent=2)+"\n")
            print(json.dumps({"case":Path(job["case"]).name,"arm":job["arm"],"returncode":code}),flush=True)
            if not execution or execution.get("guard_stop_reason"):
                raise RuntimeError("pilot stopped on a resource/orchestration failure; no implicit retry")
    (output / "completed.json").write_text(json.dumps({"executed":len(rows),"qualified":False})+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=("declare","run"))
    parser.add_argument("preparation",type=Path)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.preparation.absolute(),args.output.absolute())
    else:
        run(args.output.absolute())
