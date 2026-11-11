#!/usr/bin/env python3
"""Predeclared E subset: no-MSA source-template vs no-template proposals.

Twelve new development calls; existing template arms remain unchanged. No remote
MSA, homolog retrieval or validation-case tuning. RFD3 has no matched MSA arm.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback
from dataclasses import replace
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import audit, inventory
from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save
from scripts.refine_mac_reconstruction_smoke import run as refine

from dvbfixer.model.diffusion.contract import DiffusionRequest

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-reconstruction-comparison-2026-11-11"
OLD = PROJECT / "tmp/mac-reconstruction-2026-11-11"
RUNTIME = PROJECT / ".artifacts/mac-reconstruction-2026-11-11"


def declare(output: Path):
    case = ROOT / "modeller-attempt02/cases/1UBQ-A1-72-gap10-13"
    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    output.mkdir(parents=True, exist_ok=False)
    jobs = []
    for seed in DEVELOPMENT_SEEDS:
        for model in ("boltz", "protenix"):
            work = output / f"{model}-seed{seed}"
            work.mkdir()
            native = json.loads((OLD / f"{model}-smoke-attempt02/command.json").read_text())
            arguments = list(native[4:])
            if model == "boltz":
                sequence = request.target_sequences[0]
                save(work / "input.yaml", {"version": 1, "sequences": [{"protein": {"id": sequence.chain,
                    "sequence": sequence.sequence, "msa": "empty"}}]})
                arguments[1] = str(work / "input.yaml")
                arguments[arguments.index("--seed")+1] = str(seed)
                arguments[arguments.index("--out_dir")+1] = str(work / "run/predictions")
            else:
                current = replace(request, seeds=(seed,), candidate_count=1)
                (work / "request.json").write_text(current.to_json())
                (work / "source.pdb").write_bytes((case / "source.pdb").read_bytes())
                save(work / "input.json", [{"name": "private_sequence_evidence_control", "sequences": [{"proteinChain": {
                    "sequence": current.target_sequences[0].sequence, "count": 1}}], "modelSeeds": [seed]}])
                for name, value in (("--request", work / "request.json"), ("--input", work / "input.json"), ("--output", work / "run/predictions")):
                    arguments[arguments.index(name)+1] = str(value)
                arguments += ["--seed", str(seed)]
            module = "mps_sequence_evidence_smoke"
            command = [sys.executable, "scripts/run_mac_reconstruction_guarded.py", "--python", native[0], "--module", module,
                "--output", str(work / "run"), "--resource-policy", "pressure-v2", "--timeout", "1800", "--", *arguments]
            jobs.append({"model": model, "seed": seed, "case": str(case), "work": str(work), "command": command})
    files = [Path(__file__), PROJECT / "scripts/audit_mac_reconstruction_smoke.py", PROJECT / "scripts/refine_mac_reconstruction_smoke.py",
        PROJECT / "scripts/run_mac_reconstruction_guarded.py", PROJECT / "deploy/protenix-v1/mps_research_smoke.py"]
    files += [PROJECT / "deploy" / name / "mps_sequence_evidence_smoke.py" for name in ("boltz-2", "protenix-v1")]
    files += [p for p in output.rglob("*") if p.is_file()]
    files += [case / name for name in ("evaluation-mask.json", "source.pdb", "reference-projection-NOT-CONDITIONING.pdb", "manifest.json")]
    save(output / "plan.json", {"calls": 12, "jobs": jobs, "series": "E: no MSA in either arm, structural-template information removed",
        "protenix_semantics": "per-update fixed-source callback unchanged; template module enabled but input has no templatesPath",
        "boltz_semantics": "sequence-only domain proposal, no source observations provided to model; observed-only fit/restoration retained",
        "baseline": "reuse existing six-seed source-template/no-MSA cohorts; platform-match any necessary common repairs before final comparison",
        "MSA_or_homolog_template_cells": "not executed: no frozen local homolog/MSA assets; no remote upload authorized",
        "validation_cases_used": False, "heldout_based_selection": False, "public_capability_enabled": False, "qualification": False,
        "resource_policy": "12 serial pressure-v2 calls, 1800s/call, RSS8GiB/MPS6GiB, no new checkpoints/environments",
        "files": {str(path.resolve()): digest(path) for path in files}})
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    print("Declared 12 sequence-evidence controls; no new template-arm sampling.")


def run(output):
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    plan = json.loads((output / "plan.json").read_text())
    for name, expected in plan["files"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"declared E input/code changed: {name}")
    rows = []
    for job in plan["jobs"]:
        case, work = Path(job["case"]), Path(job["work"])
        paths = [PROJECT, PROJECT / "src", PROJECT / "scripts"]
        paths += [PROJECT / "deploy" / ("protenix-v1" if job["model"] == "protenix" else "boltz-2")]
        if job["model"] == "protenix":
            paths += [RUNTIME / "engines/protenix-v1"]
        env = {**os.environ, "PYTHONPATH": os.pathsep.join(map(str, paths)), "PROTENIX_ROOT_DIR": str(RUNTIME / "assets/protenix"), "LAYERNORM_TYPE": "torch",
            "DVBFIXER_SEQUENCE_EVIDENCE_RECEIPT": str(work / "run/sequence-evidence-mask.json")}
        with (work / "orchestration.log").open("xb") as log:
            code = subprocess.run(job["command"], env=env, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT).returncode
        execution = work / "run/execution.json"
        status = json.loads(execution.read_text()) if execution.exists() else {}
        if code or not status or status.get("guard_stop_reason"):
            save(work / "failure.json", {"returncode": code, "execution": status, "not_retried": True})
            raise RuntimeError("E technical/resource failure retained; remaining jobs not launched")
        masks = json.loads((work / "run/sequence-evidence-mask.json").read_text())
        if not masks["no_structural_template_features"]:
            raise ValueError("E final consumer did not prove empty template masks")
        predictions = list((work / "run/predictions").rglob("*_model_0.pdb" if job["model"] == "boltz" else "candidate.pdb"))
        if len(predictions) != 1:
            raise ValueError("native E proposal missing/ambiguous")
        raw_report = audit(case, predictions[0], work / "audit", seed=job["seed"])
        repair_failure = None
        try:
            refine(work / "audit", work / "common-repair", platform="CPU")
        except Exception:
            repair_failure = traceback.format_exc()
            save(work / "repair-failure.json", {"error": repair_failure})
        final = work / "common-repair/refined.pdb"
        repaired = final.is_file()
        if not repaired:
            final = work / "audit/source-restored.pdb"
        request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
        assessed = assessment_adapter(request, (case / "source.pdb").read_text())(final.read_text())
        structure = _structure(final.read_text(), "E-final.pdb")
        all_l = True
        try:
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
        except Exception:
            all_l = False
        reference, coordinates = inventory(case / "reference-projection-NOT-CONDITIONING.pdb"), inventory(final)
        rmsd = float(np.sqrt(np.mean([np.sum((reference[a][2]-coordinates[a][2])**2) for a in request.generated_atoms])))
        row = {"backend": job["model"], "seed": job["seed"], "case": str(case), "arm": "no-template-no-MSA",
            "raw_geometry_passed": raw_report["assessments"]["source-restored"]["geometry_passed"],
            "native_observed_drift": raw_report["native_observed_drift_after_observed_only_fit_angstrom"],
            "repaired": repaired, "repair_failure": repair_failure, "final_geometry_passed": assessed.geometry_passed,
            "final_admission_failures": list(assessed.admission_failures), "global_all_l": all_l,
            "heldout_heavy_rmsd": rmsd, "final_path": str(final), "final_sha256": digest(final), "qualification": False}
        rows.append(row)
        (output / "progress.json").write_text(json.dumps(rows, indent=2)+"\n")
        print(json.dumps({"backend": job["model"], "seed": job["seed"], "final_geometry_passed": assessed.geometry_passed}), flush=True)
    save(output / "completed.json", {"executed": len(rows), "qualification": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    (declare if args.mode == "declare" else run)(args.output.resolve())
