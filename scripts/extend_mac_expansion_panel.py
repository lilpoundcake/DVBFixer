#!/usr/bin/env python3
"""Complete the two private expansion masks with the five remaining dev seeds.

Thirty new calls, unchanged source/context and sampler settings; six seed7
pilots are reused, never resampled. Validation seeds remain untouched.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.extend_reference_reconstruction_panel import declare as declare_extension
from scripts.extend_reference_reconstruction_panel import run as run_extension
from scripts.mac_reconstruction_comparison import digest, save
from scripts.refine_mac_reconstruction_smoke import run as refine


def declare(pilot: Path, output: Path):
    previous = json.loads((pilot / "final-results.json").read_text())
    if len(previous) != 6:
        raise ValueError("all six seed7 outcomes must exist, including any failures")
    output.mkdir(parents=True, exist_ok=False)
    # Explicit byte-identical view for the existing extension declarer's path
    # convention. No historical directory or source plan is rewritten.
    view = output / "pilot-view/reference-pilot"
    view.mkdir(parents=True)
    (view / "plan.json").write_bytes((pilot / "runs/plan.json").read_bytes())
    declare_extension(output / "pilot-view", output / "runs")
    plan = json.loads((output / "runs/plan.json").read_text())
    if plan["calls"] != 30:
        raise ValueError("exactly thirty remaining-seed calls permitted")
    plan.update(scope="terminal A1-4 and atomic pair A10-13/A46-49 in unchanged 1UBQ A1-72",
        same_seed7_protocol=True, reference_training_disjointness="not established",
        pilot_results_sha256=digest(pilot / "final-results.json"),
        common_repair={"platform": "CPU", "threads": 1, "seed": 7, "iterations": 1000, "restarts": 0,
            "preconditioning": False, "observed_mobility": False, "all_generated_mobile": True, "clean_not_rebuilt": True},
        budget="30 serial ML calls; 1800s per call; pressure-v2 RSS8GiB/MPS6GiB; no new assets or environments")
    for name in ("extend_mac_expansion_panel.py", "mac_reconstruction_comparison.py", "refine_mac_reconstruction_smoke.py"):
        path = Path(__file__).with_name(name).resolve()
        plan["files"][str(path)] = digest(path)
    (output / "runs/plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    save(output / "declaration.json", {"plan_sha256": digest(output / "runs/plan.json"), "calls": 30,
        "seed7_reused_from": str(pilot), "reserved_validation_seeds_used": False, "qualification": False})
    print(json.dumps({"calls": 30, "plan": str(output / "runs/plan.json")}, indent=2))


def run(output: Path):
    import fcntl

    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    declared = json.loads((output / "declaration.json").read_text())
    if digest(output / "runs/plan.json") != declared["plan_sha256"]:
        raise ValueError("declared extension plan changed")
    run_extension(output / "pilot-view", output / "runs")
    plan = json.loads((output / "runs/plan.json").read_text())
    rows = []
    project = Path(__file__).resolve().parents[1]
    with (project / ".artifacts/mac-reconstruction-2026-11-11/ml.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for job in plan["jobs"]:
            case, work = Path(job["case"]), Path(job["work"])
            if not (work / "audit/audit.json").is_file():
                rows.append({"case": case.name, "backend": job["arm"], "seed": job["seed"],
                    "status": "technical failure; no auditable coordinates, no repair"})
                continue
            refine(work / "audit", work / "common-repair", platform="CPU")
            final = work / "common-repair/refined.pdb"
            repaired = final.is_file()
            if not repaired:
                final = work / "audit/source-restored.pdb"
            request = DiffusionRequest.from_json((work / "audit/evaluation-mask.json").read_text())
            assessed = assessment_adapter(request, (case / "source.pdb").read_text())(final.read_text())
            structure = _structure(final.read_text(), "extension-final.pdb")
            global_l = True
            try:
                assert_all_l(structure.pdb.topology, structure.pdb.positions)
            except Exception:
                global_l = False
            reference, coords = inventory(case / "reference-projection-NOT-CONDITIONING.pdb"), inventory(final)
            rmsd = float(np.sqrt(np.mean([np.sum((reference[a][2]-coords[a][2])**2) for a in request.generated_atoms])))
            rows.append({"case": case.name, "backend": job["arm"], "seed": job["seed"], "repaired": repaired,
                "geometry_passed": assessed.geometry_passed, "admission_failures": list(assessed.admission_failures),
                "global_all_l": global_l, "heldout_heavy_rmsd": rmsd, "final_path": str(final), "sha256": digest(final),
                "assessment": asdict(assessed), "qualification": False})
            (output / "final-progress.json").write_text(json.dumps(rows, indent=2)+"\n")
    save(output / "final-results.json", rows)
    print(json.dumps({"executed": len(plan["jobs"]), "final_outcomes": len(rows), "qualified": False}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pilot", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        if args.pilot is None:
            parser.error("declare requires --pilot")
        declare(args.pilot.resolve(), args.output.resolve())
    else:
        run(args.output.resolve())
