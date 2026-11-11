#!/usr/bin/env python3
"""Six prospective seed7 transport/science pilots, not six-seed qualification."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save
from scripts.prepare_mac_expansion_inputs import prepare
from scripts.refine_mac_reconstruction_smoke import run as refine
from scripts.run_reference_reconstruction_pilot import declare as declare_pilots
from scripts.run_reference_reconstruction_pilot import run as run_pilots


def declare(preparation: Path, output: Path):
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "inputs"
    inputs.mkdir()
    rows = json.loads((preparation / "results.json").read_text())
    admitted = []
    for row in rows:
        if row["mask_admitted"]:
            prepare(preparation / row["case"], inputs / row["case"])
            admitted.append({"case_id": row["case"], "ready": True})
    save(inputs / "preparation-results.json", admitted)
    runs = output / "runs"
    declare_pilots(inputs, runs)
    plan = json.loads((runs / "plan.json").read_text())
    if len(plan["jobs"]) != 6:
        raise ValueError("this manifest permits exactly two masks x three backends x seed7")
    plan.update(purpose="prospective N-terminal and two-gap atomic-group seed7 private transport/science pilots",
        preparation_results_sha256=digest(preparation / "results.json"),
        group_policy="both pair members come from the same model call; no cross-run mixing",
        context="full declared 1UBQ A1-72 projection; no added cropping",
        stages="native -> source restoration -> necessary-only common generated-only CPU repair",
        no_new_model_assets_or_environment_changes=True,
        total_calls=6, resource_policy="pressure-v2, serial; 1800s/call; RSS8GiB; MPS6GiB; 10GiB evidence cap",
        baseline_status="30 retrospective MODELLER results exist; new masks need separate matched baselines",
        interface_case="refused during preparation due B2 partial occupancy; not replaced or counted as sampled")
    for filename in ("run_mac_expansion_pilots.py", "prepare_mac_expansion_inputs.py", "prepare_mac_comparison_expansion.py", "mac_reconstruction_comparison.py"):
        path = Path(__file__).with_name(filename).resolve()
        plan["file_sha256"][str(path)] = digest(path)
    (runs / "plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    save(output / "declaration.json", {"runs_plan_sha256": digest(runs / "plan.json"), "calls": 6,
        "reserved_validation_seeds_used": False, "qualification": False})
    print(json.dumps({"calls": 6, "plan": str(runs / "plan.json")}, indent=2))


def run(output: Path):
    import fcntl

    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    declaration = json.loads((output / "declaration.json").read_text())
    if digest(output / "runs/plan.json") != declaration["runs_plan_sha256"]:
        raise ValueError("declared pilot plan changed")
    if sum(p.stat().st_size for p in output.rglob("*") if p.is_file()) > 10*1024**3:
        raise ValueError("pilot evidence cap exceeded")
    run_pilots(output / "runs")
    plan = json.loads((output / "runs/plan.json").read_text())
    rows = []
    project = Path(__file__).resolve().parents[1]
    with (project / ".artifacts/mac-reconstruction-2026-11-11/ml.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for job in plan["jobs"]:
            case, work = Path(job["case"]), Path(job["work"])
            if not (work / "audit/audit.json").exists():
                rows.append({"case": case.name, "backend": job["arm"], "seed": 7, "status": "no completed audit; retained technical failure"})
                continue
            refine(work / "audit", work / "common-repair", platform="CPU")
            final = work / "common-repair/refined.pdb"
            repaired = final.is_file()
            if not repaired:
                final = work / "audit/source-restored.pdb"
            request = DiffusionRequest.from_json((work / "audit/evaluation-mask.json").read_text())
            source = (case / "source.pdb").read_text()
            assessed = assessment_adapter(request, source)(final.read_text())
            structure = _structure(final.read_text(), "pilot-final.pdb")
            global_l = True
            try:
                assert_all_l(structure.pdb.topology, structure.pdb.positions)
            except Exception:
                global_l = False
            native = inventory(case / "reference-projection-NOT-CONDITIONING.pdb")
            proposed = inventory(final)
            generated = request.generated_atoms
            rmsd = float(np.sqrt(np.mean([np.sum((native[a][2]-proposed[a][2])**2) for a in generated])))
            rows.append({"case": case.name, "backend": job["arm"], "seed": 7, "repaired": repaired,
                "geometry_passed": assessed.geometry_passed, "admission_failures": list(assessed.admission_failures),
                "global_all_l": global_l, "heldout_heavy_rmsd": rmsd, "final_path": str(final), "sha256": digest(final),
                "assessment": asdict(assessed), "qualification": False})
    save(output / "final-results.json", rows)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--preparation", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        if args.preparation is None:
            parser.error("declare requires --preparation")
        declare(args.preparation.resolve(), args.output.resolve())
    else:
        run(args.output.resolve())
