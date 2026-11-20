#!/usr/bin/env python3
"""Declare and execute six saved-state controls; zero new optimization calls."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from openmm import XmlSerializer
from scripts.anchor_calibration_grid import identity, read_temporary
from scripts.anchor_refinement_guards import (
    REVISION,
    immutable_context_admission,
    returned_state_guard,
    xh_equilibria,
)
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest

CONTROLS = ("1KTZ-fixed", "6Z38-fixed", "1NON-fixed", "1NON-w1-e1-k10000",
            "2BE1-fixed", "2BE1-w1-e3-k10000")


def declare(grid, output):
    output.mkdir(parents=True, exist_ok=False)
    plan_path = grid / "plan.json"
    parent = json.loads(plan_path.read_text())
    verify_bindings(parent["source_code_and_input_sha256"], (plan_path,))
    files = [plan_path, Path(sys.executable)]
    files += list(Path("src/dvbfixer").rglob("*.py"))
    files += list(Path("scripts").glob("*.py"))
    files += [Path("tests/test_anchor_refinement_guards.py")]
    jobs = []
    for cell in CONTROLS:
        job = next(j for j in parent["jobs"] if j["case_id"] == cell)
        work = grid / f"execution/{cell}-seed-{job['seed']}"
        snapshot = Path(job["command"][job["command"].index("--snapshot")+1])
        files += [Path(job["case"]) / "request.json", work / "system.xml", work / "result.json",
                  work / "phase-4/returned-state.npz"]
        files += [snapshot / name for name in ("temporary-axis.json", "base-system.xml", "initial-state.npz")]
        jobs.append({"case_id": cell, "seed": job["seed"], "case": job["case"], "snapshot": str(snapshot),
            "saved_work": str(work), "timeout_seconds": 180,
            "command": [sys.executable, "-m", "scripts.check_saved_anchor_guards", "worker",
                        "--output", str(output / f"execution/{cell}-seed-{job['seed']}")]})
    save(output / "plan.json", {"protocol": REVISION, "count": len(jobs), "jobs": jobs,
        "python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p.resolve()): digest(p) for p in files},
        "prospective_reconstruction_cohort": False, "controls": "retrospective saved-state witnesses, not independent qualification",
        "new_optimization_calls": 0, "new_h_placement_calls": 0, "new_source_repairs": 0,
        "admission": "Python-v2 severe immutable heavy contacts; complete supplied/template connectivity; waters retained",
        "returned_state": "finite, exact massless fixation, active Cartesian component RMS <=10, FF-based covalent X-H deviation <20%, existing 0.05nm H coincidence floor, global L, generated stereo",
        "geometry_contract": "existing full geometry and movement gates remain additional requirements; guard pass is not acceptance",
        "hard_constraint_status": "blocked on massless X-H parents; no penalty substitution or release of observations",
        "resources": "serial 180s/control, existing 8GiB RSS/host pressure/10GiB evidence guards"})
    freeze_code(output / "plan.json")


def worker(output):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("guard control requires fresh executor reservation")
    plan_path = output.parent.parent / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"])
    job = next(j for j in plan["jobs"] if Path(j["command"][-1]).resolve() == output.resolve())
    snapshot, work = Path(job["snapshot"]), Path(job["saved_work"])
    topology, _, initial = read_temporary(snapshot)
    system = XmlSerializer.deserialize((work / "system.xml").read_text())
    request = DiffusionRequest.from_json((Path(job["case"]) / "request.json").read_text())
    generated = {a.index for a in topology.atoms() if identity(a) in set(request.generated_atoms)}
    if len(generated) != len(request.generated_atoms):
        raise ValueError("generated atoms missing from saved FF axis")
    before = XmlSerializer.serialize(system)
    admission = immutable_context_admission(topology, system, initial)
    constraints = xh_equilibria(topology, system)
    with np.load(work / "phase-4/returned-state.npz", allow_pickle=False) as data:
        positions, forces = data["positions_nm"].copy(), data["active_forces"].copy()
    # Read back failures too, even if the new admission gate would refuse before
    # optimization. This is a retrospective diagnosis, not an admission bypass.
    returned = returned_state_guard(topology, system, initial, positions, forces, generated)
    if before != XmlSerializer.serialize(system):
        raise ValueError("guard mutated the saved system")
    save(output / "result.json", {"case_id": job["case_id"], "entry": Path(job["case"]).name,
        "admission": admission, "returned_state": returned, "xh_constraint_inventory": constraints,
        "new_optimization_calls": 0, "Context_created": False, "system_unchanged": True,
        "historical_result": json.loads((work / "result.json").read_text())["completed"],
        "scientific_acceptance": False, "interpretation": "guards only; does not replace geometry/movement validation"})


def summarize(output):
    plan = json.loads((output / "plan.json").read_text())
    rows = []
    for job in plan["jobs"]:
        path = output / f"execution/{job['case_id']}-seed-{job['seed']}/result.json"
        if not path.is_file():
            rows.append({"case_id": job["case_id"], "state": "no_guard_result"})
            continue
        row = json.loads(path.read_text())
        rows.append({"case_id": job["case_id"], "state": "reported", "admission": row["admission"]["reasons"],
            "returned_state": row["returned_state"]["reasons"],
            "active_force_rms": row["returned_state"]["active_cartesian_component_rms"],
            "xh_bonds": len(row["xh_constraint_inventory"]),
            "massless_xh_blockers": sum(not b["openmm_hard_constraint_supported"] for b in row["xh_constraint_inventory"])})
    verify_bindings(plan["source_code_and_input_sha256"])
    save(output / "summary.json", {"plan_sha256": digest(output / "plan.json"), "controls": rows,
        "all_controls_reported": all(r["state"] == "reported" for r in rows), "new_optimization_calls": 0,
        "bound_inputs_and_code_unchanged": True, "production_changes": False,
        "constrained_refiner": "not implemented; actual OpenMM massless-particle constraint blocker"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grid", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        if args.grid is None:
            parser.error("declare requires --grid")
        declare(args.grid.resolve(), args.output.resolve())
    elif args.mode == "worker":
        worker(args.output.resolve())
    else:
        from scripts.run_mac_modeller_comparison import execute
        try:
            execute(args.output.resolve())
        finally:
            summarize(args.output.resolve())
