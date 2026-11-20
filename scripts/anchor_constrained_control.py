#!/usr/bin/env python3
"""Private bounded reduced-coordinate controls, not a new calibration grid."""
from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np
from openmm import Context, OpenMMException, Platform, VerletIntegrator, XmlSerializer
from openmm.unit import kilojoule_per_mole, nanometer
from scipy.optimize import minimize
from scripts.anchor_calibration_grid import identity, independent_coordinate_audit, read_temporary
from scripts.anchor_constrained_coordinates import REVISION, ReducedCoordinates
from scripts.anchor_refinement_guards import immutable_context_admission, returned_state_guard
from scripts.audit_anchor_calibration_grid import diagnose, evaluation_view
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

CONTROLS = ("1KTZ-fixed", "1NON-fixed", "2BE1-fixed")
PHASES = (("bonded-geometry", 1, 0., 1., 100),
          ("soft-steric-1", 1 | (1 << 30), 1., 1., 75),
          ("soft-steric-2", 1 | (1 << 30), 1., .2, 75),
          ("full-forcefield", 1 | (1 << 31), 0., .2, 250))


class Objective:
    def __init__(self, chart, context, groups):
        self.chart, self.context, self.groups = chart, context, groups
        self.calls = 0
        self.last = None
        self.failed_state = None

    def evaluate(self, q):
        if self.last is not None and np.array_equal(q, self.last[0]):
            return self.last[1], self.last[2]
        points = self.chart.positions(q)
        # Enforce equalities/fixation at EVERY evaluation, before FF evaluation.
        errors = np.abs(np.linalg.norm(points[self.chart.h]-points[self.chart.parent], axis=1)-self.chart.length)
        if (not np.array_equal(points[~self.chart.mobile], self.chart.reference[~self.chart.mobile])
                or (errors > 1.e-12).any()):
            raise ValueError("reduced-coordinate equality/fixation invariant failed")
        self.context.setPositions(points*nanometer)
        state = self.context.getState(getEnergy=True, getForces=True, groups=self.groups)
        energy = float(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
        forces = state.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer)
        self.calls += 1
        if not math.isfinite(energy) or not np.isfinite(forces).all():
            self.failed_state = (q.copy(), points, forces, energy)
            raise FloatingPointError("nonfinite reduced-coordinate energy/forces")
        gradient = self.chart.pullback(q, -forces)
        self.last = (q.copy(), energy, gradient, points, forces)
        return energy, gradient


class Trace:
    def __init__(self, objective, directory):
        self.objective, self.directory = objective, directory
        self.calls = 0

    def __call__(self, x):
        objective = self.objective
        objective.evaluate(x)
        self.calls += 1
        audit = objective.chart.audit(objective.last[3], objective.last[4])
        with (self.directory / "callbacks.jsonl").open("a") as stream:
            stream.write(json.dumps({"callback": self.calls, "energy_kj_mol": objective.last[1], **audit}, allow_nan=False)+"\n")
        if self.calls <= 5 or self.calls % 50 == 0:
            np.savez_compressed(self.directory / f"callback-{self.calls:05d}.npz",
                positions_nm=objective.last[3], raw_forces=objective.last[4], reduced_coordinates=x)


def declare(grid, output):
    output.mkdir(parents=True, exist_ok=False)
    parent_path = grid / "plan.json"
    parent = json.loads(parent_path.read_text())
    verify_bindings(parent["source_code_and_input_sha256"], (parent_path,))
    files = [parent_path, Path(sys.executable)]
    files += list(Path("scripts").glob("*.py"))+list(Path("src/dvbfixer").rglob("*.py"))
    files += [Path("tests/test_anchor_constrained_coordinates.py")]
    jobs = []
    for cell in CONTROLS:
        job = next(j for j in parent["jobs"] if j["case_id"] == cell)
        work = grid / f"execution/{cell}-seed-{job['seed']}"
        snapshot = Path(job["command"][job["command"].index("--snapshot")+1])
        files += list(Path(job["case"]).glob("*"))
        files += [work / "system.xml", work / "result.json"]
        if (work / "generated.dat").is_file():
            files.append(work / "generated.dat")
        files += [snapshot / name for name in ("temporary-axis.json", "base-system.xml", "initial-state.npz")]
        jobs.append({"case_id": cell, "case": job["case"], "seed": job["seed"], "snapshot": str(snapshot),
            "saved_work": str(work), "timeout_seconds": job["timeout_seconds"],
            "command": [sys.executable, "-m", "scripts.anchor_constrained_control", "worker", "--output",
                        str(output / f"execution/{cell}-seed-{job['seed']}")]})
    save(output / "plan.json", {"protocol": REVISION, "jobs": jobs, "count": len(jobs),
        "python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p.resolve()): digest(p) for p in files if p.is_file()},
        "scope": "two numerical implementation smokes + 2BE1 admission refusal; not matched-budget calibration or qualification",
        "phase_caps": [100, 75, 75, 250], "optimizer": "SciPy L-BFGS-B; gtol1e-6,ftol1e-12,maxls20,maxfun20*(cap+1); no restarts/rechart",
        "stationarity": "independent Cartesian tangent projected force norm / sqrt(3*Nmobile-Nxh) <=10 kJ/mol/nm",
        "xh": "exact coordinate elimination at FF covalent equilibrium; numeric equality audit1e-12nm; original temporary H directions, disclosed radial projection",
        "chart_limit": "local stereographic charts exclude their initial antipode; no automatic rechart; independent Cartesian tangent stationarity prevents chart-gradient false passes",
        "observation_policy": "fixed coordinates excluded from variables; original bytes/masks/context/ANISOU retained",
        "chemistry": "unchanged captured ff19SB/TIP3P/NoCutoff; temporary H only; no added solvent; no extra LJ or penalty forces",
        "platform": "OpenMM CPU, one thread", "observed_mobile_window": 0,
        "new_proposals": 0, "new_h_placement": 0, "restarts": 0,
        "resources": "existing serial1800s less historical shared preparation,RSS8GiB,host pressure,10GiB evidence",
        "public_acceptance": False})
    freeze_code(output / "plan.json")


def worker(output):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("control requires fresh guarded reservation")
    plan = json.loads((output.parent.parent / "plan.json").read_text())
    verify_bindings(plan["source_code_and_input_sha256"])
    job = next(j for j in plan["jobs"] if Path(j["command"][-1]).resolve() == output.resolve())
    case, snapshot, work = Path(job["case"]), Path(job["snapshot"]), Path(job["saved_work"])
    topology, _, initial = read_temporary(snapshot)
    system = XmlSerializer.deserialize((work / "system.xml").read_text())
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    axis = list(topology.atoms())
    metadata = json.loads((snapshot / "temporary-axis.json").read_text())["atoms"]
    if any(a.element.symbol == "H" and not metadata[a.index]["temporary_addition"] for a in axis):
        raise ValueError("observed H cannot be radially projected")
    row = {"entry": case.name, "case_id": job["case_id"], "public_acceptance": False,
           "Context_created": False, "phase_calls": [], "restoration_mask_passed": False}
    try:
        admission = immutable_context_admission(topology, system, initial)
        save(output / "admission.json", admission)
        if not admission["admitted"]:
            row.update(outcome_flags=admission["reasons"], optimization_executed=False)
            return
        chart = ReducedCoordinates(topology, system, initial)
        q = chart.initial.copy()
        projected_initial = chart.positions(q)
        np.savez_compressed(output / "initial-projection.npz", original_positions_nm=initial,
                            projected_positions_nm=projected_initial, reduced_coordinates=q)
        h_motion = np.linalg.norm(projected_initial[chart.h]-initial[chart.h], axis=1)
        save(output / "projection.json", {"max_temporary_h_radial_shift_angstrom": float(h_motion.max()*10),
            "heavy_coordinates_exact": np.array_equal(initial[chart.heavy], projected_initial[chart.heavy]),
            "fixed_exact": np.array_equal(initial[~chart.mobile], projected_initial[~chart.mobile]),
            "xh_inventory": chart.bonds, "source_repair": False})
        (output / "system.xml").write_text(XmlSerializer.serialize(system))
        integrator = VerletIntegrator(.001)
        context = Context(system, integrator, Platform.getPlatformByName("CPU"), {"Threads": "1"})
        row["Context_created"] = True
        for index, (name, groups, scale, alpha, cap) in enumerate(PHASES, 1):
            context.setParameter("soft_scale", scale)
            context.setParameter("soft_alpha", alpha)
            objective = Objective(chart, context, groups)
            directory = output / f"phase-{index}"
            directory.mkdir()
            callback = Trace(objective, directory)
            row["optimization_executed"] = True
            result = minimize(objective.evaluate, q, jac=True, method="L-BFGS-B", callback=callback,
                options={"maxiter": cap, "maxfun": 20*(cap+1), "maxls": 20, "gtol": 1.e-6, "ftol": 1.e-12})
            q = result.x.copy()
            objective.evaluate(q)
            _, energy, _, points, forces = objective.last
            audit = chart.audit(points, forces)
            np.savez_compressed(directory / "returned-state.npz", positions_nm=points, raw_forces=forces,
                                projected_forces=chart.project_forces(points, forces), reduced_coordinates=q)
            phase = {"phase": name, "cap": cap, "iterations": int(result.nit),
                "scipy_nfev": int(result.nfev), "actual_evaluations": objective.calls,
                "scipy_status": int(result.status), "scipy_message": str(result.message),
                "scipy_success": bool(result.success), "energy_kj_mol": energy, **audit}
            save(directory / "completed.json", phase)
            row["phase_calls"].append(phase)
            assert_all_l(topology, points*nanometer)
        projected = chart.project_forces(points, forces)
        generated = {a.index for a in axis if identity(a) in set(request.generated_atoms)}
        guard = returned_state_guard(topology, system, initial, points, projected, generated)
        save(output / "guard.json", {"v1_geometry_guard_with_projected_force_input": guard,
                                     "independent_reduced_dof_audit": audit})
        flags = list(guard["reasons"])
        if not audit["stationarity_passed"]:
            flags.append("nonstationary_constrained_tangent_state")
        if not audit["fixed_exact"] or not audit["xh_equalities_satisfied"]:
            raise ValueError("returned reduced-coordinate invariants failed")
        proposal = (case / "proposal.pdb").read_text()
        updates = {identity(a): points[a.index]*10 for a in axis if identity(a) in set(request.generated_atoms)}
        candidate = rewrite_generated_coordinates(proposal, updates)
        independent_coordinate_audit(proposal, candidate, set(request.generated_atoms))
        (output / "candidate.pdb").write_text(candidate)
        view = evaluation_view(request, (case / "source.pdb").read_text(), output / "candidate.pdb", set())
        (output / "protein-evaluation-view.pdb").write_text(view)
        assessment = assessment_adapter(request, (case / "source.pdb").read_text())(view)
        save(output / "assessment.json", asdict(assessment))
        diagnose(output / "candidate.pdb", output / "full-context-cli.json")
        if not assessment.geometry_passed or assessment.admission_failures:
            flags.append("restoration_geometry_or_admission_failure")
        # Reuse canonical ownership only, never add temporary H/OXT to the .dat.
        record = DatRecord.load(work / "generated.dat")
        owned = {(a["chain"], a["resid"], a["icode"], a["atom"]) for a in record.added_atoms}
        expected = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
        if owned != expected or len(record.added_atoms) != len(expected):
            raise ValueError("historical generated-only ownership mismatch")
        record.save(output / "generated.dat")
        row.update(outcome_flags=flags, optimization_executed=True, native_fixed_exact=True,
            candidate_sha256=digest(output / "candidate.pdb"),
            restoration_mask_passed=not flags, independent_tangent_audit=audit,
            diagnostic_scope="mask geometry + full-context CLI; deferred OXT remains separate")
    except (FloatingPointError, OpenMMException, ChiralityError) as exc:
        row.update(outcome_flags=["numerical_or_chirality_failure"], reason=str(exc),
                   optimization_executed=row.get("optimization_executed", False))
        if "objective" in locals() and objective.failed_state is not None:
            failed_q, failed_points, failed_forces, failed_energy = objective.failed_state
            np.savez_compressed(output / "failed-evaluation.npz", reduced_coordinates=failed_q,
                positions_nm=failed_points, raw_forces=failed_forces, potential_energy_kj_mol=failed_energy)
        (output / "failure.txt").write_text(traceback.format_exc())
    finally:
        # Input/code revalidation applies to failures and admission refusals too.
        verify_bindings(plan["source_code_and_input_sha256"])
        save(output / "result.json", row)


def summarize(output):
    plan = json.loads((output / "plan.json").read_text())
    rows = []
    for job in plan["jobs"]:
        directory = output / f"execution/{job['case_id']}-seed-{job['seed']}"
        path = directory / "result.json"
        rows.append(json.loads(path.read_text()) if path.exists() else {"case_id": job["case_id"], "result_available": False})
    save(output / "summary.json", {"plan_sha256": digest(output / "plan.json"), "cases": rows,
        "all_results_available": all("outcome_flags" in r for r in rows), "production_defaults_changed": False,
        "qualification": False, "scope": "bounded implementation controls; shorter full phase, not paired calibration"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "worker", "run"))
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
