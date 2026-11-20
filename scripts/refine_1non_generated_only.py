#!/usr/bin/env python3
"""One source-bound CPU refinement of the retained full-tetramer cc94 proposal.

Only generated heavy atoms and disposable auxiliary H are variables. No external
MODELLER, PDBFixer heavy rebuild, source repair, water or observed-anchor release.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
import traceback
from dataclasses import asdict
from io import StringIO
from pathlib import Path

import numpy as np
from openmm import Context, OpenMMException, Platform, VerletIntegrator, XmlSerializer
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile
from openmm.unit import dalton, nanometer
from scipy.optimize import minimize
from scripts.anchor_calibration_grid import (
    build_cell,
    cell_by_id,
    identity,
    independent_coordinate_audit,
)
from scripts.anchor_constrained_control import Objective, Trace
from scripts.anchor_constrained_coordinates import ReducedCoordinates
from scripts.anchor_refinement_guards import immutable_context_admission, returned_state_guard
from scripts.audit_anchor_calibration_grid import diagnose
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l, repair_misplaced_hydrogens
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

FORCEFIELD = ("amber19/protein.ff19SB.xml", "amber19/tip3p.xml")
PHASES = (("bonded-geometry", 1, 0., 1., 100),
          ("soft-steric-1", 1 | (1 << 30), 1., 1., 75),
          ("soft-steric-2", 1 | (1 << 30), 1., .2, 75),
          ("full-forcefield", 1 | (1 << 31), 0., .2, 1000))


def validate_inputs(request, source, proposal, record):
    fixed = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.fixed_atoms}
    generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    if len(fixed) != 5040 or len(generated) != 556 or set(source) != fixed or set(proposal) != fixed | generated:
        raise ValueError("this control requires the complete retained 1NON fixed/generated axis")
    if any(proposal[a]["line"] != source[a]["line"] for a in fixed):
        raise ValueError("raw diffusion candidate changed an observed protein record")
    if any(r["line"].startswith("HETATM") or r["element"] in {"H", "D"} for r in proposal.values()):
        raise ValueError("this control requires canonical protein-only heavy input")
    owned = {(a["chain"], a["resid"], a["icode"], a["atom"]) for a in record.added_atoms}
    if owned != generated or len(record.added_atoms) != 556:
        raise ValueError("canonical .dat does not match diffusion ownership")
    axis = target_axis(request)
    for target in request.target_sequences:
        for residue, letter in zip(axis[target.chain], target.sequence, strict=True):
            expected = {(residue.chain, residue.residue_number, residue.insertion_code, name)
                        for name in CANONICAL_HEAVY_ATOMS[_CANONICAL_ONE_TO_THREE[letter]]}
            if not expected <= proposal.keys():
                raise ValueError("input needs heavy-atom completion; this control never rebuilds heavy atoms")


def prepare_hydrogens(request, source_text, proposal_text, directory):
    directory.mkdir()
    pdb = PDBFile(StringIO(proposal_text))
    assert_all_l(pdb.topology, pdb.positions)
    original = {identity(a): np.asarray(pdb.positions[a.index].value_in_unit(nanometer)).copy()
                for a in pdb.topology.atoms()}
    np.savez_compressed(directory / "input-heavy-state.npz", positions_nm=np.array(list(original.values())))
    forcefield = ForceField(*FORCEFIELD)
    modeller = Modeller(pdb.topology, pdb.positions)
    state = random.getstate()
    random.seed(23)
    try:
        # OpenMM Modeller places auxiliary H only. This is NOT MODELLER modeling.
        modeller.addHydrogens(forcefield, pH=7.0, platform=Platform.getPlatformByName("CPU"))
    finally:
        random.setstate(state)
    repairs = repair_misplaced_hydrogens(modeller.topology, modeller.positions)
    heavy = {identity(a): np.asarray(modeller.positions[a.index].value_in_unit(nanometer))
             for a in modeller.topology.atoms() if a.element.symbol != "H"}
    if heavy.keys() != original.keys() or any(not np.array_equal(heavy[a], xyz) for a, xyz in original.items()):
        raise ValueError("auxiliary-H preparation moved/changed an input heavy atom")
    assert_all_l(modeller.topology, modeller.positions)
    base = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
    (directory / "base-system.xml").write_text(XmlSerializer.serialize(base))
    positions = np.asarray(modeller.positions.value_in_unit(nanometer)).copy()
    axis = list(modeller.topology.atoms())
    save(directory / "temporary-axis.json", {"atoms": [{"index": a.index, "identity": asdict(identity(a)),
        "resname": a.residue.name, "chain_index": a.residue.chain.index, "residue_index": a.residue.index,
        "element": a.element.symbol, "mass_dalton": float(base.getParticleMass(a.index).value_in_unit(dalton)),
        "temporary_addition": identity(a) not in original} for a in axis],
        "bonds": [[a.index, b.index] for a, b in modeller.topology.bonds()]})
    np.savez_compressed(directory / "initial-state.npz", positions_nm=positions)
    build_cell(request, source_text, proposal_text, modeller.topology, base, positions, cell_by_id("fixed"))
    actual_fixed = {identity(a) for a in axis if base.getParticleMass(a.index).value_in_unit(dalton) == 0}
    generated = {a.index for a in axis if identity(a) in set(request.generated_atoms)}
    if actual_fixed != set(request.fixed_atoms) or len(generated) != 556:
        raise ValueError("final FF mass/ownership mask differs after H insertion")
    return modeller.topology, base, positions, generated, repairs


def declare(parent, output):
    import openmm.app

    receipt_path = parent / "independent-readback/receipt.json"
    receipt = json.loads(receipt_path.read_text())
    parent_plan = json.loads((parent / "plan.json").read_text())
    if (not receipt["native_fixed_passed"] or not receipt["sampling_execution_completed"]
            or receipt["plan_sha256"] != digest(parent / "plan.json")):
        raise ValueError("refinement requires verified native scaffold preservation")
    raw = parent / "independent-readback/raw-protein-only.pdb"
    if digest(raw) != receipt["candidate_sha256"]:
        raise ValueError("retained diffusion raw bytes changed")
    output.mkdir(parents=True, exist_ok=False)
    case = output / "case"
    case.mkdir()
    for source, name in ((Path(parent_plan["source"]), "source.pdb"), (raw, "proposal.pdb"),
                         (parent / "native-sampling/research-group.json", "request.json"),
                         (parent / "independent-readback/generated.dat", "generated.dat")):
        (case / name).write_bytes(source.read_bytes())
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    if ([(t.chain, len(t.sequence)) for t in request.target_sequences] != [(c, 179) for c in "ABCD"]
            or not parent_plan.get("fixed_scaffold_inpainting")
            or digest(case / "source.pdb") != receipt["trace"]["source_sha256"]):
        raise ValueError("requires the declared complete protein-only fixed-scaffold tetramer")
    validate_inputs(request, atoms(case / "source.pdb"), atoms(case / "proposal.pdb"), DatRecord.load(case / "generated.dat"))
    files = [receipt_path, parent / "plan.json", raw, Path(__file__), Path(sys.executable)]
    files += list(case.iterdir())+list(Path("scripts").glob("*.py"))+list(Path("src/dvbfixer").rglob("*.py"))
    files += sorted((Path(openmm.app.__file__).parent / "data").rglob("*.xml"))
    files += list(Path("tests").glob("test*1non*py"))
    work = output / "execution/1NON-generated-only-seed-23"
    job = {"case_id": "1NON-generated-only", "seed": 23, "case": str(case), "snapshot": str(work / "prepared"),
        "timeout_seconds": 1800, "command": [sys.executable, "-m", "scripts.refine_1non_generated_only", "worker",
                                            "--output", str(work)]}
    save(output / "plan.json", {"protocol": "private-1non-diffusion-generated-only-reduced-xh-v1",
        "python": sys.executable, "python_sha256": digest(Path(sys.executable)), "jobs": [job], "count": 1,
        "source_code_and_input_sha256": {str(p.resolve()): digest(p) for p in files},
        "parent_diffusion_manifest_sha256": digest(parent / "plan.json"), "new_diffusion_calls": 0,
        "external_MODELLER_calls": 0, "heavy_rebuild_calls": 0, "fixed_heavy_atoms": 5040,
        "mobile_generated_heavy_atoms": 556, "auxiliary_hydrogens": "mobile, disposable, excluded from output/.dat",
        "preparation": "OpenMM Modeller seeded auxiliary H atpH7/CPU; H-only placement may optimize H, every heavy coordinate must remain exact; no heavy repair",
        "forcefield": FORCEFIELD, "nonbonded": "NoCutoff; no added solvent; no retained waters", "constraints": "FF-equilibrium covalent X-H eliminated by reduced coordinates",
        "phase_caps": [100, 75, 75, 1000], "optimizer": "same reduced-coordinate L-BFGS-B; no restart/rechart; maxls20;gtol1e-6;ftol1e-12",
        "stationarity": "independent Cartesian tangent DOF-normalized RMS <=10; not solver success alone",
        "acceptance": "native+serialized fixed custody, global L, H/steric/stereo guards, unchanged protein mask geometry and ownership",
        "resources": "serial1800s including H preparation, CPU1thread,RSS8GiB,host-pressure,10GiB evidence",
        "observed_anchor_release": False, "public_acceptance": False})
    freeze_code(output / "plan.json")


def worker(output, *, prepare=prepare_hydrogens):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("requires a fresh executor reservation")
    root = output.parent.parent
    plan = json.loads((root / "plan.json").read_text())
    verify_bindings(plan["source_code_and_input_sha256"])
    job = plan["jobs"][0]
    if Path(job["command"][-1]) != output:
        raise ValueError("worker path differs from declaration")
    case = Path(job["case"])
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    source, proposal = [(case / name).read_text() for name in ("source.pdb", "proposal.pdb")]
    row = {"case_id": job["case_id"], "new_diffusion_calls": 0, "public_acceptance": False,
           "optimization_executed": False, "Context_created": False,
           "phase_calls": [], "restoration_mask_passed": False}
    started = time.monotonic()
    try:
        topology, system, initial, generated, repairs = prepare(request, source, proposal, output / "prepared")
        row["preparation_seconds"] = time.monotonic()-started
        row["temporary_h_repairs"] = repairs
        row["heavy_preparation_exact"] = True
        (output / "system.xml").write_text(XmlSerializer.serialize(system))
        admission = immutable_context_admission(topology, system, initial)
        save(output / "admission.json", admission)
        if not admission["admitted"]:
            row["outcome_flags"] = admission["reasons"]
            return
        chart = ReducedCoordinates(topology, system, initial)
        q = chart.initial.copy()
        np.savez_compressed(output / "xh-projection.npz", before=initial, after=chart.positions(q), reduced_coordinates=q)
        integrator = VerletIntegrator(.001)
        context = Context(system, integrator, Platform.getPlatformByName("CPU"), {"Threads": "1"})
        row["Context_created"] = True
        for index, (name, groups, scale, alpha, cap) in enumerate(PHASES, 1):
            context.setParameter("soft_scale", scale)
            context.setParameter("soft_alpha", alpha)
            objective = Objective(chart, context, groups)
            directory = output / f"phase-{index}"
            directory.mkdir()
            row["optimization_executed"] = True
            result = minimize(objective.evaluate, q, jac=True, method="L-BFGS-B", callback=Trace(objective, directory),
                options={"maxiter": cap, "maxfun": 20*(cap+1), "maxls": 20, "gtol": 1.e-6, "ftol": 1.e-12})
            q = result.x.copy()
            objective.evaluate(q)
            _, energy, _, points, forces = objective.last
            audit = chart.audit(points, forces)
            np.savez_compressed(directory / "returned-state.npz", positions_nm=points, raw_forces=forces,
                                projected_forces=chart.project_forces(points, forces), reduced_coordinates=q)
            phase = {"phase": name, "cap": cap, "iterations": int(result.nit), "actual_evaluations": objective.calls,
                "scipy_status": int(result.status), "scipy_message": str(result.message), "energy_kj_mol": energy, **audit}
            save(directory / "completed.json", phase)
            row["phase_calls"].append(phase)
            assert_all_l(topology, points*nanometer)
        guard = returned_state_guard(topology, system, initial, points, chart.project_forces(points, forces), generated)
        save(output / "guard.json", {"geometry_guard_with_projected_forces": guard, "independent_tangent_audit": audit})
        flags = list(guard["reasons"])
        if not audit["stationarity_passed"]:
            flags.append("nonstationary_constrained_tangent_state")
        if not audit["fixed_exact"] or not audit["xh_equalities_satisfied"]:
            raise ValueError("native fixed/X-H equality invariant failed")
        updates = {identity(a): points[a.index]*10 for a in topology.atoms() if a.index in generated}
        candidate = rewrite_generated_coordinates(proposal, updates)
        independent_coordinate_audit(proposal, candidate, set(request.generated_atoms))
        (output / "candidate.pdb").write_text(candidate)
        assessment = assessment_adapter(request, source)(candidate)
        save(output / "assessment.json", asdict(assessment))
        diagnose(output / "candidate.pdb", output / "full-cli.json")
        if not assessment.geometry_passed or assessment.admission_failures:
            flags.append("restoration_geometry_or_admission_failure")
        DatRecord.load(case / "generated.dat").save(output / "generated.dat")
        row.update(outcome_flags=flags, restoration_mask_passed=not flags, native_fixed_exact=True,
                   candidate_sha256=digest(output / "candidate.pdb"), tangent_audit=audit)
    except (FloatingPointError, OpenMMException, ChiralityError) as exc:
        row.update(outcome_flags=["numerical_or_chirality_failure"], reason=str(exc))
        (output / "failure.txt").write_text(traceback.format_exc())
        if "objective" in locals() and objective.failed_state is not None:
            fq, fp, ff, fe = objective.failed_state
            np.savez_compressed(output / "failed-evaluation.npz", reduced_coordinates=fq, positions_nm=fp,
                                raw_forces=ff, energy_kj_mol=fe)
    except Exception as exc:
        row.update(outcome_flags=["implementation_or_input_error"], reason=str(exc))
        (output / "failure.txt").write_text(traceback.format_exc())
        raise
    finally:
        row["elapsed_seconds"] = time.monotonic()-started
        verify_bindings(plan["source_code_and_input_sha256"])
        save(output / "result.json", row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "worker", "run"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parent", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        if args.parent is None:
            parser.error("declare requires --parent")
        declare(args.parent.resolve(), args.output.resolve())
    elif args.mode == "worker":
        worker(args.output.resolve())
    else:
        from scripts.run_mac_modeller_comparison import execute
        execute(args.output.resolve())
