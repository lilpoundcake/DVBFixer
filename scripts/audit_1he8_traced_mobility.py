#!/usr/bin/env python3
"""Independent readback: geometry, bounds, rounding and actual optimizer traces."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_1he8_joint_relaxation import read_restraints
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits, check_anchor_movements
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.geometry_conditioning import BOUNDARY_WALL_K
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def verify_wall_system(path, particles, restraints, wall_limits):
    from openmm import CustomCompoundBondForce, CustomExternalForce, XmlSerializer

    system = XmlSerializer.deserialize(path.read_text())
    axis = {p["index"]: AtomIdentity(p["chain"], p["residue_number"], p["insertion_code"], p["atom_name"]) for p in particles}
    expected = {r.atom: r for r in restraints}
    particle_terms, backbone_terms = 0, 0
    for force in system.getForces():
        if isinstance(force, CustomExternalForce):
            names = [force.getPerParticleParameterName(i) for i in range(force.getNumPerParticleParameters())]
            if names != ["k", "x0", "y0", "z0", "radius"]:
                continue
            for i in range(force.getNumParticles()):
                index, values = force.getParticleParameters(i)
                atom = axis[index]
                radius = wall_limits.heavy_displacement_angstrom_max
                if atom.atom_name == "CA":
                    radius = min(radius, wall_limits.ca_displacement_angstrom_max)
                if atom not in expected or not np.allclose(values,
                    [BOUNDARY_WALL_K, *(x/10 for x in expected[atom].reference_angstrom), radius/10], rtol=0, atol=1e-12):
                    raise ValueError("actual particle wall/reference differs from declared activation limits")
                particle_terms += 1
        elif isinstance(force, CustomCompoundBondForce):
            names = [force.getPerBondParameterName(i) for i in range(force.getNumPerBondParameters())]
            if "radius" not in names:
                continue
            for i in range(force.getNumBonds()):
                indices, values = force.getBondParameters(i)
                atoms = [axis[index] for index in indices]
                if (len(atoms) != 4 or {a.atom_name for a in atoms} != {"N", "CA", "C", "O"}
                        or len({(a.chain, a.residue_number, a.insertion_code) for a in atoms}) != 1
                        or any(a not in expected for a in atoms)):
                    raise ValueError("backbone wall targets wrong particles")
                target = [BOUNDARY_WALL_K, wall_limits.backbone_rmsd_angstrom_max/10,
                          *(x/10 for a in atoms for x in expected[a].reference_angstrom)]
                if not np.allclose(values, target, rtol=0, atol=1e-12):
                    raise ValueError("actual backbone wall/reference differs from declared activation limits")
                backbone_terms += 1
    residue_count = len({(r.atom.chain, r.atom.residue_number, r.atom.insertion_code) for r in restraints})
    if particle_terms != len(restraints) or backbone_terms != residue_count:
        raise ValueError("actual wall inventory does not cover declared observed mobility")
    return {"particle_wall_terms": particle_terms, "backbone_wall_terms": backbone_terms,
            "actual_force_parameters_and_original_references_verified": True}


def audit(root, output):
    plan = json.loads((root / "plan.json").read_text())
    completed = json.loads((root / "completed.json").read_text())
    if completed["completed"] != len(plan["jobs"]) or completed["executed"] != len(plan["jobs"]):
        raise ValueError("traced mobility cohort incomplete")
    for name, expected in plan["source_code_and_input_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"declared traced-mobility code/input differs: {name}")
    rows = []
    for job in plan["jobs"]:
        work = root / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        execution = json.loads((work / "execution.json").read_text())
        if execution["returncode"] or execution["guard_stop_reason"]:
            raise ValueError("technical/resource failure is not complete scientific readback")
        result = json.loads((work / "result.json").read_text())
        if "final_sha256" not in result:
            rows.append({"case": result["case"], "arm": result["arm"], "accepted": False,
                         "refinement_failure": result.get("numerical_or_refinement_failure")})
            continue
        final_path = work / "refined.pdb"
        if digest(final_path) != result["final_sha256"]:
            raise ValueError("traced final structure changed")
        request = DiffusionRequest.from_json((Path(job["case"]) / "request.json").read_text())
        source = (Path(job["case"]) / "source.pdb").read_text()
        if digest(Path(job["case"]) / "source.pdb") != request.normalized_pdb.sha256:
            raise ValueError("original observation reference changed")
        mobility = json.loads((work / "mobility.json").read_text())
        restraints = read_restraints(mobility)
        allowed = frozenset(r.atom for r in restraints)
        limits = AnchorMovementLimits(**mobility["limits"])
        final = final_path.read_text()
        coordinates, _ = _coordinates(final)
        proposal = (work / "input.pdb").read_text()
        mobile = set(request.generated_atoms) | set(allowed)
        if final != rewrite_generated_coordinates(proposal, {a: coordinates[a] for a in mobile}):
            raise ValueError("nonmobile source records changed")
        assessment = assessment_adapter(request, source)(final, allowed)
        movements = check_anchor_movements(restraints, coordinates, limits)
        structure = _structure(final, "traced-mobility-readback.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        accepted = assessment.geometry_passed and not assessment.admission_failures and all(m.within_limits for m in movements)
        if (accepted != result["accepted"] or assessment.geometry_passed != result["geometry_passed"]
                or list(assessment.admission_failures) != result["admission_failures"]):
            raise ValueError("independent traced acceptance differs")
        particles = json.loads((work / "temporary-particle-axis.json").read_text())
        wall_limits = limits
        policy_path = work / "wall-activation-policy.json"
        if policy_path.exists():
            policy = json.loads(policy_path.read_text())
            if policy["validation_limits_unchanged"] != asdict(limits):
                raise ValueError("buffered-wall trial changed independent movement acceptance")
            wall_limits = AnchorMovementLimits(**policy["wall_activation_limits"])
            if any(v >= asdict(limits)[k] for k, v in asdict(wall_limits).items()):
                raise ValueError("buffered wall must activate strictly inside original validator")
        wall_evidence = verify_wall_system(work / "temporary-system.xml", particles, restraints, wall_limits)
        atom_indices = {AtomIdentity(p["chain"], p["residue_number"], p["insertion_code"], p["atom_name"]): p["index"] for p in particles}
        phases = []
        for index, call in enumerate(result["minimizer_calls"]):
            phase = work / f"minimize-{index+1:02d}"
            callbacks_path = phase / "callbacks.jsonl"
            callbacks = [json.loads(line) for line in callbacks_path.read_text().splitlines()] if callbacks_path.exists() else []
            if len(callbacks) != call["callback_count"] or any(c["requested_termination"] for c in callbacks):
                raise ValueError("passive callback trace incomplete or requested early termination")
            if callbacks and callbacks[-1]["optimizer_iteration_index"] != call["last_optimizer_iteration_index"]:
                raise ValueError("last actual optimizer iteration differs")
            phases.append({"phase_call": index+1, "callback_count": len(callbacks),
                "configured_max_iterations": call["maxIterations"],
                "last_iteration_index": call["last_optimizer_iteration_index"],
                "callback_count_reached_configured_cap": len(callbacks) >= call["maxIterations"],
                "last_reported_objective_gradient_component_rms": callbacks[-1]["objective_gradient_component_rms"] if callbacks else None,
                "termination_reason_exposed": False})
        last_phase = work / f"minimize-{len(phases):02d}" / "returned-state.npz"
        with np.load(last_phase, allow_pickle=False) as archive:
            positions = np.asarray(archive["positions_nm"])
            raw = {a: positions[i]*10 for a, i in atom_indices.items()}
        raw_movements = check_anchor_movements(restraints, raw, limits)
        raw_by_residue = {m.residue: m for m in raw_movements}
        comparison = []
        for m in movements:
            original = raw_by_residue[m.residue]
            excess_nm = max(0., original.backbone_rmsd_angstrom-wall_limits.backbone_rmsd_angstrom_max)/10
            comparison.append({"residue": asdict(m.residue), "rounded": asdict(m), "solver_precision": asdict(original),
                "rounding_only_violation": not m.within_limits and original.within_limits,
                "BB_wall_penalty_kj_mol_at_solver_state": .5*BOUNDARY_WALL_K*excess_nm**2,
                "wall_is_finite_penalty_not_exact_constraint": True})
        rows.append({"case": result["case"], "arm": result["arm"], "accepted": accepted,
            "geometry_passed": assessment.geometry_passed, "assessment": asdict(assessment),
            "movement_violations_rounded": sum(not m.within_limits for m in movements),
            "movement_violations_solver_precision": sum(not m.within_limits for m in raw_movements),
            "anchor_residues_mobile": len(movements), "generated_atoms_mobile": len(request.generated_atoms),
            "movement_comparison": comparison, "phase_traces": phases,
            "nonmobile_records_exact": True, "global_all_l": True, "final_sha256": digest(final_path),
            "final_path": str(final_path), "original_validation_limits": asdict(limits),
            "actual_wall_activation_limits": asdict(wall_limits), "wall_force_readback": wall_evidence,
            "RMSD_reference_used": False, "not_whole_entry_cofactor_qualification": True})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "independent-readback.json", {"rows": rows, "auditor_sha256": digest(Path(__file__)),
        "plan_sha256": digest(root / "plan.json"), "all_trials_retained": True, "policy_loosened": False})
    print(json.dumps({"independently_read": len(rows), "accepted": sum(r["accepted"] for r in rows),
        "geometry_passes": sum(r.get("geometry_passed", False) for r in rows),
        "all_nonmobile_records_and_L_checked": all(r.get("nonmobile_records_exact") and r.get("global_all_l") for r in rows)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve())
