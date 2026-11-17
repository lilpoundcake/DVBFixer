#!/usr/bin/env python3
"""Eight frozen CPU mobility controls on two retained 1HE8 proposals, no ML."""
from __future__ import annotations

import argparse
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_1he8_joint_relaxation import CASES, PROJECT
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.run_mac_modeller_comparison import execute
from scripts.trace_private_minimization import PassiveReporter

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion import boundary_refinement as boundary
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.anchor_policy import (
    AnchorMovementLimits,
    build_anchor_restraints,
    check_anchor_movements,
)
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.spatial_relaxation import spatial_mobility
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

ARMS = ("fixed", "direct", "sequence2", "spatial5-padding1")


def anchor_residues(request, source, proposal, arm):
    if arm == "fixed":
        return ()
    direct = {r for gap in request.gaps for r in (gap.left_anchor, gap.right_anchor) if r is not None}
    if arm == "direct":
        return tuple(sorted(direct))
    if arm == "spatial5-padding1":
        return spatial_mobility(request, source, proposal).observed_residues
    if arm != "sequence2":
        raise ValueError("undeclared mobility arm")
    observed = {r for placement in request.sequence_placements for r in placement.observed_residues}
    axis = target_axis(request)
    selected = set(direct)
    for residue in direct:
        chain = axis[residue.chain]
        index = chain.index(residue)
        selected.update(r for r in chain[max(0, index-2):index+3] if r in observed)
    return tuple(sorted(selected))


def worker(index, arm, output):
    from openmm import XmlSerializer
    from openmm.unit import dalton, kilojoule_per_mole, nanometer

    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("mobility trial requires fresh executor reservation")
    case = CASES[index]
    root = PROJECT / case["root"]
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    source = (root / "source.pdb").read_text()
    proposal = (root / case["input"]).read_text()
    if digest(root / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("retained 1HE8 source changed")
    restraints = build_anchor_restraints(source, anchor_residues(request, source, proposal, arm),
        observed_atoms=request.fixed_atoms, base_k_kj_mol_nm2=10000., bfactor_source="disabled")
    allowed = frozenset(r.atom for r in restraints)
    limits = AnchorMovementLimits()
    save(output / "mobility.json", {"arm": arm, "limits": asdict(limits), "restraints": [asdict(r) for r in restraints],
        "all_generated_mobile": len(request.generated_atoms), "uniform_restraint_strength": 10000.,
        "source_reference_sha256": request.normalized_pdb.sha256, "new_neural_sampling": False})
    (output / "input.pdb").write_text(proposal)
    assess = assessment_adapter(request, source)
    save(output / "initial-assessment.json", asdict(assess(proposal)))
    calls = []
    original = boundary.Simulation

    class TracedSimulation(original):
        def minimizeEnergy(self, *args, **kwargs):
            if kwargs.get("reporter") is not None:
                raise ValueError("do not replace an existing minimizer reporter")
            phase_dir = output / f"minimize-{len(calls)+1:02d}"
            reporter = PassiveReporter(phase_dir)
            record = {"call": len(calls)+1, "integration_force_groups": self.integrator.getIntegrationForceGroups(),
                "parameters": dict(self.context.getParameters()), "maxIterations": kwargs.get("maxIterations"),
                "tolerance": str(kwargs.get("tolerance")), "passive_callback": True}
            calls.append(record)
            save(phase_dir / "declaration.json", record)
            if len(calls) == 1:
                (output / "temporary-system.xml").write_text(XmlSerializer.serialize(self.system))
                save(output / "temporary-particle-axis.json", [{"index": a.index, "chain": a.residue.chain.id,
                    "residue_number": str(a.residue.id), "insertion_code": a.residue.insertionCode.strip(),
                    "atom_name": a.name, "mass_dalton": float(self.system.getParticleMass(a.index).value_in_unit(dalton))}
                    for a in self.topology.atoms()])
            super().minimizeEnergy(*args, **{**kwargs, "reporter": reporter})
            state = self.context.getState(getPositions=True, getForces=True, getEnergy=True,
                groups=self.integrator.getIntegrationForceGroups())
            force = state.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer)
            with (phase_dir / "returned-state.npz").open("xb") as stream:
                np.savez_compressed(stream, positions_nm=state.getPositions(asNumpy=True).value_in_unit(nanometer),
                                    active_force_kj_mol_nm=force)
            record.update(callback_count=reporter.calls, last_optimizer_iteration_index=reporter.last_iteration,
                returned_energy_kj_mol=float(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)),
                engine_termination_reason="not exposed by API; do not infer convergence from return alone")
            save(phase_dir / "completed.json", record)

    boundary.Simulation = TracedSimulation
    result_row = {"case": case["name"], "arm": arm, "new_neural_sampling": False,
                  "input_sha256": digest(root / case["input"]), "scope": "retained protein-only source/proposal, not whole-entry cofactor qualification"}
    try:
        result = boundary.refine_generated_region(proposal,
            generated_residues=tuple(r for gap in request.gaps for r in gap.generated_residues),
            generated_atoms=request.generated_atoms, max_iterations=1000, restart_count=0, random_seed=7,
            platform_name="CPU", geometry_preconditioning=True,
            observed_anchor_restraints=restraints, anchor_movement_limits=limits if restraints else None)
        if set(result.coordinates_angstrom) != set(request.generated_atoms) | set(allowed):
            raise ValueError("returned movable atom mask differs")
        final = rewrite_generated_coordinates(proposal, result.coordinates_angstrom)
        (output / "refined.pdb").write_text(final)
        assessed = assess(final, allowed)
        coords, _ = _coordinates(final)
        movements = check_anchor_movements(restraints, coords, limits)
        structure = _structure(final, "traced-1HE8.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        result_row.update(geometry_passed=assessed.geometry_passed, admission_failures=list(assessed.admission_failures),
            assessment=asdict(assessed), movements=[asdict(m) for m in movements], global_all_l=True,
            accepted=assessed.geometry_passed and not assessed.admission_failures and all(m.within_limits for m in movements),
            initial_energy_kj_mol=result.initial_energy_kj_mol, final_energy_kj_mol=result.final_energy_kj_mol,
            final_sha256=digest(output / "refined.pdb"), search_diagnostics=result.search_diagnostics)
    except boundary.BoundaryRefinementError as exc:
        result_row.update(accepted=False, numerical_or_refinement_failure=str(exc))
        (output / "refinement-failure.txt").write_text(traceback.format_exc())
    finally:
        boundary.Simulation = original
    save(output / "result.json", {**result_row, "minimizer_calls": calls, "RMSD_reference_used": False,
        "new_donor_manually_approved": False, "published_as_accepted_bundle": False,
        "passive_reporter_reference_control_only": True})


def declare(output):
    output.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), PROJECT / "scripts/trace_private_minimization.py", PROJECT / "scripts/audit_1he8_joint_relaxation.py",
        PROJECT / "scripts/mac_reconstruction_comparison.py", PROJECT / "scripts/run_mac_modeller_comparison.py",
        PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    for case in CASES:
        root = PROJECT / case["root"]
        files += [root / name for name in ("request.json", "source.pdb", case["input"])]
    jobs = []
    for index, case in enumerate(CASES):
        for arm in ARMS:
            case_id = f"{case['name']}-{arm}"
            jobs.append({"case_id": case_id, "case": str(PROJECT / case["root"]), "seed": 7,
                "command": [sys.executable, "-m", "scripts.run_1he8_traced_mobility", "worker", "--index", str(index),
                    "--arm", arm, "--output", str(output / "execution" / f"{case_id}-seed-7")]})
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)), "jobs": jobs,
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "count": len(jobs), "ML_calls": 0,
        "optimizer_seed": 7, "max_iterations_full_phase": 1000, "restarts": 0, "preconditioning": True,
        "arms": ARMS, "bounds": asdict(AnchorMovementLimits()), "uniform_anchor_k": 10000.,
        "all_generated_mobile_every_arm": 1712, "original_observation_reference_never_ratchets": True,
        "scope": "prospective private eight-arm repair comparison on retained protein-only1HE8 inputs; no crystal scoring",
        "no_publication_or_manual_donor_approval_implied": True,
        "resource_guard": "serial CPU executor,1800s/job,RSS8GiB,host pressure,10GiB evidence cap"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--index", type=int, choices=(0, 1))
    parser.add_argument("--arm", choices=ARMS)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        worker(args.index, args.arm, args.output.resolve())
