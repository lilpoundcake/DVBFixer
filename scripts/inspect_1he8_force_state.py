#!/usr/bin/env python3
"""Reconstruct temporary FF states; abort before heavy-atom minimization.

Disposable hydrogen placement includes OpenMM's own internal H optimization.
This is NOT the historical H state, historical gradient or a repaired proposal.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_1he8_joint_relaxation import CASES, PROJECT, read_restraints
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion import boundary_refinement as boundary
from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest


class ConstructionCaptured(RuntimeError):
    pass


def worker(index, output):
    from openmm import XmlSerializer
    from openmm.app import PDBFile
    from openmm.unit import dalton, kilojoule_per_mole, nanometer

    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("force inspector requires a fresh executor-owned reservation")
    case = CASES[index]
    root = PROJECT / case["root"]
    declaration = json.loads((root / case["declaration"]).read_text())
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    source = root / "source.pdb"
    if digest(source) != request.normalized_pdb.sha256:
        raise ValueError("force inspection source changed")
    text = (root / case["final"]).read_text()
    restraints = read_restraints(declaration)
    generated, anchors = set(request.generated_atoms), {r.atom for r in restraints}
    limits = AnchorMovementLimits(**declaration[case["limits"]])
    captured = {}
    original_simulation = boundary.Simulation

    class CaptureSimulation(original_simulation):
        def minimizeEnergy(self, *args, **kwargs):
            if captured:
                raise ValueError("unexpected second heavy-atom minimization call")
            atoms = list(self.topology.atoms())
            axes = [AtomIdentity(a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.name) for a in atoms]
            lookup = dict(zip(axes, atoms, strict=True))
            if len(lookup) != len(axes) or not generated | anchors <= lookup.keys():
                raise ValueError("final temporary topology lost or duplicated movable identities")
            original_groups = [f.getForceGroup() for f in self.system.getForces()]
            if len(original_groups) > 32:
                raise ValueError("force decomposition exceeds OpenMM group limit")
            original_xml = XmlSerializer.serialize(self.system)
            (output / "temporary-system-original-groups.xml").write_text(original_xml)
            parameters = dict(self.context.getParameters())
            for name, value in (("soft_scale", 0.0), ("soft_alpha", 0.2)):
                if name in parameters:
                    self.context.setParameter(name, value)
            for group, force in enumerate(self.system.getForces()):
                force.setForceGroup(group)
            self.context.reinitialize(preserveState=True)
            state = self.context.getState(getPositions=True, getParameters=True)
            (output / "temporary-state.xml").write_text(XmlSerializer.serialize(state))
            with (output / "temporary-H-topology-NOT-A-CANDIDATE.pdb").open("x") as stream:
                PDBFile.writeFile(self.topology, state.getPositions(), stream, keepIds=True)
            mass_rows = []
            for identity, atom in lookup.items():
                hydrogen = atom.element is not None and atom.element.symbol == "H"
                role = ("generated-heavy" if identity in generated else "released-observed-anchor" if identity in anchors
                        else "temporary-hydrogen" if hydrogen else "fixed-observed-heavy" if identity in request.fixed_atoms
                        else "temporary-heavy-completion")
                mass = float(self.system.getParticleMass(atom.index).value_in_unit(dalton))
                if role in {"generated-heavy", "released-observed-anchor"} and mass <= 0:
                    raise ValueError("requested mobile atom is massless in final FF topology")
                if role == "fixed-observed-heavy" and mass != 0:
                    raise ValueError("unreleased observed atom is mobile in final FF topology")
                mass_rows.append({"atom": asdict(identity), "index": atom.index, "role": role, "mass_dalton": mass})
            matrices, terms = [], []
            for index, force in enumerate(self.system.getForces()):
                component = self.context.getState(getEnergy=True, getForces=True, groups=1 << index)
                energy = float(component.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
                forces = np.asarray(component.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer))
                if not np.isfinite(forces).all() or not np.isfinite(energy):
                    raise ValueError("nonfinite reconstructed force component")
                matrices.append(forces)
                terms.append({"index": index, "class": type(force).__name__, "name": force.getName(),
                              "original_group": original_groups[index], "analysis_group": index, "energy_kj_mol": energy})
            force_matrix = np.asarray(matrices)
            total = force_matrix.sum(axis=0)
            norms = np.linalg.norm(total, axis=1)
            targets = []
            for identity, atom in lookup.items():
                if identity.chain != "A" or identity.residue_number not in {"321", "353", "581", "895", "901", "902"}:
                    continue
                if atom.element is not None and atom.element.symbol == "H":
                    continue
                targets.append({"atom": asdict(identity), "index": atom.index,
                    "role": mass_rows[atom.index]["role"], "total_force_kj_mol_nm": total[atom.index].tolist(),
                    "norm_kj_mol_nm": float(norms[atom.index]),
                    "component_forces_kj_mol_nm": force_matrix[:, atom.index, :].tolist()})
            np.savez_compressed(output / "temporary-force-components.npz", forces_kj_mol_nm=force_matrix,
                               positions_nm=state.getPositions(asNumpy=True).value_in_unit(nanometer))
            captured.update({"case": case["name"], "mass_rows": mass_rows, "force_terms": terms, "target_atoms": targets,
                "generated_heavy_atoms_mobile": sum(r["role"] == "generated-heavy" for r in mass_rows),
                "observed_anchor_atoms_mobile": sum(r["role"] == "released-observed-anchor" for r in mass_rows),
                "observed_nonanchor_atoms_fixed": sum(r["role"] == "fixed-observed-heavy" for r in mass_rows),
                "temporary_forcefield": list(boundary.BOUNDARY_REFINEMENT_FORCEFIELD),
                "default_tolerance_kj_mol_nm": boundary.BOUNDARY_REFINEMENT_TOLERANCE_KJ_MOL_NM,
                "main_minimization_executed": False, "temporary_H_placement_internal_optimization_executed": True,
                "historical_hydrogens_or_gradients_reproduced": False, "historical_convergence_established": False,
                "force_state": "reconstructed temporary H/FF around retained final heavy atoms, full nonbonded/soft_scale0",
                "groups_reassigned_for_readonly_decomposition": True, "source_output_rewritten": False})
            raise ConstructionCaptured()

    boundary.Simulation = CaptureSimulation
    try:
        boundary.refine_generated_region(text,
            generated_residues=tuple(r for gap in request.gaps for r in gap.generated_residues),
            generated_atoms=request.generated_atoms, platform_name="CPU", random_seed=7,
            max_iterations=1, restart_count=0, observed_anchor_restraints=restraints,
            anchor_movement_limits=limits, geometry_preconditioning=True)
    except ConstructionCaptured:
        if not captured:
            raise ValueError("FF interception did not produce evidence") from None
    finally:
        boundary.Simulation = original_simulation
    save(output / "result.json", {**captured, "input_final_sha256": digest(root / case["final"]),
        "declaration_sha256": digest(root / case["declaration"]), "inspector_sha256": digest(Path(__file__)),
        "refiner_sha256": digest(Path(boundary.__file__)), "not_a_candidate_or_historical_replay": True})


def declare(output):
    output.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), PROJECT / "scripts/audit_1he8_joint_relaxation.py"]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    for case in CASES:
        root = PROJECT / case["root"]
        files += [root / name for name in ("request.json", "source.pdb", case["declaration"], case["final"])]
    jobs = [{"case_id": case["name"], "case": str(PROJECT / case["root"]), "seed": 7,
             "command": [sys.executable, "-m", "scripts.inspect_1he8_force_state", "worker", "--index", str(index),
                         "--output", str(output / "execution" / f"{case['name']}-seed-7")]}
            for index, case in enumerate(CASES)]
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)), "jobs": jobs,
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "count": 2,
        "purpose": "temporary-system mass/force readback, not new heavy minimization or historical trajectory reproduction",
        "temporary_H_placement_has_internal_optimization": True, "ML_calls": 0,
        "resource_guard": "existing serial CPU executor, 1800s/job, RSS8GiB, host pressure and 10GiB output cap"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=(0, 1))
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        worker(args.index, args.output.resolve())
