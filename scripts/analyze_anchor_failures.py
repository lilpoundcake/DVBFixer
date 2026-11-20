#!/usr/bin/env python3
"""Saved-state refusal investigation. No Context, optimization, repair or sampling.

The topology-only steric control uses the same radii, 0.6/0.9 Å thresholds,
water inclusion, same-residue and 1-3/1-4 exclusions as the recorded study.
Only proximity-inferred extra covalent exclusions are withheld in this control.
This is a private diagnostic comparison, not a production behavior change.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import numpy as np
from openmm import CustomCompoundBondForce, NonbondedForce, XmlSerializer
from openmm.unit import dalton, elementary_charge, kilojoule_per_mole, nanometer
from scipy.spatial import cKDTree
from scripts.anchor_calibration_grid import read_temporary
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import verify_bindings

from dvbfixer.diagnose import steric
from dvbfixer.diagnose.chemistry import run_all
from dvbfixer.ffutils.geometry import find_d_residues
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.quality import check_ramachandran, check_sidechain_chi12

COULOMB_KJ_NM_MOL_E2 = 138.935456
WITNESSES = {"6Z38-fixed", "1NON-fixed", "1NON-w2-e2-k100000",
             "2BE1-w1-e3-k10000", "1KTZ-fixed"}


def historical_covalent_cutoff(first, second):
    """Frozen v1 exclusion rule (cd95dc1), NOT current diagnostic chemistry."""
    if "H" in (first.upper(), second.upper()):
        return 1.30
    return 2.25 if first.upper() == second.upper() == "S" else 1.90


def identity(row):
    return AtomIdentity(**row["identity"])


def neighbor_graph(bonds, included):
    graph = {i: set() for i in included}
    for first, second in bonds:
        if first in included and second in included:
            graph[first].add(second)
            graph[second].add(first)
    return graph


def topology_only_contacts(axis, bonds, indices, coordinates_nm, generated, mobile):
    """Expose severe contacts hidden by adding distance-derived fake bonds."""
    indices = sorted(indices)
    actual = neighbor_graph(bonds, set(indices))
    actual_exclusions = steric._expand_exclusions(actual)
    augmented = {i: set(js) for i, js in actual.items()}
    pairs = []
    for first, second in cKDTree(coordinates_nm[indices]*10).query_pairs(steric._NEIGHBOR_CUTOFF_A, output_type="ndarray"):
        i, j = indices[first], indices[second]
        distance = float(np.linalg.norm(coordinates_nm[i]-coordinates_nm[j])*10)
        pairs.append((i, j, distance))
        if distance <= historical_covalent_cutoff(axis[i]["element"], axis[j]["element"]):
            augmented[i].add(j)
            augmented[j].add(i)
    inferred_exclusions = steric._expand_exclusions(augmented)
    severe = []
    for i, j, distance in pairs:
        pair = frozenset((i, j))
        if pair in actual_exclusions or axis[i]["residue_index"] == axis[j]["residue_index"]:
            continue
        ei, ej = axis[i]["element"], axis[j]["element"]
        if (ei in steric._HBOND_ACCEPTORS and ej in steric._HBOND_ACCEPTORS
                and steric._HBOND_HEAVY_MIN_A <= distance <= steric._HBOND_HEAVY_MAX_A):
            continue
        overlap = steric._VDW_A.get(ei, 1.7)+steric._VDW_A.get(ej, 1.7)-distance
        if overlap < .9-1.e-12:
            continue
        severe.append({"atoms": [axis[i]["identity"], axis[j]["identity"]],
            "resnames": [axis[i]["resname"], axis[j]["resname"]], "indices": [i, j],
            "distance_angstrom": distance, "overlap_angstrom": overlap,
            "hidden_by_distance_inferred_exclusions": pair in inferred_exclusions,
            "touches_generated": i in generated or j in generated,
            "touches_mobile": i in mobile or j in mobile,
            "both_fixed": i not in mobile and j not in mobile})
    return sorted(severe, key=lambda r: r["distance_angstrom"])


def pair_terms(index, xyz, charges, sigmas, epsilons, exceptions, axis):
    """Analytic NoCutoff pair terms from recorded NB parameters, not a replay."""
    delta = xyz[index]-xyz
    radius = np.linalg.norm(delta, axis=1)
    radius[index] = np.inf
    charge_product = charges[index]*charges
    sigma = (sigmas[index]+sigmas)/2
    epsilon = np.sqrt(epsilons[index]*epsilons)
    exceptional = set()
    for pair, parameters in exceptions.items():
        if index in pair:
            other = pair[1] if pair[0] == index else pair[0]
            charge_product[other], sigma[other], epsilon[other] = parameters
            exceptional.add(other)
    radius[(charge_product == 0) & (epsilon == 0)] = np.inf
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        coulomb = COULOMB_KJ_NM_MOL_E2*charge_product/radius
        lj = 4*epsilon*((sigma/radius)**12-(sigma/radius)**6)
        coefficient = COULOMB_KJ_NM_MOL_E2*charge_product/radius**3 + 24*epsilon*(2*(sigma/radius)**12-(sigma/radius)**6)/radius**2
        pair_forces = delta*coefficient[:, None]
    order = sorted((j for j in range(len(axis)) if j != index and charge_product[j] < 0), key=lambda j: coulomb[j])[:3]
    def finite(value):
        return float(value) if math.isfinite(float(value)) else None

    return [{"partner": axis[j]["identity"], "partner_resname": axis[j]["resname"],
             "partner_temporary": axis[j]["temporary_addition"], "distance_angstrom": float(radius[j]*10),
             "charge_product_e2": float(charge_product[j]), "pair_epsilon_kj_mol": float(epsilon[j]),
             "coulomb_energy_kj_mol": finite(coulomb[j]), "lj_energy_kj_mol": finite(lj[j]),
             "pair_force_norm_kj_mol_nm": finite(np.linalg.norm(pair_forces[j])),
             "nonbonded_exception": j in exceptional} for j in order]


def force_summary(forces, masses, generated):
    mobile = masses > 0
    norms = np.linalg.norm(forces, axis=1)
    return {"rms_mobile_components": float(np.sqrt(np.mean(forces[mobile]**2))),
            "mobile_gradient_rms_over_all_components": float(np.sqrt(np.sum(forces[mobile]**2)/forces.size)),
            "generated_heavy_rms_components": float(np.sqrt(np.mean(forces[sorted(generated)]**2))),
            "max_mobile_force": float(norms[mobile].max()),
            "max_generated_heavy_force": float(norms[sorted(generated)].max())}


def heavy_topology(topology, public):
    from openmm.app import Topology

    result = Topology()
    chains, residues, mapped, native_indices = {}, {}, {}, []
    for atom in topology.atoms():
        a = AtomIdentity(atom.residue.chain.id, str(atom.residue.id), atom.residue.insertionCode.strip(), atom.name)
        if a not in public:
            continue
        chain, residue = atom.residue.chain, atom.residue
        if chain.index not in chains:
            chains[chain.index] = result.addChain(chain.id)
        if residue.index not in residues:
            residues[residue.index] = result.addResidue(residue.name, chains[chain.index], str(residue.id),
                                                       insertionCode=residue.insertionCode)
        mapped[atom.index] = result.addAtom(atom.name, atom.element, residues[residue.index])
        native_indices.append(atom.index)
    for first, second in topology.bonds():
        if first.index in mapped and second.index in mapped:
            result.addBond(mapped[first.index], mapped[second.index])
    return result, native_indices


def run(grid, output):
    output.mkdir(parents=True, exist_ok=False)
    plan_path = grid / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    (output / "executed-analysis.py").write_bytes(Path(__file__).read_bytes())
    bindings = {str(plan_path): digest(plan_path), str(Path(__file__).resolve()): digest(Path(__file__))}
    bindings[str(Path(sys.executable))] = digest(Path(sys.executable))
    for path in (Path(steric.__file__), Path(__file__).with_name("anchor_calibration_grid.py")):
        bindings[str(path.resolve())] = digest(path)
    for job in plan["jobs"]:
        work = grid / f"execution/{job['case_id']}-seed-{job['seed']}"
        for path in work.rglob("*"):
            if path.is_file() and path.suffix in {".npz", ".json", ".jsonl", ".xml", ".pdb"}:
                bindings[str(path)] = digest(path)
    save(output / "declaration.json", {"grid_plan_sha256": digest(plan_path), "bindings": bindings,
        "scope": "all100 endpoint force/contact scan; 5named controls + all6D refusal phase histories",
        "optimizer_calls": 0, "Context_created": False, "thresholds": {"warning_overlap_A": .6, "error_overlap_A": .9},
        "steric_control": "same rules/thresholds; actual FF graph exclusions rather than inferred proximity edges",
        "causality_limit": "saved endpoints/checkpoints, not continuous trajectories or counterfactual optimization"})
    cache, rows, histories = {}, [], []
    for job in plan["jobs"]:
        case, cell_id = Path(job["case"]), job["case_id"]
        work = grid / f"execution/{cell_id}-seed-{job['seed']}"
        result = json.loads((work / "result.json").read_text())
        snapshot = Path(job["command"][job["command"].index("--snapshot")+1])
        if case.name not in cache:
            request = DiffusionRequest.from_json((case / "request.json").read_text())
            topology, _, initial = read_temporary(snapshot)
            metadata = json.loads((snapshot / "temporary-axis.json").read_text())
            axis = metadata["atoms"]
            public = {AtomIdentity(*a) for a in atoms(case / "proposal.pdb")}
            public_indices = {i for i, a in enumerate(axis) if identity(a) in public}
            generated = {i for i, a in enumerate(axis) if identity(a) in request.generated_atoms}
            base = XmlSerializer.deserialize((snapshot / "base-system.xml").read_text())
            nb = next(f for f in base.getForces() if isinstance(f, NonbondedForce))
            parameters = [nb.getParticleParameters(i) for i in range(base.getNumParticles())]
            charges = np.array([q.value_in_unit(elementary_charge) for q, _, _ in parameters])
            sigmas = np.array([s.value_in_unit(nanometer) for _, s, _ in parameters])
            epsilons = np.array([e.value_in_unit(kilojoule_per_mole) for _, _, e in parameters])
            exceptions = {}
            for k in range(nb.getNumExceptions()):
                i, j, q, s, e = nb.getExceptionParameters(k)
                exceptions[tuple(sorted((i, j)))] = (q.value_in_unit(elementary_charge**2), s.value_in_unit(nanometer), e.value_in_unit(kilojoule_per_mole))
            heavy, heavy_indices = heavy_topology(topology, public)
            cache[case.name] = (request, topology, initial, axis, metadata["bonds"], public_indices,
                generated, charges, sigmas, epsilons, exceptions, heavy, heavy_indices, XmlSerializer.serialize(nb))
        request, topology, initial, axis, bonds, public_indices, generated, charges, sigmas, epsilons, exceptions, heavy, heavy_indices, nb_definition = cache[case.name]
        system = XmlSerializer.deserialize((work / "system.xml").read_text())
        nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
        nb.setForceGroup(0)
        if XmlSerializer.serialize(nb) != nb_definition:
            raise ValueError("cell nonbonded definition differs from shared chemistry")
        masses = np.array([system.getParticleMass(i).value_in_unit(dalton) for i in range(len(axis))])
        mobile = {i for i in public_indices if masses[i] > 0}
        with np.load(work / "phase-4/returned-state.npz", allow_pickle=False) as data:
            xyz, forces = data["positions_nm"].copy(), data["active_forces"].copy()
        summary = force_summary(forces, masses, generated)
        contacts = topology_only_contacts(axis, bonds, public_indices, xyz, generated, mobile)
        hidden = [c for c in contacts if c["hidden_by_distance_inferred_exclusions"]]
        ranked = np.argsort(np.where(masses > 0, np.linalg.norm(forces, axis=1), -1))[-3:][::-1]
        top = [{"atom": axis[int(i)], "force_norm_kj_mol_nm": float(np.linalg.norm(forces[i])),
                "attractive_NB_partners": pair_terms(int(i), xyz, charges, sigmas, epsilons, exceptions, axis)} for i in ranked]
        row = {"entry": case.name, "cell_id": cell_id, "original_completed": result["completed"],
            "original_outcome_flags": result["outcome_flags"], "force_state": summary, "top_mobile_forces": top,
            "topology_only_severe_contacts": len(contacts), "hidden_severe_contacts": len(hidden),
            "hidden_touching_generated": sum(c["touches_generated"] for c in hidden),
            "hidden_both_fixed": sum(c["both_fixed"] for c in hidden), "hidden_contacts": hidden,
            "callback_count": sum(1 for _ in (work / "phase-4/callbacks.jsonl").open()),
            "main_cap": 1000, "termination_reason_exposed": False}
        rows.append(row)
        if cell_id in WITNESSES or not result["completed"]:
            centres = []
            for force in system.getForces():
                if isinstance(force, CustomCompoundBondForce) and "vref" in force.getEnergyFunction():
                    centres = [force.getBondParameters(i) for i in range(force.getNumBonds())]
            states = [("initial", snapshot / "initial-state.npz")]
            states += [(p.parent.name, p) for p in sorted(work.glob("phase-*/returned-state.npz"))]
            phases = []
            for label, path in states:
                with np.load(path, allow_pickle=False) as data:
                    points = data["positions_nm"].copy()
                    force = data["full_force_kj_mol_nm" if label == "initial" else "active_forces"].copy()
                diagnostics = run_all(heavy, points[heavy_indices]*nanometer)
                diagnostics += check_ramachandran(heavy, points[heavy_indices]*nanometer)
                diagnostics += check_sidechain_chi12(heavy, points[heavy_indices]*nanometer)
                protected = []
                for particles, parameters in centres:
                    n, ca, c, cb = particles
                    signed = float(np.dot(np.cross(points[n]-points[ca], points[c]-points[ca]), points[cb]-points[ca]))
                    if signed < 0:
                        protected.append({"CA": axis[ca]["identity"], "signed_nm3": signed,
                                          "reference_nm3": parameters[0], "q": signed/parameters[0]})
                moving_findings = [asdict(f) for f in diagnostics if f.severity.value != "INFO"]
                local_contacts = topology_only_contacts(axis, bonds, public_indices, points, generated, mobile)
                dominant = int(np.argmax(np.where(masses > 0, np.linalg.norm(force, axis=1), -1)))
                phases.append({"phase": label, "state_sha256": digest(path), "force_state": force_summary(force, masses, generated),
                    "D_ca": find_d_residues(topology, points*nanometer), "negative_protected_centres": protected,
                    "geometry_findings": moving_findings,
                    "geometry_category_counts": dict(Counter(f["category"] for f in moving_findings)),
                    "severe_contacts": local_contacts,
                    "dominant_mobile_atom": axis[dominant],
                    "dominant_attractive_partners": pair_terms(dominant, points, charges, sigmas, epsilons, exceptions, axis)})
            history = {"entry": case.name, "cell_id": cell_id, "phases": phases,
                       "topology_basis": "recorded actual FF heavy-atom graph; exact saved native states"}
            save(output / f"{cell_id}-history.json", history)
            histories.append(cell_id)
        print(cell_id, "read-only scan complete", flush=True)
    save(output / "endpoint-scan.json", {"rows": rows, "optimizer_calls": 0, "histories": histories})
    summary = []
    for entry in plan["admitted"]:
        selected = [r for r in rows if r["entry"] == entry]
        summary.append({"entry": entry, "cells": len(selected),
            "returned_with_masked_force_rms_over_tolerance": sum(r["original_completed"] and r["force_state"]["mobile_gradient_rms_over_all_components"] > 10 for r in selected),
            "cells_with_hidden_severe_generated_contacts": sum(r["hidden_touching_generated"] > 0 for r in selected),
            "cells_with_hidden_severe_fixed_contacts": sum(r["hidden_both_fixed"] > 0 for r in selected),
            "hidden_severe_contact_range": [min(r["hidden_severe_contacts"] for r in selected), max(r["hidden_severe_contacts"] for r in selected)]})
    save(output / "summary.json", {"cases": summary, "histories": histories,
        "optimizer_calls": 0, "repair_calls": 0, "Context_created": False,
        "analysis_sha256": digest(Path(__file__)), "declaration_sha256": digest(output / "declaration.json"),
        "scope_note": "private topology-grounded diagnostic control, not a silent rewrite of prior acceptance",
        "original_results_changed": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.grid.resolve(), args.output.resolve())
