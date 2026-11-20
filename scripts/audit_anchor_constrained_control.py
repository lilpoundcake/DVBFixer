#!/usr/bin/env python3
"""Independent sparse-Jacobian/custody audit of private reduced-coordinate controls."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from openmm import XmlSerializer
from openmm.app import PDBFile
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scripts.anchor_calibration_grid import read_temporary
from scripts.anchor_refinement_guards import xh_equilibria
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import verify_bindings

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.contract import DiffusionRequest


def independent_projection(system, bonds, positions, forces):
    from openmm.unit import dalton
    masses = np.array([system.getParticleMass(i).value_in_unit(dalton) for i in range(system.getNumParticles())])
    free = np.flatnonzero(np.repeat(masses > 0, 3))
    columns = np.full(3*len(masses), -1, dtype=int)
    columns[free] = np.arange(len(free))
    rows, cols, values = [], [], []
    for row, bond in enumerate(bonds):
        h, parent = bond["hydrogen_index"], bond["parent_index"]
        normal = positions[h]-positions[parent]
        normal /= np.linalg.norm(normal)
        for index, sign in ((h, 1), (parent, -1)):
            if masses[index] == 0:
                continue
            for component in range(3):
                rows.append(row)
                cols.append(columns[3*index+component])
                values.append(sign*normal[component])
    j = coo_matrix((values, (rows, cols)), shape=(len(bonds), len(free))).tocsr()
    f = forces.ravel()[free]
    projected = f-j.T@spsolve((j@j.T).tocsc(), j@f)
    dof = len(free)-len(bonds)
    return float(np.linalg.norm(projected)/np.sqrt(dof)), float(np.max(np.abs(j@projected)))


def run(root, output):
    output.mkdir(parents=True, exist_ok=False)
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"])
    results = []
    for job in plan["jobs"]:
        work = root / f"execution/{job['case_id']}-seed-{job['seed']}"
        result = json.loads((work / "result.json").read_text())
        if not result.get("optimization_executed"):
            if not (work / "admission.json").is_file():
                results.append({"case_id": job["case_id"], "scope": "preparation failed; no optimization",
                                "outcome_flags": result["outcome_flags"]})
                continue
            admission = json.loads((work / "admission.json").read_text())
            if result["Context_created"] or admission["admitted"]:
                raise ValueError("invalid retained admission refusal")
            results.append({"case_id": job["case_id"], "scope": "refused before Context",
                            "refusal": admission["reasons"]})
            continue
        if "candidate_sha256" not in result:
            results.append({"case_id": job["case_id"], "scope": "failed; no serialized candidate",
                            "outcome_flags": result["outcome_flags"]})
            continue
        topology, _, initial = read_temporary(Path(job["snapshot"]))
        system = XmlSerializer.deserialize((work / "system.xml").read_text())
        bonds = xh_equilibria(topology, system)
        from openmm.unit import dalton
        fixed = np.array([system.getParticleMass(i).value_in_unit(dalton) == 0 for i in range(len(initial))])
        state_rows = []
        for path in sorted(work.glob("phase-*/*.npz")):
            with np.load(path, allow_pickle=False) as data:
                positions = data["positions_nm"].copy()
                forces = data["raw_forces"].copy()
            if not np.array_equal(positions[fixed], initial[fixed]):
                raise ValueError("independent native fixation failure")
            maximum = max(abs(float(np.linalg.norm(positions[b["hydrogen_index"]]-positions[b["parent_index"]]))-b["equilibrium_nm"]) for b in bonds)
            if maximum > 1.e-12:
                raise ValueError("independent FF X-H equality failure")
            if path.name == "returned-state.npz":
                rms, tangency = independent_projection(system, bonds, positions, forces)
                completed = json.loads((path.parent / "completed.json").read_text())
                if not np.isclose(rms, completed["projected_cartesian_dof_rms_kj_mol_nm"], rtol=1.e-10, atol=1.e-9):
                    raise ValueError("independent sparse vs block projection disagreement")
                state_rows.append({"phase": path.parent.name, "projected_dof_rms": rms,
                                   "J_projected_force_max": tangency, "xh_max_error_nm": maximum})
        case = Path(job["case"])
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        before, after = atoms(case / "proposal.pdb"), atoms(work / "candidate.pdb")
        generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
        if set(before) != set(after) or any(before[a]["line"] != after[a]["line"] for a in set(before)-generated):
            raise ValueError("independent serialized observation/inventory failure")
        first_lines = (case / "proposal.pdb").read_text().splitlines(keepends=True)
        second_lines = (work / "candidate.pdb").read_text().splitlines(keepends=True)
        if len(first_lines) != len(second_lines):
            raise ValueError("physical record count changed")
        for first, second in zip(first_lines, second_lines, strict=True):
            key = (first[21], first[22:26].strip(), first[26].strip(), first[12:16].strip()) if first.startswith(("ATOM  ", "HETATM")) else None
            if key in generated:
                if first[:30] != second[:30] or first[54:] != second[54:]:
                    raise ValueError("generated-atom noncoordinate metadata changed")
            elif first != second:
                raise ValueError("fixed/context/ANISOU/annotation record changed")
        # Native-to-serialized identity audit, with no source-coordinate restore.
        axis = json.loads((Path(job["snapshot"]) / "temporary-axis.json").read_text())["atoms"]
        with np.load(work / "phase-4/returned-state.npz", allow_pickle=False) as data:
            native = data["positions_nm"]*10
        for index, row in enumerate(axis):
            a = row["identity"]
            key = (a["chain"], a["residue_number"], a["insertion_code"], a["atom_name"])
            if key in generated and any(f"{value:8.3f}" != after[key]["line"][start:start+8]
                                        for value, start in zip(native[index], (30, 38, 46), strict=True)):
                raise ValueError("native generated coordinate serialization mismatch")
        pdb = PDBFile(str(work / "candidate.pdb"))
        assert_all_l(pdb.topology, pdb.positions)
        dat = DatRecord.load(work / "generated.dat")
        if len(dat.added_atoms) != len(generated) or {(a["chain"], a["resid"], a["icode"], a["atom"]) for a in dat.added_atoms} != generated:
            raise ValueError("independent generated-only ownership failure")
        if digest(work / "candidate.pdb") != result["candidate_sha256"]:
            raise ValueError("candidate bytes changed")
        results.append({"case_id": job["case_id"], "native_and_serialized_fixed_exact": True,
                        "generated_dat_exact": True, "serialized_global_l": True, "phase_checks": state_rows,
                        "worker_outcome_flags": result["outcome_flags"],
                        "restoration_mask_passed": result["restoration_mask_passed"]})
    verify_bindings(plan["source_code_and_input_sha256"])
    save(output / "receipt.json", {"plan_sha256": digest(plan_path), "audit_code_sha256": digest(Path(__file__)),
        "cases": results, "optimizer_calls": 0, "Context_created": False,
        "scope": "independent sparse-Jacobian tangent projection and native/serialized custody, not new optimization"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
