"""Read-only independent preparation, reduced-state and geometry audit for 1NON."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from openmm.app import PDBFile
from openmm.unit import nanometer
from scripts.audit_anchor_constrained_control import run as audit_states
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def audit_prepared(case, prepared):
    raw = PDBFile(str(case / "proposal.pdb"))
    expected = {(a.residue.chain.id, a.residue.id, a.residue.insertionCode.strip(), a.name):
                np.asarray(raw.positions[a.index].value_in_unit(nanometer)) for a in raw.topology.atoms()}
    axis = json.loads((prepared / "temporary-axis.json").read_text())["atoms"]
    with np.load(prepared / "initial-state.npz", allow_pickle=False) as data:
        points = data["positions_nm"].copy()
    heavy = {}
    for row in axis:
        if row["element"] == "H":
            if not row["temporary_addition"]:
                raise ValueError("auxiliary H is not marked disposable")
            continue
        a = row["identity"]
        key = (a["chain"], a["residue_number"], a["insertion_code"], a["atom_name"])
        if key in heavy or row["temporary_addition"]:
            raise ValueError("duplicate or rebuilt heavy atom in preparation")
        heavy[key] = points[row["index"]]
    if heavy.keys() != expected.keys() or any(not np.array_equal(heavy[a], p) for a, p in expected.items()):
        raise ValueError("independent preparation heavy-coordinate custody failure")
    return {"heavy_count": len(heavy), "auxiliary_h_count": len(axis)-len(heavy),
            "raw_to_prepared_heavy_bits_exact": True}


def run(root, output):
    plan = json.loads((root / "plan.json").read_text())
    verify_bindings(plan["source_code_and_input_sha256"])
    # This auditor computes sparse Cartesian tangent projection independently of
    # the optimizer's block projection and checks native/serialized atom custody.
    audit_states(root, output)
    job = plan["jobs"][0]
    case, prepared = Path(job["case"]), Path(job["snapshot"])
    work = prepared.parent
    result = json.loads((work / "result.json").read_text())
    report = {"plan_sha256": digest(root / "plan.json"), "audit_sha256": digest(Path(__file__)),
              "optimizer_calls": 0, "new_diffusion_calls": 0, "public_acceptance": False}
    if (prepared / "initial-state.npz").is_file():
        report["preparation"] = audit_prepared(case, prepared)
    if (work / "candidate.pdb").is_file():
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        assessed = assessment_adapter(request, (case / "source.pdb").read_text())((work / "candidate.pdb").read_text())
        save(output / "independent-assessment.json", asdict(assessed))
        first, second = atoms(case / "proposal.pdb"), atoms(work / "candidate.pdb")
        generated = [(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms]
        movement = np.array([np.linalg.norm(np.asarray(first[a]["xyz"])-second[a]["xyz"]) for a in generated])
        report.update(geometry_passed=assessed.geometry_passed, admission_failures=assessed.admission_failures,
                      generated_movement_rms_angstrom=float(np.sqrt(np.mean(movement**2))),
                      generated_movement_max_angstrom=float(movement.max()))
    report["worker_outcome_flags"] = result["outcome_flags"]
    verify_bindings(plan["source_code_and_input_sha256"])
    save(output / "refinement-receipt.json", report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
