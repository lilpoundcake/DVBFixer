#!/usr/bin/env python3
"""Per-gap motion/junction readback, not crystal scoring or topology qualification."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from openmm import HarmonicBondForce, XmlSerializer
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.diagnose.chemistry import _dihedral_deg
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates


def motion_summary(displacements):
    values = np.asarray(displacements, dtype=float)
    if values.size == 0 or not np.isfinite(values).all():
        raise ValueError("motion summary requires finite generated displacements")
    return {"atoms": len(values), "rms_angstrom": float(np.sqrt(np.mean(values**2))),
            "mean_angstrom": float(np.mean(values)), "p95_angstrom": float(np.quantile(values, .95)),
            "max_angstrom": float(np.max(values))}


def junction(coordinates, left, right):
    atoms = [AtomIdentity(r.chain, r.residue_number, r.insertion_code, name)
             for r, name in ((left, "CA"), (left, "C"), (right, "N"), (right, "CA"))]
    points = [coordinates[a] for a in atoms]
    return {"left": asdict(left), "right": asdict(right),
            "C_N_distance_angstrom": float(np.linalg.norm(points[1]-points[2])),
            "omega_degrees": _dihedral_deg(*(tuple(p) for p in points))}


def analyze(audits, output):
    bindings, rows = {}, []
    for audit in audits:
        path = audit / "independent-readback.json"
        bindings[str(path)] = digest(path)
        for checked in json.loads(path.read_text())["rows"]:
            final_path = Path(checked["final_path"])
            if digest(final_path) != checked["final_sha256"]:
                raise ValueError("independently audited endpoint changed")
            work = final_path.parent
            receipt = json.loads((work / "execution.json").read_text())
            request_path = Path(receipt["case"]) / "request.json"
            request = DiffusionRequest.from_json(request_path.read_text())
            proposal_path = work / "input.pdb"
            before, _ = _coordinates(proposal_path.read_text())
            after, _ = _coordinates(final_path.read_text())
            particles_path = work / "temporary-particle-axis.json"
            particles = json.loads(particles_path.read_text())
            axis = {AtomIdentity(p["chain"], p["residue_number"], p["insertion_code"], p["atom_name"]): p["index"]
                    for p in particles}
            system_path = work / "temporary-system.xml"
            system = XmlSerializer.deserialize(system_path.read_text())
            bonds = set()
            for force in system.getForces():
                if isinstance(force, HarmonicBondForce):
                    for i in range(force.getNumBonds()):
                        first, second, *_ = force.getBondParameters(i)
                        bonds.add(tuple(sorted((first, second))))
            for p in (request_path, proposal_path, final_path, particles_path, system_path):
                bindings[str(p)] = digest(p)
            generated = set(request.generated_atoms)
            gap_rows, covered = [], set()
            for index, gap in enumerate(request.gaps):
                residue_keys = {(r.chain, r.residue_number, r.insertion_code) for r in gap.generated_residues}
                atoms = {a for a in generated if (a.chain, a.residue_number, a.insertion_code) in residue_keys}
                if atoms & covered:
                    raise ValueError("generated gap ownership overlaps")
                covered.update(atoms)
                motions = [np.linalg.norm(after[a]-before[a]) for a in atoms]
                seams = []
                for left, right in ((gap.left_anchor, gap.generated_residues[0]),
                                    (gap.generated_residues[-1], gap.right_anchor)):
                    if left is None or right is None:
                        continue
                    c = AtomIdentity(left.chain, left.residue_number, left.insertion_code, "C")
                    n = AtomIdentity(right.chain, right.residue_number, right.insertion_code, "N")
                    if tuple(sorted((axis[c], axis[n]))) not in bonds:
                        raise ValueError("junction absent from actual temporary force-field bond graph")
                    seams.append({"before": junction(before, left, right), "after": junction(after, left, right),
                                  "actual_FF_C_N_bond_verified": True})
                gap_rows.append({"gap_index": index, "chain": gap.chain,
                    "target_interval": asdict(gap.target_interval), "residues": len(gap.generated_residues),
                    "motion_from_own_initial_proposal": motion_summary(motions), "outer_junctions": seams})
            if covered != generated:
                raise ValueError("generated ownership not completely accounted")
            rows.append({"case": checked["case"], "arm": checked["arm"], "audit_directory": str(audit),
                "accepted_under_original_outer_limits": checked["accepted"], "geometry_passed": checked["geometry_passed"],
                "all_generated_motion_from_own_initial_proposal": motion_summary([
                    np.linalg.norm(after[a]-before[a]) for a in generated]), "gaps": gap_rows,
                "not_crystal_RMSD_or_donor_ranking": True, "no_new_optimization_or_sampling": True,
                "no_continuous_trajectory_or_topology_guarantee": True})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "gap-endpoint-readback.json", {"rows": rows, "bindings": bindings,
        "analyzer_sha256": digest(Path(__file__)), "comparison_is_within_each_saved_proposal_only": True})
    for row in rows:
        if row["accepted_under_original_outer_limits"]:
            print(json.dumps({"case": row["case"], "arm": row["arm"],
                "generated_motion": row["all_generated_motion_from_own_initial_proposal"],
                "per_gap_rms_angstrom": [g["motion_from_own_initial_proposal"]["rms_angstrom"] for g in row["gaps"]]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, action="append", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    analyze([p.resolve() for p in args.audit], args.output.resolve())
