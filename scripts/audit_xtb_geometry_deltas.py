"""Measure atom-identity-preserving bond/angle/stereovolume changes from audit XYZs.

Relative changes are not ideal-geometry scores. Signed-volume flips are potential
stereochemical changes, not absolute CIP assignments (symmetric centres included).
Never writes a candidate PDB or accepts an unconverged structure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from openmm.app import PDBFile
from openmm.unit import nanometer

from dvbfixer.diagnose.geometry import angle_degrees


def geometry_deltas(topology, before, after):
    atoms = list(topology.atoms())
    before, after = np.asarray(before), np.asarray(after)
    if before.shape != (len(atoms), 3) or after.shape != before.shape:
        raise ValueError("Coordinate shape does not match topology")
    if not np.isfinite(before).all() or not np.isfinite(after).all():
        raise ValueError("Nonfinite coordinates")
    neighbors = {a.index: set() for a in atoms}
    bonds, angles, volumes = [], [], []
    for a, b in topology.bonds():
        neighbors[a.index].add(b.index)
        neighbors[b.index].add(a.index)
        first = float(np.linalg.norm(before[a.index] - before[b.index]))
        second = float(np.linalg.norm(after[a.index] - after[b.index]))
        bonds.append({"atoms": [a.index, b.index], "before_angstrom": first,
                      "after_angstrom": second, "delta_angstrom": second - first})
    for center, adjacent in neighbors.items():
        for a, b in combinations(sorted(adjacent), 2):
            first = angle_degrees(before[a], before[center], before[b])
            second = angle_degrees(after[a], after[center], after[b])
            if first is not None and second is not None:
                angles.append({"atoms": [a, center, b], "before_degrees": first,
                               "after_degrees": second, "delta_degrees": second - first})
        if atoms[center].element.symbol != "C" or len(adjacent) != 4:
            continue
        if sum(atoms[i].element.symbol != "H" for i in adjacent) < 3:
            continue
        order = sorted(adjacent)

        def volume(coords, _order=order):
            p = coords[_order]
            return float(np.linalg.det((p[:3] - p[3]).T))

        first, second = volume(before), volume(after)
        volumes.append({"center": center, "neighbors": order, "before_angstrom_cubed": first,
                        "after_angstrom_cubed": second, "potential_flip": first * second < 0,
                        "near_planar": min(abs(first), abs(second)) < .01})
    return {
        "bond_count": len(bonds), "angle_count": len(angles), "tetrahedral_candidate_count": len(volumes),
        "bond_delta_rms_angstrom": float(np.sqrt(np.mean([b["delta_angstrom"]**2 for b in bonds]))) if bonds else None,
        "angle_delta_rms_degrees": float(np.sqrt(np.mean([a["delta_degrees"]**2 for a in angles]))) if angles else None,
        "bonds_changed_over_0_1_angstrom": [b for b in bonds if abs(b["delta_angstrom"]) > .1],
        "angles_changed_over_10_degrees": [a for a in angles if abs(a["delta_degrees"]) > 10],
        "potential_stereovolume_flips": [v for v in volumes if v["potential_flip"]],
        "near_planar_tetrahedral_candidates": [v for v in volumes if v["near_planar"]],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True, help="Directory containing summary.json and extraction metadata")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--glycam-bonds", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output must not exist")
    pdb = PDBFile(str(args.baseline))
    if args.glycam_bonds:
        from dvbfixer.ffutils import build_glycam_system

        build_glycam_system(["amber14-all.xml", "amber14/GLYCAM_06j-1.xml", "amber14/tip3pfb.xml"],
                            pdb.topology, pdb.positions)
    before = np.asarray(pdb.positions.value_in_unit(nanometer)) * 10
    digest = hashlib.sha256(args.baseline.read_bytes()).hexdigest()
    report = {"baseline_sha256": digest, "limitations": __doc__, "runs": [],
              "atom_identities": [{"index": a.index, "chain": a.residue.chain.id, "resid": a.residue.id,
                                   "icode": a.residue.insertionCode.strip(), "atom": a.name,
                                   "resname": a.residue.name} for a in pdb.topology.atoms()]}
    for run in json.loads((args.run / "summary.json").read_text())["runs"]:
        metadata_path = args.run / f"{run['case']}-extraction.json"
        if not metadata_path.is_file():
            continue
        metadata = json.loads(metadata_path.read_text())
        if metadata["input_sha256"] != digest:
            raise ValueError("Baseline digest does not match the frozen extraction")
        work = args.run / f"{run['case']}-{run['engine']}"
        xyz = work / ("xtbopt.xyz" if run.get("final_xyz_present") else "last-trajectory-frame.xyz")
        if not xyz.exists():
            continue
        rows = xyz.read_text().splitlines()[2:]
        expected = metadata["subsystem_atom_identities"]
        if len(rows) != len(expected) or any(row.split()[0] != atom["element"] for row, atom in zip(rows, expected)):
            raise ValueError("XYZ identity/count mismatch")
        coordinates = np.array([[float(x) for x in row.split()[1:4]] for row in rows])
        if not np.isfinite(coordinates).all():
            raise ValueError("Nonfinite XYZ")
        after = before.copy()
        frozen = set(run["frozen_indices_zero_based"])
        for full, sub in metadata["full_to_sub_indices"].items():
            if sub not in frozen:
                after[int(full)] = coordinates[sub]
        report["runs"].append({"engine": run["engine"], "converged": run.get("converged", False),
                               "coordinate_source": run.get("coordinate_source"),
                               "restored_geometry": geometry_deltas(pdb.topology, before, after)})
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
