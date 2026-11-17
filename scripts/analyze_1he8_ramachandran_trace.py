#!/usr/bin/env python3
"""Read retained A581 frames; no new optimization or exact transition claims."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from openmm import HarmonicBondForce, XmlSerializer
from scripts.mac_reconstruction_comparison import digest, save, target_axis

from dvbfixer.diagnose.chemistry import _dihedral_deg
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest, ResidueIdentity
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.quality import _ramachandran_reference_populated


def angles(points):
    """Canonical dihedral semantics; points are converted from nm to Angstrom."""
    phi = _dihedral_deg(*(tuple(p) for p in points[:4]))
    psi = _dihedral_deg(*(tuple(p) for p in points[1:]))
    defined = phi is not None and psi is not None
    return {"phi_degrees": phi, "psi_degrees": psi, "defined": defined,
            "gross_outlier": defined and not _ramachandran_reference_populated(phi, psi)}


def analyze(root, output):
    plan = json.loads((root / "plan.json").read_text())
    completion = json.loads((root / "completed.json").read_text())
    if completion["completed"] != len(plan["jobs"]):
        raise ValueError("trace cohort incomplete")
    rows = []
    target = ResidueIdentity("A", "581")
    for job in plan["jobs"]:
        work = root / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        request_path = Path(job["case"]) / "request.json"
        if digest(request_path) != plan["source_code_and_input_sha256"][str(request_path)]:
            raise ValueError("request changed")
        request = DiffusionRequest.from_json(request_path.read_text())
        axis = target_axis(request)[target.chain]
        index = axis.index(target)
        previous, following = axis[index-1], axis[index+1]
        identities = [AtomIdentity(r.chain, r.residue_number, r.insertion_code, name)
                      for r, name in ((previous, "C"), (target, "N"), (target, "CA"),
                                      (target, "C"), (following, "N"))]
        particles = json.loads((work / "temporary-particle-axis.json").read_text())
        lookup = {AtomIdentity(p["chain"], p["residue_number"], p["insertion_code"], p["atom_name"]): p["index"]
                  for p in particles}
        indices = [lookup[a] for a in identities]
        system = XmlSerializer.deserialize((work / "temporary-system.xml").read_text())
        bonds = set()
        for force in system.getForces():
            if isinstance(force, HarmonicBondForce):
                for bond in range(force.getNumBonds()):
                    first, second, *_ = force.getBondParameters(bond)
                    bonds.add(tuple(sorted((first, second))))
        if any(tuple(sorted(pair)) not in bonds for pair in zip(indices, indices[1:])):
            raise ValueError("actual FF graph lacks required phi/psi bonds")
        result = json.loads((work / "result.json").read_text())
        frames, bindings = [], {}
        for filename, label in (("input.pdb", "proposal"), ("refined.pdb", "final-rounded")):
            path = work / filename
            coords, _ = _coordinates(path.read_text())
            frames.append({"frame": label, **angles(np.asarray([coords[a] for a in identities]))})
            bindings[str(path)] = digest(path)
        initial, final = frames
        frames = [initial]
        for phase in sorted(work.glob("minimize-*")):
            for snapshot in [*sorted(phase.glob("callback-*.npz")), phase / "returned-state.npz"]:
                with np.load(snapshot, allow_pickle=False) as archive:
                    positions = np.asarray(archive["positions_nm"])
                    if positions.shape != (len(particles), 3) or not np.isfinite(positions).all():
                        raise ValueError("invalid retained particle frame")
                    points = positions[indices] * 10
                frames.append({"frame": str(snapshot.relative_to(work)), **angles(points)})
                bindings[str(snapshot)] = digest(snapshot)
        frames.append(final)
        first_outlier = next((f["frame"] for f in frames if f["gross_outlier"]), None)
        rows.append({"case": result["case"], "arm": result["arm"],
            "target_residue": "A/GLU581", "actual_target_neighbors": [str(previous), str(following)],
            "actual_FF_phi_psi_bonds_verified": True, "first_retained_outlier_frame": first_outlier,
            "frames": frames, "frame_sha256": bindings,
            "limitation": "sparse retained frames; neither exact transition timing nor no-passage proof",
            "no_crystal_reference_or_new_optimization": True})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "ramachandran-trace.json", {"rows": rows, "analyzer_sha256": digest(Path(__file__)),
        "plan_sha256": digest(root / "plan.json")})
    for row in rows:
        print(json.dumps({"case": row["case"], "arm": row["arm"],
            "first_retained_outlier": row["first_retained_outlier_frame"],
            "proposal": row["frames"][0], "final": row["frames"][-1]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    analyze(args.root.resolve(), args.output.resolve())
