#!/usr/bin/env python3
"""Independent source-frame movement/custody/geometry readback for all grid cells."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def audit(root: Path, output: Path):
    rows = json.loads((root / "progress.json").read_text())
    plan = json.loads((root / "plan.json").read_text())
    completion = json.loads((root / "completed.json").read_text())
    if len(rows) != plan["count"] or completion["completed"] != plan["count"]:
        raise ValueError("cannot summarize an incomplete calibration grid")
    checked, groups = [], defaultdict(list)
    for row in rows:
        result, case = row["result"], Path(row["case"])
        final = Path(result["final_path"])
        if digest(final) != result["final_sha256"]:
            raise ValueError("calibration endpoint changed")
        request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
        source = (case / "source.pdb").read_text()
        if digest(case / "source.pdb") != request.normalized_pdb.sha256:
            raise ValueError("calibration source changed")
        declaration = json.loads((final.parent / "declaration.json").read_text())
        allowed_residues = {(r["chain"], r["residue_number"], r["insertion_code"]) for r in declaration["anchors"]}
        source_atoms, final_atoms = inventory(case / "source.pdb"), inventory(final)
        allowed = frozenset(a for a in source_atoms if (a.chain, a.residue_number, a.insertion_code) in allowed_residues)
        def lines(text):
            return {(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()): line
                    for line in text.splitlines(keepends=True) if line.startswith("ATOM  ")}
        original_lines, final_lines = lines(source), lines(final.read_text())
        for atom, line in original_lines.items():
            changed = final_lines[atom]
            if atom[:3] not in allowed_residues and changed != line:
                raise ValueError("observed atom outside the declared anchor mask changed")
            if changed[:30]+changed[54:] != line[:30]+line[54:]:
                raise ValueError("anchor movement changed a non-coordinate PDB field")
        movements = []
        limits = declaration["limits"]
        for residue in sorted(allowed_residues):
            axis = [a for a in allowed if (a.chain, a.residue_number, a.insertion_code) == residue]
            distances = {a: float(np.linalg.norm(source_atoms[a][2]-final_atoms[a][2])) for a in axis}
            bb = [d for a, d in distances.items() if a.atom_name in {"N", "CA", "C", "O"}]
            if len(bb) != 4:
                raise ValueError("anchor backbone incomplete in independent readback")
            rms = float(np.sqrt(np.mean(np.square(bb))))
            ca = next(d for a, d in distances.items() if a.atom_name == "CA")
            maximum = max(distances.values())
            partners = [a for a in source_atoms if (a.chain, a.residue_number, a.insertion_code) != residue]
            ref_contacts = np.linalg.norm(np.array([source_atoms[a][2] for a in axis])[:, None]-np.array([source_atoms[a][2] for a in partners])[None, :], axis=2) <= 4.
            end_contacts = np.linalg.norm(np.array([final_atoms[a][2] for a in axis])[:, None]-np.array([final_atoms[a][2] for a in partners])[None, :], axis=2) <= 4.
            movements.append({"residue": residue, "backbone_rmsd": rms, "ca_displacement": ca, "maximum_heavy": maximum,
                "within_limits": rms <= limits["backbone_rmsd_angstrom_max"] and ca <= limits["ca_displacement_angstrom_max"] and maximum <= limits["heavy_displacement_angstrom_max"],
                "source_observed_contact_pairs_4A": int(ref_contacts.sum()), "contact_pairs_lost": int((ref_contacts & ~end_contacts).sum()),
                "contact_pairs_gained": int((end_contacts & ~ref_contacts).sum())})
        assessment = assessment_adapter(request, source)(final.read_text(), allowed)
        structure = _structure(final.read_text(), "independent-anchor-final.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        accepted = assessment.geometry_passed and not assessment.admission_failures and all(m["within_limits"] for m in movements)
        expected = result["geometry_passed"] and not result["admission_failures"] and result["movement_limits_passed"] and result["global_all_l"]
        if accepted != expected:
            raise ValueError("independent movement/geometry acceptance disagrees")
        checked.append({"case": result["case"], "proposal_seed": result["proposal_seed"], "window": result["window"],
            "strength": result["strength"], "envelope": result["envelope"], "movements": movements,
            "source_outside_mask_byte_identical": True, "anchor_noncoordinate_fields_identical": True,
            "global_all_l": True, "geometry_passed": assessment.geometry_passed, "within_limits": all(m["within_limits"] for m in movements),
            "jointly_passed": accepted, "heldout_heavy_rmsd": result["heldout_heavy_rmsd"],
            "final_path": str(final), "final_sha256": digest(final)})
        groups[result["case"], result["proposal_seed"]].append(checked[-1])
    summary = [{"case": key[0], "proposal_seed": key[1], "attempts": len(items),
        "geometry_passes": sum(r["geometry_passed"] for r in items), "within_limits": sum(r["within_limits"] for r in items),
        "joint_passes": sum(r["jointly_passed"] for r in items),
        "unique_passing_structures": len({r["final_sha256"] for r in items if r["jointly_passed"]})} for key, items in groups.items()]
    output.mkdir(parents=True, exist_ok=False)
    save(output / "independent-results.json", {"rows": checked, "summary": summary, "qualification": False,
        "calibrated_physical_uncertainty_bounds": False, "positive_production_limits_approved": False,
        "contact_counts_are_source_only_cutoff_diagnostics_not_interaction_energy": True})
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve())
