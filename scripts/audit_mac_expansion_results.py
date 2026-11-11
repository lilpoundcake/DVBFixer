#!/usr/bin/env python3
"""Independent readback of completed private terminal/pair pilot final results."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.gap_neighbors import source_gap_components
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def audit(root: Path, output: Path):
    rows = json.loads((root / "final-results.json").read_text())
    plan = json.loads((root / "runs/plan.json").read_text())
    declared = {(Path(j["case"]).name, j["arm"], j["seed"]): j for j in plan["jobs"]}
    if len(rows) != len(declared):
        raise ValueError("incomplete pilot outcomes")
    checked, seen = [], set()
    for row in rows:
        key = row["case"], row["backend"], row["seed"]
        if key not in declared or key in seen:
            raise ValueError("undeclared/repeated outcome")
        seen.add(key)
        job = declared[key]
        case, work = Path(job["case"]), Path(job["work"])
        final = Path(row["final_path"])
        if digest(final) != row["sha256"]:
            raise ValueError("final bytes changed")
        manifest = json.loads((case / "manifest.json").read_text())
        for name, expected in manifest["files"].items():
            if digest(case / name) != expected:
                raise ValueError("source/reference/mask bytes changed")
        request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
        source, final_text = (case / "source.pdb").read_text(), final.read_text()
        def atom_lines(text):
            result = {}
            for line in text.splitlines(keepends=True):
                if not line.startswith("ATOM  "):
                    continue
                key = line[21], line[22:26], line[26], line[12:16]
                if key in result:
                    raise ValueError("duplicate atom in readback")
                result[key] = line
            return result
        originals, finals = atom_lines(source), atom_lines(final_text)
        if any(finals.get(atom) != line for atom, line in originals.items()):
            raise ValueError("observed ATOM bytes changed")
        assessment = assessment_adapter(request, source)(final_text)
        structure = _structure(final_text, "independent-final.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        coords = inventory(final)
        reference = inventory(case / "reference-projection-NOT-CONDITIONING.pdb")
        raw = inventory(work / "audit/source-restored.pdb")
        if set(coords) != set(request.fixed_atoms) | set(request.generated_atoms):
            raise ValueError("final atom inventory changed")
        def rms(axis, second=reference, current=coords):
            return float(np.sqrt(np.mean([np.sum((current[a][2]-second[a][2])**2) for a in axis])))
        generated = request.generated_atoms
        heavy_rmsd = rms(generated)
        if abs(heavy_rmsd-row["heldout_heavy_rmsd"]) > 1e-10:
            raise ValueError("independent RMSD disagrees")
        raw_audit = json.loads((work / "audit/audit.json").read_text())
        components, edges = source_gap_components(request, source)
        checked.append({"case": key[0], "backend": key[1], "seed": key[2], "final_path": str(final), "final_sha256": digest(final),
            "observed_atom_records_byte_identical": True, "global_all_l": True,
            "raw_restored_geometry_passed": raw_audit["assessments"]["source-restored"]["geometry_passed"],
            "native_observed_drift": raw_audit["native_observed_drift_after_observed_only_fit_angstrom"],
            "repaired": row["repaired"], "final_geometry_passed": assessment.geometry_passed,
            "final_admission_failures": list(assessment.admission_failures),
            "final_heldout_heavy_rmsd": heavy_rmsd,
            "final_heldout_backbone_rmsd": rms([a for a in generated if a.atom_name in {"N", "CA", "C", "O"}]),
            "generated_repair_rms_movement": rms(generated, raw),
            "source_only_gap_components": components, "source_only_gap_edges": edges,
            "per_gap_heavy_rmsd": [rms([a for a in generated if (a.chain, a.residue_number, a.insertion_code) in
                {(r.chain, r.residue_number, r.insertion_code) for r in gap.generated_residues}]) for gap in request.gaps],
            "assessment": asdict(assessment), "qualification": False})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "independent-results.json", {"cohort": "two private masks; one seed per backend; not qualification", "rows": checked,
        "result_sha256": digest(root / "final-results.json"), "plan_sha256": digest(root / "runs/plan.json"),
        "audit_script_sha256": digest(Path(__file__))})
    lines = ["# Terminal / joint-pair seed7 pilot readback", "", "One seed per cell, not a backend ranking or broad qualification.", "",
        "| Case | Backend | Restored raw geometry | Repair | Final geometry/custody/L | Heavy RMSD Å | Backbone RMSD Å |",
        "|---|---|---|---|---|---:|---:|"]
    for row in checked:
        passed = row["final_geometry_passed"] and not row["final_admission_failures"] and row["global_all_l"]
        lines.append(f"| {row['case']} | {row['backend']} | {row['raw_restored_geometry_passed']} | {row['repaired']} | "
            f"{passed} | {row['final_heldout_heavy_rmsd']:.3f} | {row['final_heldout_backbone_rmsd']:.3f} |")
    (output / "README.md").write_text("\n".join(lines)+"\n")
    print(json.dumps([{k: v for k, v in row.items() if k not in {"assessment", "final_path"}} for row in checked], indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve())
