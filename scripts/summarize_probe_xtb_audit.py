"""Create a compact scientific index without confusing execution with acceptance."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from audit_xtb_controls import geometry_snapshot, read_xyz


def compare_components(audit, directories, baselines):
    from openmm.app import PDBFile
    from openmm.unit import nanometer

    from dvbfixer.ffutils import build_glycam_system

    baseline = baselines / "hinge-openmm.pdb"
    pdb = PDBFile(str(baseline))
    build_glycam_system(["amber14-all.xml", "amber14/GLYCAM_06j-1.xml", "amber14/tip3pfb.xml"],
                        pdb.topology, pdb.positions)
    before = np.asarray(pdb.positions.value_in_unit(nanometer)) * 10
    separate, joint = before.copy(), before.copy()
    changed = set()
    joint_indices = set()
    baseline_digest = hashlib.sha256(baseline.read_bytes()).hexdigest()
    for directory, target in [(audit / "fixed-glycan-extended", joint), *[(p, separate) for p in directories]]:
        metadata = json.loads(next(directory.glob("*extraction.json")).read_text())
        if metadata["input_sha256"] != baseline_digest:
            raise ValueError("Component baseline digest mismatch")
        run = next(r for r in json.loads((directory / "summary.json").read_text())["runs"]
                   if r["case"].startswith("pdb-") and r["engine"] == "default")
        if not run.get("converged") or run.get("returncode") != 0:
            return {"comparison_available": False, "reason": "Unconverged component"}
        _, coordinates = read_xyz(directory / f"{run['case']}-default/xtbopt.xyz")
        frozen = set(run["frozen_indices_zero_based"])
        for full, sub in metadata["full_to_sub_indices"].items():
            if sub not in frozen:
                target[int(full)] = coordinates[sub]
                if target is separate:
                    if int(full) in changed:
                        raise ValueError("Component runs overlap")
                    changed.add(int(full))
                else:
                    joint_indices.add(int(full))
    if changed != joint_indices:
        raise ValueError("Separate components do not cover joint mobile atoms")
    delta = np.linalg.norm(separate[sorted(changed)] - joint[sorted(changed)], axis=1)
    return {"comparison_available": True, "mobile_atoms": len(changed),
            "joint_vs_separate_rms_angstrom": float(np.sqrt(np.mean(delta**2))),
            "joint_vs_separate_max_angstrom": float(delta.max()),
            "combined_separate_geometry": geometry_snapshot(pdb.topology, separate),
            "note": "Both omit the protein environment; separating components is not a pocket model"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fixed-component-run", type=Path, action="append", default=[])
    parser.add_argument("--baselines", type=Path, help="Baseline directory for frozen-coordinate replays")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output must not exist")
    report = {"execution": json.loads((args.audit / "report.json").read_text()),
              "scientific_acceptance": "Not inferred by this index; review contract and geometry gates", "xtb": {}, "probe": {}}
    for path in sorted(args.audit.glob("*/summary.json")):
        data = json.loads(path.read_text())
        if path.parent.name.startswith("probe-"):
            report["probe"][path.parent.name] = [
                {"input": Path(r["source_input"]).name, "python_findings": len(r["python_findings"]),
                 "probe_overlap_pairs": r["direct-self"]["unique_overlap_pairs_ge_0.4"],
                 "python_seconds": r["python_seconds"], "probe_seconds": r["direct-self"]["seconds"],
                 "bonded_control": r.get("bonded_control"), "directional_control": r.get("directional_control")}
                for r in data["runs"]]
            continue
        rows = []
        for row in data["runs"]:
            selected = {k: row.get(k) for k in ["case", "engine", "returncode", "converged", "seconds",
                "last_cycle_banner", "raw_frozen_drift", "frozen_drift_pass", "frozen_mobile_bonds",
                "whole_structure_geometry", "fragmented_diagonalization_observed", "validation_error", "timeout"]}
            selected["pocket"] = [{k: v for k, v in p.items() if not k.endswith("contacts")} | {
                k: len(v) for k, v in p.items() if k.endswith("contacts")}
                for p in row.get("heterogen_pocket_contacts", [])]
            rows.append(selected)
        report["xtb"][path.parent.name] = rows
    if len(args.fixed_component_run) == 2:
        report["fixed_joint_vs_separate"] = compare_components(
            args.audit, args.fixed_component_run, args.baselines or args.audit / "baselines")
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
