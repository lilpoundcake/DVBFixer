#!/usr/bin/env python3
"""Independent readback of a completed frozen cohort; never selects candidates."""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def atom_records(text):
    result = {}
    for line in text.splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = line[21], line[22:26], line[26], line[12:16]
        if key in result:
            raise ValueError("duplicate atom in independent readback")
        result[key] = line
    return result


def check_records(source, final):
    originals, finals = atom_records(source), atom_records(final)
    if not originals or any(finals.get(key) != line for key, line in originals.items()):
        raise ValueError("observed atom records changed, including terminal OXT")
    return sum(line[12:16].strip() == "OXT" for line in originals.values())


def check_keys(rows, expected):
    keys = [(r["case"], r["backend"], r["seed"]) for r in rows]
    if len(keys) != len(set(keys)) or set(keys) != set(expected):
        raise ValueError("incomplete, repeated or undeclared validation outcomes")


def audit(root, output):
    complete = json.loads((root / "completed.json").read_text())
    ml_complete = json.loads((root / "ML-completed.json").read_text())
    cpu_complete = json.loads((root / "cpu/completed.json").read_text())
    if (complete["ML_executed"] != 90 or complete["MODELLER_executed"] != 30
            or ml_complete["executed"] != 90 or cpu_complete["completed"] != 30
            or cpu_complete["executed"] != 30):
        raise ValueError("validation completion receipts do not cover the frozen budget")
    plan = json.loads((root / "plan.json").read_text())
    if digest(root / "policy.json") != plan["policy_sha256"]:
        raise ValueError("validation policy changed")
    policy = json.loads((root / "policy.json").read_text())
    if digest(root / "jobs.json") != policy["jobs_sha256"]:
        raise ValueError("validation matrix changed")
    for name, expected in plan["code_and_input_files"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"frozen input/code changed: {name}")
    jobs = json.loads((root / "jobs.json").read_text())
    cpu_plan = json.loads((root / "cpu/plan.json").read_text())
    if digest(root / "cpu/plan.json") != plan["cpu_plan_sha256"]:
        raise ValueError("CPU matrix changed")
    jobs += [{**job, "model": "modeller", "work": str(root / "cpu/execution" / f"{job['case_id']}-seed-{job['seed']}")}
             for job in cpu_plan["jobs"]]
    keys = [(Path(j["case"]).name, j["model"], j["seed"]) for j in jobs]
    if len(keys) != 120 or len(set(keys)) != 120:
        raise ValueError("declared matrix is not 120 distinct outcomes")
    rows = [json.loads((Path(j["work"]) / "result.json").read_text()) for j in jobs]
    check_keys(rows, keys)
    by_key = {(r["case"], r["backend"], r["seed"]): r for r in rows}
    checked = []
    for job, key in zip(jobs, keys, strict=True):
        row = by_key[key]
        case, work, final = Path(job["case"]), Path(job["work"]), Path(row["final_path"])
        if final.resolve() not in {(work / "common-repair/refined.pdb").resolve(), (work / "audit/source-restored.pdb").resolve()}:
            raise ValueError("final output is outside its declared candidate")
        if digest(final) != row["final_sha256"]:
            raise ValueError("final output changed")
        execution = json.loads((work / "execution.json").read_text())
        if execution.get("guard_stop_reason") or execution.get("returncode") != 0:
            raise ValueError("execution receipt is not successful")
        request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
        source, text = (case / "source.pdb").read_text(), final.read_text()
        oxt_count = check_records(source, text)
        positions, reference = inventory(final), inventory(case / "reference-projection-NOT-CONDITIONING.pdb")
        transport = json.loads((work / "audit/transport.json").read_text())
        if transport != row["native_transport"] or transport["source_sha256"] != digest(case / "source.pdb"):
            raise ValueError("native transport/source custody changed")
        if job["model"] == "modeller":
            natives = [work / "proposal/mapped-native.pdb"]
        else:
            pattern = {"rfd3": "*.cif.gz", "boltz": "*_model_0.pdb", "protenix": "candidate.pdb"}[job["model"]]
            natives = list((work / "run/predictions").rglob(pattern))
        if len(natives) != 1 or digest(natives[0]) != transport["prediction_sha256"]:
            raise ValueError("native prediction changed or is ambiguous")
        source_positions = inventory(case / "source.pdb")
        axis = sorted(set(request.fixed_atoms) | set(request.generated_atoms))
        with np.load(work / "audit/aligned-arrays.npz", allow_pickle=False) as arrays:
            coordinates, represented = arrays["coordinates"], arrays["native_represented"]
            if coordinates.shape != (len(axis), 3) or represented.shape != (len(axis),):
                raise ValueError("native aligned-array inventory changed")
            fixed_indices = [i for i, atom in enumerate(axis) if atom in source_positions and represented[i]]
            drift = np.array([np.linalg.norm(coordinates[i]-source_positions[axis[i]][2]) for i in fixed_indices])
        if (len(fixed_indices) != transport["fixed_atoms_represented_in_native"]
                or abs(float(np.sqrt(np.mean(drift**2)))-transport["native_observed_rms_drift"]) > 1e-10
                or abs(float(drift.max())-transport["native_observed_max_drift"]) > 1e-10):
            raise ValueError("independent native observed drift disagrees")
        if set(positions) != set(request.fixed_atoms) | set(request.generated_atoms):
            raise ValueError("final requested atom inventory changed")
        if any(not np.isfinite(value[2]).all() for value in positions.values()):
            raise ValueError("nonfinite final coordinates")
        assessment = assessment_adapter(request, source)(text)
        raw = assessment_adapter(request, source)((work / "audit/source-restored.pdb").read_text())
        all_l = True
        structure = _structure(text, "independent-fixed-validation.pdb")
        try:
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
        except ChiralityError:
            all_l = False
        accepted = (row["repair_error"] is None and assessment.geometry_passed
                    and not assessment.admission_failures and all_l)
        flags = {"final_geometry_passed": assessment.geometry_passed, "global_all_l": all_l,
                 "raw_restored_geometry_passed": raw.geometry_passed, "final_accepted": accepted}
        if any(row[name] != value for name, value in flags.items()):
            raise ValueError("independent acceptance disagrees with candidate result")
        if list(assessment.admission_failures) != row["final_admission_failures"]:
            raise ValueError("independent admission findings disagree")
        errors = {a: float(np.sum((positions[a][2]-reference[a][2])**2)) for a in request.generated_atoms}
        heavy = float(np.sqrt(np.mean(list(errors.values()))))
        backbone = float(np.sqrt(np.mean([e for a, e in errors.items() if a.atom_name in {"N", "CA", "C", "O"}])))
        if abs(heavy-row["heldout_heavy_rmsd"]) > 1e-10 or abs(backbone-row["heldout_backbone_rmsd"]) > 1e-10:
            raise ValueError("independent held-out RMSD disagrees")
        checked.append({**row, **flags, "observed_OXT_records_preserved": oxt_count,
                        "independent_heavy_rmsd": heavy, "independent_backbone_rmsd": backbone})
    groups = defaultdict(list)
    for row in checked:
        groups[row["case"], row["backend"]].append(row)
    summaries = []
    for (case, backend), members in sorted(groups.items()):
        if Counter(r["seed"] for r in members) != Counter(policy["seeds"]):
            raise ValueError("cell does not retain exactly five declared seeds")
        summaries.append({"case": case, "backend": backend, "count": len(members),
            "accepted": sum(r["final_accepted"] for r in members),
            "raw_passed": sum(r["raw_restored_geometry_passed"] for r in members),
            "repaired": sum(r["repaired"] for r in members),
            "all_outcome_heavy_median": float(np.median([r["independent_heavy_rmsd"] for r in members])),
            "all_outcome_backbone_median": float(np.median([r["independent_backbone_rmsd"] for r in members])),
            "distinct_final_hashes": len({r["final_sha256"] for r in members})})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "independent-results.json", {"rows": checked, "cells": summaries,
        "policy_sha256": plan["policy_sha256"], "auditor_sha256": digest(Path(__file__)),
        "qualification": False, "selection_performed": False, "statistical_units": "two entries, not 120 proteins"})
    lines = ["# Fixed-source validation independent readback", "",
        "All outcomes retained, including rejections. Two entries; not public qualification.", "",
        "| Case | Backend | Raw passes /5 | Final passes /5 | Repairs | Heavy median Å | Backbone median Å |",
        "|---|---|---:|---:|---:|---:|---:|"]
    for cell in summaries:
        lines.append(f"| {cell['case']} | {cell['backend']} | {cell['raw_passed']} | {cell['accepted']} | {cell['repaired']} | "
                     f"{cell['all_outcome_heavy_median']:.3f} | {cell['all_outcome_backbone_median']:.3f} |")
    (output / "README.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({"independently_read": len(checked), "accepted": sum(r["final_accepted"] for r in checked), "qualification": False}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve())
