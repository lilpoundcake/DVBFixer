#!/usr/bin/env python3
"""Read back both E arms after execution, retaining all six seeds per cell."""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_fixed_validation import check_records
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save

from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def check_masks(receipt):
    if not receipt["no_structural_template_features"] or not receipt["masks"]:
        raise ValueError("no-template consumer proof missing")
    if any(row["nonzero"] != 0 for row in receipt["masks"].values()):
        raise ValueError("no-template arm consumed populated masks")


def audit(experimental, baseline, output):
    if json.loads((experimental / "completed.json").read_text())["executed"] != 12:
        raise ValueError("E budget incomplete")
    if json.loads((baseline / "completed.json").read_text())["completed"] != 12:
        raise ValueError("matched baseline incomplete")
    experiment_plan = json.loads((experimental / "plan.json").read_text())
    baseline_plan = json.loads((baseline / "plan.json").read_text())
    for bindings in (experiment_plan["files"], baseline_plan["source_code_and_input_sha256"]):
        for name, expected in bindings.items():
            if digest(Path(name)) != expected:
                raise ValueError(f"declared E/baseline file changed: {name}")
    prospective = json.loads((experimental / "progress.json").read_text())
    if len(prospective) != 12 or len({(r["backend"], r["seed"]) for r in prospective}) != 12:
        raise ValueError("E outcomes missing or repeated")
    no_template = {(r["backend"], r["seed"]): r for r in prospective}
    jobs = [("no-template-no-MSA", j) for j in experiment_plan["jobs"]]
    jobs += [("source-template-no-MSA", j) for j in baseline_plan["jobs"]]
    checked = []
    for arm, job in jobs:
        backend = job.get("model", job.get("backend"))
        if arm == "no-template-no-MSA":
            work, case = Path(job["work"]), Path(job["case"])
            row = no_template[backend, job["seed"]]
            source_audit = work / "audit"
            check_masks(json.loads((work / "run/sequence-evidence-mask.json").read_text()))
            execution = json.loads((work / "run/execution.json").read_text())
            failure = row["repair_failure"]
        else:
            work = baseline / "execution" / f"{job['case_id']}-seed-{job['seed']}"
            row = json.loads((work / "result.json").read_text())
            case = Path(experiment_plan["jobs"][0]["case"])
            source_audit = Path(job["case"])
            execution = json.loads((work / "execution.json").read_text())
            failure = row["repair_failure"]
            if digest(source_audit / "source-restored.pdb") != row["input_sha256"]:
                raise ValueError("retained baseline input changed")
        if execution.get("guard_stop_reason") or execution["returncode"] != 0:
            raise ValueError("E/baseline execution unsuccessful")
        final = Path(row["final_path"])
        if digest(final) != row["final_sha256"]:
            raise ValueError("E/baseline output changed")
        source, text = (case / "source.pdb").read_text(), final.read_text()
        if (source_audit / "source.pdb").read_bytes() != (case / "source.pdb").read_bytes():
            raise ValueError("E arms have different observations")
        request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
        audit_request = DiffusionRequest.from_json((source_audit / "evaluation-mask.json").read_text())
        if set(request.fixed_atoms) != set(audit_request.fixed_atoms) or set(request.generated_atoms) != set(audit_request.generated_atoms):
            raise ValueError("E arms have different ownership masks")
        check_records(source, text)
        positions, reference = inventory(final), inventory(case / "reference-projection-NOT-CONDITIONING.pdb")
        if set(positions) != set(request.fixed_atoms) | set(request.generated_atoms):
            raise ValueError("E final inventory differs")
        if any(not np.isfinite(value[2]).all() for value in positions.values()):
            raise ValueError("nonfinite E output")
        assess = assessment_adapter(request, source)
        final_assessment = assess(text)
        raw_assessment = assess((source_audit / "source-restored.pdb").read_text())
        structure = _structure(text, "sequence-evidence-independent.pdb")
        all_l = True
        try:
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
        except ChiralityError:
            all_l = False
        errors = {a: float(np.sum((positions[a][2]-reference[a][2])**2)) for a in request.generated_atoms}
        heavy = float(np.sqrt(np.mean(list(errors.values()))))
        backbone = float(np.sqrt(np.mean([v for a, v in errors.items() if a.atom_name in {"N", "CA", "C", "O"}])))
        raw_report = json.loads((source_audit / "audit.json").read_text())
        checked.append({"arm": arm, "backend": backend, "seed": job["seed"], "final_path": str(final),
            "final_sha256": digest(final), "raw_geometry_passed": raw_assessment.geometry_passed,
            "repaired": row["repaired"], "repair_failure": failure, "final_geometry_passed": final_assessment.geometry_passed,
            "final_admission_failures": list(final_assessment.admission_failures), "global_all_l": all_l,
            "observations_byte_identical": True, "heldout_heavy_rmsd": heavy, "heldout_backbone_rmsd": backbone,
            "native_observed_drift": raw_report["native_observed_drift_after_observed_only_fit_angstrom"],
            "accepted": failure is None and all_l and final_assessment.geometry_passed and not final_assessment.admission_failures})
    groups = defaultdict(list)
    for row in checked:
        groups[row["backend"], row["arm"]].append(row)
    if len(groups) != 4:
        raise ValueError("E matrix incomplete")
    cells = []
    for (backend, arm), members in sorted(groups.items()):
        if Counter(r["seed"] for r in members) != Counter(DEVELOPMENT_SEEDS):
            raise ValueError("E cell lacks exactly six declared seeds")
        cells.append({"backend": backend, "arm": arm, "count": 6, "raw_passes": sum(r["raw_geometry_passed"] for r in members),
            "final_passes": sum(r["accepted"] for r in members), "repairs": sum(r["repaired"] for r in members),
            "all_outcome_heavy_median": float(np.median([r["heldout_heavy_rmsd"] for r in members])),
            "all_outcome_backbone_median": float(np.median([r["heldout_backbone_rmsd"] for r in members]))})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "independent-results.json", {"rows": checked, "cells": cells, "qualification": False,
        "new_ML_calls": 12, "baseline_resampled": False, "comparison": "one retrospective 1UBQ mask, no MSA in either arm",
        "auditor_sha256": digest(Path(__file__)), "reference_based_selection": False})
    lines = ["# Sequence-evidence E subset: independent readback", "",
        "One retrospective 1UBQ mask; all outcomes retained; no public qualification.", "",
        "| Backend | Arm | Raw /6 | Final /6 | Repairs | Heavy median Å | Backbone median Å |",
        "|---|---|---:|---:|---:|---:|---:|"]
    for cell in cells:
        lines.append(f"| {cell['backend']} | {cell['arm']} | {cell['raw_passes']} | {cell['final_passes']} | {cell['repairs']} | "
                     f"{cell['all_outcome_heavy_median']:.3f} | {cell['all_outcome_backbone_median']:.3f} |")
    (output / "README.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(cells, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experimental", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    audit(args.experimental.resolve(), args.baseline.resolve(), args.output.resolve())
