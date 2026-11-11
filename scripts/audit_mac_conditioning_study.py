#!/usr/bin/env python3
"""Independent A/C/D assembly readback and platform-matched baseline derivatives."""
from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_conditioning_contracts import materialize, native_coordinates
from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save
from scripts.refine_mac_reconstruction_smoke import run as refine

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-reconstruction-comparison-2026-11-11"
OLD = PROJECT / "tmp/mac-reconstruction-2026-11-11"


def baseline(output):
    import fcntl
    output.mkdir(parents=True, exist_ok=False)
    save(output / "declaration.json", {"purpose": "same-CPU derivative of old full-heavy raw proposals; original Reference results unchanged",
        "seeds": DEVELOPMENT_SEEDS, "sampling_calls": 0, "max_repairs": 2, "platform": "CPU", "optimizer_seed": 7,
        "iterations": 1000, "restarts": 0, "clean_not_rebuilt": True, "no_reference_to_refiner": True})
    rows = []
    with (PROJECT / ".artifacts/mac-reconstruction-2026-11-11/ml.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for seed in DEVELOPMENT_SEEDS:
            audit = OLD / f"rfd3-explicit-seed{seed}/audit"
            work = output / f"seed-{seed}"
            refine(audit, work, platform="CPU")
            final = work / "refined.pdb"
            repaired = final.is_file()
            if not repaired:
                final = audit / "source-restored.pdb"
            rows.append({"arm": "AC-full-heavy-baseline", "seed": seed,
                "case": str(ROOT / "modeller-attempt02/cases/1UBQ-A1-72-gap10-13"),
                "raw_path": str(audit / "source-restored.pdb"), "repaired": repaired,
                "final_path": str(final), "final_sha256": digest(final)})
    if sum(row["repaired"] for row in rows) > 2:
        raise ValueError("baseline repair count exceeded the separately declared two derivatives")
    save(output / "results.json", rows)


def check(row):
    case, final = Path(row["case"]), Path(row["final_path"])
    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    if digest(final) != row["final_sha256"] or digest(case / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("final/source custody changed")
    source, text = (case / "source.pdb").read_text(), final.read_text()
    originals = [line for line in source.splitlines(keepends=True) if line.startswith("ATOM  ")]
    if any(line not in text.splitlines(keepends=True) for line in originals):
        raise ValueError("original observed ATOM bytes changed")
    raw = Path(row["raw_path"]).read_text()
    assessor = assessment_adapter(request, source)
    before, after = assessor(raw), assessor(text)
    structure = _structure(text, "independent-conditioning-final.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    reference, coords = inventory(case / "reference-projection-NOT-CONDITIONING.pdb"), inventory(final)
    def rms(axis):
        return float(np.sqrt(np.mean([np.sum((coords[a][2]-reference[a][2])**2) for a in axis])))
    return {**row, "raw_geometry_passed": before.geometry_passed, "final_geometry_passed": after.geometry_passed,
        "final_admission_failures": list(after.admission_failures), "source_atom_bytes_identical": True, "global_all_l": True,
        "heldout_heavy_rmsd": rms(request.generated_atoms),
        "heldout_backbone_rmsd": rms([a for a in request.generated_atoms if a.atom_name in {"N", "CA", "C", "O"}])}


def audit(study, baseline_root, output):
    completion = json.loads((study / "completed.json").read_text())
    if completion["calls"] != 30 or completion["assembled_outcomes"] != 24:
        raise ValueError("incomplete A/C/D cohort")
    plan = json.loads((study / "plan.json").read_text())
    jobs = {(j["arm"], j["seed"]): j for j in plan["jobs"]}
    rows = []
    for row in json.loads((study / "final-results.json").read_text()):
        arm, seed = row["arm"], row["seed"]
        names = ("D-independent-gap0", "D-independent-gap1") if arm == "D-independent-pair" else (arm,)
        groups = [native_coordinates(Path(jobs[name, seed]["inputs"]), Path(jobs[name, seed]["work"]) / "run")[0] for name in names]
        reproduced = materialize(Path(row["case"]), groups)
        raw_path = study / "assembled" / f"{arm}-seed{seed}" / "source-restored.pdb"
        if reproduced != raw_path.read_text():
            raise ValueError("independent native-to-assembly replay differs")
        checked = check({**row, "raw_path": str(raw_path)})
        if abs(checked["heldout_heavy_rmsd"]-row["heldout_heavy_rmsd"]) > 1e-10:
            raise ValueError("independent score differs")
        rows.append(checked)
    rows.extend(check(row) for row in json.loads((baseline_root / "results.json").read_text()))
    combined = json.loads((ROOT / "expansion-six-seed-summary/summary.json").read_text())
    for row in combined["rows"]:
        if row["backend"] != "rfd3" or "pair" not in row["case"]:
            continue
        seed = row["seed"]
        root = ROOT / ("expansion-seed7-pilots-attempt01" if seed == 7 else "expansion-development-extension-attempt01")
        joint_plan = json.loads((root / "runs/plan.json").read_text())
        job = next(j for j in joint_plan["jobs"] if j["arm"] == "rfd3" and j["seed"] == seed and "pair" in j["case"])
        rows.append(check({"arm": "D-joint-baseline", "seed": seed, "case": job["case"], "final_path": row["final_path"],
            "final_sha256": row["final_sha256"], "raw_path": str(Path(job["work"]) / "audit/source-restored.pdb"), "repaired": row["repaired"]}))
    groups = defaultdict(list)
    for row in rows:
        groups[row["arm"]].append(row)
    summary = []
    for arm, items in sorted(groups.items()):
        if sorted(row["seed"] for row in items) != sorted(DEVELOPMENT_SEEDS):
            raise ValueError("missing/repeated seeds in a conditioning arm")
        summary.append({"arm": arm, "raw_passes": sum(row["raw_geometry_passed"] for row in items),
            "final_passes": sum(row["final_geometry_passed"] and not row["final_admission_failures"] for row in items),
            "median_heavy_rmsd": statistics.median(row["heldout_heavy_rmsd"] for row in items),
            "median_backbone_rmsd": statistics.median(row["heldout_backbone_rmsd"] for row in items)})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "independent-results.json", {"summary": summary, "rows": rows, "qualification": False,
        "new_sample_count": 30, "baseline_resampled": False, "baseline_repair_platform_matched": "CPU",
        "context_and_entity_layout_not_separable": True, "joint_vs_independent_task_sizes_not_equal": True})
    lines = ["# RFD3 conditioning A/C/D: independent readback", "", "One protein; no universal ranking or public qualification.", "",
        "| Arm | Raw passes | Final passes | Median heavy RMSD Å | Median backbone RMSD Å |", "|---|---:|---:|---:|---:|"]
    for group in summary:
        lines.append(f"| {group['arm']} | {group['raw_passes']}/6 | {group['final_passes']}/6 | {group['median_heavy_rmsd']:.3f} | {group['median_backbone_rmsd']:.3f} |")
    lines += ["", "All original observations are byte-identical. Native nuisance sidechains in C are not newly generated ownership.",
        "Spatial context changes entity cuts/layout as well as observations; D changes coupling and complementary-gap task size.",
        "Old full-heavy Reference repairs remain retained; separate CPU derivatives match the downstream platform without resampling."]
    (output / "README.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("baseline", "audit"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--study", type=Path)
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    if args.mode == "baseline":
        baseline(args.output.resolve())
    else:
        audit(args.study.resolve(), args.baseline.resolve(), args.output.resolve())
