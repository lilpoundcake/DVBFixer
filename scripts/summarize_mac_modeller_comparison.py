#!/usr/bin/env python3
"""Independent final readback and descriptive summary; no sampling/selection."""
from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def summarize(root: Path, output: Path):
    plan = json.loads((root / "plan.json").read_text())
    completion = json.loads((root / "completed.json").read_text())
    rows = json.loads((root / "progress.json").read_text())
    if completion["completed"] != len(plan["jobs"]) or len(rows) != len(plan["jobs"]):
        raise ValueError("cannot summarize an incomplete cohort as completed")
    groups = defaultdict(list)
    readbacks = []
    for row in rows:
        if not row["completed"] or row["returncode"]:
            raise ValueError("technical failures must not disappear from summary")
        case, result = Path(row["case"]), row["result"]
        source_path, final_path = case / "source.pdb", Path(result["final_path"])
        for path in (source_path, case / "evaluation-mask.json", case / "reference-projection-NOT-CONDITIONING.pdb"):
            if digest(path) != plan["source_code_and_input_sha256"][str(path)]:
                raise ValueError("declared case bytes changed")
        if digest(final_path) != result["final_sha256"]:
            raise ValueError("final bytes changed")
        source, text = source_path.read_text(), final_path.read_text()
        def fixed_lines(value):
            return {(line[21], line[22:26].strip(), line[26], line[12:16]): line
                    for line in value.splitlines(keepends=True) if line.startswith("ATOM  ")}
        original, final = fixed_lines(source), fixed_lines(text)
        if any(final.get(key) != line for key, line in original.items()):
            raise ValueError("observed ATOM records are not byte-identical")
        request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
        assessed = assessment_adapter(request, source)(text)
        geometry = assessed.geometry_passed and not assessed.admission_failures
        structure = _structure(text, "independent-final.pdb")
        all_l = True
        try:
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
        except Exception:
            all_l = False
        if (geometry and all_l) != result["final_source_geometry_and_l_passed"]:
            raise ValueError("independent final acceptance disagrees with worker")
        readbacks.append({"case": row["case_id"], "seed": row["seed"], "source_atom_records_byte_identical": True,
            "final_sha256": digest(final_path), "global_all_l": all_l, "geometry_passed": geometry})
        groups[row["case_id"]].append(row)
    summary = []
    for case, values in groups.items():
        if sorted(row["seed"] for row in values) != sorted(DEVELOPMENT_SEEDS):
            raise ValueError("missing/repeated declared seeds")
        results = [row["result"] for row in values]
        def stats(key, entries=results):
            numbers = [result[key] for result in entries]
            return {"median": statistics.median(numbers), "min": min(numbers), "max": max(numbers)}
        summary.append({"case": case, "n": len(values), "raw_restored_geometry_passes": sum(r["raw_restored_geometry_passed"] for r in results),
            "final_source_geometry_and_l_passes": sum(r["final_source_geometry_and_l_passed"] for r in results),
            "repaired": sum(r["repaired"] for r in results), "all_outcome_heavy_rmsd": stats("final_heldout_heavy_rmsd"),
            "all_outcome_backbone_rmsd": stats("final_heldout_backbone_rmsd"),
            "max_native_observed_drift": max(r["native_observed_drift"]["max"] for r in results),
            "seconds_total": sum(r["elapsed_seconds"] for r in values), "peak_rss_bytes": max(r["peak_rss_bytes"] for r in values),
            "seed7_heavy_rmsd": next(r["final_heldout_heavy_rmsd"] for r in results if r["seed"] == 7),
            "distinct_final_outputs": len({r["final_sha256"] for r in results})})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "summary.json", {"stage": "retrospective private seeded one-loop MODELLER comparator",
        "declared_outcomes": len(plan["jobs"]), "all_declared_outcomes_included_in_rmsd_summaries": True,
        "rejected_structures_not_accepted_for_accuracy": True,
        "groups": summary, "independent_readback": readbacks, "plan_sha256": digest(root / "plan.json"), "qualification": False})
    lines = ["# MODELLER comparison: retained results", "", "Private one-loop seeded comparator, not production-default MODELLER qualification.",
        "All six outcomes, including geometry rejections, enter descriptive RMSD medians. No reference-based candidate selection.", "",
        "| Mask | Restored raw passes | Final passes | Median heavy RMSD, Å | Median backbone RMSD, Å |",
        "|---|---:|---:|---:|---:|"]
    for group in summary:
        lines.append(f"| {group['case']} | {group['raw_restored_geometry_passes']}/6 | {group['final_source_geometry_and_l_passes']}/6 | "
            f"{group['all_outcome_heavy_rmsd']['median']:.3f} | {group['all_outcome_backbone_rmsd']['median']:.3f} |")
    lines += ["", "Every final structure has independent source-record custody and global-L readback; geometry rejections remain rejected.",
        "Native MODELLER drift is measured separately and not excused by final source restoration.",
        "ML cells with only seed7 must not be compared with six-seed MODELLER medians as matched distributions."]
    if plan.get("technical_amendment") is not None:
        lines.append("The original 1UBQ gap10-13 seed7 was recovered from its unchanged native file after metadata serialization failed; it was not resampled.")
    (output / "README.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    summarize(args.root.resolve(), args.output.resolve())
