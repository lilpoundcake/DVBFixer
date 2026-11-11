#!/usr/bin/env python3
"""Stage-aware six-seed ML/MODELLER expansion summary, with no new selection."""
from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save


def aggregate(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["case"], row["backend"]].append(row)
    summary = []
    for (case, backend), items in sorted(groups.items()):
        if sorted(row["seed"] for row in items) != sorted(DEVELOPMENT_SEEDS):
            raise ValueError("each admitted cell must contain every development seed exactly once")
        heavy = [row["final_heldout_heavy_rmsd"] for row in items]
        backbone = [row["final_heldout_backbone_rmsd"] for row in items]
        summary.append({"case": case, "backend": backend, "n": 6,
            "raw_restored_geometry_passes": sum(row["raw_restored_geometry_passed"] for row in items),
            "final_source_geometry_and_l_passes": sum(row["final_geometry_passed"] and not row["final_admission_failures"]
                and row["global_all_l"] and row["observed_atom_records_byte_identical"] for row in items),
            "repairs": sum(row["repaired"] for row in items),
            "heavy_rmsd": {"median": statistics.median(heavy), "min": min(heavy), "max": max(heavy)},
            "backbone_rmsd": {"median": statistics.median(backbone), "min": min(backbone), "max": max(backbone)},
            "distinct_final_outputs": len({row["final_sha256"] for row in items}),
            "native_observed_max_drift": max(row["native_observed_drift"]["max"] for row in items)})
    return summary


def summarize(root: Path, output: Path):
    inputs = [root / "expansion-seed7-independent-readback-attempt02/independent-results.json",
        root / "expansion-development-independent-readback/independent-results.json",
        root / "expansion-modeller-independent-summary-attempt02/summary.json",
        root / "expansion-modeller-attempt01/progress.json"]
    rows = []
    for path in inputs[:2]:
        rows.extend(json.loads(path.read_text())["rows"])
    independent = json.loads(inputs[2].read_text())["independent_readback"]
    audits = {(r["case"], r["seed"]): r for r in independent}
    for execution in json.loads(inputs[3].read_text()):
        result = execution["result"]
        checked = audits[result["case"], result["seed"]]
        if result["final_sha256"] != checked["final_sha256"]:
            raise ValueError("MODELLER independent custody differs")
        rows.append({"case": result["case"], "seed": result["seed"], "backend": "modeller",
            "raw_restored_geometry_passed": result["raw_restored_geometry_passed"],
            "final_geometry_passed": checked["geometry_passed"], "final_admission_failures": result["final_admission_failures"],
            "global_all_l": checked["global_all_l"], "observed_atom_records_byte_identical": checked["source_atom_records_byte_identical"],
            "repaired": result["repaired"], "final_sha256": result["final_sha256"], "final_path": result["final_path"],
            "native_observed_drift": result["native_observed_drift"],
            "final_heldout_heavy_rmsd": result["final_heldout_heavy_rmsd"],
            "final_heldout_backbone_rmsd": result["final_heldout_backbone_rmsd"]})
    if len(rows) != 48 or len({(r["case"], r["backend"], r["seed"]) for r in rows}) != 48:
        raise ValueError("expected two masks x four backends x six unique seeds")
    for row in rows:
        if digest(Path(row["final_path"])) != row["final_sha256"]:
            raise ValueError("an independently audited final structure changed")
    noise = defaultdict(set)
    for directory in (root / "expansion-seed7-pilots-attempt01/runs", root / "expansion-development-extension-attempt01/runs"):
        plan = json.loads((directory / "plan.json").read_text())
        for job in plan["jobs"]:
            if job["arm"] != "rfd3":
                continue
            consumer = Path(job["work"]) / "run/transport/consumer.json"
            receipt = json.loads(consumer.read_text())
            noise[Path(job["case"]).name].add(receipt["first_noisy_sha256"])
            inputs.append(consumer)
    if len(noise) != 2 or any(len(hashes) != 6 for hashes in noise.values()):
        raise ValueError("RFD3 lacks six distinct first-noise arrays per mask")
    summary = aggregate(rows)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "summary.json", {"scope": "two masks of one cropped protein; descriptive development comparison",
        "rows": rows, "groups": summary, "rfd3_distinct_first_noise_arrays": {case: len(v) for case, v in noise.items()},
        "input_sha256": {str(path): digest(path) for path in inputs},
        "reserved_validation_seeds_used": False, "heldout_based_selection": False, "qualification": False})
    lines = ["# Six-seed terminal/pair comparison", "", "Two masks of one protein, not independent structure-level qualification.",
        "Raw counts mean source-restored geometry, NOT native source-preserving acceptance (especially Boltz/MODELLER).", "",
        "| Case | Backend | Raw restored passes | Final passes | Median heavy RMSD Å | Range Å | Median backbone RMSD Å |",
        "|---|---|---:|---:|---:|---|---:|"]
    for group in summary:
        h, b = group["heavy_rmsd"], group["backbone_rmsd"]
        lines.append(f"| {group['case']} | {group['backend']} | {group['raw_restored_geometry_passes']}/6 | "
            f"{group['final_source_geometry_and_l_passes']}/6 | {h['median']:.3f} | {h['min']:.3f}–{h['max']:.3f} | {b['median']:.3f} |")
    lines += ["", "Seed7 reused without resampling. No reference-based best-of-six selection. Same integer seeds do not mean matched model noise.",
        "MODELLER is the private seeded one-fast-loop protocol, not the production-default two-loop pipeline.",
        "All 48 outputs retain exact observations and pass final geometry/global-L readback; this is not topology or trajectory qualification."]
    (output / "README.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    summarize(args.root.resolve(), args.output.resolve())
