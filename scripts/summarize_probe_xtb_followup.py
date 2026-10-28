"""Index the bounded follow-up without converting optimizer exits to acceptance."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def physical_clash_changes(new_findings, lost_findings):
    # Raw diagnostic keys retain severity/location metadata. A threshold
    # crossing on the same atom pair is not a newly appearing physical pair.
    new = {tuple(tuple(atom) for atom in finding[2]) for finding in new_findings}
    lost = {tuple(tuple(atom) for atom in finding[2]) for finding in lost_findings}
    return len(new - lost), len(lost - new)


def summarize_context(directory):
    parent = json.loads((directory / "report.json").read_text())
    rows = []
    for job in parent["jobs"]:
        report = directory / job["name"] / "report.json"
        data = json.loads(report.read_text()) if report.exists() else {}
        delta = data.get("geometry_deltas", {})
        new_pairs, lost_pairs = physical_clash_changes(data.get("new_clash_pairs", []), data.get("lost_clash_pairs", []))
        evaluations = directory / job["name"] / "evaluations.jsonl"
        history = [json.loads(line) for line in evaluations.read_text().splitlines()] if evaluations.exists() else []
        rows.append({"name": job["name"], "returncode": job.get("returncode"), "timeout": job.get("timeout", False),
                     "seconds": job["seconds"], "atoms": data.get("atoms"), "mobile_atoms": data.get("mobile_atoms"),
                     "baseline_sha256": data.get("baseline_sha256"), "threads": data.get("threads"),
                     "library_sha256": parent["library_sha256"], "numeric_libraries": data.get("numeric_libraries"),
                     "converged": data.get("converged", False), "scipy_success": data.get("scipy_success"),
                     "iterations": data.get("iterations"), "optimizer_message": data.get("optimizer_message"),
                     "gradient_max_hartree_angstrom": data.get("projected_gradient_max_hartree_angstrom"),
                     "screens_pass": data.get("screens_pass", False), "gates": data.get("gates"),
                     "raw_frozen_drift_max_angstrom": data.get("raw_frozen_drift_max_angstrom"),
                     "clashes_before": data.get("clash_pairs_before"), "clashes_after": data.get("clash_pairs_after"),
                     "new_clashes": new_pairs if "new_clash_pairs" in data else None,
                     "lost_clashes": lost_pairs if "lost_clash_pairs" in data else None,
                     "completed_evaluations": len(history), "initial_evaluation": history[0] if history else None,
                     "last_completed_evaluation": history[-1] if history else None,
                     "relative_geometry_changes": {k: len(delta[k]) for k in ["bonds_changed_over_0_1_angstrom",
                         "angles_changed_over_10_degrees", "potential_stereovolume_flips"]} if delta else None,
                     "pockets": [{k: len(v) if k.endswith("contacts") else v for k, v in pocket.items()}
                                 for pocket in data["pockets"]] if "pockets" in data else None})
    invocation = parent.get("invocation", [])
    expected = (1 if "--baseline" in invocation else 3) * (1 if "--single-arm" in invocation else 3)
    return {"directory": str(directory), "completed_jobs": len(parent["jobs"]), "expected_jobs": expected,
            "execution_complete": len(parent["jobs"]) == expected, "runs": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context", type=Path, action="append", required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--api-control", type=Path, action="append", default=[])
    parser.add_argument("--native-sp", type=Path)
    parser.add_argument("--restraint-comparison", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output must not exist")
    probe = json.loads((args.probe / "report.json").read_text())
    report = {"context": [summarize_context(d) for d in args.context],
              "probe": {"valid_input_runs": len(probe["runs"]), "all_complete": all(r["status"] == "complete" for r in probe["runs"]),
                        "format_control": probe["format_control"],
                        "failure_controls": [{k: row.get(k) for k in ["case", "status", "actual_engine", "fallback_findings_count", "failure"]}
                                             for row in probe["failures"]],
                        "runs": [{"input": Path(r["source"]).name, "aligned": r["aligned"], "nuclear": r["nuclear"],
                                  "overlap_pairs": len(r["findings"]), "ambiguous_hb_overlap_pairs": r.get("ambiguous_hb_overlap_pairs")}
                                 for r in probe["runs"]]},
              "note": "Comparative screens and bounded failures; neither a general chemistry qualification nor an affinity estimate"}
    report["api_controls"] = {str(p): json.loads((p / "validation.json").read_text()) for p in args.api_control}
    if args.native_sp:
        native = json.loads((args.native_sp / "report.json").read_text())
        raw = (args.native_sp / "stdout.raw").read_text(errors="replace")
        energies = re.findall(r"TOTAL ENERGY\s+([-+\d.Ee]+)", raw)
        report["native_full_singlepoint"] = {**native, "energy_hartree": float(energies[-1]) if energies else None,
                                             "warnings": [line.strip() for line in raw.splitlines() if "WARNING" in line]}
    if args.restraint_comparison:
        report["matched_restraint_ablation"] = json.loads(args.restraint_comparison.read_text())
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
