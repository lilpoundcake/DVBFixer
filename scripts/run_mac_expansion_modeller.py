#!/usr/bin/env python3
"""The unchanged six-seed MODELLER comparator on the two new admitted masks."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save
from scripts.run_mac_modeller_comparison import execute


def declare(pilot: Path, baseline: Path, output: Path):
    original = json.loads((baseline / "plan.json").read_text())
    for name, expected in original["source_code_and_input_sha256"].items():
        if "/scripts/" in name or "/src/" in name:
            if digest(Path(name)) != expected:
                raise ValueError(f"baseline method changed: {name}")
    admitted = json.loads((pilot / "inputs/preparation-results.json").read_text())
    if len(admitted) != 2 or any(not row["ready"] for row in admitted):
        raise ValueError("exactly the two declared expansion masks are required")
    output.mkdir(parents=True, exist_ok=False)
    jobs, bindings = [], {}
    for row in admitted:
        source_case = pilot / "inputs" / row["case_id"]
        case = output / "cases" / row["case_id"]
        case.mkdir(parents=True)
        manifest = json.loads((source_case / "manifest.json").read_text())
        for name, expected in manifest["files"].items():
            if digest(source_case / name) != expected:
                raise ValueError("declared expansion source/reference/mask changed")
        for name in ("source.pdb", "reference-projection-NOT-CONDITIONING.pdb", "evaluation-mask.json", "manifest.json"):
            (case / name).write_bytes((source_case / name).read_bytes())
            bindings[str(case / name)] = digest(case / name)
        jobs.extend({"case": str(case), "case_id": case.name, "seed": seed} for seed in DEVELOPMENT_SEEDS)
    bindings.update({name: expected for name, expected in original["source_code_and_input_sha256"].items()
                     if "/scripts/" in name or "/src/" in name})
    bindings[str(Path(__file__).resolve())] = digest(Path(__file__))
    plan = {**original, "phase": "same comparator on prospectively frozen terminal/pair masks; not validation",
        "jobs": jobs, "count": 12, "technical_amendment": None,
        "source_code_and_input_sha256": bindings, "previous_baseline_plan_sha256": digest(baseline / "plan.json"),
        "scope": "two fixed expansion masks x six development seeds, all settings unchanged, no heldout-based tuning"}
    save(output / "plan.json", plan)
    print(json.dumps({"jobs": 12, "plan": str(output / "plan.json")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        if args.pilot is None or args.baseline is None:
            parser.error("declare requires --pilot and --baseline")
        declare(args.pilot.resolve(), args.baseline.resolve(), args.output.resolve())
    else:
        execute(args.output.resolve())
