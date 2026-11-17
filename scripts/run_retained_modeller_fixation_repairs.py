#!/usr/bin/env python3
"""Frozen common CPU repair for saved seed101 fixation-arm proposals only."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from scripts import refine_retained_sampler_returns as common
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

CONTROLS = common.ROOT / "retained-modeller-fixation-seed101-attempt01"
READBACK = common.ROOT / "retained-modeller-fixation-seed101-readback-attempt01/independent-readback.json"


def declare(output):
    evidence = json.loads(READBACK.read_text())
    if (len(evidence["rows"]) != 2 or any(not row["technical_completion"]
            or not row["native_initial_generated_selection_verified"] or row["source_restored_passed"]
            or not row["restored_global_all_l"] for row in evidence["rows"])):
        raise ValueError("require audited globally-L source-restored geometry failures, not arbitrary rebuilding")
    for name, expected in evidence["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError("native fixation-arm evidence changed")
    files = [Path(__file__), Path(common.__file__), READBACK,
        common.PROJECT / "scripts/run_mac_modeller_comparison.py", common.PROJECT / "scripts/run_1he8_traced_mobility.py",
        common.PROJECT / "scripts/run_1he8_buffered_walls.py", common.PROJECT / "scripts/audit_retained_sampler_return_repairs.py",
        common.PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += [Path(name) for name in evidence["bindings"]]
    files += list((common.PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    output.mkdir(parents=True, exist_ok=False)
    jobs = [{"case_id": case, "case": str(CONTROLS / "cases" / case), "seed": 11,
        "command": [sys.executable, "-m", "scripts.run_retained_modeller_fixation_repairs", "worker", "--case", case,
            "--output", str(output / "execution" / f"{case}-seed-11")]} for case in ("1NON", "2BE1")]
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": 2,
        "native_backend": "MODELLER", "native_protocol": "initial-generated-selection", "native_seed": 101, "optimizer_seed": 11,
        "same_common_repair_policy_as_previous_two_arms": True, "new_proposal_calls": 0,
        "generated_only_then_necessary_sequence2": True, "max_iterations_per_phase": 1000, "restarts": 0,
        "source_reference_unchanged": True, "wall_buffer_angstrom": .05,
        "outer_limits_angstrom": [.5, .5, 1.25], "solvent": False,
        "resource_guard": "serial CPU1800s/job,RSS8GiB,host-pressure,10GiB evidence",
        "native_fixation_not_qualified_by_repair_success": True, "manual_approval_or_publication": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--case", choices=("1NON", "2BE1"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        case = CONTROLS / "cases" / args.case
        common.worker(args.case, args.output.resolve(), request_path=case / "request.json", source_path=case / "source.pdb",
            raw_path=READBACK.parent / args.case / "source-restored.pdb", expected_native_seed=101, native_backend="MODELLER")
