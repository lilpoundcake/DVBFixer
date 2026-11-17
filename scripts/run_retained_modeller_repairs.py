#!/usr/bin/env python3
"""Same bounded common repair, correct restored MODELLER bytes, no new proposals."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from scripts import refine_retained_sampler_returns as common
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

ROOT = common.ROOT
CONTROLS = ROOT / "retained-loop-modeller-controls-attempt02"
READBACK = ROOT / "retained-loop-modeller-controls-readback-attempt02/independent-readback.json"


def declare(output):
    readback = json.loads(READBACK.read_text())
    if not readback["restoration_recovery"] or any(not row["technical_completion"] for row in readback["rows"]):
        raise ValueError("require correctly restored, independently audited native proposals")
    for name, expected in readback["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError("native transport/readback changed")
    files = [Path(__file__), Path(common.__file__), READBACK,
             common.PROJECT / "scripts/audit_retained_sampler_return_repairs.py",
             common.PROJECT / "scripts/run_mac_modeller_comparison.py",
             common.PROJECT / "scripts/run_1he8_traced_mobility.py",
             common.PROJECT / "scripts/run_1he8_buffered_walls.py",
             common.PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += [Path(name) for name in readback["bindings"]]
    files += list((common.PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    output.mkdir(parents=True, exist_ok=False)
    jobs = []
    for row in readback["rows"]:
        entry = row["case"]
        case = CONTROLS / "cases" / entry
        jobs.append({"case_id": entry, "case": str(case), "seed": 11,
            "command": [sys.executable, "-m", "scripts.run_retained_modeller_repairs", "worker", "--case", str(case),
                "--output", str(output / "execution" / f"{entry}-seed-11")]})
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": 2,
        "native_backend": "MODELLER", "native_seed": 42, "optimizer_seed": 11,
        "same_common_repair_as_Protenix_seed11": True, "input": "corrected independent source-restored derivative, never invalid worker restoration",
        "all_generated_atoms_jointly_mobile": True, "generated_only_first_then_necessary_sequence2": True,
        "max_iterations_per_phase": 1000, "restarts": 0, "wall_buffer_angstrom": .05,
        "outer_movement_limits_unchanged": {"backbone": .5, "CA": .5, "heavy": 1.25},
        "new_proposal_calls": 0, "solvent": False, "resource_guard": "serial CPU1800s/job,RSS8GiB,host-pressure,10GiB evidence",
        "no_method_superiority_or_donor_approval": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--case", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        case = args.case.resolve()
        common.worker(case.name, args.output.resolve(), request_path=case / "request.json",
            source_path=case / "source.pdb", raw_path=READBACK.parent / case.name / "source-restored.pdb",
            expected_native_seed=42, native_backend="MODELLER")
