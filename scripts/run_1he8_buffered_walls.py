#!/usr/bin/env python3
"""Prospective finite-wall guard band; original movement acceptance is unchanged."""
from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

from scripts import run_1he8_traced_mobility as traced
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion import boundary_refinement as boundary
from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits

BUFFER_ANGSTROM = .05
CASES = ((0, "direct"), (1, "direct"), (1, "sequence2"), (1, "spatial5-padding1"))


def buffered_limits(validation):
    values = {name: value-BUFFER_ANGSTROM for name, value in asdict(validation).items()}
    if any(value <= 0 for value in values.values()):
        raise ValueError("wall guard band must remain inside positive validation bounds")
    return AnchorMovementLimits(**values)


def worker(index, arm, output):
    original = boundary.refine_generated_region
    calls = []

    def guarded_refine(*args, **kwargs):
        validation = kwargs.get("anchor_movement_limits")
        if validation is None or calls:
            raise ValueError("buffered trial requires exactly one observed-anchor refinement call")
        inner = buffered_limits(validation)
        evidence = {"wall_buffer_angstrom": BUFFER_ANGSTROM, "wall_activation_limits": asdict(inner),
            "validation_limits_unchanged": asdict(validation), "source_references_unchanged": True,
            "finite_penalty_not_a_feasibility_guarantee": True,
            "only_modified_refiner_argument": "anchor_movement_limits for wall construction",
            "post_refinement_movement_check": "original traced worker uses its untouched original limits"}
        calls.append(evidence)
        save(output / "wall-activation-policy.json", evidence)
        return original(*args, **{**kwargs, "anchor_movement_limits": inner})

    boundary.refine_generated_region = guarded_refine
    try:
        # This worker still owns the original 0.5/0.5/1.25 independent validator.
        traced.worker(index, arm, output)
    finally:
        boundary.refine_generated_region = original


def declare(output):
    output.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), Path(traced.__file__), traced.PROJECT / "scripts/trace_private_minimization.py",
        traced.PROJECT / "scripts/audit_1he8_joint_relaxation.py", traced.PROJECT / "scripts/mac_reconstruction_comparison.py",
        traced.PROJECT / "scripts/run_mac_modeller_comparison.py", traced.PROJECT / "scripts/audit_1he8_traced_mobility.py",
        traced.PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += list((traced.PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    jobs = []
    for index, arm in CASES:
        case = traced.CASES[index]
        root = traced.PROJECT / case["root"]
        files += [root / name for name in ("request.json", "source.pdb", case["input"])]
        case_id = f"{case['name']}-{arm}-wall-buffer005"
        jobs.append({"case_id": case_id, "case": str(root), "seed": 7,
            "command": [sys.executable, "-m", "scripts.run_1he8_buffered_walls", "worker", "--index", str(index),
                "--arm", arm, "--output", str(output / "execution" / f"{case_id}-seed-7")]})
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": len(jobs),
        "scope": "post-result private engineering control; no new NN samples, no crystal scoring or public qualification",
        "motivation": "three richer-context arms passed geometry but finite wall minima exceeded validation bounds",
        "wall_buffer_angstrom": BUFFER_ANGSTROM, "wall_activation_limits": asdict(buffered_limits(AnchorMovementLimits())),
        "validation_limits_unchanged": asdict(AnchorMovementLimits()), "all_generated_atoms_mobile": 1712,
        "same_starting_proposals": True, "main_full_phase_max_iterations": 1000, "restarts": 0,
        "negative_control": "original full-patchwork/direct retains independent generated Ramachandran issue",
        "no_old_outcome_reclassified": True, "donor_approval_or_publication_implied": False,
        "resource_guard": "serial CPU executor,1800s/job,RSS8GiB,host pressure,10GiB evidence cap"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--index", type=int, choices=(0, 1))
    parser.add_argument("--arm", choices=("direct", "sequence2", "spatial5-padding1"))
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        if (args.index, args.arm) not in CASES:
            parser.error("undeclared buffered trial")
        worker(args.index, args.arm, args.output.resolve())
