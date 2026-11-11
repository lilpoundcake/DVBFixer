#!/usr/bin/env python3
"""Private all-generated + direct observed anchor refinement of retained patchwork.

No sampling, no solvent, no source repair, no selected-error mobility mask.
Every generated atom participates in ONE shared full-complex refinement problem.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path


def main():
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.anchor_policy import (
        AnchorMovementLimits,
        build_anchor_restraints,
        check_anchor_movements,
    )
    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.gap_portfolio import _coordinates
    from dvbfixer.model.diffusion.portfolio_evidence import json_bytes
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--patchwork", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--max-iterations", type=int, default=500)
    parser.add_argument("--restart-count", type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()

    def write(name, value):
        data = value if isinstance(value, bytes) else json_bytes(value)
        with (args.output / name).open("xb") as stream:
            stream.write(data)

    request_bytes, source, proposal = args.request.read_bytes(), args.source.read_bytes(), args.patchwork.read_bytes()
    write("request.json", request_bytes)
    write("source.pdb", source)
    write("patchwork-input.pdb", proposal)
    write("run-script.py", Path(__file__).read_bytes())
    request = DiffusionRequest.from_json(request_bytes.decode())
    if hashlib.sha256(source).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("patchwork source digest differs")
    source_text, proposal_text = source.decode(), proposal.decode()
    anchors = tuple(sorted({r for g in request.gaps for r in (g.left_anchor, g.right_anchor) if r is not None}))
    restraints = build_anchor_restraints(source_text, anchors, observed_atoms=request.fixed_atoms, bfactor_source="auto")
    allowed = frozenset(r.atom for r in restraints)
    generated = set(request.generated_atoms)
    mobile = generated | set(allowed)
    xyz, rows = _coordinates(proposal_text)
    if not mobile <= xyz.keys():
        raise ValueError("patchwork omits a generated/anchor atom")
    assess = assessment_adapter(request, source_text)
    initial = assess(proposal_text)
    if initial.admission_failures:
        raise ValueError("patchwork changed observed source records")
    limits = AnchorMovementLimits()
    write("protocol.json", {"scope": __doc__, "new_sampling": False, "solvent": False,
        "simultaneously_mobile_generated_atoms": len(generated), "gap_count": len(request.gaps),
        "anchors": [asdict(r) for r in anchors], "restraints": [asdict(r) for r in restraints],
        "generated_atoms": [asdict(a) for a in request.generated_atoms], "limits": asdict(limits),
        "seed": args.seed, "platform": "CPU", "max_iterations": args.max_iterations,
        "restart_count": args.restart_count, "request_sha256": hashlib.sha256(request_bytes).hexdigest(),
        "source_sha256": hashlib.sha256(source).hexdigest(), "patchwork_sha256": hashlib.sha256(proposal).hexdigest()})
    write("initial-assessment.json", asdict(initial))
    trials = []

    def score(coordinates):
        candidate = rewrite_generated_coordinates(proposal_text, coordinates)
        assessment = assess(candidate, allowed)
        movements = check_anchor_movements(restraints, _coordinates(candidate)[0], limits)
        write(f"trial-{len(trials):03d}.pdb", candidate.encode())
        trials.append({"assessment": asdict(assessment), "movements": [asdict(m) for m in movements]})
        return len(assessment.admission_failures), sum(not m.within_limits for m in movements), sum(i.actionable for i in assessment.issues)

    try:
        result = refine_generated_region(proposal_text,
            generated_residues=tuple(r for g in request.gaps for r in g.generated_residues),
            generated_atoms=request.generated_atoms, max_iterations=args.max_iterations,
            restart_count=args.restart_count, random_seed=args.seed, perturbation_angstrom=0.25,
            platform_name="CPU", geometry_preconditioning=True, candidate_score=score,
            observed_anchor_restraints=restraints, anchor_movement_limits=limits)
        final = rewrite_generated_coordinates(proposal_text, result.coordinates_angstrom)
        final_xyz, final_rows = _coordinates(final)
        if set(final_rows) != set(rows) or final != rewrite_generated_coordinates(
                proposal_text, {a: final_xyz[a] for a in mobile}):
            raise ValueError("patchwork refinement changed records outside generated atoms and direct anchors")
        structure = _structure(final, "patchwork-refined.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        assessment = assess(final, allowed)
        movements = check_anchor_movements(restraints, final_xyz, limits)
        passed = assessment.geometry_passed and not assessment.admission_failures and all(m.within_limits for m in movements)
        write("patchwork-refined.pdb", final.encode())
        write("final-assessment.json", asdict(assessment))
        write("movements.json", [asdict(m) for m in movements])
        summary = {"passed": passed, "new_sampling": False, "gap_count": len(request.gaps),
            "simultaneously_mobile_generated_atoms": len(generated), "observed_anchor_residues": len(anchors),
            "observed_anchor_atoms": len(allowed), "outside_mask_exact": True, "global_all_l": True,
            "observed_movement_violations": sum(not m.within_limits for m in movements),
            "actionable_categories": dict(Counter(i.category for i in assessment.issues if i.actionable)),
            "initial_actionable_categories": dict(Counter(i.category for i in initial.issues if i.actionable)),
            "final_sha256": hashlib.sha256(final.encode()).hexdigest(),
            "wall_seconds": time.monotonic() - started, "scientific_qualification": False,
            "published": False, "search_diagnostics": result.search_diagnostics}
        write("summary.json", summary)
        print(json.dumps(summary), flush=True)
    except BaseException as exc:
        write("failure.json", {"error_type": type(exc).__name__, "error": str(exc),
            "wall_seconds": time.monotonic() - started, "new_sampling": False})
        raise
    finally:
        write("trials.json", trials)


if __name__ == "__main__":
    main()
