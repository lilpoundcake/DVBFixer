#!/usr/bin/env python3
"""Bounded necessity-gated reference-blind anchor experiment, not policy approval."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion.contract import DiffusionRequest, ResidueIdentity

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-reconstruction-comparison-2026-11-11"
ENVELOPES = ((.10, .10, .50), (.25, .25, .75), (.50, .50, 1.25))  # BB, CA, maximum heavy
STRENGTHS = (10000., 100000.)


def affected_anchors(request, assessment):
    selected = set()
    for issue in assessment["issues"]:
        if not issue["actionable"] or issue["category"] not in {"non_planar_amide", "junction-peptide-connectivity", "bond_angle", "bond_length"}:
            continue
        identities = {ResidueIdentity(**r) for r in issue["residues"]}
        for gap in request.gaps:
            if identities & set(gap.generated_residues):
                selected.update(r for r in (gap.left_anchor, gap.right_anchor) if r is not None and r in identities)
    return tuple(sorted(selected))


def declare(output):
    baseline = ROOT / "modeller-attempt02"
    rows = json.loads((baseline / "progress.json").read_text())
    output.mkdir(parents=True, exist_ok=False)
    jobs, excluded = [], []
    python = Path(json.loads((baseline / "plan.json").read_text())["python"])
    for row in rows:
        result, case = row["result"], Path(row["case"])
        request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
        anchors = affected_anchors(request, result["final_assessment"])
        if result["final_source_geometry_and_l_passed"] or result["final_admission_failures"] or not anchors:
            excluded.append({"case": row["case_id"], "seed": row["seed"], "eligible": False,
                "reason": "already passes, admission failure, or no residual junction defect implicating an observed direct flank"})
            continue
        for window in ("direct", "expanded2"):
            for strength in STRENGTHS:
                for index, limits in enumerate(ENVELOPES):
                    case_id = f"{row['case_id']}-{window}-k{int(strength)}-envelope{index}"
                    work = output / "execution" / f"{case_id}-seed-{row['seed']}"
                    command = [str(python), "-m", "scripts.run_mac_anchor_calibration", "worker", "--output", str(work),
                        "--case", str(case), "--input", result["final_path"], "--proposal-seed", str(row["seed"]),
                        "--window", window, "--strength", str(strength), "--envelope", str(index)]
                    jobs.append({"case": str(case), "case_id": case_id, "seed": row["seed"], "command": command,
                        "input": result["final_path"], "input_sha256": result["final_sha256"]})
    if len(jobs) > 36:
        raise ValueError("declared calibration slice exceeds 36 CPU jobs; do not expand automatically")
    files = [PROJECT / "scripts" / name for name in ("run_mac_anchor_calibration.py", "run_mac_modeller_comparison.py",
        "mac_reconstruction_comparison.py", "audit_mac_reconstruction_smoke.py", "run_mac_reconstruction_guarded.py")]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    files += [Path(job["input"]) for job in jobs]
    files += [Path(job["case"]) / name for job in jobs for name in ("source.pdb", "evaluation-mask.json", "reference-projection-NOT-CONDITIONING.pdb")]
    save(output / "plan.json", {"protocol": "private-necessary-anchor-grid-v1", "python": str(python), "python_sha256": digest(python),
        "jobs": jobs, "count": len(jobs), "source_code_and_input_sha256": {str(path): digest(path) for path in files},
        "excluded_zero_motion_controls": excluded, "envelopes_BB_CA_maxheavy_A": ENVELOPES, "strengths_kj_mol_nm2": STRENGTHS,
        "windows": ["affected direct observed flank only", "affected flank plus at most two further observed residues on that side"],
        "selection": "all eligible residual observed-junction cases from 30 fixed-repair baseline outcomes, no hidden-RMSD selection",
        "reference_in_refinement": False, "bfactor_weighting": False, "optimizer_seed": 7, "max_iterations": 1000,
        "geometry_preconditioning": False, "restarts": 0, "all_generated_atoms_mobile": True,
        "resource_policy": "serial CPU pressure-v2, 1800s/job, 8GiB process RSS, shared study lock",
        "no_validation_seeds_or_new_neural_samples": True, "production_limits_changed": False, "qualification": False})
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps({"jobs": len(jobs), "eligible_proposals": len(jobs)//12, "plan": str(output / "plan.json")}, indent=2))


def worker(case, input_path, output, proposal_seed, window, strength, envelope):
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.anchor_policy import (
        AnchorMovementLimits,
        build_anchor_restraints,
        check_anchor_movements,
    )
    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
    from dvbfixer.model.diffusion.gap_portfolio import _coordinates
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    source, text = (case / "source.pdb").read_text(), input_path.read_text()
    if digest(case / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("source changed")
    assess = assessment_adapter(request, source)
    before = assess(text)
    anchors = set(affected_anchors(request, asdict(before)))
    if before.geometry_passed or before.admission_failures or not anchors:
        raise ValueError("anchor mobility requires an unresolved observed-junction defect")
    if window == "expanded2":
        mapping = target_axis(request)
        observed = {r for p in request.sequence_placements for r in p.observed_residues}
        for anchor in tuple(anchors):
            ordered = mapping[anchor.chain]
            index = ordered.index(anchor)
            anchors.update(r for r in ordered[max(0, index-2):index+3] if r in observed)
    restraints = build_anchor_restraints(source, tuple(sorted(anchors)), observed_atoms=request.fixed_atoms,
        base_k_kj_mol_nm2=strength, bfactor_source="disabled")
    allowed = frozenset(r.atom for r in restraints)
    limits = AnchorMovementLimits(backbone_rmsd_angstrom_max=ENVELOPES[envelope][0],
        ca_displacement_angstrom_max=ENVELOPES[envelope][1], heavy_displacement_angstrom_max=ENVELOPES[envelope][2])
    save(output / "declaration.json", {"source_sha256": digest(case / "source.pdb"), "input_sha256": digest(input_path),
        "proposal_seed": proposal_seed, "optimizer_seed": 7, "window": window, "anchors": [asdict(a) for a in sorted(anchors)],
        "limits": asdict(limits), "strength": strength, "all_generated_mobile": len(request.generated_atoms), "reference_to_optimizer": False})
    (output / "input.pdb").write_text(text)
    def score(coordinates):
        candidate = rewrite_generated_coordinates(text, coordinates)
        assessed = assess(candidate, allowed)
        movements = check_anchor_movements(restraints, _coordinates(candidate)[0], limits)
        return len(assessed.admission_failures), sum(not m.within_limits for m in movements), sum(i.actionable for i in assessed.issues)
    result = refine_generated_region(text, generated_residues=tuple(r for gap in request.gaps for r in gap.generated_residues),
        generated_atoms=request.generated_atoms, platform_name="CPU", random_seed=7, max_iterations=1000, restart_count=0,
        geometry_preconditioning=False, observed_anchor_restraints=restraints, anchor_movement_limits=limits, candidate_score=score)
    if set(result.coordinates_angstrom) != set(request.generated_atoms) | set(allowed):
        raise ValueError("calibration movable mask changed")
    final = rewrite_generated_coordinates(text, result.coordinates_angstrom)
    (output / "refined.pdb").write_text(final)
    assessed = assess(final, allowed)
    positions, _ = _coordinates(final)
    movements = check_anchor_movements(restraints, positions, limits)
    structure = _structure(final, "anchor-grid-final.pdb")
    global_l = True
    try:
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
    except Exception:
        global_l = False
    # First reference-coordinate use: independent post-refinement scoring.
    reference, final_inventory = inventory(case / "reference-projection-NOT-CONDITIONING.pdb"), inventory(output / "refined.pdb")
    rmsd = float(np.sqrt(np.mean([np.sum((reference[a][2]-final_inventory[a][2])**2) for a in request.generated_atoms])))
    save(output / "result.json", {"case": case.name, "proposal_seed": proposal_seed, "window": window, "strength": strength,
        "envelope": envelope, "movements": [asdict(m) for m in movements], "global_all_l": global_l,
        "geometry_passed": assessed.geometry_passed, "admission_failures": list(assessed.admission_failures),
        "movement_limits_passed": all(m.within_limits for m in movements), "assessment": asdict(assessed),
        "heldout_heavy_rmsd": rmsd, "final_path": str(output / "refined.pdb"), "final_sha256": digest(output / "refined.pdb"),
        "calibrated_positive_production_limits": False, "qualification": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case", type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--proposal-seed", type=int)
    parser.add_argument("--window", choices=("direct", "expanded2"))
    parser.add_argument("--strength", type=float, choices=STRENGTHS)
    parser.add_argument("--envelope", type=int, choices=range(len(ENVELOPES)))
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        worker(args.case.resolve(), args.input.resolve(), args.output.resolve(), args.proposal_seed, args.window, args.strength, args.envelope)
