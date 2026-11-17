#!/usr/bin/env python3
"""Bounded shared CPU repair of saved raw returns; no neural sampling/solvent."""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_1he8_buffered_walls import buffered_limits
from scripts.run_1he8_traced_mobility import anchor_residues
from scripts.run_mac_modeller_comparison import execute

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
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"
STUDY = ROOT / "sampler-return-study-seed11-attempt01"
READBACK = ROOT / "sampler-return-study-seed11-readback-attempt01/raw-readback.json"


def release_anchors(assessment):
    if assessment.admission_failures:
        raise ValueError("observation/admission failure cannot be repaired by releasing anchors")
    return not assessment.geometry_passed


def worker(case, output, *, request_path=None, source_path=None, raw_path=None,
           expected_native_seed=11, native_backend="Protenix"):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("repair requires fresh guarded reservation")
    started = time.monotonic()
    request_path = request_path if request_path is not None else STUDY / "inputs" / case / "context-request.json"
    request = DiffusionRequest.from_json(request_path.read_text())
    source_path = source_path if source_path is not None else request_path.parent / request.normalized_pdb.path
    raw_path = raw_path if raw_path is not None else STUDY / "execution" / case / "sampler-return/raw-materialized.pdb"
    source, raw = source_path.read_text(), raw_path.read_text()
    if digest(source_path) != request.normalized_pdb.sha256 or request.seeds != (expected_native_seed,):
        raise ValueError("repair source/seed mismatch")
    assess = assessment_adapter(request, source)
    initial = assess(raw)
    release_anchors(initial)
    structure = _structure(raw, "raw.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    for name, text in (("raw.pdb", raw), ("source.pdb", source), ("request.json", request.to_json())):
        with (output / name).open("x") as stream:
            stream.write(text)
    save(output / "initial-assessment.json", asdict(initial))
    limits, rows, current = AnchorMovementLimits(), [], raw
    for stage in ("generated-only", "generated-plus-sequence2"):
        if stage == "generated-plus-sequence2" and not release_anchors(assess(current)):
            break
        selected = () if stage == "generated-only" else anchor_residues(request, source, raw, "sequence2")
        restraints = build_anchor_restraints(source, selected, observed_atoms=request.fixed_atoms,
            base_k_kj_mol_nm2=10000., bfactor_source="disabled") if selected else ()
        allowed = frozenset(restraint.atom for restraint in restraints)
        stage_output = output / stage
        stage_output.mkdir(exist_ok=False)
        save(stage_output / "protocol.json", {"generated_atoms_mobile": len(request.generated_atoms),
            "native_backend": native_backend, "native_seed": expected_native_seed,
            "gap_count": len(request.gaps), "anchor_residues": [asdict(r) for r in selected],
            "restraints": [asdict(r) for r in restraints], "original_limits": asdict(limits),
            "inner_wall_limits": asdict(buffered_limits(limits)), "seed": 11, "restarts": 0,
            "max_iterations": 1000, "platform": "CPU", "solvent": False,
            "source_reference_never_ratcheted": True, "all_generated_jointly_mobile": True})
        trials = []

        def score(points, start_text=current, allowed_atoms=allowed, anchor_refs=restraints, recorded=trials):
            text = rewrite_generated_coordinates(start_text, points)
            assessment = assess(text, allowed_atoms)
            movement = check_anchor_movements(anchor_refs, _coordinates(text)[0], limits)
            recorded.append({"assessment": asdict(assessment), "movements": [asdict(m) for m in movement]})
            return len(assessment.admission_failures), sum(not m.within_limits for m in movement), sum(i.actionable for i in assessment.issues)

        try:
            result = refine_generated_region(current,
                generated_residues=tuple(r for gap in request.gaps for r in gap.generated_residues),
                generated_atoms=request.generated_atoms, max_iterations=1000, restart_count=0,
                random_seed=11, platform_name="CPU", geometry_preconditioning=True,
                candidate_score=score, observed_anchor_restraints=restraints,
                anchor_movement_limits=buffered_limits(limits))
            final = rewrite_generated_coordinates(current, result.coordinates_angstrom)
            coordinates, _ = _coordinates(final)
            global_structure = _structure(final, "refined.pdb")
            assert_all_l(global_structure.pdb.topology, global_structure.pdb.positions)
            assessment = assess(final, allowed)
            movement = check_anchor_movements(restraints, coordinates, limits)
            passed = assessment.geometry_passed and not assessment.admission_failures and all(m.within_limits for m in movement)
            with (stage_output / "refined.pdb").open("x") as stream:
                stream.write(final)
            row = {"stage": stage, "passed": passed, "assessment": asdict(assessment),
                "movements": [asdict(m) for m in movement], "global_all_l": True,
                "final_sha256": digest(stage_output / "refined.pdb"),
                "generated_atoms_jointly_mobile": len(request.generated_atoms),
                "observed_atoms_released": len(allowed), "search_diagnostics": result.search_diagnostics}
            save(stage_output / "result.json", row)
            rows.append(row)
            if stage == "generated-only":
                current = final
            if passed:
                break
        finally:
            save(stage_output / "trials.json", trials)
    summary = {"case": case, "stages": rows, "accepted": bool(rows and rows[-1]["passed"]),
        "native_backend": native_backend, "native_seed": expected_native_seed, "optimizer_seed": 11,
        "wall_seconds": time.monotonic()-started, "new_sampling_calls": 0,
        "raw_native_failures_preserved": True, "public_or_donor_approval": False}
    save(output / "result.json", summary)
    print(json.dumps({"case": case, "accepted": summary["accepted"], "stages": len(rows)}), flush=True)


def declare(output):
    readback = json.loads(READBACK.read_text())
    rows = {row["entry_id"]: row for row in readback["rows"]}
    if set(rows) != {"1NON", "2BE1"} or any(not row["raw_endpoint_available"] or not row["global_all_l"] for row in rows.values()):
        raise ValueError("repair requires independently audited, globally L raw endpoints")
    files = [Path(__file__), READBACK, PROJECT / "scripts/run_mac_modeller_comparison.py",
             PROJECT / "scripts/run_1he8_traced_mobility.py", PROJECT / "scripts/run_1he8_buffered_walls.py",
             PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += [Path(name) for name in readback["bindings"]]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    output.mkdir(parents=True, exist_ok=False)
    jobs = [{"case_id": case, "case": str(STUDY / "inputs" / case), "seed": 11,
             "command": [sys.executable, "-m", "scripts.refine_retained_sampler_returns", "worker", "--case", case,
                         "--output", str(output / "execution" / f"{case}-seed-11")]} for case in ("1NON", "2BE1")]
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": 2,
        "scope": "two saved raw endpoints; first generated-only, sequence2 release only if geometry still fails",
        "optimizer_seed": 11, "max_iterations_per_phase": 1000, "restarts": 0,
        "all_generated_atoms_jointly_mobile": True, "outer_limits_unchanged": asdict(AnchorMovementLimits()),
        "wall_buffer_angstrom": .05, "new_neural_calls": 0, "solvent": False,
        "resource_guard": "existing serial CPU1800s/job,RSS8GiB,host pressure,10GiB evidence",
        "no_old_outcome_reclassified": True, "no_publication_or_donor_approval": True})


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
        worker(args.case, args.output.resolve())
