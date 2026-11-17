#!/usr/bin/env python3
"""One private, source-bound 6Z38 MODELLER attempt; no public capability."""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute
from scripts.run_retained_loop_modeller_controls import multichain_alignment
from scripts.run_retained_modeller_fixation_controls import ROOT, SEED

from dvbfixer.model.diffusion.contract import DiffusionRequest


def observation_audit(request, source, native):
    """Measure actual native coordinates, before any fit or observation copying."""
    if set(source) != set(request.fixed_atoms):
        raise ValueError("source observation inventory differs from request")
    if set(native) != set(request.fixed_atoms) | set(request.generated_atoms):
        raise ValueError("native atom ownership differs from request")
    changed = []
    for atom in request.fixed_atoms:
        if source[atom][:2] != native[atom][:2]:
            raise ValueError("native observed atom chemistry changed")
        distance = float(np.linalg.norm(source[atom][2] - native[atom][2]))
        if distance != 0:
            changed.append({"atom": asdict(atom), "distance_angstrom": distance})
    return {"fixed_atoms": len(source), "changed_atoms": changed,
            "max_motion_angstrom": max((r["distance_angstrom"] for r in changed), default=0.0),
            "native_observations_exact_at_PDB_precision": not changed,
            "global_fit_or_source_coordinate_copy_applied": False}


def signature(finding):
    return json.dumps(asdict(finding), sort_keys=True, ensure_ascii=False)


def compare_findings(baseline, candidate):
    remaining = Counter(signature(f) for f in baseline)
    unchanged, new = [], []
    for finding in candidate:
        key = signature(finding)
        if remaining[key]:
            remaining[key] -= 1
            unchanged.append(asdict(finding))
        else:
            new.append(asdict(finding))
    return {"unchanged_source_findings": unchanged, "new_or_changed_findings": new}


def full_findings(path):
    from openmm.app import PDBFile

    from dvbfixer.diagnose.chemistry import run_all
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.diagnose.steric import clashes_python
    from dvbfixer.model.diffusion.quality import check_ramachandran, check_sidechain_chi12

    pdb = PDBFile(str(path))
    findings = run_all(pdb.topology, pdb.positions)
    findings += check_ramachandran(pdb.topology, pdb.positions)
    findings += check_sidechain_chi12(pdb.topology, pdb.positions)
    findings += clashes_python(pdb.topology, pdb.positions, include_water=False,
                              clash_warn_a=.6, clash_error_a=.9)
    return [f for f in findings if f.severity != Severity.INFO]


def declare(output, *, seed=SEED, previous=None):
    review = ROOT / "retained-source-review-matrix-attempt01/execution/retained-source-review-seed-0/report/entries/6Z38.json"
    row = json.loads(review.read_text())
    source, original_request = Path(row["source_path"]), Path(row["request_path"])
    if digest(source) != row["source_sha256"] or digest(original_request) != row["request_sha256"]:
        raise ValueError("retained source/request changed")
    request = DiffusionRequest.from_json(original_request.read_text())
    if (len(request.target_sequences) != 1 or request.target_sequences[0].chain != "A"
            or len(request.target_sequences[0].sequence) != 432
            or [g.target_interval.stop - g.target_interval.start for g in request.gaps] != [3, 9, 6]
            or not request.target_sequences[0].sequence.endswith("HHHHHH")
            or request.retained_explicit_links or request.heterogen_contexts
            or row["omitted_heterogen_atom_counts"]):
        raise ValueError("outside the declared complete-context 6Z38 task")
    parent = original_request.parent
    evidence = [review, source, original_request, parent / "full-assembly-context.pdb",
                parent / "protein-projection.pdb", parent / "target.fasta", parent / "result.json"]
    evidence += sorted((parent / "dvbfixer-input-completion-4phbtyou").rglob("*"))
    evidence += [ROOT / name / "entries/6Z38.json" for name in (
        "retained-context-connections-attempt02", "retained-construct-annotation-audit-attempt01",
        "retained-construct-difference-audit-attempt02")]
    files = [p for p in evidence if p.is_file()]
    project = Path(__file__).resolve().parents[1]
    files += [project / "scripts" / name for name in (
        "run_retained_6z38_restoration.py", "run_mac_modeller_proposal.py",
        "run_retained_modeller_fixation_controls.py", "run_retained_loop_modeller_controls.py",
        "mac_reconstruction_comparison.py", "run_mac_modeller_comparison.py", "run_mac_reconstruction_guarded.py")]
    files += list((project / "src/dvbfixer/model/diffusion").glob("*.py"))
    files += list((project / "src/dvbfixer/diagnose").glob("*.py"))
    files += [project / "src/dvbfixer/model/modeller_run.py", project / "src/dvbfixer/ffutils/geometry.py"]
    if seed not in (101, 102):
        raise ValueError("undeclared development seed")
    if seed == 102:
        if previous is None:
            raise ValueError("seed102 amendment requires retained seed101 failure")
        previous_plan = json.loads((previous / "plan.json").read_text())
        receipt_path = previous / "execution/6Z38-seed-101/execution.json"
        receipt = json.loads(receipt_path.read_text())
        if previous_plan["seed"] != 101 or receipt["completed"] or receipt["returncode"] != 1:
            raise ValueError("unexpected previous execution state")
        files += [previous / "plan.json", receipt_path,
                  previous / "execution/6Z38-seed-101/output.log"]
        files += list((previous / "implementation-before-terminal-selection-fix").iterdir())
    output.mkdir(parents=True, exist_ok=False)
    case = output / "cases/6Z38"
    case.mkdir(parents=True)
    shutil.copyfile(source, case / "source.pdb")
    (case / "request.json").write_text(replace(request, seeds=(seed,), protein_context=None).to_json())
    (case / "alignment.pir").write_text(multichain_alignment(request))
    files += list(case.iterdir())
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "count": 1,
        "jobs": [{"case_id": "6Z38", "case": str(case), "seed": seed,
            "command": [sys.executable, "-m", "scripts.run_retained_6z38_restoration", "worker",
                        "--case", str(case), "--seed", str(seed), "--output", str(output / f"execution/6Z38-seed-{seed}")]}],
        "backend": "MODELLER CPU; documented generated-only initial selection", "seed": seed,
        "previous_attempt": None if previous is None else str(previous),
        "cumulative_model_make_budget": 1 if previous is None else 2,
        "engine_terminal_OXT": "immobile, outside request mask; excluded by existing audited normalizer; never published",
        "actual_model_make_budget": 1, "initial_models": 1, "fast_loop_derivatives": 1,
        "neural_calls": 0, "repair_calls": 0, "automatic_retries": 0, "solvent": False,
        "resource_guard": "existing serial CPU1800s/job,RSS8GiB,host-pressure,10GiB evidence",
        "construct": "deposited TodX deltaS2S3 including six C-terminal expression-tag HIS; no database-deletion restoration",
        "fixed_reference": "retained prepared source, including two previously completed heavy atoms; original bytes/evidence separately bound",
        "native_fixation_gate": "exact coordinates at PDB precision before fit/copy; no observed-anchor motion authorized",
        "full_diagnostics": "chemistry including sidechain stereocentres and S-S, Ramachandran, chi12, sterics .6/.9 A; same source/candidate settings",
        "terminal_chemistry": "no additions outside frozen mask; no H, caps, charges or protonation",
        "new_downloads": 0, "public_capability_or_qualification": False})


def worker(case, output, *, seed=SEED):
    from scripts.audit_mac_reconstruction_smoke import inventory
    from scripts.run_retained_modeller_fixation_controls import worker as native_worker

    request = DiffusionRequest.from_json((case / "request.json").read_text())
    native_worker(case, output, seed=seed)
    native_result = json.loads((output / "result.json").read_text())
    custody = observation_audit(request, inventory(case / "source.pdb"),
                                inventory(output / "proposal/mapped-native.pdb"))
    baseline = full_findings(case / "source.pdb")
    stages = {name: compare_findings(baseline, full_findings(path)) for name, path in (
        ("native", output / "proposal/mapped-native.pdb"),
        ("source_restored_diagnostic_only", output / "source-restored.pdb"))}
    save(output / "strict-readback.json", {"native_observation_audit": custody,
        "source_findings": [asdict(f) for f in baseline], "stages": stages,
        "restoration_accepted": custody["native_observations_exact_at_PDB_precision"]
            and native_result["source_restored_geometry_and_L_passed"]
            and native_result["native_global_all_l"]
            and not stages["source_restored_diagnostic_only"]["new_or_changed_findings"],
        "source_copy_does_not_repair_native_fixation_failure": True,
        "public_output_published": False, "no_refinement_or_resampling": True,
        "bindings": {str(p): digest(p) for p in (case / "source.pdb", case / "request.json",
            output / "proposal/native.pdb", output / "proposal/mapped-native.pdb", output / "source-restored.pdb")}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--case", type=Path)
    parser.add_argument("--seed", type=int, choices=(101, 102), default=SEED)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve(), seed=args.seed, previous=args.previous.resolve() if args.previous else None)
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        if args.case is None:
            parser.error("worker requires --case")
        worker(args.case.resolve(), args.output.resolve(), seed=args.seed)
