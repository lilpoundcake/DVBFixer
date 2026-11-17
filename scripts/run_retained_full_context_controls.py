#!/usr/bin/env python3
"""Bounded, serial full-context controls on explicitly retained local inputs."""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import traceback
from dataclasses import replace
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import freeze_code
from scripts.retained_modeller_full_context import alignment, merge_projection, read_atoms, run
from scripts.run_mac_modeller_comparison import execute
from scripts.run_retained_modeller_fixation_controls import ROOT

from dvbfixer.model.diffusion.contract import DiffusionRequest

PROJECT = Path(__file__).resolve().parents[1]
SEED = 103
ENTRIES = ("1CC6", "9M3K")
ALLOWED_ENTRIES = (*ENTRIES, "1X0X", "5HJI", "3C7D", "6ZH0")


def declare(output, entries=ENTRIES, seed=SEED):
    if not 1 <= len(entries) <= 2 or len(set(entries)) != len(entries) or any(entry not in ALLOWED_ENTRIES for entry in entries):
        raise ValueError("declare one or two distinct retained cases")
    if seed not in (103, 104):
        raise ValueError("undeclared development seed")
    output.mkdir(parents=True, exist_ok=False)
    files, jobs = [], []
    for entry in entries:
        review = ROOT / f"retained-source-review-matrix-attempt01/execution/retained-source-review-seed-0/report/entries/{entry}.json"
        row = json.loads(review.read_text())
        original_request, original_source = Path(row["request_path"]), Path(row["source_path"])
        if digest(original_source) != row["source_sha256"] or digest(original_request) != row["request_sha256"]:
            raise ValueError("retained source/request changed")
        request = DiffusionRequest.from_json(original_request.read_text())
        parent = original_request.parent
        context_path = parent / "full-assembly-context.pdb"
        projection_path = parent / "projection.json"
        projection = json.loads(projection_path.read_text())
        if digest(context_path) != projection["full_context_sha256"]:
            raise ValueError("original full context changed")
        source = original_source.read_text()
        context = merge_projection(source, context_path.read_text())
        if (set(read_atoms(source)) != set(request.fixed_atoms) or request.retained_explicit_links
                or request.heterogen_contexts or len(request.target_sequences) != 1):
            raise ValueError("outside the declared unlinked single-protein task")
        case = output / "cases" / entry
        case.mkdir(parents=True)
        shutil.copyfile(original_source, case / "source.pdb")
        shutil.copyfile(context_path, case / "original-full-context.pdb")
        (case / "prepared-full-context.pdb").write_text(context)
        (case / "request.json").write_text(replace(request, seeds=(seed,), protein_context=None).to_json())
        (case / "alignment.pir").write_text(alignment(request, context))
        annotation = ROOT / f"retained-construct-annotation-audit-attempt01/entries/{entry}.json"
        original_cif = json.loads(annotation.read_text())
        source_cif = Path(original_cif["original_source_path"])
        if digest(source_cif) != original_cif["original_source_sha256"]:
            raise ValueError("deposited source CIF changed")
        files += [review, original_request, original_source, context_path, projection_path, annotation, source_cif,
                  parent / "result.json", parent / "target.fasta"]
        files += [ROOT / name / f"entries/{entry}.json" for name in (
            "retained-context-connections-attempt02", "retained-construct-difference-audit-attempt02")]
        if row["prepared_source_derivative"]:
            files += [p for p in original_source.parent.rglob("*") if p.is_file()]
        save(case / "lineage.json", {"entry_id": entry, "title": row["title"],
            "original_request": str(original_request), "original_request_sha256": digest(original_request),
            "original_source": str(source_cif), "original_source_sha256": digest(source_cif),
            "assembly_id": "1", "all_retained_context_in_native_input": True,
            "prepared_source_derivative": row["prepared_source_derivative"],
            "fixed_reference": "exact retained prepared protein plus unchanged deposited context",
            "target_sequences_unchanged": True, "declared_gap_lengths": [g["length"] for g in row["gaps"]],
            "generated_atoms": len(request.generated_atoms), "fixed_protein_atoms": len(request.fixed_atoms),
            "fixed_full_context_atoms": len(read_atoms(context, allow_fractional_heterogens=True)),
            "context_atom_counts": row["omitted_heterogen_atom_counts"],
            "fractional_context_policy": "retain unique nonalternate HETATM records with original0<occupancy<=1 as fixed conservative context; never normalize"})
        files += list(case.iterdir())
        jobs.append({"case_id": entry, "case": str(case), "seed": seed,
            "command": [sys.executable, "-m", "scripts.run_retained_full_context_controls", "worker",
                        "--case", str(case), "--seed", str(seed), "--output", str(output / "execution" / f"{entry}-seed-{seed}")]})
    files += [PROJECT / "scripts" / name for name in ("run_retained_full_context_controls.py",
        "retained_modeller_full_context.py", "run_mac_modeller_comparison.py", "run_mac_reconstruction_guarded.py",
        "mac_reconstruction_comparison.py", "run_retained_6z38_restoration.py", "retained_evidence.py",
        "finish_retained_context_batch.py", "recover_retained_full_context_analysis.py",
        "audit_retained_full_context_controls.py", "export_retained_context_review.py")]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    files += list((PROJECT / "src/dvbfixer/diagnose").glob("*.py"))
    files += [PROJECT / "src/dvbfixer/model/modeller_run.py", PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(path): digest(path) for path in files}, "jobs": jobs, "count": len(jobs),
        "backend": "MODELLER 10.8 CPU; full native protein/ligand/cofactor/water context", "seed": seed,
        "model_make_budget": len(jobs), "initial_models_per_case": 1, "fast_loop_models_per_case": 1,
        "generated_only_CG_continuations_per_case": 1, "CG_max_iterations": 1000,
        "CG_movement_limit_from_raw_angstrom": 3., "observed_anchor_motion": False,
        "IUPAC_policy": "skip documented full-model to_iupac; instrument native method calls; no post-hoc atom swaps",
        "new_neural_calls": 0, "new_downloads": 0, "solvent_box": False, "automatic_retries": 0,
        "resource_guard": "existing serial1800s/job,RSS8GiB,host-pressure,10GiB evidence",
        "chemistry_scope": "unchanged deposited construct sequences; no protonation, H/caps/charges or external OXT additions",
        "fractional_context_policy": "original nonzero HETATM occupancy retained; selected protein requires occupancy1; alternate/zero-occupancy records refused",
        "full_diagnostics": "same chemistry/Ramachandran/chi12/sterics for source and results; include waters; .6/.9 A thresholds",
        "scientific_qualification_or_public_capability": False})
    freeze_code(output / "plan.json")


def worker(case, output, seed=SEED):
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    if request.seeds != (seed,):
        raise ValueError("undeclared seed")
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("worker requires a fresh guard reservation")
    try:
        run(request, (case / "source.pdb").read_text(), (case / "prepared-full-context.pdb").read_text(), output / "native", seed)
        result = json.loads((output / "native/stages.json").read_text())
        save(output / "result.json", {"entry_id": case.name, "seed": seed, "native": result,
                                      "public_output_published": False})
    except Exception:
        save(output / "failure.json", {"traceback": traceback.format_exc(), "automatic_retry": False})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--case", type=Path)
    parser.add_argument("--entries", nargs="+", choices=ALLOWED_ENTRIES, default=ENTRIES)
    parser.add_argument("--seed", type=int, choices=(103, 104), default=SEED)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve(), entries=args.entries, seed=args.seed)
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        if args.case is None:
            parser.error("worker requires --case")
        worker(args.case.resolve(), args.output.resolve(), seed=args.seed)
