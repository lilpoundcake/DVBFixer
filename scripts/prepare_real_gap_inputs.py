#!/usr/bin/env python3
"""Prepare eligible canonical source projections; backend admission is separate."""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

PROJECT = Path(__file__).resolve().parents[1]


def worker(entry_id, selected, output):
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.ffutils import sanitize_protein_hetatm
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.input_geometry import assess_input_geometry
    from dvbfixer.model.diffusion.input_preparation import prepare_input
    from dvbfixer.model.diffusion.scope import assess_diffusion_scope
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request
    from dvbfixer.model.pipeline import _strip_hetatm_lines

    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("source preparer requires a fresh executor-owned reservation")
    row = {"entry_id": entry_id, "status": "source_preparation_refused", "backend_admission_pending": True,
           "checkpoint_sampling_started": False, "public_capability_enabled": False}
    stage = "explicit-protein-projection"
    try:
        text = selected.read_text()
        full = output / "full-assembly-context.pdb"
        full.write_text(text)
        sanitized = Path(sanitize_protein_hetatm(str(full)))
        sanitized_text = sanitized.read_text()
        omitted = Counter(line[17:20].strip() for line in sanitized_text.splitlines() if line.startswith("HETATM"))
        model_input = output / "protein-projection.pdb"
        model_input.write_text("".join(_strip_hetatm_lines(sanitized_text.splitlines(keepends=True), keep_water=False)))
        save(output / "projection.json", {"original_selected_sha256": digest(selected), "full_context_sha256": digest(full),
            "protein_projection_sha256": digest(model_input), "omitted_heterogen_atom_counts": dict(omitted),
            "not_whole_entry_chemistry": True, "full_context_retained": True,
            "final_whole_context_obstacle_audit_still_required": True})
        stage = "shared-sequence-placement-and-request"
        protein_chains = {line[21] for line in model_input.read_text().splitlines() if line.startswith("ATOM  ")}
        seqres_chains = {line[11] for line in model_input.read_text().splitlines() if line.startswith("SEQRES") and len(line) > 11}
        if not protein_chains <= seqres_chains:
            raise ValueError("assembly lost deposited SEQRES for an observed protein chain; no observed-only fallback")
        targets = _target_sequences(model_input, None)
        if not targets:
            raise ValueError("assembly has no deposited target sequences; no observed-only fallback")
        (output / "target.fasta").write_text("".join(f">{chain}\n{sequence}\n" for chain, sequence in targets.items()))
        request = build_cli_diffusion_request(model_input, targets, seeds=(7,),
            profile="research-natural-gap-source-v1", no_terminal=False)
        (output / "request-before-preparation.json").write_text(request.to_json())
        stage = "shared-observed-input-preparation"
        prepared = prepare_input(request, model_input.read_bytes(), output)
        if prepared is not None:
            request, model_input = prepared
        stage = "independent-source-geometry"
        findings = assess_input_geometry(request, model_input.read_bytes())
        save(output / "source-findings.json", [asdict(f) for f in findings])
        errors = [f for f in findings if f.severity is Severity.ERROR]
        if errors:
            raise ValueError(f"source projection retains {len(errors)} severe geometry findings")
        structure = _structure(model_input.read_text(), "source-prepared.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        # This is explicitly private input readiness, NOT a public profile bypass.
        public_scope = assess_diffusion_scope(request, model_input.read_bytes())
        save(output / "existing-public-scope.json", asdict(public_scope))
        (output / "request-source-prepared.json").write_text(request.to_json())
        row.update(status="private_source_projection_prepared", source_input=str(model_input),
            source_input_sha256=digest(model_input), prepared_source_derivative=prepared is not None,
            request_path=str(output / "request-source-prepared.json"), request_sha256=digest(output / "request-source-prepared.json"),
            target_residues=sum(len(t.sequence) for t in request.target_sequences),
            gap_count=len(request.gaps), gap_lengths=[g.target_interval.stop-g.target_interval.start for g in request.gaps],
            fixed_atoms=len(request.fixed_atoms), generated_atoms=len(request.generated_atoms), global_all_l=True,
            severe_source_geometry_findings=0, existing_public_scope=asdict(public_scope),
            original_context_and_source_construct_review_pending=True)
    except (ImportError, FileNotFoundError):
        (output / "environment-failure.txt").write_text(traceback.format_exc())
        raise
    except (ValueError, RuntimeError) as exc:
        row.update(stage=stage, error_type=type(exc).__name__, error=str(exc))
        (output / "source-preparation-failure.txt").write_text(traceback.format_exc())
    except Exception:
        (output / "technical-failure.txt").write_text(traceback.format_exc())
        raise
    save(output / "result.json", row)


def declare(assemblies, output):
    complete = json.loads((assemblies / "completed.json").read_text())
    plan = json.loads((assemblies / "plan.json").read_text())
    if complete["completed"] != len(plan["jobs"]):
        raise ValueError("assembly inventory incomplete")
    selected = []
    for result_path in sorted((assemblies / "execution").glob("*/result.json")):
        for row in json.loads(result_path.read_text())["rows"]:
            if row["status"] == "assembly_inventory_complete" and not row["strict_protein_route_obstacles"]:
                path = result_path.parent / "cases" / row["entry_id"] / "selected.pdb"
                if digest(path) != row["selected_sha256"]:
                    raise ValueError("selected source assembly changed")
                selected.append((row["entry_id"], path))
    output.mkdir(parents=True, exist_ok=False)
    jobs = [{"case_id": entry_id, "case": str(path), "seed": 0,
        "command": [sys.executable, "-m", "scripts.prepare_real_gap_inputs", "worker", "--entry", entry_id,
            "--selected", str(path), "--output", str(output / "execution" / f"{entry_id}-seed-0")]}
        for entry_id, path in selected]
    files = [Path(__file__), PROJECT / "src/dvbfixer/model/diffusion_cli.py", PROJECT / "src/dvbfixer/model/pipeline.py"]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    files += [path for _, path in selected]
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": len(jobs),
        "scope": "canonical no-inventory-obstacle subset only, not completed200-entry benchmark",
        "private_profile": "research-natural-gap-source-v1, never registered as public capability",
        "observed_heavy_completion_and_bounded_source_repair": "existing prepare_input; evidence/new source ownership retained",
        "seed_0": "executor reservation only", "request_seed_7": "placeholder not consumed by a neural sampler",
        "ML_calls": 0, "source_cases": len(jobs), "resource_guard": "serial CPU,1800s/case,RSS8GiB,10GiB evidence cap",
        "heterogen_projection": "explicit, full context retained; not whole-entry conditioning/FF qualification",
        "remaining_assembly_context_reviews_not_dropped_from_cohort_ledger": True})
    print(json.dumps({"source_preparation_cases_declared": len(jobs), "ML_calls": 0}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--assemblies", type=Path)
    parser.add_argument("--selected", type=Path)
    parser.add_argument("--entry")
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.assemblies.resolve(), args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        worker(args.entry, args.selected.resolve(), args.output.resolve())
