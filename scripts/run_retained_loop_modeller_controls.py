#!/usr/bin/env python3
"""Two source-matched MODELLER controls, no hidden reference or native fallback."""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import digest, modeller_alignment, save
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion.contract import DiffusionRequest

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"
ORIGINAL = ROOT / "retained-loop-protenix-encoder-attempt01"
SEED = 42


def multichain_alignment(request):
    text = modeller_alignment(request)
    first, last = request.target_sequences[0].chain, request.target_sequences[-1].chain
    if any(len(chain) != 1 or chain in {" ", ":", "@"} for chain in (first, last)):
        raise ValueError("private comparator needs explicit supported endpoint chain IDs")
    old = "structureX:source:FIRST:@:LAST:@::::"
    if text.count(old) != 1:
        raise ValueError("unexpected original PIR header")
    return text.replace(old, f"structureX:source:FIRST:{first}:LAST:{last}::::", 1)


def declare(output, previous=None):
    preflight_path = ROOT / "modeller-multichain-alignment-preflight-attempt01/result.json"
    preflight = {row["entry_id"]: row for row in json.loads(preflight_path.read_text())["rows"]}
    output.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), preflight_path, PROJECT / "scripts/run_mac_modeller_proposal.py",
        PROJECT / "scripts/run_mac_modeller_comparison.py", PROJECT / "scripts/mac_reconstruction_comparison.py",
        PROJECT / "src/dvbfixer/model/modeller_run.py", PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    jobs = []
    for entry in ("1NON", "2BE1"):
        original = ORIGINAL / entry / "context-request.json"
        request = DiffusionRequest.from_json(original.read_text())
        # MODELLER uses the whole observed assembly, unlike the local NN view.
        view = replace(request, protein_context=None, seeds=(SEED,))
        case = output / "cases" / entry
        case.mkdir(parents=True, exist_ok=False)
        source = ORIGINAL / entry / request.normalized_pdb.path
        shutil.copyfile(source, case / "source.pdb")
        (case / "request.json").write_text(view.to_json())
        (case / "alignment.pir").write_text(multichain_alignment(view))
        if (not preflight[entry]["native_alignment_check_passed"]
                or preflight[entry]["PIR_sha256"] != digest(case / "alignment.pir")
                or preflight[entry]["source_sha256"] != digest(source)):
            raise ValueError("native PIR/source preflight differs from this declaration")
        save(case / "lineage.json", {"original_request_path": str(original), "original_request_sha256": digest(original),
            "method_view_changes": ["protein_context=None; whole-source template", f"MODELLER seed{SEED}", "explicit first/last source chain PIR header"],
            "fixed_generated_gaps_target_axes_unchanged": True, "not_context_matched_to_local_Protenix": True})
        files += [original, source] + list(case.iterdir())
        jobs.append({"case_id": entry, "case": str(case), "seed": SEED,
            "command": [sys.executable, "-m", "scripts.run_retained_loop_modeller_controls", "worker",
                "--case", str(case), "--output", str(output / "execution" / f"{entry}-seed-{SEED}")]})
    failure = None
    if previous is not None:
        path = previous / "alignment-failure-readback.json"
        failure = json.loads(path.read_text())
        if failure["native_coordinate_proposal_available"] or failure["second_case_launched"]:
            raise ValueError("unexpected previous comparator execution state")
        for name, expected in failure["bindings"].items():
            if digest(Path(name)) != expected:
                raise ValueError("previous alignment failure evidence changed")
        files += [path] + [Path(name) for name in failure["bindings"]]
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": 2,
        "actual_MODELLER_proposal_budget": 2, "seed": SEED, "one_initial_and_one_fast_loop_derivative": True,
        "previous_failure": failure, "old_seed7_not_replayed": True,
        "cumulative_model_make_attempts_including_previous_failure": 3 if failure else 2,
        "gap_selection": "all missing residues together, exact existing PIR; no initial-model fallback",
        "no_hidden_reference_or_crystal_RMSD": True, "new_neural_calls": 0, "repair_calls": 0,
        "native_drift_before_source_restoration": True, "native_D_not_reflected": True,
        "restoration_policy": "generated-only source copy",
        "comparison_scope": "same source/gap task; different native context/seed/method, not superiority evidence",
        "resource_guard": "existing serial CPU1800s/job,RSS8GiB,host-pressure,10GiB evidence",
        "new_sources_or_production_limits": False, "publication_or_donor_approval": False})


def worker(case, output):
    from scripts import run_mac_modeller_proposal as proposal
    from scripts.audit_mac_reconstruction_smoke import inventory

    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _assessment, _structure
    from dvbfixer.model.diffusion.geometry import weighted_kabsch
    from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("MODELLER control needs fresh executor reservation")
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    source_path = case / "source.pdb"
    if request.seeds != (SEED,) or digest(source_path) != request.normalized_pdb.sha256 or request.protein_context is not None:
        raise ValueError("MODELLER method-view source changed")
    old_alignment = proposal.modeller_alignment
    proposal.modeller_alignment = multichain_alignment
    try:
        proposal.run(case / "request.json", source_path, output / "proposal", SEED)
    finally:
        proposal.modeller_alignment = old_alignment
    source, native = inventory(source_path), inventory(output / "proposal/mapped-native.pdb")
    expected = set(request.fixed_atoms) | set(request.generated_atoms)
    if set(native) != expected or set(source) != set(request.fixed_atoms):
        raise ValueError("MODELLER atom ownership changed")
    axis, fixed = tuple(sorted(native)), tuple(sorted(source))
    transform = weighted_kabsch(np.array([native[a][2] for a in fixed]), np.array([source[a][2] for a in fixed]))
    aligned = transform.apply(np.array([native[a][2] for a in axis]))
    coordinates = dict(zip(axis, aligned, strict=True))
    drift = np.array([np.linalg.norm(coordinates[a]-source[a][2]) for a in fixed])
    mapped = (output / "proposal/mapped-native.pdb").read_text()
    native_text = rewrite_generated_coordinates(mapped, coordinates)
    source_text = source_path.read_text()
    generated_axis = tuple(request.generated_atoms)
    restored = materialize_candidate_pdb(source_text, generated_axis, [coordinates[a] for a in generated_axis],
        tuple(native[a][0] for a in generated_axis), request, elements=tuple(native[a][1] for a in generated_axis))
    for name, text in (("aligned-native.pdb", native_text), ("source-restored.pdb", restored)):
        with (output / name).open("x") as stream:
            stream.write(text)
    native_structure, source_structure = _structure(native_text, "native.pdb"), _structure(source_text, "source.pdb")
    native_assessment, _ = _assessment(request, source_structure, native_structure)
    native_l, native_chirality_failure = True, None
    try:
        assert_all_l(native_structure.pdb.topology, native_structure.pdb.positions)
    except ChiralityError as exc:
        native_l, native_chirality_failure = False, str(exc)
    restored_assessment = assessment_adapter(request, source_text)(restored)
    restored_structure = _structure(restored, "restored.pdb")
    restored_l, restored_chirality_failure = True, None
    try:
        assert_all_l(restored_structure.pdb.topology, restored_structure.pdb.positions)
    except ChiralityError as exc:
        restored_l, restored_chirality_failure = False, str(exc)
    with (output / "native-transport.npz").open("xb") as stream:
        np.savez_compressed(stream, native=np.array([native[a][2] for a in axis]), aligned=aligned,
                            rotation=transform.rotation, translation=transform.translation)
    save(output / "result.json", {"case": case.name, "backend": "MODELLER", "seed": SEED,
        "native_assessment": asdict(native_assessment), "native_global_all_l": native_l,
        "native_chirality_failure": native_chirality_failure,
        "source_restored_assessment": asdict(restored_assessment), "source_restored_global_all_l": restored_l,
        "source_restored_chirality_failure": restored_chirality_failure,
        "source_restored_geometry_and_L_passed": restored_assessment.geometry_passed and not restored_assessment.admission_failures and restored_l,
        "native_observed_drift_after_global_observed_fit": {"rms_angstrom": float(np.sqrt(np.mean(drift**2))), "max_angstrom": float(drift.max())},
        "full_atom_identity_axis": [asdict(a) for a in axis], "request_sha256": digest(case / "request.json"),
        "native_sha256": digest(output / "proposal/native.pdb"), "restored_sha256": digest(output / "source-restored.pdb"),
        "neural_calls": 0, "repair_calls": 0, "manual_approval_or_publication": False,
        "no_reference_coordinates_or_RMSD": True, "not_context_or_seed_matched_to_local_Protenix": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--case", type=Path)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve(), args.previous.resolve() if args.previous else None)
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        worker(args.case.resolve(), args.output.resolve())
