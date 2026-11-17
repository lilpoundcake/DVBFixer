#!/usr/bin/env python3
"""Fresh seed101 native generated-only initial optimization; no post-hoc clamp."""
from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import asdict, replace
from pathlib import Path

from scripts import run_retained_loop_modeller_controls as baseline
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest

SEED = 101
ROOT = baseline.ROOT


def selected_generated_initial(model, request):
    axis = target_axis(request)
    expected_residues = tuple(r for target in request.target_sequences for r in axis[target.chain])
    native_residues = tuple(model.residues)
    if len(native_residues) != len(expected_residues):
        raise ValueError("MODELLER selection residue axis changed")
    mapping = {residue.index: source for residue, source in zip(native_residues, expected_residues, strict=True)}
    selection = model.select_loop_atoms()
    selected = tuple(AtomIdentity(mapping[atom.residue.index].chain,
        mapping[atom.residue.index].residue_number, mapping[atom.residue.index].insertion_code, atom.name)
        for atom in selection)
    if len(set(selected)) != len(selected):
        raise ValueError("native initial optimization mask is not exactly all generated atoms")
    # MODELLER supplies a terminal OXT even when the frozen completion mask
    # excludes terminal chemistry. Keep that engine-only atom immobile; the
    # existing native normalizer separately records and excludes it on export.
    terminal = expected_residues[-1]
    terminal_oxt = AtomIdentity(terminal.chain, terminal.residue_number, terminal.insertion_code, "OXT")
    if set(selected) - set(request.generated_atoms) == {terminal_oxt}:
        from modeller import Selection

        selection = Selection(*(atom for atom, identity in zip(selection, selected, strict=True)
                                if identity != terminal_oxt))
        selected = tuple(identity for identity in selected if identity != terminal_oxt)
    if len(set(selected)) != len(selected) or set(selected) != set(request.generated_atoms):
        raise ValueError("native initial optimization mask is not exactly all generated atoms")
    return selection, selected


def declare(output):
    output.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), Path(baseline.__file__), baseline.PROJECT / "scripts/run_mac_modeller_proposal.py",
        baseline.PROJECT / "scripts/run_mac_modeller_comparison.py", baseline.PROJECT / "scripts/mac_reconstruction_comparison.py",
        baseline.PROJECT / "src/dvbfixer/model/modeller_run.py", baseline.PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += list((baseline.PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    jobs = []
    for entry in ("1NON", "2BE1"):
        original_case = ROOT / "retained-loop-modeller-controls-attempt02/cases" / entry
        original = original_case / "request.json"
        request = DiffusionRequest.from_json(original.read_text())
        view = replace(request, seeds=(SEED,))
        case = output / "cases" / entry
        case.mkdir(parents=True)
        shutil.copyfile(original_case / "source.pdb", case / "source.pdb")
        (case / "request.json").write_text(view.to_json())
        (case / "alignment.pir").write_text(baseline.multichain_alignment(view))
        save(case / "lineage.json", {"original_request_path": str(original), "original_request_sha256": digest(original),
            "source_mask_context_unchanged": True, "fresh_seed": SEED,
            "method_change": "AutoModel.select_atoms delegates to the same exact gap-only loop selection"})
        files += list(case.iterdir())+[original]
        jobs.append({"case_id": entry, "case": str(case), "seed": SEED,
            "command": [sys.executable, "-m", "scripts.run_retained_modeller_fixation_controls", "worker",
                "--case", str(case), "--output", str(output / "execution" / f"{entry}-seed-{SEED}")]})
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": 2,
        "seed": SEED, "actual_MODELLER_proposal_budget": 2, "repair_calls": 0, "new_neural_calls": 0,
        "native_interface": "documented AutoModel.select_atoms; existing LoopModel gap selection",
        "all_generated_atoms_selected_jointly": True, "native_observation_fixation_not_assumed_until_measured": True,
        "one_initial_and_one_fast_loop_derivative": True, "restoration_policy": "generated-only source copy",
        "full_observed_template_context_unchanged": True, "automatic_retries": 0,
        "resource_guard": "serial CPU1800s/job,RSS8GiB,host-pressure,10GiB evidence",
        "not_final_evaluation_or_public_capability": True})


def worker(case, output, *, seed=SEED):
    from modeller.automodel import LoopModel

    request = DiffusionRequest.from_json((case / "request.json").read_text())
    if request.seeds != (seed,):
        raise ValueError("undeclared fixation-control seed")
    original_select, original_seed = LoopModel.select_atoms, baseline.SEED
    records = []

    def select(model):
        selected, atoms = selected_generated_initial(model, request)
        row = {"native_interface": "AutoModel.select_atoms", "atoms": [asdict(a) for a in atoms],
               "exact_generated_mask": True, "source_observation_motion_measured_downstream": True}
        records.append(row)
        save(output / f"initial-selection-{len(records):03d}.json", row)
        return selected

    LoopModel.select_atoms, baseline.SEED = select, seed
    try:
        baseline.worker(case, output)
        if not records:
            raise ValueError("native initial-selection interface was not invoked")
        save(output / "native-fixation-evidence.json", {"selection_calls": len(records),
            "generated_atoms_each": len(request.generated_atoms), "no_posthoc_native_coordinate_clamp": True,
            "all_generated_selected_jointly": True, "seed": seed})
    finally:
        LoopModel.select_atoms, baseline.SEED = original_select, original_seed


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
        worker(args.case.resolve(), args.output.resolve())
