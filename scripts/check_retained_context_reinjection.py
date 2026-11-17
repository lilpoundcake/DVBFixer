#!/usr/bin/env python3
"""Synthetic MPS callback conformance, never a diffusion proposal or repair."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import digest, save

PROJECT = Path(__file__).resolve().parents[1]
RUNTIME = PROJECT / ".artifacts/mac-reconstruction-2026-11-11"


def represented_fixed(request, atom_axis, coordinates):
    plan = request.protein_context
    if plan is None:
        raise ValueError("context callback control requires an explicit plan")
    selected = {item.atom for item in plan.fixed_conditioning}
    preserved = set(plan.preserved_only)
    if selected & preserved or selected | preserved != set(request.fixed_atoms):
        raise ValueError("represented and preserved observations must partition the original fixed mask")
    if not selected <= set(atom_axis) or preserved & set(atom_axis):
        raise ValueError("callback must not condition on preserved-only observations")
    if not set(request.fixed_atoms) <= coordinates.keys():
        raise ValueError("original observation coordinates are incomplete")
    return {a: coordinates[a] for a in selected}


def worker(case, output):
    os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "0"
    sys.path.insert(0, str(PROJECT / "deploy/protenix-v1"))
    import torch
    from mps_research_smoke import MPSReinjector

    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
    from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE
    from dvbfixer.model.diffusion.protein_context import context_evidence, source_coordinates

    if not torch.backends.mps.is_available():
        raise RuntimeError("MPS control requires real MPS; no CPU inference fallback")
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("callback control requires fresh guarded reservation")
    request_path = case / "context-request.json"
    request = DiffusionRequest.from_json(request_path.read_text())
    source = case / request.normalized_pdb.path
    if digest(source) != request.normalized_pdb.sha256:
        raise ValueError("retained source changed")
    coordinates = source_coordinates(request, source.read_bytes())
    axis = tuple(sorted({item.atom for item in request.protein_context.fixed_conditioning} | set(request.generated_atoms)))
    selected = represented_fixed(request, axis, coordinates)
    try:
        MPSReinjector(axis, coordinates)
    except ValueError as exc:
        if "absent from sampled atoms" not in str(exc):
            raise
        legacy_refusal = str(exc)
    else:
        raise ValueError("negative control must expose omitted observations in legacy callback wiring")
    evidence = context_evidence(request, axis)
    callback = MPSReinjector(axis, selected, alignment_atoms=evidence.alignment_atoms)
    center = np.mean(list(selected.values()), axis=0)
    # Generated positions are synthetic arithmetic fixtures, not predictions.
    points = np.array([selected[a] if a in selected else center+np.array([i*.001, .2, -.3])
                       for i, a in enumerate(axis)], dtype=np.float32)
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]], dtype=np.float32)
    transformed = points @ rotation + np.array([3., -2., 1.], dtype=np.float32)
    corrected = callback(torch.tensor(transformed, device="mps"), 0, 2)
    selected_indices = [i for i, a in enumerate(axis) if a in selected]
    generated_indices = [i for i, a in enumerate(axis) if a in request.generated_atoms]
    target = np.array([selected[axis[i]] for i in selected_indices], dtype=np.float32)
    actual = corrected.detach().cpu().numpy()
    if not np.array_equal(actual[selected_indices], target):
        raise ValueError("represented fixed observations not restored exactly at float32 precision")
    max_generated_rigid_error = float(np.max(np.linalg.norm(actual[generated_indices]-points[generated_indices], axis=1)))
    if max_generated_rigid_error > .0005:
        raise ValueError("synthetic global-frame control failed")
    perturbed = corrected.clone()
    perturbed[selected_indices[0]] += .25
    corrected = callback(perturbed, 1, 2)
    actual = corrected.detach().cpu().numpy()
    if not np.array_equal(actual[selected_indices], target):
        raise ValueError("represented fixed observation perturbation survived projection")
    targets = {t.chain: t.sequence for t in request.target_sequences}
    names = {(r.chain, r.residue_number, r.insertion_code): _CANONICAL_ONE_TO_THREE[targets[placement.chain][i]]
             for placement in request.sequence_placements
             for i, r in zip(placement.observed_target_indices, placement.observed_residues, strict=True)}
    for gap in request.gaps:
        sequence = next(t.sequence for t in request.target_sequences if t.chain == gap.chain)
        for i, r in zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True):
            names[r.chain, r.residue_number, r.insertion_code] = _CANONICAL_ONE_TO_THREE[sequence[i]]
    text = materialize_candidate_pdb(source.read_text(), axis, actual,
        tuple(names[a.chain, a.residue_number, a.insertion_code] for a in axis), request,
        elements=tuple(a.atom_name[0] for a in axis))

    def fixed_records(pdb):
        from dvbfixer.model.diffusion.contract import AtomIdentity

        fixed = set(request.fixed_atoms)
        return {AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()): line
                for line in pdb.splitlines(keepends=True) if line.startswith(("ATOM  ", "HETATM"))
                and AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()) in fixed}

    if fixed_records(text) != fixed_records(source.read_text()):
        raise ValueError("synthetic materialization changed original fixed records")
    save(output / "result.json", {"case": case.name, "legacy_all_observations_refused": legacy_refusal,
        "represented_fixed_atoms": len(selected), "preserved_only_atoms": len(request.protein_context.preserved_only),
        "alignment_atoms": len(evidence.alignment_atoms), "synthetic_MPS_callback_calls": len(callback.steps),
        "generated_global_rigid_error_angstrom": max_generated_rigid_error,
        "all_source_fixed_records_preserved_in_memory": True,
        "driver_bytes_after_controls": torch.mps.driver_allocated_memory(),
        "request_sha256": digest(request_path), "checkpoint_loaded": False, "diffusion_steps": 0,
        "synthetic_coordinates_not_saved_as_candidate": True, "production_or_frozen_runner_changed": False,
        "not_native_sampler_conformance_or_scientific_geometry_acceptance": True})


def declare(encoded, output):
    output.mkdir(parents=True, exist_ok=False)
    native_python = RUNTIME / "envs/protenix/bin/python"
    files = [Path(__file__), PROJECT / "deploy/protenix-v1/mps_research_smoke.py",
             PROJECT / "deploy/protenix-v1/reinjection.py", PROJECT / "deploy/protenix-v1/checkpoint_gap_smoke.py"]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    jobs = []
    for entry in ("1NON", "2BE1"):
        case = encoded / entry
        files += [p for p in case.rglob("*") if p.is_file()]
        jobs.append({"case_id": entry+"-callback", "case": str(case), "seed": 0,
            "command": [str(native_python), "-m", "scripts.check_retained_context_reinjection", "worker",
                "--case", str(case), "--output", str(output / "execution" / f"{entry}-callback-seed-0")]})
    save(output / "plan.json", {"python": str(native_python), "python_sha256": digest(native_python),
        "jobs": jobs, "count": len(jobs), "source_code_and_input_sha256": {str(p): digest(p) for p in files},
        "purpose": "synthetic callback/frame/materialization controls, not sampler steps or model inference",
        "new_checkpoints_or_sampling": 0, "limits_changed": False,
        "resource_guard": "existing shared serial CPU/RSS8GiB/1800s/host-pressure guard; tiny explicit MPS coordinate operations"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--encoded", type=Path)
    parser.add_argument("--case", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.encoded.resolve(), args.output.resolve())
    elif args.mode == "run":
        from scripts.run_mac_modeller_comparison import execute

        execute(args.output.resolve())
    else:
        worker(args.case.resolve(), args.output.resolve())
