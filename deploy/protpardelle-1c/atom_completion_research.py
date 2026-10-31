#!/usr/bin/env python3
"""Private cc94 completion experiment; not a production runner or public profile.

Retains original inputs, native coordinates, masks, trace and rejected outputs.
Never emits a production RunnerResult or publishes a bundle.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

from dvbfixer.model.diffusion import cc94
from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
from dvbfixer.model.diffusion.atom_completion_trace import completion_conditioning
from dvbfixer.model.diffusion.private_artifact import read_private_artifact


def _load(name):
    spec = importlib.util.spec_from_file_location(f"completion_{name}", Path(__file__).with_name(name + ".py"))
    if spec is None or spec.loader is None:
        raise RuntimeError("pinned cc94 adapter is unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def conditioning_config(fixed_positions):
    """Frozen cc94 native conditioning, no cc89 phase or reinjection changes."""
    return completion_conditioning(tuple(fixed_positions))


def refine_worker(directory, profile):
    from dvbfixer.model.diffusion.atom_completion_refinement import (
        CompletionRefinementError,
        refine_completion_candidate,
    )

    if any((directory / name).exists() for name in ("refined.pdb", "raw-validation.json", "refined-validation.json", "refinement.json")):
        raise FileExistsError("completion refinement never overwrites a retained attempt")
    request = AtomCompletionRequest.from_json((directory / "request.json").read_text())
    source, raw = (directory / "source.pdb").read_bytes(), (directory / "raw.pdb").read_bytes()
    try:
        result = refine_completion_candidate(source, raw, request, seed=request.seeds[0],
                                             platform_name="CPU" if profile == cc94.MPS_PROFILE else "Reference")
    except CompletionRefinementError as exc:
        _json(directory / "raw-validation.json", asdict(exc.raw_validation))
        _json(directory / "refinement.json", {"status": "failed", "error": str(exc), "published": False})
        return
    (directory / "refined.pdb").write_bytes(result.refined_candidate)
    _json(directory / "raw-validation.json", asdict(result.raw_validation))
    _json(directory / "refined-validation.json", asdict(result.refined_validation))
    _json(directory / "refinement.json", {
        "status": "accepted-private" if result.refined_validation.summary.passed else "rejected",
        "revision": result.revision, "kernel_revision": result.kernel_revision,
        "platform": result.platform, "seed": result.seed,
        "initial_energy_kj_mol": result.initial_energy_kj_mol, "final_energy_kj_mol": result.final_energy_kj_mol,
        "source_sha256": result.refined_validation.source_sha256,
        "raw_sha256": result.raw_validation.candidate_sha256,
        "refined_sha256": result.refined_validation.candidate_sha256,
        "published": False, "scientific_qualification": False,
    })


def run(args):
    from dvbfixer.model.diffusion.atom_completion_materialize import materialize_atom_completion
    from dvbfixer.model.diffusion.atom_completion_sampling import (
        assess_completion_sampling,
        completion_engine_view,
    )
    from dvbfixer.model.diffusion.atom_completion_validation import _coordinates
    from dvbfixer.model.diffusion.geometry import weighted_kabsch

    if not args.acknowledge_experimental_completion:
        raise ValueError("explicit --acknowledge-experimental-completion is required; no public capability is enabled")
    if args.request.is_symlink() or not args.request.is_file() or args.request.stat().st_size > 2_000_000:
        raise ValueError("completion request must be a bounded regular file")
    request_path = args.request.resolve(strict=True)
    request_bytes = request_path.read_bytes()
    request = AtomCompletionRequest.from_json(request_bytes.decode("utf-8"))
    source = read_private_artifact(request_path.parent, request.normalized_pdb,
                                   label="completion source", max_bytes=2_000_000, encoding="utf-8").encode("utf-8")
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    (directory / "source.pdb").write_bytes(source)
    (directory / "request.json").write_bytes(request_bytes)
    trace = {"schema_version": 1, "kind": "private-cc94-atom-completion-trace", "status": "admission",
             "request_sha256": hashlib.sha256(request_bytes).hexdigest(), "source_sha256": hashlib.sha256(source).hexdigest(),
             "adapter_sha256": _digest(Path(__file__)), "checkpoint_sampling_executed": False,
             "public_profile_enabled": False, "published": False, "scientific_qualification": False,
             "seed": request.seeds[0], "ode_updates": 0, "chain_phase": "native", "step_scale": 1.0, "s_churn": 0.0}
    _json(directory / "trace.json", trace)
    started = time.monotonic()
    probe = None
    try:
        admission = assess_completion_sampling(request, source)
        _json(directory / "admission.json", asdict(admission))
        if not admission.supported:
            raise ValueError("unsupported completion source: " + ", ".join(admission.reasons))
        launcher = _load("cc94_runner")
        preflight = launcher.preflight_report(args.profile, args.checkpoint, args.refinement_python)
        (directory / "preflight.json").write_text(preflight.to_json())
        if not preflight.passed:
            raise RuntimeError("; ".join(issue.message for issue in preflight.issues))
        view = completion_engine_view(request, source)
        _json(directory / "engine-view.json", asdict(view))
        motif = directory / "motif.pdb"
        motif.write_text(view.motif_text)
        # Upstream resolves auxiliary model assets at module import time.
        os.environ["PROTPARDELLE_MODEL_PARAMS"] = str(args.checkpoint.parent.parent)
        import numpy as np
        import torch
        from protpardelle.common import residue_constants
        from protpardelle.core.models import load_model
        from protpardelle.data.pdb_io import read_pdb
        from protpardelle.utils import seed_everything

        protein, _, _ = read_pdb(motif)
        row_for_index = {int(index) - 1: row for row, index in enumerate(protein.residue_index)}
        fixed = set(view.represented_fixed_atoms)
        sequence = request.target_sequences[0].sequence
        authoritative = _coordinates(source.decode().splitlines())
        expected_positions = set(view.fixed_positions)
        if set(row_for_index) != expected_positions:
            raise ValueError("pinned reader changed the observed target residue axis")
        decode_errors = []
        for index, row in row_for_index.items():
            if protein.aatype[row] != residue_constants.restype_order[sequence[index]]:
                raise ValueError("pinned reader changed an observed residue's chemistry")
        for atom, index in zip(view.atom_axis, view.atom_target_indices, strict=True):
            present = bool(protein.atom_mask[row_for_index[index], residue_constants.atom_order[atom.atom_name]]) if index in row_for_index else False
            if present != (atom in fixed):
                raise ValueError("pinned reader changed a fixed/partial/gap atom mask")
            if present:
                decode_errors.append(float(np.linalg.norm(
                    protein.atom_positions[row_for_index[index], residue_constants.atom_order[atom.atom_name]] - authoritative[atom])))
        if max(decode_errors, default=0.0) > 1e-5:
            raise ValueError("pinned reader changed observed motif coordinates beyond decoding precision")
        trace["reader_max_decode_error_angstrom"] = max(decode_errors, default=0.0)
        np.savez_compressed(directory / "reader.npz", atom_mask=protein.atom_mask, residue_index=protein.residue_index,
                            atom_positions=protein.atom_positions)
        adapter = _load("checkpoint_gap_smoke")
        device = adapter._resolve_device(torch, "mps" if args.profile == cc94.MPS_PROFILE else "cuda:0")
        model = load_model(args.checkpoint.parent.parent / "configs/cc94.yaml", args.checkpoint, device=str(device))
        state = torch.load(args.checkpoint, map_location="cpu", weights_only=False)["model_state_dict"]
        model.load_state_dict(state, strict=True)
        del state
        adapter._assert_model_device(model, device)
        if model.task != "ai-allatom-nomask" or model.config.model.struct_model.uvit.position_embedding_type != "relative":
            raise ValueError("cc94 task/relative encoding mismatch")
        seed_everything(request.seeds[0])
        if device.type == "mps":
            torch.mps.manual_seed(request.seeds[0])
        seq_mask, _, _ = model.make_seq_mask_for_sampling(prot_lens_per_chain=torch.tensor([[len(sequence)]]))
        seq_mask = seq_mask.to(device)
        residue_index = torch.tensor([view.residue_index], dtype=torch.long, device=device)
        chain_index = torch.tensor([view.chain_index], dtype=torch.long, device=device)
        aatype = torch.tensor([[residue_constants.restype_order[letter] for letter in sequence]], dtype=torch.long, device=device)
        probe = adapter._RelativePositionProbe(view.residue_index, view.chain_index, device)
        trace.update(status="sampling", checkpoint_sampling_started=True, motif_sha256=_digest(motif), checkpoint_sha256=cc94.CHECKPOINT_SHA256,
                     config_sha256=cc94.CONFIG_SHA256, engine_revision=cc94.ENGINE_REVISION, device=str(device),
                     conditioning=conditioning_config(view.fixed_positions))
        _json(directory / "trace.json", trace)
        with adapter._probe_hook(model.struct_model, probe):
            sampled = model.sample(
                seq_mask=seq_mask, residue_index=residue_index, chain_index=chain_index, gt_aatype=aatype,
                num_steps=500, step_scale=1.0, s_churn=0.0, motif_file_path=str(motif),
                motif_placements_full=[view.motif_placement], dummy_fill_mode=model.config.data.dummy_fill_mode,
                partial_diffusion={"enabled": False, "pdb_file_path": None, "num_steps": 100},
                conditional_cfg=conditioning_config(view.fixed_positions), sidechain_mode=False,
                skip_mpnn_proportion=1.0, use_fullmpnn=False, use_fullmpnn_for_final=False, jump_steps=False,
                uniform_steps=True, **({"record_trajectory": False} if device.type == "mps" else {}),
            )
        trace.update(checkpoint_sampling_executed=True, ode_updates=probe.calls)
        _json(directory / "trace.json", trace)
        if probe.calls != 500 or not torch.equal(sampled["s"], aatype):
            raise RuntimeError("cc94 did not preserve its declared sequence and 500-update schedule")
        adapter._assert_tensor_device("completion output", sampled["x"], device)
        native = adapter._coordinates_to_numpy(torch, sampled["x"][0])
        np.savez_compressed(directory / "native.npz", atom37=native)
        if native.shape != (len(sequence), 37, 3) or not np.isfinite(native).all():
            raise ValueError("invalid native cc94 atom37 coordinates")
        coordinates = np.asarray([native[index, residue_constants.atom_order[atom.atom_name]]
                                  for atom, index in zip(view.atom_axis, view.atom_target_indices, strict=True)], dtype=np.float64)
        by_atom = dict(zip(view.atom_axis, coordinates, strict=True))
        transform = weighted_kabsch(np.asarray([by_atom[atom] for atom in view.represented_fixed_atoms]),
                                    np.asarray([authoritative[atom] for atom in view.represented_fixed_atoms]))
        fitted = transform.apply(coordinates)
        fitted_by_atom = dict(zip(view.atom_axis, fitted, strict=True))
        errors = [float(np.linalg.norm(fitted_by_atom[atom] - authoritative[atom])) for atom in view.represented_fixed_atoms]
        trace.update(status="sampled", conditioning_fit_rmsd_angstrom=float(np.sqrt(np.mean(np.square(errors)))),
                     conditioning_fit_max_angstrom=max(errors), represented_fixed_atoms=len(fixed),
                     preserved_only_atoms=len(view.preserved_only_atoms), native_sha256=_digest(directory / "native.npz"))
        raw = materialize_atom_completion(source, request, request.generated_atoms,
                                          [tuple(float(value) for value in fitted_by_atom[atom]) for atom in request.generated_atoms])
        (directory / "raw.pdb").write_bytes(raw)
        trace["raw_sha256"] = _digest(directory / "raw.pdb")
        _json(directory / "trace.json", trace)
        with (directory / "refinement.log").open("x") as log:
            subprocess.run([str(args.refinement_python), str(Path(__file__).resolve()), "--refine-worker",
                            "--output", str(directory), "--profile", args.profile],
                           env=launcher.refinement_environment(args.refinement_python), stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        trace.update(status=json.loads((directory / "refinement.json").read_text())["status"])
    except BaseException as exc:
        trace.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        if probe is not None:
            trace.update(ode_updates=probe.calls, checkpoint_sampling_executed=probe.calls > 0)
        trace["wall_seconds"] = time.monotonic() - started
        _json(directory / "trace.json", trace)
    if trace["status"] == "accepted-private":
        from dvbfixer.model.diffusion.atom_completion_evidence import attest_completion_run

        attest_completion_run(directory, Path(__file__).read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", choices=sorted(cc94.PROFILES), required=True)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--refinement-python", type=Path)
    parser.add_argument("--acknowledge-experimental-completion", action="store_true")
    parser.add_argument("--refine-worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.refine_worker:
        refine_worker(args.output.resolve(), args.profile)
    else:
        if args.request is None or args.checkpoint is None or args.refinement_python is None:
            parser.error("sampling requires --request, --checkpoint and --refinement-python")
        run(args)


if __name__ == "__main__":
    main()
