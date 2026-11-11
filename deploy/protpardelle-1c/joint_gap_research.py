#!/usr/bin/env python3
"""Private joint-gap MPS pilot, not a public RunnerResult or profile override.

One model.sample call holds all declared gaps in one shared coordinate tensor.
Original observed coordinates alone condition it. Saved predictions are used
only AFTER sampling for a labelled full-structure graft, never as a motif.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path


def load(name):
    spec = importlib.util.spec_from_file_location("joint_" + name, Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main():
    from dvbfixer.model.diffusion import cc94
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.joint_gap_research import (
        COMPLETE_REVISION,
        MOTIF_COMPONENT_REVISION,
        REVISION,
        STRUCTURE_REVISION,
        build_joint_gap_view,
        build_joint_structure_view,
        complete_context_admission,
        motif_component_admission,
        plan_complete_context,
        plan_motif_component,
    )
    from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-request", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True, help="Output graft frame, NOT sampling conditioning")
    parser.add_argument("--gap-indices", nargs="+", type=int, required=True)
    parser.add_argument("--whole-structure", action="store_true",
        help="Select every gap; does not merge separate neighbor components")
    parser.add_argument("--complete-observed-context", action="store_true",
        help="Explicit historical v3 full-complex experiment; NOT the default motif policy")
    parser.add_argument("--legacy-local-union", action="store_true",
        help="Explicit historical v1/v2 replay only; permits incomplete observed context")
    parser.add_argument("--prepare-only", action="store_true",
        help="Write complete inputs and admission without loading weights or sampling, including over-budget plans")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--refinement-python", type=Path, required=True, help="Pinned runtime preflight only")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--acknowledge-private-joint", action="store_true", required=True)
    args = parser.parse_args()
    if args.complete_observed_context and args.legacy_local_union:
        parser.error("complete context and historical local union are mutually exclusive")
    args.output.mkdir(parents=True, exist_ok=False)
    directory = args.output.resolve(strict=True)

    def write(name, value):
        data = value if isinstance(value, bytes) else (json.dumps(value, sort_keys=True, ensure_ascii=False,
            allow_nan=False, indent=2) + "\n").encode()
        with (directory / name).open("xb") as stream:
            stream.write(data)

    def read(path):
        return read_private_file_bytes(path.parent.resolve(strict=True), path.name, max_bytes=16_000_000)

    trace = {"kind": ((STRUCTURE_REVISION if args.whole_structure else REVISION)
                      if args.legacy_local_union else
                      COMPLETE_REVISION if args.complete_observed_context else MOTIF_COMPONENT_REVISION),
        "schema_version": 1, "public_profile_enabled": False,
        "scientific_qualification": False, "published": False, "status": "admission",
        "seed": args.seed, "gap_indices": args.gap_indices, "sample_invocations": 0, "ode_updates": 0,
        "step_scale": 1.0, "s_churn": 0.0, "chain_phase": "native", "fallback_disabled": True}

    def update():
        # A trace is diagnostics, never a completion index. Retain it on failures.
        (directory / "trace.json").write_text(json.dumps(trace, ensure_ascii=False, indent=2) + "\n")

    update()
    started = time.monotonic()
    probe = None
    try:
        parent_data, source, baseline = read(args.parent_request), read(args.source), read(args.baseline)
        parent = DiffusionRequest.from_json(parent_data.decode())
        if args.whole_structure:
            if args.gap_indices != list(range(len(parent.gaps))):
                raise ValueError("whole-structure mode must declare every gap exactly once")
        if args.complete_observed_context:
            view = plan_complete_context(parent, source, tuple(args.gap_indices), seed=args.seed)
            admission = complete_context_admission(view)
        elif not args.legacy_local_union:
            view = plan_motif_component(parent, source, tuple(args.gap_indices), seed=args.seed)
            admission = motif_component_admission(view)
        elif args.whole_structure:
            view = build_joint_structure_view(parent, source, seed=args.seed)
            admission = {"supported": True, "historical_incomplete_context": True}
        else:
            view = build_joint_gap_view(parent, source, tuple(args.gap_indices), seed=args.seed)
            admission = {"supported": True, "historical_incomplete_context": True}
        trace["kind"] = view.request.protein_context.policy.name
        write("parent-request.json", parent_data)
        write("source.pdb", source)
        write("baseline-output-only.pdb", baseline)
        write("research-group.json", view.request.to_json().encode())
        write("engine-view.json", asdict(view))
        write("admission.json", admission)
        write("adapter-source.py", Path(__file__).read_bytes())
        import dvbfixer.model.diffusion.joint_gap_research as transport

        write("transport-source.py", Path(transport.__file__).read_bytes())
        adapter = load("checkpoint_gap_smoke")
        engine = view.engine_view.request
        source_to_engine = {s: e for e, s in view.engine_view.atom_to_source}
        engine_to_source = dict(view.engine_view.atom_to_source)
        motif = directory / "motif.pdb"
        adapter._write_motif_input(directory / "source.pdb", motif, engine, source_to_engine)
        placement, fixed_positions = adapter._motif_placement(engine)
        trace.update(tokens=len(view.positions), generated_gaps=len(engine.gaps),
            generated_atoms=len(view.request.generated_atoms), parent_request_sha256=view.parent_request_sha256,
            source_sha256=hashlib.sha256(source).hexdigest(), motif_placement=placement,
            physical_positions=view.positions, physical_chains=view.chains, fixed_positions=fixed_positions,
            execution_admission=admission)
        if args.prepare_only:
            trace.update(status="input-prepared" if admission["supported"] else "input-prepared-budget-blocked",
                         model_weights_loaded=False)
            update()
            print(json.dumps(trace), flush=True)
            return
        if not admission["supported"]:
            raise ValueError("input execution refused: " + "; ".join(admission["reasons"]))
        launcher = load("cc94_runner")
        preflight = launcher.preflight_report(cc94.MPS_PROFILE, args.checkpoint, args.refinement_python)
        write("preflight.json", preflight.to_json().encode())
        if not preflight.passed:
            raise ValueError("pinned cc94 preflight failed: " + "; ".join(i.message for i in preflight.issues))
        os.environ["PROTPARDELLE_MODEL_PARAMS"] = str(args.checkpoint.parent.parent.resolve())
        import numpy as np
        import torch
        from protpardelle.common import residue_constants
        from protpardelle.core.models import load_model
        from protpardelle.data.pdb_io import read_pdb
        from protpardelle.utils import seed_everything

        from dvbfixer.model.diffusion.atom_completion_trace import completion_conditioning
        from dvbfixer.model.diffusion.geometry import weighted_kabsch

        protein, _, labels = read_pdb(motif)
        rows = {(int(chain), int(index)): row for row, (chain, index) in enumerate(
            zip(protein.chain_index, protein.residue_index, strict=True))}
        expected_atoms = {}
        for line in motif.read_text().splitlines():
            if line.startswith("ATOM  "):
                expected_atoms[(labels[line[21]], int(line[22:26]), line[12:16].strip())] = np.asarray(
                    [float(line[i:i + 8]) for i in (30, 38, 46)])
        decoded = set()
        max_error = 0.0
        for (chain, index), row in rows.items():
            for name, column in residue_constants.atom_order.items():
                if protein.atom_mask[row, column]:
                    key = (chain, index, name)
                    if key not in expected_atoms:
                        raise ValueError("pinned reader added a motif atom")
                    decoded.add(key)
                    max_error = max(max_error, float(np.linalg.norm(protein.atom_positions[row, column] - expected_atoms[key])))
        if decoded != set(expected_atoms) or max_error > 1e-5:
            raise ValueError("pinned reader changed the joint conditioning mask/coordinates")
        trace["reader_max_decode_error_angstrom"] = max_error
        write("reader-axis.json", {"chain_labels": labels, "residue_index": protein.residue_index.tolist(),
                                  "chain_index": protein.chain_index.tolist(), "fixed_atoms": len(decoded)})
        device = adapter._resolve_device(torch, "mps")
        model = load_model(args.checkpoint.parent.parent / "configs/cc94.yaml", args.checkpoint, device=str(device))
        state = torch.load(args.checkpoint, map_location="cpu", weights_only=False)["model_state_dict"]
        model.load_state_dict(state, strict=True)
        del state
        adapter._assert_model_device(model, device)
        if model.task != "ai-allatom-nomask" or model.config.model.struct_model.uvit.position_embedding_type != "relative":
            raise ValueError("cc94 task/position encoding mismatch")
        seed_everything(args.seed)
        torch.mps.manual_seed(args.seed)
        sequence = "".join(t.sequence for t in engine.target_sequences)
        seq_mask, _, _ = model.make_seq_mask_for_sampling(
            prot_lens_per_chain=torch.tensor([[len(t.sequence) for t in engine.target_sequences]]))
        seq_mask = seq_mask.to(device)
        residue_index = torch.tensor([view.positions], dtype=torch.long, device=device)
        chain_index = torch.tensor([view.chains], dtype=torch.long, device=device)
        aatype = torch.tensor([[residue_constants.restype_order[a] for a in sequence]], dtype=torch.long, device=device)
        class GuardedProbe(adapter._RelativePositionProbe):
            peak_rss = 0
            peak_mps = 0

            def __call__(self, module, hook_args, kwargs):
                super().__call__(module, hook_args, kwargs)
                if self.calls == 1 or self.calls % 10 == 0:
                    import psutil

                    self.peak_rss = max(self.peak_rss, psutil.Process().memory_info().rss)
                    self.peak_mps = max(self.peak_mps, torch.mps.driver_allocated_memory())
                    if self.peak_rss > 8 * 1024**3 or self.peak_mps > 6 * 1024**3:
                        raise RuntimeError("joint research live memory budget exceeded; no fallback")

        probe = GuardedProbe(view.positions, view.chains, device)
        trace.update(status="sampling", sample_invocations=1, tokens=len(sequence),
            generated_gaps=len(engine.gaps), generated_atoms=len(view.request.generated_atoms),
            checkpoint_sha256=cc94.CHECKPOINT_SHA256, engine_revision=cc94.ENGINE_REVISION,
            config_sha256=cc94.CONFIG_SHA256, device=str(device), parent_request_sha256=view.parent_request_sha256,
            source_sha256=hashlib.sha256(source).hexdigest(), baseline_output_only_sha256=hashlib.sha256(baseline).hexdigest(),
            motif_placement=placement, fixed_positions=fixed_positions,
            physical_positions=view.positions, physical_chains=view.chains)
        update()
        with adapter._probe_hook(model.struct_model, probe):
            sampled = model.sample(seq_mask=seq_mask, residue_index=residue_index, chain_index=chain_index,
                gt_aatype=aatype, num_steps=500, step_scale=1.0, s_churn=0.0, motif_file_path=str(motif),
                motif_placements_full=[placement], dummy_fill_mode=model.config.data.dummy_fill_mode,
                partial_diffusion={"enabled": False, "pdb_file_path": None, "num_steps": 100},
                conditional_cfg=completion_conditioning(tuple(fixed_positions)), sidechain_mode=False,
                skip_mpnn_proportion=1.0, use_fullmpnn=False, use_fullmpnn_for_final=False,
                jump_steps=False, uniform_steps=True, record_trajectory=False)
        trace.update(ode_updates=probe.calls, status="sampled")
        update()
        if probe.calls != 500 or not torch.equal(sampled["s"], aatype):
            raise RuntimeError("joint sampling sequence/schedule mismatch")
        adapter._assert_tensor_device("joint coordinates", sampled["x"], device)
        native = adapter._coordinates_to_numpy(torch, sampled["x"][0])
        np.savez_compressed(directory / "native.npz", atom37=native)
        identities, xyz, names = adapter._model_atoms(engine, native)
        by_atom = dict(zip(identities, xyz, strict=True))
        fixed = adapter._read_fixed_coordinates(directory / "source.pdb", view.request.fixed_atoms)
        pairs = [(e, s) for e, s in view.engine_view.atom_to_source if e in engine.fixed_atoms and e in by_atom]
        excluded_alignment = [(e, s) for e, s in view.engine_view.atom_to_source if e in engine.fixed_atoms and e not in by_atom]
        if any(e.atom_name != "OXT" for e, _ in excluded_alignment):
            raise ValueError("native model omitted an unsupported fixed alignment atom")
        trace["fixed_atoms_excluded_from_alignment_only"] = [asdict(s) for _, s in excluded_alignment]
        transform = weighted_kabsch(np.asarray([by_atom[e] for e, _ in pairs]), np.asarray([fixed[s] for _, s in pairs]))
        fitted = transform.apply(xyz)
        errors = [float(np.linalg.norm(transform.apply(np.asarray([by_atom[e]]))[0] - fixed[s])) for e, s in pairs]
        generated = set(view.request.generated_atoms)
        proposed = {engine_to_source[e]: coordinate for e, coordinate in zip(identities, fitted, strict=True)
            if e in engine_to_source and engine_to_source[e] in generated}
        if set(proposed) != generated:
            raise ValueError("joint output omitted generated atoms")
        # Materialization is in the CPU process: retain scientific evidence first.
        np.savez_compressed(directory / "generated.npz", coordinates=np.asarray([proposed[a] for a in view.request.generated_atoms]))
        write("generated-axis.json", [asdict(a) for a in view.request.generated_atoms])
        trace.update(status="native-proposal-retained", native_sha256=hashlib.sha256((directory / "native.npz").read_bytes()).hexdigest(),
            conditioning_fit_rmsd_angstrom=float(np.sqrt(np.mean(np.square(errors)))),
            conditioning_fit_max_angstrom=max(errors), generated_axis_sha256=hashlib.sha256((directory / "generated-axis.json").read_bytes()).hexdigest())
        print(json.dumps(trace), flush=True)
    except BaseException as exc:
        trace.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        if probe is not None:
            trace["ode_updates"] = probe.calls
            trace["peak_process_rss_bytes"] = probe.peak_rss
            trace["peak_mps_driver_bytes"] = probe.peak_mps
        trace["wall_seconds"] = time.monotonic() - started
        update()


if __name__ == "__main__":
    main()
