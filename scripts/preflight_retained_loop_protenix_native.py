#!/usr/bin/env python3
"""Pinned native feature readback at existing sizes; no model or sampling calls."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import digest, save

PROJECT = Path(__file__).resolve().parents[1]
RUNTIME = PROJECT / ".artifacts/mac-reconstruction-2026-11-11"
CLONE = RUNTIME / "engines/protenix-v1"
NATIVE_PYTHON = RUNTIME / "envs/protenix/bin/python"
ASSETS = RUNTIME / "assets/protenix"
KALIGN = RUNTIME / "tools/bin/kalign"


def check_template_mask(mask, generated_tokens, token_count):
    array = np.asarray(mask)
    hidden = np.asarray(generated_tokens, dtype=bool)
    if hidden.shape != (token_count,) or not hidden.any() or hidden.all():
        raise ValueError("template check needs explicit generated and observed native tokens")
    if array.ndim >= 3 and array.shape[-2:] == (token_count, token_count):
        if array[..., hidden, :].any() or array[..., :, hidden].any():
            raise ValueError("generated native token entered a valid template pair mask")
        coverage = array.astype(bool).reshape(-1, token_count, token_count).any(axis=(0, 2))
    elif array.ndim in (1, 2) and array.shape[-1] == token_count:
        if array[..., hidden].any():
            raise ValueError("generated native token entered a valid template mask")
        coverage = array.astype(bool).reshape(-1, token_count).any(axis=0)
    else:
        raise ValueError(f"unsupported native template mask shape: {array.shape}")
    if not coverage[~hidden].all():
        raise ValueError("native template masks lost observed context token coverage")
    return {"shape": list(array.shape), "generated_tokens_zero": True,
            "all_observed_context_tokens_covered": True}


def token_ownership(axis, atom_to_token, generated_atoms, observed_atoms=()):
    mapping = np.asarray(atom_to_token)
    if mapping.shape != (len(axis),) or not np.issubdtype(mapping.dtype, np.integer) or mapping.min() < 0:
        raise ValueError("native atom-to-token mapping invalid")
    count = int(mapping.max())+1
    generated_residues = {(a.chain, a.residue_number, a.insertion_code) for a in generated_atoms}
    if generated_residues & {(a.chain, a.residue_number, a.insertion_code) for a in observed_atoms}:
        raise ValueError("this whole-gap preflight cannot reinterpret partial-residue ownership")
    residues, generated = {}, {}
    for atom, index in zip(axis, mapping, strict=True):
        index = int(index)
        residue = (atom.chain, atom.residue_number, atom.insertion_code)
        # Native tokenization can add disposable terminal atoms; a generated
        # residue's extra OXT must not turn its token into an observed position.
        role = residue in generated_residues
        if (index in residues and residues[index] != residue) or (index in generated and generated[index] != role):
            raise ValueError("native token mixes source residues or generated/observed ownership")
        residues[index], generated[index] = residue, role
    if set(residues) != set(range(count)) or len(set(residues.values())) != count:
        raise ValueError("native token residue axis is incomplete or duplicated")
    return np.array([generated[i] for i in range(count)], dtype=bool)


def worker(case, output):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("native preflight needs fresh guarded reservation")
    os.environ["PROTENIX_ROOT_DIR"] = str(ASSETS)
    os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "0"
    sys.path.insert(0, str(CLONE))
    sys.path.insert(0, str(PROJECT / "deploy/protenix-v1"))
    import torch

    # Guards are installed before importing upstream data modules, not just
    # before iterating the loader: import-time asset access must also fail closed.
    def forbidden(*args, **kwargs):
        raise RuntimeError("native data-only preflight forbids checkpoint/model/network execution")

    import socket

    torch.load = forbidden
    socket.create_connection = forbidden
    socket.socket.connect = forbidden
    socket.socket.connect_ex = forbidden
    handle = torch.nn.modules.module.register_module_forward_pre_hook(forbidden)

    from checkpoint_gap_smoke import _model_atom_axis
    from configs.configs_base import configs as configs_base
    from configs.configs_data import data_configs
    from configs.configs_inference import inference_configs
    from configs.configs_model_type import model_configs
    from ml_collections import ConfigDict
    from protenix.config.config import parse_configs
    from protenix.data.inference.infer_dataloader import get_inference_dataloader

    from dvbfixer.model.diffusion.contract import DiffusionRequest

    request_path = case / "context-request.json"
    request = DiffusionRequest.from_json(request_path.read_text())
    configs = parse_configs(configs={**configs_base, "data": data_configs, **inference_configs}, fill_required_with_null=True)
    configs.update(ConfigDict(model_configs["protenix_base_default_v1.0.0"]))
    configs.model_name = "protenix_base_default_v1.0.0"
    configs.seeds, configs.use_seeds_in_json = [0], False
    configs.model.N_cycle = 1
    configs.sample_diffusion.N_sample, configs.sample_diffusion.N_step = 1, 200
    configs.dtype, configs.use_msa, configs.use_template, configs.use_rna_msa = "fp32", False, True, False
    configs.triangle_multiplicative = configs.triangle_attention = "torch"
    configs.enable_diffusion_shared_vars_cache = configs.enable_efficient_fusion = configs.enable_tf32 = False
    configs.sample_diffusion.guidance.enable = False
    configs.data.template.kalign_binary_path = str(KALIGN)
    configs.input_json_path = str(case / "encoded/input.json")
    configs.dump_dir, configs.load_checkpoint_dir = str(output / "native-data"), str(ASSETS / "checkpoint")
    configs.num_workers = 0
    try:
        data, atom_array, error = next(iter(get_inference_dataloader(configs=configs)))[0]
        if error:
            raise ValueError(error)
        axis, _, _ = _model_atom_axis(atom_array, request)
        features = data["input_feature_dict"]
        hidden = token_ownership(axis, features["atom_to_token_idx"].cpu().numpy(), set(request.generated_atoms),
                                 (item.atom for item in request.protein_context.fixed_conditioning))
        count = len(hidden)
        masks = {name: check_template_mask(features[name].cpu().numpy(), hidden, count)
                 for name in ("template_pseudo_beta_mask", "template_backbone_frame_mask")}
        expected_tokens = sum(s.stop-s.start for s in request.protein_context.spans)
        if count != expected_tokens:
            raise ValueError("native features changed selected context size")
        np.savez_compressed(output / "native-axis.npz", atom_to_token=features["atom_to_token_idx"].cpu().numpy(),
                            generated_tokens=hidden)
        save(output / "result.json", {"case": case.name, "native_tokens": count, "native_atoms": len(axis),
            "generated_tokens": int(hidden.sum()), "template_masks": masks,
            "request_sha256": digest(request_path), "checkpoint_loaded": False,
            "model_forward_executed": False, "sampling": False, "limits_changed": False,
            "preprocessing_on_CPU_is_not_an_inference_fallback": True})
    finally:
        handle.remove()


def declare(encoded, output):
    manifest = json.loads((encoded / "results.json").read_text())
    files = [Path(__file__), encoded / "results.json", encoded / "declaration.json", KALIGN,
             PROJECT / "deploy/protenix-v1/checkpoint_gap_smoke.py"]
    files += list((CLONE / "configs").glob("*.py")) + list((CLONE / "protenix/data").rglob("*.py"))
    files += [p for p in (ASSETS / "common").rglob("*") if p.is_file()]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    jobs = []
    for row in manifest["rows"]:
        case = encoded / row["entry_id"]
        if digest(case / "encoded/input.json") != row["input_json_sha256"]:
            raise ValueError("retained encoded input changed")
        files += [p for p in case.rglob("*") if p.is_file()]
        case_id = row["entry_id"]+"-native-data"
        jobs.append({"case_id": case_id, "case": str(case), "seed": 0,
            "command": [str(NATIVE_PYTHON), "-m", "scripts.preflight_retained_loop_protenix_native", "worker",
                "--case", str(case), "--output", str(output / "execution" / f"{case_id}-seed-0")]})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "plan.json", {"python": str(NATIVE_PYTHON), "python_sha256": digest(NATIVE_PYTHON),
        "jobs": jobs, "count": len(jobs), "source_code_and_input_sha256": {str(p): digest(p) for p in files},
        "purpose": "native data/mask/identity correctness only, current354/365-token inputs, not size calibration",
        "model_weights_or_sampling_calls": 0, "no_downloads_or_model_forward": True,
        "resource_guard": "shared serial CPU executor,1800s/job,RSS8GiB,pressure/headroom,10GiB evidence cap",
        "seed0_is_data_only_not_sampler_seed": True})


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
