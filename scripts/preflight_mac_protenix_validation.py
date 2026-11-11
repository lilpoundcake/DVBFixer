#!/usr/bin/env python3
"""Real Protenix feature/atom-mask preflight without checkpoint loading/sampling."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def run(case: Path, assets: Path, kalign: Path, output: Path):
    from checkpoint_gap_smoke import _model_atom_axis, _requested_atom_indices
    from configs.configs_base import configs as configs_base
    from configs.configs_data import data_configs
    from configs.configs_inference import inference_configs
    from configs.configs_model_type import model_configs
    from ml_collections import ConfigDict
    from protenix.config.config import parse_configs
    from protenix.data.inference.infer_dataloader import get_inference_dataloader

    from dvbfixer.model.diffusion.contract import DiffusionRequest

    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    configs = parse_configs(configs={**configs_base, "data": data_configs, **inference_configs}, fill_required_with_null=True)
    configs.update(ConfigDict(model_configs["protenix_base_default_v1.0.0"]))
    configs.model_name = "protenix_base_default_v1.0.0"
    configs.seeds, configs.use_seeds_in_json = [0], False  # data only, no reserved sampler seed
    configs.model.N_cycle = 1
    configs.sample_diffusion.N_sample, configs.sample_diffusion.N_step = 1, 200
    configs.dtype, configs.use_msa, configs.use_template, configs.use_rna_msa = "fp32", False, True, False
    configs.triangle_multiplicative = configs.triangle_attention = "torch"
    configs.enable_diffusion_shared_vars_cache = configs.enable_efficient_fusion = configs.enable_tf32 = False
    configs.sample_diffusion.guidance.enable = False
    configs.data.template.kalign_binary_path = str(kalign)
    configs.input_json_path = str(case / "protenix-v1/input.json")
    configs.dump_dir, configs.load_checkpoint_dir = str(output), str(assets / "checkpoint")
    configs.num_workers = 0
    output.mkdir(parents=True, exist_ok=False)
    data, atom_array, error = next(iter(get_inference_dataloader(configs=configs)))[0]
    if error:
        raise ValueError(error)
    axis, _, _ = _model_atom_axis(atom_array, request)
    selected = _requested_atom_indices(axis, request)
    if {axis[i] for i in selected} != set(request.fixed_atoms) | set(request.generated_atoms):
        raise ValueError("Protenix native axis omits requested observations/terminal OXT")
    hidden = [i for gap in request.gaps for i in range(gap.target_interval.start, gap.target_interval.stop)]
    features, masks = data["input_feature_dict"], {}
    for name in ("template_pseudo_beta_mask", "template_backbone_frame_mask"):
        mask = features[name]
        if bool(mask[..., hidden, :].any()) or bool(mask[..., :, hidden].any()):
            raise ValueError("a hidden target position entered a valid structural template mask")
        masks[name] = {"shape": list(mask.shape), "hidden_entries_zero": True}
    (output / "preflight.json").write_text(json.dumps({"model_axis_atoms": len(axis), "requested_atoms": len(selected),
        "request_sha256": hashlib.sha256((case / "evaluation-mask.json").read_bytes()).hexdigest(),
        "all_requested_observations_including_OXT_represented": True, "template_masks": masks,
        "checkpoint_loaded": False, "sampling": False, "reserved_seeds_sampled": False}, indent=2)+"\n")
    print(json.dumps({"case": case.name, "requested_atoms": len(selected), "checkpoint_loaded": False}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "assets", "kalign", "output"):
        parser.add_argument("--"+name, required=True, type=Path)
    args = parser.parse_args()
    run(args.case.resolve(), args.assets.resolve(), args.kalign.resolve(), args.output.resolve())
