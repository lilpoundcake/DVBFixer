"""Private Protenix v1 MPS smoke; never changes frozen CUDA integration.

Only the 3x3 rigid-fit covariance SVD runs explicitly on CPU: Torch 2.7.1's
MPS SVD emits a built-in fallback warning even with generic fallback disabled.
All model computation and coordinate transforms remain explicitly on MPS.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from checkpoint_gap_smoke import _model_atom_axis, _read_fixed_coordinates, _requested_atom_indices
from reinjection import FixedAtomReinjector

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb


class MPSReinjector(FixedAtomReinjector):
    """Explicit CPU 3x3 SVD, MPS coordinate transform and exact observation copy."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not np.all(self._weights == 1):
            raise ValueError("this private smoke supports uniform rigid-fit weights only")
        self.steps = []
        self.last_pre_overwrite = None

    def __call__(self, coordinates, step_index, step_count):
        if coordinates.device.type != "mps" or coordinates.dtype != torch.float32:
            raise ValueError("private MPS route requires MPS float32 coordinates")
        if coordinates.shape[-2:] != (self._atom_count, 3) or not 0 <= step_index < step_count:
            raise ValueError("invalid sampling coordinates/step")
        indices = torch.as_tensor(self._fixed_indices, device=coordinates.device)
        target = torch.as_tensor(self._fixed_coordinates, device=coordinates.device, dtype=coordinates.dtype)
        mobile = coordinates.index_select(-2, indices)
        center = mobile.mean(dim=-2)
        target_center = target.mean(dim=-2)
        covariance = torch.einsum("...ki,kj->...ij", mobile - center[..., None, :], target - target_center)
        # Explicit and bounded host calculation, NOT generic unsupported-op fallback.
        left, _, right = np.linalg.svd(covariance.detach().cpu().numpy().astype(np.float64))
        correction = np.broadcast_to(np.eye(3), left.shape).copy()
        correction[..., 2, 2] = np.where(np.linalg.det(right.swapaxes(-2, -1) @ left.swapaxes(-2, -1)) > 0, 1, -1)
        rotation = torch.as_tensor(right.swapaxes(-2, -1) @ correction @ left.swapaxes(-2, -1),
                                   device=coordinates.device, dtype=coordinates.dtype)
        synchronized = (coordinates - center[..., None, :]) @ rotation.transpose(-2, -1) + target_center
        overwrite = torch.as_tensor(self._overwrite_indices, device=coordinates.device)
        authoritative = torch.as_tensor(self._overwrite_coordinates, device=coordinates.device, dtype=coordinates.dtype)
        pre = float(torch.linalg.vector_norm(synchronized[..., overwrite, :] - authoritative, dim=-1).max().item())
        self.last_pre_overwrite = synchronized.detach().cpu().numpy().copy()
        synchronized[..., overwrite, :] = authoritative
        post = float((synchronized[..., overwrite, :] - authoritative).abs().max().item())
        self.steps.append({"step_index": step_index, "step_count": step_count,
                           "pre_overwrite_max_angstrom": pre, "post_overwrite_max_angstrom": post})
        return synchronized


def run(args, *, allowed_seeds=(7, 42, 101, 102, 103, 104)):
    from configs.configs_base import configs as configs_base
    from configs.configs_data import data_configs
    from configs.configs_inference import inference_configs
    from configs.configs_model_type import model_configs
    from ml_collections import ConfigDict
    from protenix.config.config import parse_configs
    from protenix.data.inference.infer_dataloader import get_inference_dataloader
    from protenix.utils.distributed import DIST_WRAPPER
    from protenix.utils.seed import seed_everything
    from runner.inference import InferenceRunner, update_inference_configs

    class MPSRunner(InferenceRunner):
        def init_env(self):
            if not torch.backends.mps.is_available() or DIST_WRAPPER.world_size != 1:
                raise RuntimeError("single native MPS process required")
            self.device = torch.device("mps")
            self.use_cuda = False

    request = DiffusionRequest.from_json(args.request.read_text())
    if args.seed not in allowed_seeds or request.seeds != (args.seed,):
        raise ValueError("private cohort seed must match its declared request")
    lock = json.loads(Path(__file__).with_name("profile-lock.json").read_text())
    for name, digest in [("checkpoint/" + lock["checkpoint"]["name"], lock["checkpoint"]["sha256"]),
                         *[(name, lock["auxiliary_artifact_sha256"][name]) for name in
                           ("common/components.cif", "common/components.cif.rdkit_mol.pkl")]]:
        if hashlib.file_digest((args.assets / name).open("rb"), "sha256").hexdigest() != digest:
            raise ValueError(f"asset digest mismatch: {name}")
    configs = parse_configs(configs={**configs_base, "data": data_configs, **inference_configs}, fill_required_with_null=True)
    configs.update(ConfigDict(model_configs["protenix_base_default_v1.0.0"]))
    configs.model_name = "protenix_base_default_v1.0.0"
    configs.seeds = [args.seed]
    configs.use_seeds_in_json = False
    configs.model.N_cycle = 1
    configs.sample_diffusion.N_sample = 1
    configs.sample_diffusion.N_step = 200
    configs.dtype = "fp32"
    configs.use_msa = False
    configs.use_template = True
    configs.use_rna_msa = False
    configs.triangle_multiplicative = "torch"
    configs.triangle_attention = "torch"
    configs.enable_diffusion_shared_vars_cache = False
    configs.enable_efficient_fusion = False
    configs.enable_tf32 = False
    configs.sample_diffusion.guidance.enable = False
    configs.data.template.kalign_binary_path = str(args.kalign.absolute())
    configs.input_json_path = str(args.input.absolute())
    configs.dump_dir = str(args.output.absolute())
    configs.load_checkpoint_dir = str((args.assets / "checkpoint").absolute())
    configs.num_workers = 0
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "config.json").write_text(configs.to_json(indent=2))
    runner = MPSRunner(configs)
    dataloader = get_inference_dataloader(configs=runner.configs)
    data, atom_array, error = next(iter(dataloader))[0]
    if error:
        raise RuntimeError(error)
    axis, names, elements = _model_atom_axis(atom_array, request)
    source_path = args.request.absolute().parent / request.normalized_pdb.path
    callback = MPSReinjector(axis, _read_fixed_coordinates(source_path, request.fixed_atoms))
    runner.model.dvbfixer_step_callback = callback
    runner.update_model_configs(update_inference_configs(runner.configs, data["N_token"].item()))
    seed_everything(args.seed, deterministic=False)
    try:
        prediction = runner.predict(data)
    finally:
        (args.output / "reinjection-steps.json").write_text(json.dumps(callback.steps, indent=2) + "\n")
        if callback.last_pre_overwrite is not None:
            np.savez_compressed(args.output / "last-pre-overwrite.npz", coordinates=callback.last_pre_overwrite)
    xyz = prediction["coordinate"]
    if xyz.device.type != "mps" or len(callback.steps) != 200 or any(s["post_overwrite_max_angstrom"] for s in callback.steps):
        raise RuntimeError("MPS/per-step reinjection conformance failed")
    native = xyz[0, :len(axis)].detach().cpu().numpy().astype(np.float64)
    np.savez_compressed(args.output / "native-arrays.npz", coordinates=native,
                        atom_name=np.array([k.atom_name for k in axis]), residue_name=np.array(names))
    source = source_path.read_text()
    selected = _requested_atom_indices(axis, request)
    text = materialize_candidate_pdb(source, tuple(axis[i] for i in selected), native[list(selected)],
                                     tuple(names[i] for i in selected), request,
                                     elements=tuple(elements[i] for i in selected))
    (args.output / "candidate.pdb").write_text(text)
    (args.output / "reinjection-steps.json").write_text(json.dumps(callback.steps, indent=2) + "\n")
    (args.output / "execution-evidence.json").write_text(json.dumps({"device": str(xyz.device), "dtype": str(xyz.dtype),
        "steps": len(callback.steps), "cpu_alignment_only": "explicit NumPy SVD on 3x3 covariance each step",
        "upstream_terminal_regularization_not_called": True, "scientific_qualification": False,
        "unrequested_atoms_excluded_only_at_materialization": len(axis) - len(selected)}, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("request", "input", "assets", "kalign", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--seed", type=int, default=7)
    run(parser.parse_args())
