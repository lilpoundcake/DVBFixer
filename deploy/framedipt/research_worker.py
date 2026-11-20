"""Pinned native FrameDiPT CPU core; no TCR pipeline or sequence design imports."""
from __future__ import annotations

import copy
import hashlib
import json
import sys
import time
from pathlib import Path


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def project_fixed_frames(candidate, reference, fixed):
    """Explicit research sampler intervention; never changes generated frames."""
    import numpy as np

    candidate, reference, fixed = np.asarray(candidate), np.asarray(reference), np.asarray(fixed)
    if candidate.shape[-2:] != reference.shape or reference.shape != (len(fixed), 7) or fixed.dtype != bool:
        raise ValueError("invalid fixed-frame projection axis")
    if not np.isfinite(candidate).all() or not np.isfinite(reference).all():
        raise ValueError("nonfinite fixed-frame projection")
    result = candidate.copy()
    result[..., fixed, :] = reference[fixed]
    if not np.array_equal(result[..., ~fixed, :], candidate[..., ~fixed, :]):
        raise ValueError("projection changed a generated frame")
    return result


def main(root):
    import numpy as np
    import torch
    from experiments import utils as exp_utils
    from framedipt.diffusion.se3_diffuser import SE3Diffuser
    from framedipt.model.score_network import ScoreNetwork
    from framedipt.protein import all_atom
    from framedipt.protein import residue_constants as rc
    from omegaconf import OmegaConf
    from openfold.data import data_transforms
    from openfold.utils import rigid_utils

    plan = json.loads((root / "plan.json").read_text())
    transport = json.loads((root / "transport.json").read_text())
    inputs = dict(np.load(root / "transport.npz", allow_pickle=False))
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.manual_seed(plan["seed"])
    weights = Path(plan["weights"])
    if hashlib.sha256(weights.read_bytes()).hexdigest() != plan["weights_sha256"]:
        raise ValueError("checkpoint bytes changed")
    # Official pinned checkpoint includes its OmegaConf, not just tensor weights.
    checkpoint = torch.load(weights, map_location="cpu", weights_only=False)
    cfg = OmegaConf.load(Path(plan["upstream"]) / "config/base.yaml")
    cfg.model = OmegaConf.merge(cfg.model, checkpoint["conf"].model)
    cfg.diffuser.r3 = checkpoint["conf"].diffuser.r3
    cfg.diffuser.r3.seed = plan["seed"]
    cfg.diffuser.so3.seed = plan["seed"]
    cfg.diffuser.so3.cache_dir = str(root / "so3-cache")
    if not cfg.model.input_aatype or "".join(rc.restypes) != "ARNDCQEGHILKMFPSTWYV":
        raise ValueError("checkpoint is not the known-sequence inpainting model")
    if tuple(rc.atom_types) != ("N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD", "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3", "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT"):
        raise ValueError("native atom37 order differs from transport")
    diffuser = SE3Diffuser(cfg.diffuser)
    model = ScoreNetwork(cfg.model, diffuser, inpainting=True).cpu().eval()
    state = {k.removeprefix("module."): v for k, v in checkpoint["model"].items()}
    model.load_state_dict(state, strict=True)
    del checkpoint, state
    save(root / "runtime.json", {"torch": torch.__version__, "numpy": np.__version__, "python": sys.version,
         "parameters": sum(p.numel() for p in model.parameters()), "device": "CPU", "checkpoint_loaded_strictly": True,
         "upstream_revision": plan["upstream_revision"], "weights_sha256": plan["weights_sha256"],
         "borrowed_base_prefix": sys.base_prefix, "isolated_overlay_prefix": sys.prefix})
    fixed = inputs["fixed_mask"].astype(bool)
    n = len(fixed)
    # One whole-scaffold centering, never independent per-chain fitting.
    center = inputs["atom_positions"][fixed, 1].sum(axis=0) / (fixed.sum() + 1e-5)
    positions = (inputs["atom_positions"] - center) * inputs["atom_mask"][..., None]
    feats = {"aatype": torch.tensor(inputs["aatype"]).long(),
             "all_atom_positions": torch.tensor(positions).double(),
             "all_atom_mask": torch.tensor(inputs["atom_mask"]).double()}
    feats = data_transforms.atom37_to_frames(feats)
    feats = data_transforms.atom37_to_torsion_angles()(feats)
    frames = feats["rigidgroups_gt_frames"][:, 0]
    rotations = torch.eye(3, dtype=frames.dtype).repeat(n, 1, 1)
    translations = torch.zeros((n, 3), dtype=frames.dtype)
    rotations[fixed] = frames[fixed, :3, :3]
    translations[fixed] = frames[fixed, :3, 3]
    rigids = rigid_utils.Rigid(rigid_utils.Rotation(rot_mats=rotations), translations)
    reference_frames = rigids.to_tensor_7()
    torsions = feats["torsion_angles_sin_cos"]
    torsions[~fixed] = 0
    initial = diffuser.sample_ref(n_samples=n, chain_index=inputs["chain_idx"], impute=rigids,
                                  diffuse_mask=~fixed, as_tensor_7=True)["rigids_t"]
    sample = {"aatype": torch.tensor(inputs["aatype"])[None].long(),
              "seq_idx": torch.tensor(inputs["seq_idx"])[None], "chain_idx": torch.tensor(inputs["chain_idx"])[None],
              "fixed_mask": torch.tensor(inputs["fixed_mask"])[None], "res_mask": torch.ones(1, n),
              "torsion_angles_sin_cos": torsions[None], "rigids_0": reference_frames[None],
              "sc_ca_t": torch.zeros(1, n, 3), "rigids_t": initial[None]}
    np.savez_compressed(root / "native-initial.npz", rigids_0=reference_frames.numpy(),
                        rigids_t=initial.numpy(), torsions=torsions.numpy(), center=center,
                        **inputs)
    projection_count = 0
    if plan["frame_policy"] == "reference-projection":
        # Separate, declared arm. Preserve the pinned native denoiser/reverse
        # calculation, then replace only the fixed tensor7 entries to prevent
        # rotation<->rotvec conversion drift. Do not call this unmodified native
        # sampling or conceal its failed predecessor.
        native_reverse = diffuser.reverse
        (root / "projection-states").mkdir()

        def projected_reverse(*args, **kwargs):
            nonlocal projection_count
            result = native_reverse(*args, **kwargs)
            before = result.to_tensor_7().detach().cpu().numpy()
            after = project_fixed_frames(before, reference_frames.numpy(), fixed)
            projection_count += 1
            state_path = root / "projection-states" / f"{projection_count:03d}.npz"
            np.savez_compressed(state_path, before=before, after=after)
            with (root / "projection-trace.jsonl").open("a") as stream:
                stream.write(json.dumps({"step": projection_count, "generated_frames_exact": True,
                                        "state_sha256": hashlib.sha256(state_path.read_bytes()).hexdigest(),
                                        "fixed_tensor7_exact": np.array_equal(after[..., fixed, :], reference_frames.numpy()[None, fixed, :])}) + "\n")
            return rigid_utils.Rigid.from_tensor_7(torch.tensor(after), normalize_quats=False)

        diffuser.reverse = projected_reverse
    trace = []
    expected_consumer = {name: sample[name].clone() for name in ("aatype", "seq_idx", "fixed_mask", "res_mask", "chain_idx")}

    def frame_error(frame):
        candidate = rigid_utils.Rigid.from_tensor_7(frame[0])
        trans = candidate.get_trans().detach().numpy()
        rot = candidate.get_rots().get_rot_mats().detach().numpy()
        dt = float(np.linalg.norm(trans[fixed] - translations.numpy()[fixed], axis=-1).max())
        dr = float(np.max(np.abs(rot[fixed] - rotations.numpy()[fixed])))
        if not np.isfinite([dt, dr]).all() or dt > plan["frame_translation_tolerance_A"] or dr > plan["frame_rotation_tolerance"]:
            raise ValueError(f"native fixed frame drift: {dt} A, rotation {dr}")
        return dt, dr

    def consumer_hook(module, args):
        batch = args[0]
        for name in ("aatype", "seq_idx", "fixed_mask", "res_mask", "chain_idx"):
            if not torch.equal(batch[name], expected_consumer[name]):
                raise ValueError(f"native model consumer changed {name}")
        dt, dr = frame_error(batch["rigids_t"])
        trace.append({"call": len(trace)+1, "t": float(batch["t"][0]), "target_residues": n,
                      "fixed_translation_max_A": dt, "fixed_rotation_max_component": dr})
        with (root / "consumer-trace.jsonl").open("a") as stream:
            stream.write(json.dumps(trace[-1]) + "\n")

    hook = model.register_forward_pre_hook(consumer_hook)
    reference_bb = all_atom.compute_backbone(rigids[None], torsions[None, :, 2], sample["aatype"])[0][0].numpy()
    source_bb_drift = {}
    for atom_name in ("N", "CA", "C", "O", "CB"):
        j = rc.atom_order[atom_name]
        present = fixed & (inputs["atom_mask"][:, j] > 0)
        delta = reference_bb[present, j] + center - inputs["atom_positions"][present, j]
        source_bb_drift[atom_name] = {"count": int(present.sum()),
                                    "max_A": float(np.linalg.norm(delta, axis=-1).max()) if present.any() else None,
                                    "rms_A": float(np.sqrt(np.mean(np.sum(delta**2, axis=-1)))) if present.any() else None}
    save(root / "representation-audit.json", {"observed_backbone_roundtrip": source_bb_drift,
         "center": center.tolist(), "source_sidechains_encoded": False,
         "source_OXT_native_output": False, "native_atomic_preservation_passed": False,
         "note": "ideal residue frames are not exact experimental heavy-atom coordinates"})
    exp_utils._set_t_feats(sample, 1., torch.ones(1), diffuser)
    if plan["mode"] == "preflight":
        with torch.no_grad():
            predicted = model(sample)
        frame_error(predicted["rigids"])
        if not torch.isfinite(predicted["atom37"]).all():
            raise ValueError("nonfinite native backbone forward")
        save(root / "preflight.json", {"consumer_passed": True, "native_model_calls": len(trace),
             "target_residues": n, "fixed_residues": int(fixed.sum()), "generated_residues": int((~fixed).sum()),
             "sampling_completed": False, "frame_contract_passed": True, "atomic_contract_passed": False})
    else:
        # The actual native sampler is called unchanged. It retains its native
        # SciPy/NumPy CPU reverse diffusion, final denoiser step and mask policy.
        start = time.monotonic()
        result = exp_utils.inference_fn(model=model, diffuser=diffuser, data_init=copy.deepcopy(sample),
            num_t=plan["steps"], min_t=plan["min_t"], noise_scale=plan["noise_scale"],
            aux_traj=True, embed_self_conditioning=cfg.model.embed.embed_self_conditioning,
            inpainting=True, input_aatype=True)
        frame_error(torch.tensor(result["rigid_traj"][0]))
        backbone = result["prot_traj"][0, 0].copy()
        backbone[:, :5] += center
        if backbone.shape != (n, 37, 3) or not np.isfinite(backbone).all():
            raise ValueError("native endpoint is incomplete or nonfinite")
        # Only the five emitted backbone/C-beta slots contain predicted atoms.
        np.savez_compressed(root / "native-endpoint.npz", backbone=backbone,
                            final_rigids=result["rigid_traj"][0], center=center,
                            aatype=inputs["aatype"], fixed_mask=inputs["fixed_mask"])
        lines = []
        serial = 0
        for i, residue in enumerate(transport["residues"]):
            for atom_name in ("N", "CA", "C", "O", "CB"):
                if atom_name == "CB" and residue["resname"] == "GLY":
                    continue
                serial += 1
                xyz = backbone[i, rc.atom_order[atom_name]]
                lines.append(f'ATOM  {serial:5d} {atom_name:^4s} {residue["resname"]} {residue["chain"]}{int(residue["residue_number"]):4d}{residue["insertion_code"] or " "}   '
                             f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           {atom_name[0]}  \n')
            if i == n-1 or transport["residues"][i+1]["chain"] != residue["chain"]:
                lines.append("TER\n")
        (root / "native-backbone.pdb").write_text("".join(lines) + "END\n")
        save(root / "sampling.json", {"sampling_completed": True, "native_model_calls": len(trace),
             "configured_steps": plan["steps"], "seed": plan["seed"], "wall_seconds": time.monotonic()-start,
             "target_residues": n, "frame_contract_passed": True, "atomic_contract_passed": False,
             "frame_policy": plan["frame_policy"], "fixed_frame_projection_steps": projection_count,
             "full_atom_completion": False, "public_output_published": False})
    hook.remove()


if __name__ == "__main__":
    directory = Path(sys.argv[1]).resolve()
    try:
        main(directory)
    except Exception as exc:
        import traceback
        save(directory / "native-failure.json", {"error": str(exc), "traceback": traceback.format_exc()})
        raise
