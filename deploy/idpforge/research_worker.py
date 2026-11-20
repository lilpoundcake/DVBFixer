"""Private CPU native IDPForge sampler; retain raw atom14 before any restoration."""
from __future__ import annotations

import importlib.util
import json
import random
import sys
import time
from pathlib import Path


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def source_secondary_structure(md, assign_rama, source, gap_length):
    """Keep MDTraj angles in radians until the native classifier consumes them."""
    import numpy as np

    traj = md.load(str(source))
    dssp = md.compute_dssp(traj, simplified=True)[0]
    phi = np.r_[-np.pi, md.compute_phi(traj)[1][0]]
    psi = np.r_[md.compute_psi(traj)[1][0], np.pi]
    rama = assign_rama(np.stack((phi, psi), axis=-1))
    return "C" * gap_length + "".join(d if d in ("H", "E") else r for d, r in zip(dssp, rama, strict=True))


def main(root):
    import ml_collections
    import numpy as np
    import torch
    import yaml
    from idpforge.model import IDPForge
    from idpforge.utils.diff_utils import Denoiser, Diffuser
    from openfold.np import residue_constants as rc

    plan = json.loads((root / "plan.json").read_text())
    transport = json.loads((root / "transport.json").read_text())
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    np.random.seed(plan["seed"])
    random.seed(plan["seed"])
    torch.manual_seed(plan["seed"])
    upstream = Path(plan["upstream"])
    spec = importlib.util.spec_from_file_location("native_mk_template", upstream / "mk_ldr_template.py")
    prep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prep)
    template_data = prep.main(str(root / "source.pdb"), [], 1,
                              prepend_seq=transport["sequence"][:transport["gap_length"]])
    if plan.get("template_ss_policy", "stock") == "radian-consistent":
        import mdtraj as md
        from idpforge.utils.np_utils import assign_rama

        corrected = source_secondary_structure(md, assign_rama, root / "source.pdb", transport["gap_length"])
        save(root / "ss-preparation-control.json", {"stock": str(template_data["sec"]), "corrected": corrected,
             "policy": "radian-consistent", "source_coordinates_changed": False,
             "sampling_algorithm_changed": False})
        template_data["sec"] = corrected
    np.savez_compressed(root / "template.npz", **template_data)
    if (str(template_data["seq"]) != transport["sequence"]
            or template_data["mask"].tolist() != [not row["generated"] for row in transport["residues"]]):
        raise ValueError("native preparation changed sequence or mask")
    settings = yaml.safe_load((upstream / "configs/sample.yml").read_text())
    model = IDPForge(settings["diffuse"]["n_tsteps"], settings["diffuse"]["inference_steps"],
                     ml_collections.ConfigDict(settings["model"])).cpu().eval()
    checkpoint = torch.load(plan["weights"], map_location="cpu", weights_only=True)
    old = {f"trunk.structure_module.ipa.linear_{name}_points.{kind}"
           for name in ("q", "kv") for kind in ("weight", "bias")}
    # Same explicit OpenFold-2 key migration as official sample_ldr.py.
    sd = {k.replace("points.", "points.linear.") if k in old else k: v for k, v in checkpoint["ema"]["params"].items()}
    model.load_state_dict(sd, strict=True)
    save(root / "runtime.json", {"torch": torch.__version__, "numpy": np.__version__, "python": sys.version,
         "parameters": sum(p.numel() for p in model.parameters()), "device": "CPU", "strict_weights_loaded": True,
         "key_migration": sorted(old), "runtime_prefix": sys.prefix, "borrowed_base_prefix": sys.base_prefix})
    cfg = settings["diffuse"]
    diffuser = Diffuser(cfg["n_tsteps"], euclid_b0=cfg["euclid_b0"], euclid_bT=cfg["euclid_bT"],
                        tor_b0=cfg["torsion_b0"], tor_bT=cfg["torsion_bT"])
    denoiser = Denoiser(cfg["inference_steps"], diffuser)
    xt, alpha = denoiser.init_samples([transport["sequence"]])
    # Unchanged stock translation policy, with the single prepared offset.
    offset = template_data["coord_offset"][0]
    template = {"mask": torch.tensor(template_data["mask"][None], dtype=torch.long),
                "torsion": torch.tensor(template_data["torsion"][None], dtype=torch.float32),
                "coord": torch.tensor(template_data["coord"][None] - offset[None, None, None], dtype=torch.float32)}
    ss = str(template_data["sec"])
    n = len(transport["sequence"])
    expected_aa = torch.tensor([[rc.restype_order[a] for a in transport["sequence"]]])
    fixed = template["mask"].bool()
    np.savez_compressed(root / "native-initial.npz", xt=xt[0].numpy(), alpha=alpha[0].numpy(),
                        template_coord=template["coord"].numpy(), template_torsion=template["torsion"].numpy(),
                        template_mask=template["mask"].numpy())
    if plan["mode"] == "sample":
        previous = dict(np.load(Path(plan["preflight"]) / "native-initial.npz", allow_pickle=False))
        current = dict(np.load(root / "native-initial.npz", allow_pickle=False))
        if set(current) != set(previous) or not all(np.array_equal(current[k], previous[k]) for k in current):
            raise ValueError("sample initialization differs from its preflight")
        save(root / "paired-initial-state.json", {"all_arrays_identical": True})
    trace = []
    (root / "consumer-states").mkdir()
    native_forward = model.forward

    def consumer(t, tor, xyz, batch, **kwargs):
        if (not torch.equal(batch["sequence"], expected_aa) or not torch.all(batch["mask"] == 1)
                or batch["sequence"].shape != (1, n) or not torch.all(t[fixed] == 0)
                or not torch.equal(xyz[fixed], template["coord"][fixed][:, :5, :])
                or not torch.equal(tor[fixed], template["torsion"][fixed])):
            raise ValueError("native model consumer changed target/observation conditioning")
        call = len(trace) + 1
        np.savez_compressed(root / "consumer-states" / f"{call:03d}.npz", t=t.numpy(), xyz=xyz.numpy(),
                            tor=tor.numpy(), **{k: v.numpy() for k, v in batch.items()})
        row = {"call": call, "t_generated": int(t[0, 0]), "target_residues": n,
               "device": str(xyz.device), "fixed_input_exact_to_float32_template": True}
        trace.append(row)
        with (root / "consumer-trace.jsonl").open("a") as stream:
            stream.write(json.dumps(row) + "\n")
        output = native_forward(t, tor, xyz, batch, **kwargs)
        if not torch.isfinite(output["positions"]).all():
            raise ValueError("nonfinite native atom14")
        if plan["mode"] == "preflight":
            np.savez_compressed(root / "preflight-forward.npz", positions=output["positions"][-1].numpy())
            raise PreflightComplete
        return output

    model.forward = consumer
    started = time.monotonic()
    try:
        output = model.sample(denoiser, [transport["sequence"]], [ss], alpha, xt, template_cfgs=template)
    except PreflightComplete:
        save(root / "preflight.json", {"consumer_passed": True, "model_calls": len(trace),
             "sampling_completed": False, "strict_weights_loaded": True})
        return
    if plan["mode"] != "sample" or len(trace) != 40:
        raise ValueError("unexpected native sampling count")
    np.savez_compressed(root / "native-endpoint.npz", **{k: v.numpy() for k, v in output.items()})
    save(root / "sampling.json", {"sampling_completed": True, "model_calls": len(trace),
         "seed": plan["seed"], "wall_seconds": time.monotonic()-started, "source_restored": False,
         "relaxed": False, "repaired": False, "public_output_published": False})


class PreflightComplete(Exception):
    """Stop after one actual forward, before any reverse diffusion."""


if __name__ == "__main__":
    directory = Path(sys.argv[1]).resolve()
    try:
        main(directory)
    except Exception as exc:
        import traceback
        save(directory / "native-failure.json", {"error": str(exc), "traceback": traceback.format_exc()})
        raise
