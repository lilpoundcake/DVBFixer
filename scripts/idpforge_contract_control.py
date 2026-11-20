"""Private IDPForge sampler-control probe, NOT checkpoint inference.

Execute only the exact pinned upstream ``recon`` method with synthetic forward
and denoiser consumers. This tests orchestration, not trained-model quality,
runtime compatibility, valid torsion/SS preparation or chemical acceptance.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import subprocess
import typing
from pathlib import Path
from types import SimpleNamespace

import numpy as np

REVISION = "e8da98924e956fd307a6df903432e5a271637076"


def extract_recon(source):
    tree = ast.parse(source)
    models = [n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "IDPForge"]
    methods = [n for model in models for n in model.body if isinstance(n, ast.FunctionDef) and n.name == "recon"]
    if len(models) != 1 or len(methods) != 1:
        raise ValueError("require exactly one upstream IDPForge.recon")
    return methods[0]


def probe(source, torch):
    """Use real tensor indexing and the unchanged upstream orchestration body."""
    function = extract_recon(source)
    namespace = {"torch": torch, "T": typing,
                 "get_chi_mask": lambda aa: torch.ones((*aa.shape, 4), dtype=torch.bool)}
    module = ast.Module(body=[function], type_ignores=[])
    exec(compile(module, "<pinned-idpforge-recon>", "exec"), namespace)
    coords = torch.arange(6 * 5 * 3, dtype=torch.float32).reshape(1, 6, 5, 3)
    torsions = torch.zeros(1, 6, 4, 2)
    fixed = torch.tensor([[False, False, True, True, True, True]])
    template = {"mask": fixed, "coord": coords.clone(), "torsion": torsions.clone()}
    batch = {"sequence": torch.zeros(1, 6, dtype=torch.long), "mask": torch.ones(1, 6)}
    calls = []

    def forward(tstep, alpha, xyz, consumed, **kwargs):
        calls.append({"t": int(tstep[0, 0]),
                      "fixed_times_zero": bool(torch.all(tstep[fixed] == 0)),
                      "template_coordinates_exact": torch.equal(xyz[fixed], coords[fixed]),
                      "template_torsions_exact": torch.equal(alpha[fixed], torsions[fixed]),
                      "all_residues_retained": bool(torch.all(consumed["mask"] == 1))})
        # Deliberately move ALL predicted atoms. This is a surrogate, not a
        # claim that the trained model produces this motion or this geometry.
        return {"positions": (xyz + .25)[None], "angles": torch.zeros(1, 1, 6, 7, 2)}

    reverse_calls = []

    def next_pose(xyz, predicted, alpha, predicted_alpha, t, chi_mask, aa_mask, motiff_mask):
        reverse_calls.append({"t": t, "motif_mask_exact": np.array_equal(motiff_mask, fixed.numpy())})
        # Perturb even the motif to prove the next forward receives a freshly
        # restored template rather than an accidentally unchanged mock state.
        return xyz + 10., alpha + 10.

    surrogate = SimpleNamespace(n_tsteps=200, inf_tsteps=40, end_tsteps=0,
                                device=torch.device("cpu"), cfg=SimpleNamespace(self_condition=True),
                                forward=forward)
    output = namespace["recon"](surrogate, SimpleNamespace(get_next_pose=next_pose),
                                 torch.zeros_like(torsions), torch.zeros_like(coords), batch,
                                 template_cfg=template)
    returned = output["positions"][-1]
    shift = returned[fixed] - coords[fixed]
    passed = (len(calls) == 40 and len(reverse_calls) == 40
              and all(all(row[k] for k in ("fixed_times_zero", "template_coordinates_exact",
                      "template_torsions_exact", "all_residues_retained")) for row in calls)
              and all(row["motif_mask_exact"] for row in reverse_calls)
              and bool(torch.all(shift == .25)))
    if not passed:
        raise ValueError("pinned recon orchestration differs from the declared control")
    return {"protocol": "private-idpforge-synthetic-recon-control-v1", "control_passed": True,
            "checkpoint_loaded": False, "native_model_calls": 0, "diffusion_samples": 0,
            "surrogate_forward_calls": len(calls), "surrogate_reverse_calls": len(reverse_calls),
            "first_t": calls[0]["t"], "last_t": calls[-1]["t"],
            "template_restored_before_every_forward": True,
            "returned_prediction_restored_to_template": False,
            "deliberate_surrogate_shift_per_coordinate_A": .25,
            "actual_checkpoint_observed_drift_A": None,
            "chemistry_or_model_quality_claimed": False, "calls": calls,
            "note": "synthetic torsions/coordinates; not a scientific proposal or model preflight"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.upstream.resolve()
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root, text=True).strip()
    if revision != REVISION or dirty:
        raise ValueError("require clean pinned upstream")
    source = (root / "idpforge/model.py").read_bytes()
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "upstream-model.py").write_bytes(source)
    (args.output / "executed-control.py").write_bytes(Path(__file__).read_bytes())
    plan = {"upstream_revision": revision, "model_sha256": hashlib.sha256(source).hexdigest(),
            "control_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "checkpoint_inference_requested": False, "synthetic_control_only": True}
    (args.output / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    import torch

    result = probe(source.decode(), torch)
    result["torch"] = torch.__version__
    (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "calls"}))


if __name__ == "__main__":
    main()
