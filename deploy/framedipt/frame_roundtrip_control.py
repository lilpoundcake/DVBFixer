"""No model/sampling: replay the native fixed-frame representation conversions."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch
from framedipt.diffusion.se3_diffuser import _assemble_rigid, _extract_trans_rots
from openfold.utils.rigid_utils import Rigid


def main(initial_path, output):
    output.mkdir(parents=True, exist_ok=False)
    arrays = dict(np.load(initial_path, allow_pickle=False))
    fixed = arrays["fixed_mask"].astype(bool)
    reference = torch.tensor(arrays["rigids_0"])
    target = Rigid.from_tensor_7(reference)
    rot0 = target.get_rots().get_rot_mats().numpy()
    trans0 = target.get_trans().numpy()
    current = reference.clone()
    trace = []
    with (output / "plan.json").open("x") as stream:
        json.dump({"protocol": "native-fixed-frame-conversions-only-v1", "iterations": 100,
                   "initial_sha256": hashlib.sha256(initial_path.read_bytes()).hexdigest(),
                   "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "model_calls": 0, "diffusion_calls": 0, "tolerance_changed": False}, stream, indent=2)
    for step in range(100):
        # SE3Diffuser.reverse extracts these representations, applies the fixed
        # mask and reassembles, and its caller stores tensor7 for the next step.
        trans, rot = _extract_trans_rots(Rigid.from_tensor_7(current))
        current = _assemble_rigid(rot, trans).to_tensor_7()
        restored = Rigid.from_tensor_7(current)
        trace.append({"iteration": step+1,
                      "rotation_max_component": float(np.max(np.abs(restored.get_rots().get_rot_mats().numpy()[fixed] - rot0[fixed]))),
                      "translation_max_A": float(np.linalg.norm(restored.get_trans().numpy()[fixed] - trans0[fixed], axis=-1).max())})
    with (output / "result.json").open("x") as stream:
        json.dump({"trace": trace, "first_rotation_crossing_1e5": next((r["iteration"] for r in trace if r["rotation_max_component"] > 1e-5), None),
                   "maximum_rotation_component": max(r["rotation_max_component"] for r in trace),
                   "interpretation": "conversion-only control, not full-sampler equivalence or model quality"}, stream, indent=2)
    print(json.dumps(trace[-1]))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
