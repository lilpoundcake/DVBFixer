#!/usr/bin/env python3
"""Private non-mutating RFD3 final-consumer/RNG instrumentation for seed audit."""
from __future__ import annotations

import hashlib
import json
import os
import runpy
from pathlib import Path

import numpy as np
import torch
from rfd3.model.RFD3 import RFD3

ROOT = Path(os.environ["DVBFIXER_RFD3_TRANSPORT_OUTPUT"])
ROOT.mkdir(exist_ok=False)
original = RFD3.forward


def traced(self, input, *args, **kwargs):
    features = input["f"]
    report = {"feature_keys": sorted(features), "model_training": self.training,
              "cpu_rng_sha256": hashlib.sha256(torch.get_rng_state().numpy().tobytes()).hexdigest(),
              "mps_rng_sha256": hashlib.sha256(torch.mps.get_rng_state().numpy().tobytes()).hexdigest()}
    arrays = {}
    for key in ("is_motif_atom_with_fixed_coord", "is_motif_atom_with_fixed_seq", "atom_to_token_map"):
        if key in features:
            arrays[key] = features[key].detach().cpu().numpy()
            report[key] = {"shape": list(arrays[key].shape), "sum": int(arrays[key].sum())}
    expected_fixed = os.environ.get("DVBFIXER_RFD3_EXPECTED_FIXED")
    expected_generated = os.environ.get("DVBFIXER_RFD3_EXPECTED_GENERATED")
    if expected_fixed is not None and expected_generated is not None:
        fixed = arrays["is_motif_atom_with_fixed_coord"]
        if (int(fixed.sum()) != int(expected_fixed) or int((~fixed.astype(bool)).sum()) != int(expected_generated)
                or not arrays["is_motif_atom_with_fixed_seq"].all()):
            (ROOT / "consumer-admission-failure.json").write_text(json.dumps(report, indent=2)+"\n")
            raise RuntimeError("final RFD3 sampler mask differs from declared source/generated ownership")
    calls = 0
    def noisy(module, args, kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            arrays["first_noisy"] = kwargs["X_noisy_L"].detach().cpu().numpy().copy()
            report["first_noisy_sha256"] = hashlib.sha256(arrays["first_noisy"].tobytes()).hexdigest()
    handle = self.diffusion_module.register_forward_pre_hook(noisy, with_kwargs=True)
    try:
        result = original(self, input, *args, **kwargs)
        arrays["final_X_L"] = result["X_L"].detach().cpu().numpy().copy()
        return result
    finally:
        handle.remove()
        report["diffusion_module_calls"] = calls
        (ROOT / "consumer.json").write_text(json.dumps(report, indent=2) + "\n")
        np.savez_compressed(ROOT / "native-transport.npz", **arrays)


if __name__ == "__main__":
    RFD3.forward = traced
    runpy.run_module("rfd3.cli", run_name="__main__")
