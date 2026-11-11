#!/usr/bin/env python3
"""Separate explicit MPS seed-transport research arm; old cohorts stay frozen.

The first six requested seeds reached different CPU states but identical MPS RNG
and first-noisy coordinates. Set only MPS RNG immediately before the outer model
forward, after checkpoint/preprocessing, and retain the full consumer audit.
This is a new seed contract, not evidence for the untouched upstream route.
"""
from __future__ import annotations

import hashlib
import json
import os
import runpy

import torch
import trace_rfd3_seed_transport as transport


def reseeded(self, input, *args, **kwargs):
    seed = int(os.environ["DVBFIXER_RFD3_MODEL_SEED"])
    if seed not in (7, 42, 101, 102, 103, 104):
        raise ValueError("only declared development seeds allowed")
    before = hashlib.sha256(torch.mps.get_rng_state().numpy().tobytes()).hexdigest()
    torch.mps.manual_seed(seed)
    after = hashlib.sha256(torch.mps.get_rng_state().numpy().tobytes()).hexdigest()
    (transport.ROOT / "seed-amendment.json").write_text(json.dumps({"method": "private-mps-seed-at-model-v1",
        "seed": seed, "before_mps_rng_sha256": before, "after_mps_rng_sha256": after,
        "old_results_unchanged": True, "cpu_rng_reset": False, "weights_masks_sampler_unchanged": True}, indent=2) + "\n")
    return transport.traced(self, input, *args, **kwargs)


if __name__ == "__main__":
    transport.RFD3.forward = reseeded
    runpy.run_module("rfd3.cli", run_name="__main__")
