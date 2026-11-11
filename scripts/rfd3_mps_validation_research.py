#!/usr/bin/env python3
"""Frozen reserved-seed counterpart of the explicit MPS development wrapper.

Only admission differs: exact hash-bound validation policy and disjoint seeds.
Old development seed gates and all historical source files remain unchanged.
"""
from __future__ import annotations

import hashlib
import json
import os
import runpy

import torch
from mac_validation_policy import require_policy

PERMISSION = require_policy("rfd3", int(os.environ["DVBFIXER_RFD3_MODEL_SEED"]))

import trace_rfd3_seed_transport as transport  # noqa: E402 -- policy admission precedes model instrumentation


def reseeded(self, input, *args, **kwargs):
    seed = int(os.environ["DVBFIXER_RFD3_MODEL_SEED"])
    require_policy("rfd3", seed)
    before = hashlib.sha256(torch.mps.get_rng_state().numpy().tobytes()).hexdigest()
    torch.mps.manual_seed(seed)
    after = hashlib.sha256(torch.mps.get_rng_state().numpy().tobytes()).hexdigest()
    (transport.ROOT / "seed-amendment.json").write_text(json.dumps({"method": "private-mps-seed-at-model-v1",
        "cohort": "frozen independent-structure validation", "policy": PERMISSION, "seed": seed,
        "before_mps_rng_sha256": before, "after_mps_rng_sha256": after,
        "old_results_unchanged": True, "cpu_rng_reset": False, "weights_masks_sampler_unchanged": True}, indent=2)+"\n")
    return transport.traced(self, input, *args, **kwargs)


if __name__ == "__main__":
    transport.RFD3.forward = reseeded
    runpy.run_module("rfd3.cli", run_name="__main__")
