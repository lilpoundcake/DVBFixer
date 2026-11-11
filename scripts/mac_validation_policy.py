"""Explicit frozen-policy admission for private reserved-seed execution."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

SEEDS = (211, 223, 227, 229, 233)


def require_policy(model: str, seed: int, policy_path: Path | None = None):
    path = policy_path if policy_path is not None else Path(os.environ["DVBFIXER_VALIDATION_POLICY"])
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != os.environ.get("DVBFIXER_VALIDATION_POLICY_SHA256"):
        raise ValueError("reserved-seed execution requires the exact frozen policy digest")
    policy = json.loads(data)
    if (policy.get("protocol") != "private-fixed-source-validation-v1" or tuple(policy.get("seeds", ())) != SEEDS
            or seed not in SEEDS or model not in policy.get("models", ())
            or policy.get("observed_motion") != "none" or policy.get("automatic_candidate_selection") is not False
            or policy.get("public_capability_enabled") is not False or policy.get("validation_seed_sampling_authorized") is not True):
        raise ValueError("policy does not admit this fixed-source private validation request")
    return {"path": str(path), "sha256": digest, "model": model, "seed": seed, "qualification": False}
