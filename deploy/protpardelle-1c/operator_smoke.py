#!/usr/bin/env python3
"""Run a digest-pinned Protpardelle cc89 CPU or MPS operator smoke."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
import time
from contextlib import nullcontext
from pathlib import Path
from typing import Any

_CHECKPOINT_SHA256 = "dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483"
_CONFIG_SHA256 = "e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _conditional_config() -> dict[str, Any]:
    return {
        "enabled": False,
        "discontiguous_motif_assignment": {
            "enabled": False,
            "strategy": "fixed",
            "fixed_motif_pos": [],
        },
        "num_recurrence_steps": 1,
        "crop_conditional_guidance": {
            "enabled": False,
            "start": 0.0,
            "end": 2.0,
            "freq": 1,
            "freq_start": 0.0,
            "freq_end": 0.0,
            "strategy": "backbone-sidechain",
        },
        "reconstruction_guidance": {
            "enabled": False,
            "start": 0.0,
            "end": 2.0,
            "schedule": "custom",
            "max_scale": 10.0,
            "loss_weights": {"motif": 1.0},
        },
        "replacement_guidance": {"enabled": False, "start": 0.0, "end": 0.92},
    }


def _synchronize(torch: Any, device: Any) -> None:
    if device.type == "mps":
        torch.mps.synchronize()


def run(
    config: Path,
    checkpoint: Path,
    output: Path,
    *,
    device_name: str,
    length: int,
    steps: int,
    seed: int,
    mps_profile: bool,
) -> dict[str, Any]:
    import numpy as np
    import torch
    from protpardelle.common import residue_constants
    from protpardelle.core.models import load_model

    coordinate_output = output.with_name(f"{output.stem}.coordinates.npy")
    if output.exists() or coordinate_output.exists():
        raise FileExistsError(f"operator-smoke output already exists: {output}")
    if _sha256(config) != _CONFIG_SHA256 or _sha256(checkpoint) != _CHECKPOINT_SHA256:
        raise ValueError("cc89 config or checkpoint digest mismatch")
    if length <= 0 or length > 512 or steps <= 0:
        raise ValueError("length must be 1-512 and steps must be positive")
    device = torch.device(device_name)
    if device.type == "mps":
        if not torch.backends.mps.is_built() or not torch.backends.mps.is_available():
            raise RuntimeError("MPS is not built and available")
        if os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK") != "0":
            raise RuntimeError("MPS smoke requires PYTORCH_ENABLE_MPS_FALLBACK=0")
    if mps_profile and device.type != "mps":
        raise ValueError("--mps-profile requires --device mps")

    torch.manual_seed(seed)
    if device.type == "mps":
        torch.mps.manual_seed(seed)
    load_start = time.perf_counter()
    model = load_model(config, checkpoint, device=device)
    _synchronize(torch, device)
    load_seconds = time.perf_counter() - load_start
    wrong_state = [
        name
        for name, value in (*model.named_parameters(), *model.named_buffers())
        if value.device.type != device.type
    ]
    if wrong_state:
        raise RuntimeError(f"model state is on the wrong device: {wrong_state[:5]}")

    seq_mask, residue_index, chain_index = model.make_seq_mask_for_sampling(
        prot_lens_per_chain=torch.tensor([[length]])
    )
    seq_mask = seq_mask.to(device)
    residue_index = residue_index.to(device)
    chain_index = chain_index.to(device)
    gt_aatype = torch.full(
        (1, length),
        residue_constants.restype_order["A"],
        dtype=torch.long,
        device=device,
    )
    _synchronize(torch, device)
    sample_start = time.perf_counter()
    profile = torch.mps.profiler.profile() if mps_profile else nullcontext()
    with torch.no_grad(), profile:
        result = model.sample(
            seq_mask=seq_mask,
            residue_index=residue_index,
            chain_index=chain_index,
            gt_aatype=gt_aatype,
            num_steps=steps,
            step_scale=1.2,
            s_churn=200.0,
            partial_diffusion={"enabled": False, "pdb_file_path": None, "num_steps": 100},
            conditional_cfg=_conditional_config(),
            sidechain_mode=False,
            skip_mpnn_proportion=1.0,
            use_fullmpnn=False,
            use_fullmpnn_for_final=False,
            jump_steps=False,
            uniform_steps=True,
            **({"record_trajectory": False} if device.type == "mps" else {}),
        )
    coordinates = result["x"]
    if coordinates.device.type != device.type:
        raise RuntimeError(f"sample output moved to {coordinates.device}")
    _synchronize(torch, device)
    sample_seconds = time.perf_counter() - sample_start
    expected_shape = (1, length, 37, 3)
    if tuple(coordinates.shape) != expected_shape or not bool(torch.isfinite(coordinates).all()):
        raise RuntimeError("operator smoke returned an invalid coordinate tensor")
    coordinate_array = coordinates.detach().to(device="cpu", dtype=torch.float32).numpy()
    coordinate_bytes = coordinate_array.tobytes()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb",
        dir=output.parent,
        prefix=f".{coordinate_output.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        coordinate_temporary = Path(handle.name)
        np.save(handle, coordinate_array, allow_pickle=False)
    coordinate_temporary.replace(coordinate_output)
    payload = {
        "schema_version": 1,
        "status": "passed",
        "device": str(device),
        "precision": str(coordinates.dtype),
        "length": length,
        "steps": steps,
        "seed": seed,
        "shape": list(expected_shape),
        "finite": True,
        "coordinate_sha256": hashlib.sha256(coordinate_bytes).hexdigest(),
        "coordinate_artifact": coordinate_output.name,
        "coordinate_artifact_sha256": _sha256(coordinate_output),
        "config_sha256": _CONFIG_SHA256,
        "checkpoint_sha256": _CHECKPOINT_SHA256,
        "torch_version": torch.__version__,
        "model_load_seconds": load_seconds,
        "sample_seconds": sample_seconds,
        "mps_profile_enabled": mps_profile,
        "mps_current_allocated_bytes": (
            int(torch.mps.current_allocated_memory()) if device.type == "mps" else None
        ),
        "mps_driver_allocated_bytes": (
            int(torch.mps.driver_allocated_memory()) if device.type == "mps" else None
        ),
        "mps_recommended_max_memory_bytes": (
            int(torch.mps.recommended_max_memory()) if device.type == "mps" else None
        ),
    }
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=output.parent,
        prefix=f".{output.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(output)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--device", choices=("cpu", "mps"), default="mps")
    parser.add_argument("--length", default=32, type=int)
    parser.add_argument("--steps", default=1, type=int)
    parser.add_argument("--seed", default=20260929, type=int)
    parser.add_argument("--mps-profile", action="store_true")
    args = parser.parse_args()
    payload = run(
        args.config,
        args.checkpoint,
        args.output,
        device_name=args.device,
        length=args.length,
        steps=args.steps,
        seed=args.seed,
        mps_profile=args.mps_profile,
    )
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
