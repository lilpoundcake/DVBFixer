"""Private inference-only row tiling of native sinusoidal pair embeddings.

No low_memory_mode shortcut, changed weights, attention neighbourhoods or omitted
pairs. This only bounds the large sin/cos temporary tensors; the complete dense
pair embedding is still returned. Not automatically installed by any runner.
"""
from __future__ import annotations

import math

import torch


def sinusoidal_rows(module, positions, valid_mask, *, rows: int = 64):
    """Evaluate the native expression in query-row tiles, with gradients refused."""
    if torch.is_grad_enabled() or module.training:
        raise ValueError("private row tiling requires eval and no_grad inference")
    if isinstance(rows, bool) or not isinstance(rows, int) or rows <= 0:
        raise ValueError("query tile rows must be a positive integer")
    if positions.ndim not in (2, 3) or positions.shape[-1] != 3:
        raise ValueError("positions must be [L,3] or [B,L,3]")
    length = positions.shape[-2]
    if length == 0 or valid_mask.shape != (*positions.shape[:-2], length, length, 1):
        raise ValueError("mask must have the complete native pair axis")
    if positions.device != valid_mask.device or valid_mask.dtype != torch.bool:
        raise ValueError("mask must be boolean on the positions device")
    # Match native CPU frequency construction then transfer, rather than silently
    # moving exp/arange to MPS and introducing an unrelated numerical change.
    frequency = torch.exp(
        -math.log(10000.0) * torch.arange(0, module.n_freqs, dtype=torch.float32) / module.n_freqs
    ).to(positions.device)
    result = None
    for start in range(0, length, rows):
        stop = min(start + rows, length)
        delta = positions[..., start:stop, None, :] - positions[..., None, :, :]
        distance = torch.linalg.norm(delta, dim=-1)
        angles = distance.unsqueeze(-1) * frequency
        sine = torch.sin(angles)
        cosine = torch.cos(angles)
        encoded = torch.cat([sine, cosine], dim=-1)
        projected = module.output_proj(encoded)
        mask = valid_mask[..., start:stop, :, :]
        projected = projected * mask
        projected = projected + module.process_valid_mask(mask.to(projected.dtype)) * mask
        if result is None:
            result = torch.empty((*positions.shape[:-2], length, length, module.c_atompair),
                                 device=projected.device, dtype=projected.dtype)
        result[..., start:stop, :, :] = projected
        # End each tile's lifetime explicitly; do not retain LxL sin/cos arrays.
        del delta, distance, angles, sine, cosine, encoded, projected
    return result
