"""Additional private exact-expression tiling; no runner registration."""
from __future__ import annotations

import torch


def _rows(rows):
    if isinstance(rows, bool) or not isinstance(rows, int) or rows <= 0:
        raise ValueError("rows must be a positive integer")
    if torch.is_grad_enabled():
        raise ValueError("private tiling requires no_grad inference")


def reference_rows(module, positions, valid_mask, *, rows=64, accumulate_into=None):
    """Native non-frame inverse-distance embedding, over every query/key pair."""
    _rows(rows)
    if module.training or module.embed_frame:
        raise ValueError("only eval non-frame reference embeddings are supported")
    if positions.ndim not in (2, 3) or positions.shape[-1] != 3:
        raise ValueError("positions must be [L,3] or [B,L,3]")
    length = positions.shape[-2]
    if length == 0 or valid_mask.shape != (*positions.shape[:-2], length, length, 1):
        raise ValueError("mask must cover the full pair axis")
    if valid_mask.device != positions.device or valid_mask.dtype != torch.bool:
        raise ValueError("mask must be boolean and on the positions device")
    output = accumulate_into
    if output is not None and (output.shape != (*positions.shape[:-2], length, length, module.process_inverse_dist.out_features)
                               or output.device != positions.device or output.dtype != positions.dtype
                               or not output.is_contiguous() or output.requires_grad):
        raise ValueError("accumulator must be an owned contiguous inference pair tensor")
    for start in range(0, length, rows):
        stop = min(length, start + rows)
        delta = positions[..., start:stop, None, :] - positions[..., None, :, :]
        norm = torch.linalg.norm(delta, dim=-1, keepdim=True) ** 2
        norm = torch.clamp(norm, min=1e-6)
        inverse = 1 / (1 + norm)
        mask = valid_mask[..., start:stop, :, :]
        result = module.process_inverse_dist(inverse) * mask
        result = result + module.process_valid_mask(mask.to(result.dtype)) * mask
        if output is None:
            output = torch.empty((*positions.shape[:-2], length, length, result.shape[-1]),
                                 dtype=result.dtype, device=result.device)
        if accumulate_into is None:
            output[..., start:stop, :, :] = result
        else:
            output[..., start:stop, :, :].add_(result)
        del delta, norm, inverse, result
    return output


def pointwise_rows(function, values, *, rows=64, training=False):
    """Run an explicitly supplied pointwise MLP on query-row tiles, no reductions.

    The caller must establish pointwise semantics. This helper is not valid for
    attention, pooling, dropout, cross-row normalization or other coupled layers.
    """
    _rows(rows)
    if training:
        raise ValueError("private tiling requires eval inference")
    if values.ndim not in (3, 4) or values.shape[-3] == 0 or values.shape[-2] == 0:
        raise ValueError("pointwise input must be [L,M,C] or [B,L,M,C]")
    length = values.shape[-3]
    output = None
    for start in range(0, length, rows):
        stop = min(length, start + rows)
        result = function(values[..., start:stop, :, :])
        if result.shape[:-1] != values[..., start:stop, :, :].shape[:-1]:
            raise ValueError("pointwise function changed pair axes")
        if output is None:
            output = torch.empty((*values.shape[:-3], length, values.shape[-2], result.shape[-1]),
                                 dtype=result.dtype, device=result.device)
        output[..., start:stop, :, :] = result
        del result
    return output
