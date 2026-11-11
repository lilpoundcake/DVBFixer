"""Private padded-axis MPS frame fit with explicit CPU 3x3 covariance SVD."""
from __future__ import annotations

import numpy as np
import torch


class MPSFixedAtomCallback:
    def __init__(self, axis, padded_count, fixed):
        if len(set(axis)) != len(axis) or padded_count < len(axis) or len(fixed) < 3:
            raise ValueError('invalid real/padded atom axis or fixed anchors')
        lookup = {atom: i for i, atom in enumerate(axis)}
        if not set(fixed) <= lookup.keys():
            raise ValueError('fixed atom missing from real atom axis')
        ordered = tuple(sorted(fixed))
        self.indices = tuple(lookup[a] for a in ordered)
        self.target = np.asarray([fixed[a] for a in ordered], dtype=np.float64)
        if self.target.shape != (len(fixed), 3) or not np.isfinite(self.target).all():
            raise ValueError('invalid fixed coordinates')
        if np.linalg.matrix_rank(self.target - self.target.mean(axis=0)) < 2:
            raise ValueError('fixed frame anchors must not be collinear')
        self.padded_count = padded_count
        self.steps = []
        self.last_pre_overwrite = None

    def __call__(self, x, index, count):
        if x.device.type != 'mps' or x.dtype != torch.float32:
            raise ValueError('private hook requires MPS float32')
        if x.ndim != 3 or x.shape[-2:] != (self.padded_count, 3) or not 0 <= index < count:
            raise ValueError('invalid padded coordinates or sampling step')
        if index != len(self.steps) or not torch.isfinite(x).all():
            raise ValueError('nonsequential callback or nonfinite coordinates')
        indices = torch.as_tensor(self.indices, device=x.device)
        target = torch.as_tensor(self.target, device=x.device, dtype=x.dtype)
        mobile = x.index_select(-2, indices)
        center, target_center = mobile.mean(dim=-2), target.mean(dim=-2)
        covariance = torch.einsum('...ki,kj->...ij', mobile-center[..., None, :], target-target_center)
        left, _, right = np.linalg.svd(covariance.detach().cpu().numpy().astype(np.float64))
        correction = np.broadcast_to(np.eye(3), left.shape).copy()
        correction[..., 2, 2] = np.where(np.linalg.det(right.swapaxes(-2, -1) @ left.swapaxes(-2, -1)) > 0, 1, -1)
        rotation = torch.as_tensor(right.swapaxes(-2, -1) @ correction @ left.swapaxes(-2, -1), device=x.device, dtype=x.dtype)
        output = (x-center[..., None, :]) @ rotation.transpose(-2, -1) + target_center
        self.last_pre_overwrite = output.detach().cpu().numpy().copy()
        pre = float(torch.linalg.vector_norm(output[..., indices, :]-target, dim=-1).max().item())
        output[..., indices, :] = target
        post = float((output[..., indices, :]-target).abs().max().item())
        self.steps.append({'step_index': index, 'step_count': count,
                           'pre_overwrite_max_angstrom': pre, 'post_overwrite_max_angstrom': post})
        return output
