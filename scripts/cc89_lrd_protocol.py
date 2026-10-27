"""Pure helpers for the research-only cc89 local-denoising ablation."""

from __future__ import annotations

import numpy as np

from dvbfixer.model.diffusion.contract import DiffusionRequest, GapRegion, ResidueIdentity


def relaxation_residues(request: DiffusionRequest, *, include_anchors: bool = False) -> tuple[ResidueIdentity, ...]:
    """Explicit research minimization mask; never change the request's generated mask."""
    residues = []
    for gap in request.gaps:
        residues.extend(gap.generated_residues)
        if include_anchors:
            residues.extend(anchor for anchor in (gap.left_anchor, gap.right_anchor) if anchor is not None)
    return tuple(dict.fromkeys(residues))


def boundary_residues(
    request: DiffusionRequest, gap: GapRegion, *, include_observed: bool, window: int = 3,
) -> frozenset[ResidueIdentity]:
    """Use target positions, never author-number arithmetic or neighboring chains."""
    if window < 1:
        raise ValueError("boundary window must be positive")
    placement = next(item for item in request.sequence_placements if item.chain == gap.chain)
    positions = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
    generated = dict(zip(range(gap.target_interval.start, gap.target_interval.stop),
                         gap.generated_residues, strict=True))
    positions.update(generated)
    edges = []
    if gap.left_anchor is not None:
        edges.append(gap.target_interval.start)
    if gap.right_anchor is not None:
        edges.append(gap.target_interval.stop)
    if not edges:
        raise ValueError("LRD requires an observed boundary anchor")
    indices = {index for edge in edges for index in range(edge - window, edge + window)}
    return frozenset(residue for index, residue in positions.items()
                     if index in indices and (include_observed or index in generated))


def low_noise_schedule(*, steps: int = 25, sigma_start: float = 1.5,
                       sigma_min: float = 0.0064, rho: float = 7.0) -> np.ndarray:
    if (steps < 2 or not all(np.isfinite(x) and x > 0 for x in (sigma_start, sigma_min, rho))
            or sigma_start <= sigma_min):
        raise ValueError("invalid low-noise schedule")
    ramp = np.linspace(0, 1, steps)
    values = (sigma_start ** (1 / rho)
              + ramp * (sigma_min ** (1 / rho) - sigma_start ** (1 / rho))) ** rho
    return np.concatenate((values, [0.0]))


def project_static(state: np.ndarray, reference: np.ndarray, mobile: np.ndarray) -> np.ndarray:
    """Research test oracle: preserve every inactive coordinate exactly."""
    if state.shape != reference.shape or state.shape[:-1] != mobile.shape or state.shape[-1] != 3:
        raise ValueError("coordinate and mobility axes do not match")
    if mobile.dtype != np.bool_ or not np.isfinite(state).all() or not np.isfinite(reference).all():
        raise ValueError("invalid coordinates or mobility mask")
    return np.where(mobile[..., None], state, reference)


def observed_frame_center(coordinates: np.ndarray, observed: np.ndarray) -> np.ndarray:
    """One source-defined frame for BOTH mobility policies; never center by mobility."""
    if coordinates.shape[:-1] != observed.shape or coordinates.shape[-1] != 3:
        raise ValueError("observed mask does not match coordinates")
    if observed.dtype != np.bool_ or np.count_nonzero(observed) < 3:
        raise ValueError("at least three observed atoms define the frame")
    points = coordinates[observed]
    if not np.isfinite(points).all():
        raise ValueError("observed coordinates must be finite")
    return points.mean(axis=0)
