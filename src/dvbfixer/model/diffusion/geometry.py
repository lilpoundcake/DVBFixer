"""Coordinate synchronization primitives for constrained diffusion sampling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from dvbfixer.model.diffusion.contract import AtomIdentity

FloatArray = NDArray[np.float64]


class DiffusionGeometryError(ValueError):
    """Raised when coordinate synchronization inputs are invalid."""


@dataclass(frozen=True, slots=True)
class RigidTransform:
    rotation: FloatArray
    translation: FloatArray

    def apply(self, coordinates: FloatArray) -> FloatArray:
        points = _coordinates(coordinates, "coordinates")
        return points @ self.rotation.T + self.translation


def _coordinates(value: FloatArray, name: str) -> FloatArray:
    array = np.asarray(value, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 3:
        raise DiffusionGeometryError(f"{name} must have shape (N, 3)")
    if not np.isfinite(array).all():
        raise DiffusionGeometryError(f"{name} must contain only finite coordinates")
    return array


def weighted_kabsch(
    mobile: FloatArray,
    target: FloatArray,
    weights: FloatArray | None = None,
) -> RigidTransform:
    """Return the proper rigid transform that best maps mobile onto target."""
    mobile_points = _coordinates(mobile, "mobile")
    target_points = _coordinates(target, "target")
    if mobile_points.shape != target_points.shape:
        raise DiffusionGeometryError("mobile and target must have identical shape")
    if mobile_points.shape[0] < 3:
        raise DiffusionGeometryError("weighted Kabsch requires at least three points")

    if weights is None:
        fit_weights = np.ones(mobile_points.shape[0], dtype=np.float64)
    else:
        fit_weights = np.asarray(weights, dtype=np.float64)
        if fit_weights.ndim != 1 or fit_weights.shape[0] != mobile_points.shape[0]:
            raise DiffusionGeometryError("weights must have shape (N,)")
        if not np.isfinite(fit_weights).all() or np.any(fit_weights <= 0.0):
            raise DiffusionGeometryError("weights must be finite and strictly positive")

    weight_sum = float(fit_weights.sum())
    mobile_center = np.sum(mobile_points * fit_weights[:, None], axis=0) / weight_sum
    target_center = np.sum(target_points * fit_weights[:, None], axis=0) / weight_sum
    mobile_centered = mobile_points - mobile_center
    target_centered = target_points - target_center
    covariance = (mobile_centered * fit_weights[:, None]).T @ target_centered
    left, _singular_values, right_transpose = np.linalg.svd(covariance)
    determinant = float(np.linalg.det(right_transpose.T @ left.T))
    correction = np.diag([1.0, 1.0, 1.0 if determinant > 0.0 else -1.0])
    rotation = right_transpose.T @ correction @ left.T
    translation = target_center - mobile_center @ rotation.T
    return RigidTransform(rotation=rotation, translation=translation)


def synchronize_and_reinject(
    atom_identities: tuple[AtomIdentity, ...],
    coordinates: FloatArray,
    fixed_coordinates: dict[AtomIdentity, FloatArray],
    *,
    anchor_weights: dict[AtomIdentity, float] | None = None,
) -> FloatArray:
    """Synchronize a sampled frame and overwrite fixed atoms exactly.

    All fixed identities must occur exactly once in ``atom_identities``. The
    weighted Kabsch fit uses their sampled and authoritative coordinates, then
    the authoritative coordinates are reinserted byte-for-byte into the returned
    array. Non-fixed coordinates receive only the fitted rigid transform.
    """
    sampled = _coordinates(coordinates, "coordinates")
    if len(atom_identities) != sampled.shape[0]:
        raise DiffusionGeometryError("atom identity count must match coordinate count")
    if len(set(atom_identities)) != len(atom_identities):
        raise DiffusionGeometryError("atom identities must be unique")
    if len(fixed_coordinates) < 3:
        raise DiffusionGeometryError("at least three fixed anchors are required")

    index_by_identity = {identity: index for index, identity in enumerate(atom_identities)}
    missing = set(fixed_coordinates) - set(index_by_identity)
    if missing:
        raise DiffusionGeometryError("fixed coordinate identity is absent from sampled atoms")
    if anchor_weights is not None and set(anchor_weights) != set(fixed_coordinates):
        raise DiffusionGeometryError("anchor_weights must cover exactly the fixed identities")

    fixed_identities = tuple(sorted(fixed_coordinates))
    fixed_indices = [index_by_identity[identity] for identity in fixed_identities]
    mobile_anchors = sampled[fixed_indices]
    target_anchors = np.vstack(
        [_coordinate(fixed_coordinates[identity], "fixed coordinate") for identity in fixed_identities]
    )
    weights = None
    if anchor_weights is not None:
        weights = np.asarray(
            [anchor_weights[identity] for identity in fixed_identities],
            dtype=np.float64,
        )
    transform = weighted_kabsch(mobile_anchors, target_anchors, weights)
    synchronized = transform.apply(sampled)
    for identity, target_coordinate in fixed_coordinates.items():
        synchronized[index_by_identity[identity]] = _coordinate(
            target_coordinate, "fixed coordinate"
        )
    return synchronized


def _coordinate(value: FloatArray, name: str) -> FloatArray:
    coordinate = np.asarray(value, dtype=np.float64)
    if coordinate.shape != (3,):
        raise DiffusionGeometryError(f"{name} must have shape (3,)")
    if not np.isfinite(coordinate).all():
        raise DiffusionGeometryError(f"{name} must contain only finite values")
    return coordinate
