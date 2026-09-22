"""Tests for weighted Kabsch synchronization and fixed-coordinate reinjection."""

from __future__ import annotations

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.geometry import (
    DiffusionGeometryError,
    synchronize_and_reinject,
    weighted_kabsch,
)


def _rotation_z(angle: float) -> np.ndarray:
    cosine = np.cos(angle)
    sine = np.sin(angle)
    return np.asarray(
        [
            [cosine, -sine, 0.0],
            [sine, cosine, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )


def test_weighted_kabsch_recovers_known_transform() -> None:
    mobile = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ]
    )
    rotation = _rotation_z(np.pi / 3.0)
    translation = np.asarray([3.0, -2.0, 5.0])
    target = mobile @ rotation.T + translation

    transform = weighted_kabsch(mobile, target, np.asarray([5.0, 2.0, 1.0, 1.0]))

    assert np.allclose(transform.apply(mobile), target, atol=1e-12)
    assert np.allclose(transform.rotation.T @ transform.rotation, np.eye(3), atol=1e-12)
    assert np.linalg.det(transform.rotation) == pytest.approx(1.0)


def test_reinjection_restores_fixed_atoms_exactly() -> None:
    identities = tuple(
        AtomIdentity("D", str(number), "", "CA")
        for number in range(1, 6)
    )
    authoritative = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, 0.5, 0.5],
        ]
    )
    rotation = _rotation_z(-np.pi / 4.0)
    translation = np.asarray([-4.0, 3.0, 2.0])
    sampled = authoritative @ rotation.T + translation
    fixed = {
        identities[index]: authoritative[index].copy()
        for index in (0, 1, 2, 3)
    }

    synchronized = synchronize_and_reinject(identities, sampled, fixed)

    for identity, coordinate in fixed.items():
        index = identities.index(identity)
        assert np.array_equal(synchronized[index], coordinate)
    assert np.allclose(synchronized[4], authoritative[4], atol=1e-12)


def test_weighted_fit_prioritizes_reliable_anchor() -> None:
    mobile = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    target = mobile.copy()
    target[3] += np.asarray([2.0, 0.0, 0.0])

    unweighted = weighted_kabsch(mobile, target).apply(mobile)
    weighted = weighted_kabsch(
        mobile,
        target,
        np.asarray([100.0, 100.0, 100.0, 0.1]),
    ).apply(mobile)

    unweighted_reliable_error = np.linalg.norm(unweighted[:3] - target[:3])
    weighted_reliable_error = np.linalg.norm(weighted[:3] - target[:3])
    assert weighted_reliable_error < unweighted_reliable_error


def test_reinjection_preserves_case_sensitive_identity() -> None:
    upper = AtomIdentity("D", "1", "", "CA")
    lower = AtomIdentity("d", "1", "", "CA")
    third = AtomIdentity("D", "2", "A", "CA")
    identities = (upper, lower, third)
    coordinates = np.asarray(
        [
            [10.0, 0.0, 0.0],
            [10.0, 1.0, 0.0],
            [10.0, 0.0, 1.0],
        ]
    )
    fixed = {
        upper: np.asarray([0.0, 0.0, 0.0]),
        lower: np.asarray([0.0, 1.0, 0.0]),
        third: np.asarray([0.0, 0.0, 1.0]),
    }

    synchronized = synchronize_and_reinject(identities, coordinates, fixed)

    assert np.array_equal(synchronized[0], fixed[upper])
    assert np.array_equal(synchronized[1], fixed[lower])
    assert np.array_equal(synchronized[2], fixed[third])


def test_geometry_inputs_fail_closed() -> None:
    with pytest.raises(DiffusionGeometryError, match="at least three"):
        weighted_kabsch(np.zeros((2, 3)), np.zeros((2, 3)))
    with pytest.raises(DiffusionGeometryError, match="strictly positive"):
        weighted_kabsch(
            np.zeros((3, 3)),
            np.zeros((3, 3)),
            np.asarray([1.0, 0.0, 1.0]),
        )

    identities = tuple(AtomIdentity("A", str(i), "", "CA") for i in range(3))
    fixed = {identity: np.zeros(3) for identity in identities}
    with pytest.raises(DiffusionGeometryError, match="identity count"):
        synchronize_and_reinject(identities, np.zeros((4, 3)), fixed)
