"""Topology-bounded proper torsions preserve all known coordinates and handedness."""
import math

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError
from dvbfixer.model.diffusion.rfd3_terminal_torsions import (
    rotate_terminal_component,
    terminal_component,
)


def fixture():
    a, b, c, d, fixed = [AtomIdentity("A", str(i), "", "CA") for i in range(5)]
    points = {a: (1., 1., 0.), b: (1., 0., 0.), c: (0., 0., 0.), d: (0., 0., 1.), fixed: (3., 3., 3.)}
    return points, ((a, b), (b, c), (c, d), (d, fixed)), (c, d), a


def test_graph_selected_proper_rotation_preserves_distances_and_fixed_side():
    points, bonds, axis, seed = fixture()
    component = terminal_component(bonds, axis=axis, seed=seed, allowed=frozenset(points) - {axis[1], next(reversed(points))})
    after = rotate_terminal_component(points, component=component, axis=axis, degrees=90.)
    assert after[seed] == pytest.approx((-1., 1., 0.))
    for first, second in bonds:
        assert math.dist(after[first], after[second]) == pytest.approx(math.dist(points[first], points[second]))
    assert all(after[a] == points[a] for a in points.keys() - component)
    forward = np.array([after[a] for a in list(points)[:4]])
    before = np.array([points[a] for a in list(points)[:4]])
    assert np.linalg.det(forward[:3] - forward[3]) == pytest.approx(np.linalg.det(before[:3] - before[3]))


@pytest.mark.parametrize("problem", ["cycle", "fixed", "duplicate", "missing-axis", "self"])
def test_unsafe_terminal_components_are_refused(problem):
    points, bonds, axis, seed = fixture()
    allowed = frozenset(points)
    if problem == "cycle":
        bonds += ((seed, axis[1]),)
    elif problem == "fixed":
        allowed -= {seed}
    elif problem == "duplicate":
        bonds += (tuple(reversed(bonds[0])),)
    elif problem == "missing-axis":
        bonds = tuple(p for p in bonds if p != axis)
    else:
        bonds += ((seed, seed),)
    with pytest.raises(CustodyError):
        terminal_component(bonds, axis=axis, seed=seed, allowed=allowed)


@pytest.mark.parametrize("degrees", [True, float("nan"), float("inf")])
def test_invalid_rotation_angles_fail_closed(degrees):
    points, _, axis, seed = fixture()
    with pytest.raises(CustodyError):
        rotate_terminal_component(points, component=frozenset({seed}), axis=axis, degrees=degrees)
