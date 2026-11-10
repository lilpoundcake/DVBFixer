"""Numerical projection controls do not share analytic angle-force code."""
import math

import pytest
from scripts.audit_stereo_angle_readback import angle_energy, directional_force


def test_numerical_projection_has_the_independent_closed_form_sign_and_units():
    points = [[1., 0, 0], [0, 0, 0], [0, 1., 0]]
    directions = [[0, 1., 0], [0, 0, 0], [0, 0, 0]]
    assert angle_energy(points, math.pi/3, 2.) == pytest.approx(math.pi**2/36)
    assert directional_force(points, directions, math.pi/3, 2.) == pytest.approx(math.pi/3, abs=1e-8)
    opposite = [[-x for x in direction] for direction in directions]
    assert directional_force(points, opposite, math.pi/3, 2.) == pytest.approx(-math.pi/3, abs=1e-8)


def test_equilibrium_has_zero_first_derivative():
    points = [[1., 0, 0], [0, 0, 0], [0, 1., 0]]
    directions = [[0, 1., 0], [0, 0, 0], [0, 0, 0]]
    assert directional_force(points, directions, math.pi/2, 2.) == pytest.approx(0, abs=1e-8)


@pytest.mark.parametrize('points', [
    [[0, 0, 0]]*3, [[math.nan, 0, 0], [0, 0, 0], [0, 1, 0]], [[0, 0]]*3,
])
def test_invalid_numerical_energy_inputs_are_not_silently_clamped(points):
    with pytest.raises(ValueError):
        angle_energy(points, 1.9, 500.)
