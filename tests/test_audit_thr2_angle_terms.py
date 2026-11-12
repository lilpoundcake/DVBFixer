"""Independent finite-difference and native-force checks for analytic angles."""
import math

import pytest
from scripts.audit_thr2_angle_terms import angle_forces


def test_angle_forces_match_finite_differences():
    points = [[0.1, 0.3, 0.2], [0., 0., 0.], [-0.2, 0.1, -0.1]]
    theta0, strength = 1.9, 300.
    theta, forces = angle_forces(points, theta0, strength)
    assert 0 < theta < math.pi
    step = 1e-6
    for atom in range(3):
        for axis in range(3):
            energies = []
            for sign in (-1, 1):
                shifted = [p[:] for p in points]
                shifted[atom][axis] += sign*step
                angle, _ = angle_forces(shifted, theta0, strength)
                energies.append(0.5*strength*(angle-theta0)**2)
            numerical = -(energies[1]-energies[0])/(2*step)
            assert forces[atom][axis] == pytest.approx(numerical, abs=1e-6)
    for axis in range(3):
        assert sum(force[axis] for force in forces) == pytest.approx(0., abs=1e-10)


def test_degenerate_angle_is_not_silently_normalized():
    with pytest.raises(ValueError, match='singular'):
        angle_forces([[1., 0., 0.], [0., 0., 0.], [-1., 0., 0.]], 1.9, 300.)
