"""Physical finite differences must not be redefined to fit a native gradient."""
import numpy as np
import pytest
from scripts.audit_post_h_1vcu import penalty
from scripts.research_post_h_cpu_sign import compare, predictions


def test_source_sign_hypothesis_is_not_derivative_of_positive_penalty_energy():
    points = np.array([[0., 0., 0.], [0.12, 0., 0.]])
    constraints = [(0, 1, 0.1)]
    _, grad, _, _ = penalty(points, constraints, 1000.)
    left, right = points.copy(), points.copy()
    left[1, 0] -= 1e-6
    right[1, 0] += 1e-6
    derivative = (penalty(right, constraints, 1000.)[0]-penalty(left, constraints, 1000.)[0])/2e-6
    values = predictions(np.zeros((2, 3)), grad)
    assert values['physical'][1, 0] == pytest.approx(derivative, abs=1e-7)
    assert values['source_sign'][1, 0] == pytest.approx(-derivative, abs=1e-7)
    assert values['source_sign_float32'][1, 0] < 0


def test_matching_source_formula_does_not_pass_original_physical_comparison():
    policy = dict(gradient_absolute_bound=1e-3, gradient_relative_bound=2e-5)
    physical = np.array([[10., 20., 30.], [1., 2., 3.]])
    grad = np.array([[1., 0., 0.], [-1., 0., 0.]])
    values = predictions(physical, grad)
    atoms = [dict(key=['A', str(i), '', 'HOH', 'O']) for i in range(2)]
    native = values['source_sign_float32']
    assert not compare(values['physical'], native, policy, atoms, {0, 1})['passed']
    assert compare(values['source_sign_float32'], native, policy, atoms, {0, 1})['passed']
