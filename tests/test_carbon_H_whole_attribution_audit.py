"""Independent determinant and finite-difference force audit rejects drift."""
import math

import pytest
from scripts.audit_carbon_H_whole_attribution import (
    independent_gradient,
    require_close,
    scalar_volume,
)


def test_independent_nearplanar_derivatives_retain_sign_and_translation_invariance():
    points = [[0., 0., 0.], [.14, 0., 0.], [0., .14, 0.], [.05, .03, -.001]]
    assert scalar_volume(points) == pytest.approx(-.0000196)
    gradient = independent_gradient(points)
    assert gradient[3][2] == pytest.approx(.0196, abs=1e-10)
    assert all(sum(g[j] for g in gradient) == pytest.approx(0, abs=1e-9) for j in range(3))


@pytest.mark.parametrize('value', [math.nan, math.inf, 2.])
def test_nonfinite_or_wrong_force_closure_refuses(value):
    with pytest.raises(ValueError):
        require_close(value, 1., 'force closure')
