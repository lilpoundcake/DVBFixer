"""Derivative diagnostic bounds are finite and independent of native accuracy."""
import math

import pytest
from scripts.research_chi_derivatives import derivative_comparison


def test_derivative_unit_conversion_per_degree_to_per_radian():
    per_degree=.01
    per_radian=per_degree*180/math.pi
    assert per_radian*math.pi/180==pytest.approx(per_degree)
    assert derivative_comparison(per_degree,per_radian)['within_reporting_bound'] is False


def test_bounds_report_mismatch_without_altering_values():
    assert derivative_comparison(.1,.1000001)['within_reporting_bound']
    result=derivative_comparison(.1,.2)
    assert not result['within_reporting_bound']
    assert result['absolute_error']==.1


def test_nonfinite_gradient_is_not_a_pass():
    with pytest.raises(ValueError,match='nonfinite'):
        derivative_comparison(float('nan'),0.)
