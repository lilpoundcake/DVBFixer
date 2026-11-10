"""Finite sensitivity arithmetic, units and refusal without scientific runtimes."""
import math

import pytest
from scripts.research_chi_stationarity import central_slopes


def test_linear_slope_is_in_energy_per_radian_not_degree():
    h=.01
    delta=3*math.radians(h)
    result=central_slopes(2.,2.-delta,2.+delta,h)
    assert result['slope_per_radian']==pytest.approx(3.)
    assert result['curvature_per_radian_squared']==pytest.approx(0.,abs=1e-7)
    assert result['best_signed_energy_change']==pytest.approx(-delta)


def test_stationary_quadratic_does_not_imply_negative_descent():
    h=.005
    delta=math.radians(h)**2
    result=central_slopes(0.,delta,delta,h)
    assert result['slope_per_radian']==0.
    assert result['curvature_per_radian_squared']==pytest.approx(2.)
    assert result['best_signed_energy_change']>0


@pytest.mark.parametrize('step',[0.,-1.,float('nan'),float('inf')])
def test_invalid_steps_refuse(step):
    with pytest.raises(ValueError,match='finite'):
        central_slopes(1.,1.,1.,step)
