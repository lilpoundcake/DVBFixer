"""Actual endpoint-gradient stopping is bounded and independent of native error."""
import pytest
from scripts.research_chi_endpoint_guard import guard_decision


def test_actual_gradient_not_a_small_previous_step_controls_continuation():
    assert guard_decision(.01,1,3,1e-4)=='continue'
    assert guard_decision(1e-5,1,3,1e-4)=='endpoint_gradient_met'


def test_exhausted_budget_is_not_convergence():
    assert guard_decision(.01,3,3,1e-4)=='budget_exhausted'


@pytest.mark.parametrize('maximum,cycle',[(float('nan'),1),(-1.,1),(1.,0),(1.,4)])
def test_invalid_gradient_or_unbounded_cycle_refuses(maximum,cycle):
    with pytest.raises(ValueError):
        guard_decision(maximum,cycle,3,1e-4)
