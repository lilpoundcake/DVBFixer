"""Independent bounded policy refuses early stops and treats exhaustion distinctly."""
import pytest
from scripts.audit_chi_endpoint_guard import expected_decision


def test_measured_endpoint_not_native_error_controls_stopping():
    assert expected_decision(1e-5,1,3,1e-4)=='endpoint_gradient_met'
    assert expected_decision(.1,1,3,1e-4)=='continue'


def test_exhaustion_is_not_relabelled_as_endpoint_gradient_met():
    assert expected_decision(.1,3,3,1e-4)=='budget_exhausted'
    with pytest.raises(ValueError):
        expected_decision(float('nan'),3,3,1e-4)
