"""Independent derivative readback retains errors instead of fitting tolerance."""
import pytest
from scripts.audit_chi_derivatives import compare


def test_independent_reporting_bound_is_absolute_plus_relative():
    row=compare(.1,.1000001,1e-6,1e-3)
    assert row['within_reporting_bound']
    assert row['bound_per_degree']==pytest.approx(1e-6+1e-3*.1000001)


def test_mismatch_and_nonfinite_derivative_are_not_certified():
    assert not compare(.1,.2,1e-6,1e-3)['within_reporting_bound']
    with pytest.raises(ValueError,match='nonfinite'):
        compare(float('nan'),0.,1e-6,1e-3)
