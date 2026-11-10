"""Independent finite-probe accounting rejects missing/duplicate signed states."""
import math

import pytest
from scripts.audit_chi_stationarity import finite_summary


def probes(slope):
    return [dict(chi=1,step_degrees=h,sign=sign,score=2.+slope*sign*math.radians(h))
            for h in (.01,.005) for sign in (-1,1)]


def test_finite_descent_is_not_a_convergence_or_native_accuracy_measure():
    r=finite_summary(2.,probes(3.),1,[.01,.005],1e-7)
    assert r['resolved_descent']
    assert r['maximum_absolute_slope_per_radian']==pytest.approx(3.)


def test_reporting_floor_keeps_roundoff_distinct_from_resolved_descent():
    r=finite_summary(2.,probes(1e-6),1,[.01,.005],1e-7)
    assert not r['resolved_descent']


def test_missing_or_duplicate_signed_probe_refuses():
    p=probes(1.)
    for rows in (p[:3],p+[p[0]]):
        with pytest.raises(ValueError,match='inventory'):
            finite_summary(2.,rows,1,[.01,.005],1e-7)
