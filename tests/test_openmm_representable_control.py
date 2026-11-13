"""Prospective trigger algebra only; no native runs before freeze."""
import numpy as np
import pytest
from scripts.research_openmm_cpu_sign import analytic, energy_terms, finite_difference


def test_detector_interval_and_old_oversized_trigger_are_distinct():
    points = np.array([[-0.0500025, 0., 0.], [0.0500025, 0., 0.], [10., 10., 0.]])
    below, _ = analytic(points, 1e7, True, 1e6)
    trigger, _ = analytic(points, 5e7, True, 1e6)
    assert np.max(np.abs(below)) < 2**30
    assert 2**30 < np.max(np.abs(trigger)) < 2**31
    assert 1e14 > 2**31


def test_unit_scale_first_stress_step_retains_representable_trigger():
    points = np.array([[-0.0500025, 0., 0.], [0.0500025, 0., 0.], [10., 10., 0.]])
    physical, gp = analytic(points, 5e7, True, 1e6)
    # Test both declared signs; no presumed native event or line-search outcome.
    for gradient in (physical+gp, physical-gp):
        next_points = points-gradient/np.linalg.norm(gradient)
        g, _ = analytic(next_points, 5e7, True, 1e6)
        assert 2**30 < np.max(np.abs(g)) < 2**31
        assert finite_difference(next_points, 5e7, True, 1e6) == pytest.approx(
            sum(analytic(next_points, 5e7, True, 1e6)), abs=1e-5, rel=1e-6)
        assert energy_terms(next_points, 5e7, True, 1e6)[2] > 0
