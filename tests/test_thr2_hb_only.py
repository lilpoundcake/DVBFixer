"""Mask and numeric transport controls; never run experimental minimization."""
import math

import pytest
from scripts.research_thr2_hb_only import proposed_points, transport


def source_and_initial():
    source = dict(heavy={0:[0., 0., 0.], 1:[0.1, 0., 0.], 2:[0., 0.1, 0.], 3:[0., 0., 0.1]},
                  hb=dict(index=4, parent=0, neighbors=[1, 2, 3], draws=[0.2, 0.4, 0.6]))
    return source, list(source['heavy'].values())+[[0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]


def test_only_hb_changes_and_input_arrays_are_not_mutated():
    source, initial = source_and_initial()
    original = [list(p) for p in initial]
    result, selected = proposed_points(source, initial, 'HB-only-unperturbed')
    assert selected == [4] and initial == original
    assert result[:4] == initial[:4] and result[5:] == initial[5:]
    assert math.dist(result[0], result[4]) == pytest.approx(0.1)
    assert all(x < 0 for x in result[4])


def test_unchanged_is_exact_and_undeclared_arm_refuses():
    source, initial = source_and_initial()
    assert proposed_points(source, initial, 'unchanged') == (initial, [])
    with pytest.raises(ValueError, match='Undeclared'):
        proposed_points(source, initial, 'other-H')


def test_reference_transport_is_exact_not_tolerance_widened():
    assert transport([[1., 2., 3.]], [[1., 2., 3.]], 'Reference') == 0.
    with pytest.raises(ValueError, match='not exact'):
        transport([[1., 2., 3.]], [[1.+1e-15, 2., 3.]], 'Reference')


def test_captured_float32_rounding_is_bounded_before_observation():
    import numpy as np

    requested = [[1.123456789, -2.123456789, 3.123456789]]
    applied = [[float(np.float32(x)) for x in requested[0]]]
    assert transport(requested, applied, 'captured') > 0
    with pytest.raises(ValueError, match='rounding bound'):
        transport(requested, [[x+1e-4 for x in applied[0]]], 'captured')


@pytest.mark.parametrize('mode,points', [('captured', [[float('nan'), 0., 0.]]), ('Reference', [[float('inf'), 0., 0.]])])
def test_nonfinite_transport_refuses(mode, points):
    with pytest.raises(ValueError):
        transport(points, points, mode)
