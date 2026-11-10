"""Native endpoint drift is measured, not silently replaced by a callback."""
import math

import pytest
from scripts.audit_carbon_H_whole_path import finite_points, terminal_summary


def test_terminal_report_keeps_rejected_returned_coordinates_and_energy_visible():
    initial, reported, returned = [[0., 0., 0.]], [[.1, 0., 0.]], [[.11, 0., 0.]]
    callbacks = [dict(iteration=7, positions_nm=[.1, 0., 0.])]
    nodes = [dict(energy_kj_mol=e, centres=[dict(volume_angstrom3=v)])
        for e, v in ((20., 2.), (10., .02), (900., -.006))]
    result = terminal_summary(initial, returned, callbacks, nodes)
    assert result['comparison'] == 'last reported iteration'
    assert not result['returned_matches_last_reported_nm']
    assert result['returned_minus_reported_energy_kj_mol'] == 890.
    assert result['maximum_returned_vs_reported_displacement_angstrom'] == pytest.approx(.1)
    assert result['last_reported_centres'][0]['volume_angstrom3'] > 0
    assert result['returned_centres'][0]['volume_angstrom3'] < 0
    assert reported == [[.1, 0., 0.]] and returned == [[.11, 0., 0.]]


def test_zero_callback_phase_is_not_labelled_a_reported_accepted_iteration():
    points = [[0., 0., 0.]]
    nodes = [dict(energy_kj_mol=0., centres=[])]*2
    result = terminal_summary(points, points, [], nodes)
    assert result['last_reported_iteration'] is None
    assert result['comparison'] == 'initial state; no reported iteration'
    assert result['returned_matches_last_reported_nm']


@pytest.mark.parametrize('points', [[], [[0., 0.]], [[math.nan, 0., 0.]]])
def test_partial_or_nonfinite_callback_is_not_clean_geometry(points):
    with pytest.raises(ValueError):
        finite_points(points, 1)
