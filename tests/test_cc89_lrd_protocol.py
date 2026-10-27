"""Boundary policies for research LRD must not conflate source and generated atoms."""

from dataclasses import replace

import numpy as np
import pytest
from scripts.cc89_lrd_protocol import (
    boundary_residues,
    low_noise_schedule,
    observed_frame_center,
    project_static,
    relaxation_residues,
)

from tests.test_diffusion_pdb_materialize import _request


def test_generated_only_and_mobile_source_windows_are_distinct():
    request = _request()
    gap = request.gaps[0]
    assert boundary_residues(request, gap, include_observed=False) == frozenset(gap.generated_residues)
    assert boundary_residues(request, gap, include_observed=True) == frozenset(
        (*gap.generated_residues, gap.left_anchor, gap.right_anchor))


def test_anchor_minimization_is_explicit_and_does_not_mutate_generated_identity():
    request = _request()
    generated = request.generated_atoms
    gap = request.gaps[0]
    assert relaxation_residues(request) == gap.generated_residues
    assert set(relaxation_residues(request, include_anchors=True)) == {
        *gap.generated_residues, gap.left_anchor, gap.right_anchor}
    assert request.generated_atoms == generated


def test_anchor_minimization_keeps_case_and_insertion_codes():
    from dvbfixer.model.diffusion.contract import ResidueIdentity

    request = _request()
    left, right = ResidueIdentity('A', '82', 'A'), ResidueIdentity('A', '82', 'B')
    gap = replace(request.gaps[0], left_anchor=left, right_anchor=right,
                  movable_junction_residues=(left, *request.gaps[0].generated_residues, right))
    # The pure helper accepts a request-like view so unrelated contract fields
    # need not be falsified to test identity preservation.
    from types import SimpleNamespace
    residues = relaxation_residues(SimpleNamespace(gaps=(gap,)), include_anchors=True)
    assert left in residues and right in residues and len(residues) == 5


def test_boundary_window_uses_positions_not_author_numbers_or_insertion_codes():
    from dvbfixer.model.diffusion.contract import ResidueIdentity

    request = _request()
    gap = request.gaps[0]
    identities = (ResidueIdentity('A', '82', 'A'), ResidueIdentity('A', '400'), ResidueIdentity('A', '-3'))
    gap = replace(gap, generated_residues=identities,
                  movable_junction_residues=(gap.left_anchor, *identities, gap.right_anchor))
    assert boundary_residues(request, gap, include_observed=False) == frozenset(identities)
    assert all(r.chain == 'A' for r in boundary_residues(request, gap, include_observed=True))


def test_schedule_has_exactly_25_updates_and_terminates_at_zero():
    schedule = low_noise_schedule()
    assert len(schedule) == 26
    assert schedule[0] == pytest.approx(1.5)
    assert schedule[-2] == pytest.approx(0.0064)
    assert schedule[-1] == 0
    assert np.all(np.diff(schedule) < 0)


@pytest.mark.parametrize('kwargs', [{'steps': 1}, {'sigma_start': 0}, {'rho': float('nan')}, {'sigma_min': 2}])
def test_invalid_schedules_fail(kwargs):
    with pytest.raises(ValueError):
        low_noise_schedule(**kwargs)


def test_projection_changes_only_selected_atoms_without_restoring_mobile_source():
    reference = np.arange(18., dtype=float).reshape(2, 3, 3)
    state = reference + 5
    mobile = np.array([[True, False, False], [False, True, False]])
    result = project_static(state, reference, mobile)
    np.testing.assert_array_equal(result[~mobile], reference[~mobile])
    np.testing.assert_array_equal(result[mobile], state[mobile])
    with pytest.raises(ValueError):
        project_static(state, reference, mobile.astype(int))


def test_frame_depends_on_original_observations_not_generated_coordinates():
    coordinates = np.arange(18., dtype=np.float32).reshape(2, 3, 3)
    observed = np.array([[True, True, False], [True, False, False]])
    center = observed_frame_center(coordinates, observed)
    changed = coordinates.copy()
    changed[~observed] += 1000
    np.testing.assert_array_equal(observed_frame_center(changed, observed), center)
    with pytest.raises(ValueError):
        observed_frame_center(coordinates, np.zeros_like(observed))
