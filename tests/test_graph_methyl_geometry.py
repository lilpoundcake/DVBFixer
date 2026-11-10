"""Generic geometry and enrolled one-pass ownership, without sampling outcomes."""
import copy
import math

import pytest
from scripts.graph_methyl_geometry import PHASES, canonical, choose, discover, initialize


def source():
    rows = [dict(key=['a', '82', 'A', 'THR', name], element=e) for name, e in
        [('CG2', 'C'), ('CB', 'C'), ('HG21', 'H'), ('HG22', 'H'), ('HG23', 'H'), ('CA', 'C')]]
    bonds = [[0, 1], [0, 2], [0, 3], [0, 4], [1, 5]]
    points = [[0., 0., 0.], [0., 0., .15], [.01, .001, -.1], [.002, .01, -.11], [.01, .01, -.09], [.1, 0., .2]]
    group = discover(rows, bonds, {2, 3, 4}, {2, 3, 4})['groups'][0]
    return rows, bonds, points, group


def test_graph_discovery_preserves_complete_identity_and_owns_only_generated_mobile_H():
    rows, bonds, _, group = source()
    assert group['parent_key'] == ['a', '82', 'A', 'THR', 'CG2']
    assert group['frame'] == 5
    result = discover(rows, bonds, {2, 3}, {2, 3, 4})
    assert not result['groups'] and len(result['exclusions']) == 1
    assert not discover(rows, bonds, {2, 3, 4}, {2, 3})['groups']


def test_all_phases_canonical_radii_and_unowned_coordinates_independently():
    _, _, points, group = source()
    before = copy.deepcopy(points)
    for phase in PHASES:
        after = canonical(points, group, phase)
        assert all(after[i] == points[i] for i in (0, 1, 5))
        unit = []
        for i in (2, 3, 4):
            radius = math.dist(points[0], points[i])
            assert math.dist(after[0], after[i]) == pytest.approx(radius, abs=1e-12)
            vector = [x/radius for x in after[i]]
            assert vector[2] == pytest.approx(-1/3, abs=1e-12)
            unit.append(vector)
        for i, j in ((0, 1), (0, 2), (1, 2)):
            assert sum(x*y for x, y in zip(unit[i], unit[j])) == pytest.approx(-1/3, abs=1e-12)
        assert math.degrees(math.atan2(after[2][1], after[2][0]))%360 == pytest.approx(phase)
    assert points == before


def test_frame_prefers_same_residue_then_case_sensitive_identity():
    rows, bonds, _, _ = source()
    rows.append(dict(key=['A', '1', '', 'ALA', 'C'], element='C'))
    assert discover(rows, bonds+[[1, 6]], {2, 3, 4}, {2, 3, 4})['groups'][0]['frame'] == 5
    rows[5]['key'] = ['b', '2', '', 'ALA', 'C']
    assert discover(rows, bonds+[[1, 6]], {2, 3, 4}, {2, 3, 4})['groups'][0]['frame'] == 6


@pytest.mark.parametrize('fault', ['duplicate', 'bond', 'frame', 'context', 'heavy_owner', 'H_extra_bond'])
def test_invalid_native_graph_or_ownership_fails_closed(fault):
    rows, bonds, _, _ = source()
    owned = {2, 3, 4}
    if fault == 'duplicate':
        rows[5]['key'] = rows[0]['key'][:]
    elif fault == 'bond':
        bonds.append([1, 0])
    elif fault == 'frame':
        bonds.pop()
    elif fault == 'context':
        rows[0]['key'][3] = 'UNK'
    elif fault == 'heavy_owner':
        owned.add(0)
    else:
        bonds.append([2, 5])
    with pytest.raises(ValueError):
        discover(rows, bonds, owned, {2, 3, 4})


def test_one_pass_starts_all_groups_at_zero_and_does_not_revisit():
    _, _, points, group = source()
    second = dict(parent=6, axis=7, frame=11, hydrogens=[8, 9, 10],
        parent_key=['b', '1', '', 'THR', 'CG2'])
    points += [[x+1 for x in p] for p in points]
    original = copy.deepcopy(points)
    zeros, _ = initialize(points, [group, second])
    seen = []
    def energy(candidate):
        seen.append(copy.deepcopy(candidate))
        return 0.
    result, grids = initialize(points, [group, second], selected=True, energy=energy)
    assert len(seen) == 48 and len(grids) == 2
    assert result == zeros and points == original
    assert all(seen[j][i] == zeros[i] for j in range(24) for i in (8, 9, 10))
    assert all(seen[j][i] == zeros[i] for j in range(24, 48) for i in (2, 3, 4))
    assert all(g['selected_phase'] == 0 for g in grids)
    with pytest.raises(ValueError, match='order'):
        initialize(points, [second, group])


def test_phase_ties_and_nonfinite_scores():
    energies = dict.fromkeys(PHASES, 100.)
    energies[15] = energies[345] = 1.
    assert choose(energies) == 345
    energies[0] = 1.+1e-9
    assert choose(energies) == 0
    energies[0] = math.nan
    with pytest.raises(ValueError):
        choose(energies)


def test_one_pass_keeps_previous_selected_group_and_original_radii():
    _, _, points, group = source()
    second = dict(parent=6, axis=7, frame=11, hydrogens=[8, 9, 10],
        parent_key=['b', '1', '', 'THR', 'CG2'])
    points += [[x+1 for x in p] for p in points]
    groups = [group, second]
    seen = []
    def energy(candidate):
        call = len(seen)
        seen.append(copy.deepcopy(candidate))
        phase = PHASES[call%24]
        return 0. if phase == (60 if call < 24 else 345) else 1.
    result, grids = initialize(points, groups, selected=True, energy=energy)
    first_expected = canonical(points, group, 60)
    second_expected = canonical(points, second, 345)
    assert [g['selected_phase'] for g in grids] == [60, 345]
    assert len(seen) == 48
    assert all(seen[j][i] == first_expected[i] for j in range(24, 48) for i in (2, 3, 4))
    assert all(result[i] == first_expected[i] for i in (2, 3, 4))
    assert all(result[i] == second_expected[i] for i in (8, 9, 10))


@pytest.mark.parametrize('fault', ['axis', 'collinear', 'radius', 'phase', 'nonfinite'])
def test_undefined_geometry_has_no_arbitrary_fallback(fault):
    _, _, points, group = source()
    phase = 0
    if fault == 'axis':
        points[1] = points[0][:]
    elif fault == 'collinear':
        points[5] = [0., 0., .2]
    elif fault == 'radius':
        points[2] = points[0][:]
    elif fault == 'phase':
        phase = 7
    else:
        points[2][0] = math.nan
    with pytest.raises(ValueError):
        canonical(points, group, phase)
