"""Independent signed-axis-angle and fixed-grid selection controls."""
import math

import pytest
from scripts.audit_methyl_start_control import choose_independently, rotation_invariants


def source():
    names = ('CG2', 'CB', 'HG21', 'HG22', 'HG23', 'CA')
    rows = [dict(key=['A', '24', '', 'THR', n], element='H' if n.startswith('HG') else 'C') for n in names]
    points = [[0., 0., 0.], [0., 0., 1.], [1., 0., -.3], [0., 1., -.3], [-1., -1., -.3], [3., 2., 1.]]
    return rows, points


def test_independent_rotation_audit_accepts_right_handed_axis_angle():
    rows, points = source()
    rotated = [p[:] for p in points]
    for i in (2, 3, 4):
        x, y, z = points[i]
        rotated[i] = [-y, x, z]
    assert rotation_invariants(rows, points, rotated, 90) == [0, 1, 2, 3, 4]
    with pytest.raises(ValueError, match='orientation'):
        rotation_invariants(rows, points, rotated, 270)


@pytest.mark.parametrize('fault', ['other_atom', 'radius', 'axis_component', 'missing_identity'])
def test_independent_audit_refuses_nonrigid_or_undeclared_motion(fault):
    rows, points = source()
    after = [p[:] for p in points]
    if fault == 'other_atom':
        after[-1][0] += .001
    elif fault == 'radius':
        after[2][0] *= 2
    elif fault == 'axis_component':
        after[2][2] += .01
    else:
        rows[0]['key'][0] = 'a'
    with pytest.raises(ValueError):
        rotation_invariants(rows, points, after, 0)


def test_independent_selection_is_initial_energy_only_with_fixed_ties():
    grid = [dict(angle_degrees=a, energy_kj_mol=100.) for a in range(0, 360, 15)]
    grid[1]['energy_kj_mol'] = grid[-1]['energy_kj_mol'] = 10.
    assert choose_independently(grid) == 345
    grid[0]['energy_kj_mol'] = 10.+1e-9
    assert choose_independently(grid) == 0
    grid[0]['energy_kj_mol'] = math.nan
    with pytest.raises(ValueError):
        choose_independently(grid)
    with pytest.raises(ValueError):
        choose_independently(grid[1:])
