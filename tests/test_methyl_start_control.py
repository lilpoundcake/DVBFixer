"""Initial methyl rotation is generated-only rigid geometry, not native ranking."""
import math

import pytest
from scripts.research_methyl_start_control import (
    ANGLES,
    methyl_indices,
    rotate_methyl,
    selected_angle,
)


def source():
    names = ('CG2', 'CB', 'HG21', 'HG22', 'HG23', 'OTHER')
    rows = [dict(key=['A', '24', '', 'THR', n], element='H' if n.startswith('HG') else 'C') for n in names]
    points = [[1., 2., 3.], [1., 2., 4.], [1.1, 2., 2.96], [.95, 2.08, 2.96], [.95, 1.92, 2.96], [9., 8., 7.]]
    return rows, points


def test_all_declared_angles_preserve_parent_and_internal_geometry():
    rows, points = source()
    original = [p[:] for p in points]
    for angle in ANGLES:
        rotated = rotate_methyl(rows, points, angle)
        assert rotated[0] == points[0] and rotated[1] == points[1] and rotated[-1] == points[-1]
        for i in range(2, 5):
            assert math.dist(rotated[i], rotated[0]) == pytest.approx(math.dist(points[i], points[0]), abs=1e-12)
            assert math.dist(rotated[i], rotated[1]) == pytest.approx(math.dist(points[i], points[1]), abs=1e-12)
            for j in range(2, 5):
                assert math.dist(rotated[i], rotated[j]) == pytest.approx(math.dist(points[i], points[j]), abs=1e-12)
    assert rotate_methyl(rows, points, 0) == original == points


@pytest.mark.parametrize('fault', ['duplicate', 'chain_case', 'element', 'zero_axis', 'angle', 'nonfinite'])
def test_unsupported_identity_or_geometry_is_not_guessed(fault):
    rows, points = source()
    angle = 15
    if fault == 'duplicate':
        rows[-1]['key'] = rows[0]['key'][:]
    elif fault == 'chain_case':
        rows[0]['key'][0] = 'a'
    elif fault == 'element':
        rows[2]['element'] = 'C'
    elif fault == 'zero_axis':
        points[1] = points[0][:]
    elif fault == 'angle':
        angle = 17
    else:
        points[2][0] = math.nan
    with pytest.raises(ValueError):
        rotate_methyl(rows, points, angle)


def test_selection_uses_only_the_complete_fixed_initial_energy_grid():
    energies = {a: 100. for a in ANGLES}
    energies[15] = energies[345] = 10.
    assert selected_angle(energies) == 345
    energies[0] = 10.+1e-9
    assert selected_angle(energies) == 0
    energies[60] = 5.
    assert selected_angle(energies) == 60
    assert methyl_indices(source()[0]) == [0, 1, 2, 3, 4]


@pytest.mark.parametrize('fault', ['missing', 'nonfinite', 'extra'])
def test_incomplete_or_nonfinite_selection_is_refused(fault):
    energies = {a: 100. for a in ANGLES}
    if fault == 'missing':
        del energies[0]
    elif fault == 'nonfinite':
        energies[0] = math.nan
    else:
        energies[360] = 1.
    with pytest.raises(ValueError):
        selected_angle(energies)
