"""Canonical angular initialization owns only three generated H coordinates."""
import math

import pytest
from scripts.research_tetra_methyl_control import (
    ANGLES,
    ARMS,
    methyl_hch_angles,
    tetrahedral_methyl,
)


def source():
    names = ('CG2', 'CB', 'HG21', 'HG22', 'HG23', 'CA', 'OTHER')
    rows = [dict(key=['A', '24', '', 'THR', n], element='H' if n.startswith('HG') else 'C') for n in names]
    points = [[0., 0., 0.], [0., 0., 1.], [.01, 0., -.1], [.005, .006, -.1], [-.004, .002, -.1], [1., 0., 2.], [9., 8., 7.]]
    return rows, points


def unit(vector):
    norm = math.hypot(*vector)
    return [x/norm for x in vector]


def test_all_fixed_phases_are_tetrahedral_preserve_radii_and_other_atoms():
    rows, points = source()
    original = [p[:] for p in points]
    for angle in ANGLES:
        result = tetrahedral_methyl(rows, points, angle)
        assert all(result[i] == points[i] for i in (0, 1, 5, 6))
        directions = []
        for i in (2, 3, 4):
            assert math.dist(result[i], result[0]) == pytest.approx(math.dist(points[i], points[0]), abs=1e-12)
            directions.append(unit([x-y for x, y in zip(result[i], result[0])]))
        axis = unit([x-y for x, y in zip(result[1], result[0])])
        for direction in directions:
            assert sum(x*y for x, y in zip(direction, axis)) == pytest.approx(-1/3, abs=1e-12)
        for i, j in ((0, 1), (0, 2), (1, 2)):
            assert sum(x*y for x, y in zip(directions[i], directions[j])) == pytest.approx(-1/3, abs=1e-12)
        assert methyl_hch_angles(rows, result) == pytest.approx([math.degrees(math.acos(-1/3))]*3)
    assert points == original
    assert ARMS == ('unchanged', 'tetrahedral_zero', 'tetrahedral_energy_selected')


def test_frame_is_translation_and_proper_rotation_covariant():
    rows, points = source()
    def transform(point):
        x, y, z = point
        return [-y+7, x-3, z+2]

    actual = tetrahedral_methyl(rows, [transform(p) for p in points], 75)
    expected = [transform(p) for p in tetrahedral_methyl(rows, points, 75)]
    for first, second in zip(actual, expected, strict=True):
        assert first == pytest.approx(second, abs=1e-12)


@pytest.mark.parametrize('fault', ['axis', 'frame', 'missing_ca', 'zero_h', 'nonfinite', 'angle'])
def test_undefined_or_undeclared_geometry_is_refused(fault):
    rows, points = source()
    angle = 0
    if fault == 'axis':
        points[1] = points[0][:]
    elif fault == 'frame':
        points[5] = [0., 0., 2.]
    elif fault == 'missing_ca':
        rows[5]['key'][0] = 'a'
    elif fault == 'zero_h':
        points[2] = points[0][:]
    elif fault == 'nonfinite':
        points[2][0] = math.nan
    else:
        angle = 17
    with pytest.raises(ValueError):
        tetrahedral_methyl(rows, points, angle)
