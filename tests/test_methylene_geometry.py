"""Heavy-plane CH2 geometry owns generated H, never observations or phase search."""
import copy
import math

import pytest
from scripts.research_methylene_geometry import canonical, discover, initialize


def source(angle=109.47122063449069):
    rows = [dict(key=['a', '82', 'A', 'LEU', n], element=e) for n, e in
        [('CB', 'C'), ('CA', 'C'), ('CG', 'C'), ('HB2', 'H'), ('HB3', 'H'), ('OTHER', 'C')]]
    radians = math.radians(angle)
    points = [[0., 0., 0.], [.15, 0., 0.], [.15*math.cos(radians), .15*math.sin(radians), 0.], [.001, 0., -.1], [.002, .001, -.11], [9., 8., 7.]]
    bonds = [[0, 1], [0, 2], [0, 3], [0, 4]]
    group = discover(rows, bonds, {3, 4}, {3, 4})['groups'][0]
    return rows, bonds, points, group


@pytest.mark.parametrize('angle', [90., 100., 109.47122063449069, 120., 130.])
def test_each_H_heavy_angle_is_tetrahedral_but_heavy_angle_is_unchanged(angle):
    _, _, points, group = source(angle)
    original = copy.deepcopy(points)
    actual = canonical(points, group)
    assert all(actual[i] == points[i] for i in (0, 1, 2, 5))
    for h in (3, 4):
        radius = math.hypot(*points[h])
        assert math.hypot(*actual[h]) == pytest.approx(radius, abs=1e-12)
        for heavy in (1, 2):
            dot = sum(x*y for x, y in zip(actual[h], actual[heavy]))/(radius*math.hypot(*actual[heavy]))
            assert dot == pytest.approx(-1/3, abs=1e-12)
    assert actual[3][2] > 0 and actual[4][2] < 0
    assert initialize(points, [group]) == actual
    assert points == original


def test_kept_or_measured_H_is_excluded_without_rewriting():
    rows, bonds, _, _ = source()
    assert not discover(rows, bonds, {3}, {3, 4})['groups']
    assert not discover(rows, bonds, {3, 4}, {3})['groups']


def test_nonideal_heavy_angle_implies_H_H_angle_not_a_four_ideal_bond_claim():
    _, _, points, group = source(120.)
    actual = canonical(points, group)
    cosine = sum(x*y for x, y in zip(actual[3], actual[4]))/(math.hypot(*actual[3])*math.hypot(*actual[4]))
    assert cosine == pytest.approx(-1/9, abs=1e-12)
    assert cosine != pytest.approx(-1/3)


@pytest.mark.parametrize('fault', ['parallel', 'opposite', 'infeasible', 'zero_H', 'nonfinite'])
def test_undefined_or_infeasible_heavy_plane_is_refused_not_fudged(fault):
    _, _, points, group = source()
    if fault == 'parallel':
        points[2] = points[1][:]
    elif fault == 'opposite':
        points[2] = [-x for x in points[1]]
    elif fault == 'infeasible':
        _, _, points, group = source(150.)
    elif fault == 'zero_H':
        points[3] = points[0][:]
    else:
        points[3][0] = math.nan
    with pytest.raises(ValueError):
        canonical(points, group)


def test_proper_rigid_transform_covariance():
    _, _, points, group = source()
    def transform(p):
        x, y, z = p
        return [-y+3, x-2, z+1]
    expected = [transform(p) for p in canonical(points, group)]
    actual = canonical([transform(p) for p in points], group)
    for a, b in zip(actual, expected, strict=True):
        assert a == pytest.approx(b, abs=1e-12)
