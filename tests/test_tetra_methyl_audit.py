"""Independent canonical geometry/azimuth and observation ownership controls."""
import math

import pytest
from scripts.audit_tetra_methyl_control import tetrahedral_invariants


def source():
    names = ('CG2', 'CB', 'HG21', 'HG22', 'HG23', 'CA')
    rows = [dict(key=['A', '24', '', 'THR', n], element='H' if n.startswith('HG') else 'C') for n in names]
    points = [[0., 0., 0.], [0., 0., 1.], [.01, 0., -.1], [.005, .006, -.1], [-.004, .002, -.1], [1., 0., 2.]]
    return rows, points


def canonical(points, phase):
    result = [p[:] for p in points]
    for ordinal, i in enumerate((2, 3, 4)):
        radius = math.hypot(*points[i])
        theta = math.radians(phase+ordinal*120)
        result[i] = [radius*math.sqrt(8/9)*math.cos(theta), radius*math.sqrt(8/9)*math.sin(theta), -radius/3]
    return result


def test_independent_canonical_geometry_and_axis_phase():
    rows, points = source()
    for phase in range(0, 360, 15):
        assert tetrahedral_invariants(rows, points, canonical(points, phase), phase) == [0, 1, 2, 3, 4, 5]
    assert tetrahedral_invariants(rows, points, points, None) == [0, 1, 2, 3, 4, 5]


@pytest.mark.parametrize('fault', ['heavy_motion', 'radius', 'phase', 'unchanged', 'wrong_H_angle', 'frame', 'missing_CA'])
def test_noncanonical_or_unowned_changes_are_refused(fault):
    rows, points = source()
    after = canonical(points, 30)
    phase = 30
    if fault == 'heavy_motion':
        after[5][0] += .01
    elif fault == 'radius':
        after[2] = [x*2 for x in after[2]]
    elif fault == 'phase':
        phase = 45
    elif fault == 'unchanged':
        phase = None
    elif fault == 'wrong_H_angle':
        radius = math.hypot(*points[2])
        after[2] = [radius, 0., 0.]
    elif fault == 'frame':
        points[5] = after[5] = [0., 0., 2.]
    else:
        rows[5]['key'][0] = 'a'
    with pytest.raises(ValueError):
        tetrahedral_invariants(rows, points, after, phase)
