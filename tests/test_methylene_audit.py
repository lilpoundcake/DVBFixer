"""Independent H-heavy cosines, signed slots and exact nonownership checks."""
import copy
import math

import pytest
from scripts.audit_ASP49_methylene_control import check_initial


def source():
    rows = [dict(key=['A', '48', '', 'LEU', n], element=e) for n, e in
        [('CB', 'C'), ('CA', 'C'), ('CG', 'C'), ('HB2', 'H'), ('HB3', 'H')]]
    before = [[0., 0., 0.], [.15, 0., 0.], [-.05, .15*math.sqrt(8/9), 0.], [0., 0., -.1], [0., 0., -.11]]
    after = copy.deepcopy(before)
    after[3] = [-.1/3, -.1*math.sqrt(2)/3, .1*math.sqrt(2/3)]
    after[4] = [-.11/3, -.11*math.sqrt(2)/3, -.11*math.sqrt(2/3)]
    group = dict(parent=0, hydrogens=[3, 4], heavy=[1, 2])
    return rows, before, after, [group]


def test_independent_slot_geometry_and_unchanged_empty_intervention():
    rows, before, after, groups = source()
    check_initial(rows, before, after, groups)
    check_initial(rows, before, before, [])


@pytest.mark.parametrize('fault', ['heavy_motion', 'radius', 'angle', 'slot_sign', 'order'])
def test_undeclared_or_noncanonical_CH2_initial_change_is_refused(fault):
    rows, before, after, groups = source()
    if fault == 'heavy_motion':
        after[1][0] += .001
    elif fault == 'radius':
        after[3] = [2*x for x in after[3]]
    elif fault == 'angle':
        after[3] = [0., 0., .1]
    elif fault == 'slot_sign':
        after[3][2] *= -1
    else:
        groups[0]['heavy'] = [2, 1]
    with pytest.raises(ValueError):
        check_initial(rows, before, after, groups)
