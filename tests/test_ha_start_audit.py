"""Independent HA replay auditing does not use the intervention constructor."""
import math

import pytest
from scripts.audit_ha_start_control import ha_geometry, intervention, inventory, observed_motion


def rows(outward=False):
    radius = math.sqrt(8/9)
    points = [[0., 0., 0.], [radius, 0., 1/3], [-radius/2, radius*math.sqrt(3)/2, 1/3],
        [-radius/2, -radius*math.sqrt(3)/2, 1/3], [0., 0., -.1 if outward else .1]]
    return [dict(key=['A', '16', '', 'GLN', name], element=name[0], xyz=point)
            for name, point in zip(('CA', 'N', 'C', 'CB', 'HA'), points, strict=True)]


def test_independent_expected_one_ha_intervention():
    result = intervention(rows(), rows(True), 'outward_HA')
    assert result['initial_HA_movement_angstrom'] == pytest.approx(.2)
    assert result['non_HA_coordinates_exact']
    assert min(result['after']['ha_angles_degrees'].values()) > 90
    assert ha_geometry(rows(True))['ca_ha_distance_angstrom'] == pytest.approx(.1)
    assert intervention(rows(), rows(), 'unchanged')['changed_atoms'] == []


@pytest.mark.parametrize('fault', ['heavy', 'inward', 'length', 'unchanged_changed', 'element', 'missing'])
def test_independent_audit_rejects_undeclared_or_malformed_intervention(fault):
    before, after = rows(), rows(True)
    arm = 'outward_HA'
    if fault == 'heavy':
        after[0]['xyz'][0] += .001
    elif fault == 'inward':
        after = rows()
    elif fault == 'length':
        after[-1]['xyz'][-1] *= 2
    elif fault == 'unchanged_changed':
        arm = 'unchanged'
    elif fault == 'element':
        after[-1]['element'] = 'C'
    else:
        after.pop()
    with pytest.raises(ValueError):
        intervention(before, after, arm)


def test_duplicate_or_nonfinite_inventory_is_not_silently_reduced():
    data = rows()
    with pytest.raises(ValueError):
        inventory(data+data[:1])
    data[-1]['xyz'][0] = math.nan
    with pytest.raises(ValueError):
        inventory(data)


def test_observation_motion_counts_only_original_heavy_records():
    source = b'ATOM      1  CA  GLN A  16       0.000   0.000   0.000  1.00  0.00           C\n'
    data = rows()
    assert observed_motion(source, data)['exact']
    data[0]['xyz'][0] = 1.
    motion = observed_motion(source, data)
    assert motion['observed_heavy_atoms'] == 1
    assert motion['maximum_angstrom'] == 1.
    assert not motion['observations_restored']
    with pytest.raises(ValueError):
        observed_motion(source, data[1:])
