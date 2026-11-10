"""Independent H-placement ownership and descriptor controls."""
import copy
import math

import pytest
from scripts.audit_ha_placement_trace import descriptor, existing_fixed


@pytest.mark.parametrize('angles,expected', [
    (dict(N=100., C=110., CB=120.), 'all_obtuse'),
    (dict(N=60., C=70., CB=80.), 'all_acute'),
    (dict(N=90., C=110., CB=120.), 'mixed_or_boundary'),
])
def test_triad_is_a_descriptor_not_a_new_gate(angles, expected):
    assert descriptor(angles) == expected


def test_existing_atoms_are_exact_and_hydrogen_addition_is_separate():
    original = [dict(key=['A', '16', '', 'GLN', 'CA'], element='C', xyz=[0., 0., 0.])]
    before = copy.deepcopy(original)+[dict(key=['A', '16', '', 'GLN', 'HA'], element='H', xyz=[1., 0., 0.])]
    after = copy.deepcopy(before)
    after[-1]['xyz'][0] = -1.
    assert existing_fixed(original, before, after)['added_atoms'] == 1
    after[0]['xyz'][0] = 1e-10
    with pytest.raises(ValueError, match='existing observation'):
        existing_fixed(original, before, after)


def test_missing_existing_observation_is_not_ignored():
    original = [dict(key=['A', '16', '', 'GLN', 'CA'], element='C', xyz=[0., 0., 0.])]
    with pytest.raises(ValueError, match='inventory'):
        existing_fixed(original, [], [])


def test_nonfinite_triad_is_not_assigned_a_hemisphere():
    with pytest.raises(ValueError):
        descriptor(dict(N=math.nan, C=110., CB=120.))
