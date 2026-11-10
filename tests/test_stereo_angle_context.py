"""The context HA descriptor is geometric, not a chemical/chirality gate."""
import math

import pytest
from scripts.audit_stereo_angle_context import digest, triad_descriptor


@pytest.mark.parametrize('angles,expected', [
    (dict(N=64., C=78., CB=73.), 'all_three_acute'),
    (dict(N=109., C=112., CB=107.), 'all_three_obtuse'),
    (dict(N=113., C=125., CB=83.), 'mixed_or_boundary'),
    (dict(N=90., C=90., CB=90.), 'mixed_or_boundary'),
])
def test_complete_triad_descriptor(angles, expected):
    assert triad_descriptor(angles) == expected


@pytest.mark.parametrize('angles', [dict(N=109., C=109.), dict(N=math.nan, C=109., CB=109.)])
def test_incomplete_or_nonfinite_triad_is_not_assigned(angles):
    with pytest.raises(ValueError):
        triad_descriptor(angles)


def test_digest_ignores_dictionary_insertion_order_but_not_coordinates():
    assert digest(dict(a=1, b=2)) == digest(dict(b=2, a=1))
    assert digest(dict(a=1, b=2)) != digest(dict(a=1, b=3))
