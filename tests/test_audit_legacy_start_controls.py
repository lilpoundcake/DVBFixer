"""All-atom starting-state comparisons cannot silently discard H or identities."""
import pytest
from scripts.audit_legacy_start_controls import coordinate_delta


def atom(name,element,x):
    return dict(key=['A','2','','THR',name],element=element,xyz=(x,0.,0.))


def test_h_motion_is_distinct_from_heavy_context_motion():
    result=coordinate_delta([atom('CB','C',0.),atom('HB','H',1.)],[atom('CB','C',0.),atom('HB','H',2.)])
    assert result['heavy']['exact'] and not result['hydrogen']['exact']
    assert result['hydrogen']['maximum_angstrom']==1.


def test_missing_atom_does_not_turn_into_an_intersection_comparison():
    with pytest.raises(ValueError,match='identity'):
        coordinate_delta([atom('CB','C',0.)],[])


def test_changed_element_is_not_silently_accepted():
    with pytest.raises(ValueError,match='element'):
        coordinate_delta([atom('CB','C',0.)],[atom('CB','N',0.)])


def test_small_numerical_motion_is_reported_not_zeroed_or_called_exact():
    result=coordinate_delta([atom('CB','C',0.)],[atom('CB','C',1e-6)])
    assert result['heavy']['maximum_angstrom']==1e-6
    assert not result['heavy']['exact']


def test_identical_all_atom_state_requires_hydrogens_to_match_too():
    before=[atom('CB','C',0.),atom('HB','H',1.)]
    after=[atom('CB','C',0.),atom('HB','H',1.000001)]
    result=coordinate_delta(before,after)
    assert result['heavy']['exact'] and not result['hydrogen']['exact']
