"""Continuous focal volume and H-angle readback preserves full residue identity."""
import pytest
from scripts.audit_methylene_whole_control import site_geometry


def source():
    target = ['a','82','A','ALA']
    rows = [dict(key=[*target,n],element=e) for n,e in [('CA','C'),('N','N'),('C','C'),('CB','C'),('HA','H')]]
    points = [[0.,0.,0.],[.1,0.,0.],[0.,.1,0.],[0.,0.,.1],[-.06,-.06,-.06]]
    return rows,points,target


def test_focal_shape_is_continuous_and_uses_nm_to_Angstrom_conversion():
    rows,points,target = source()
    result = site_geometry(rows,points,target)
    assert result['residue'] == ['a','82','A','ALA']
    assert result['CA_volume_angstrom3'] == pytest.approx(1.)
    assert result['outward_cosine'] == pytest.approx(1.)
    assert result['CA_H_angstrom'] == pytest.approx(3**.5*.6)
    assert all(angle > 90 for angle in result['HA_angles_degrees'].values())


def test_wrong_chain_case_or_missing_insertion_code_is_not_a_focal_match():
    rows,points,target = source()
    with pytest.raises(ValueError):
        site_geometry(rows,points,['A',*target[1:]])
    with pytest.raises(ValueError):
        site_geometry(rows,points,['a','82','','ALA'])
