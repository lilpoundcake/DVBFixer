"""Mobility ownership and failure-on-motion, without PyRosetta/minimization."""
import copy

import pytest
from scripts.research_generated_chi_relax import is_mobile, movement

TARGET=['A',29,' ','VAL']
GENERATED=[TARGET+[n] for n in ('CB','CG1','CG2')]


def test_generated_cb_stays_protected_while_distal_atoms_can_move():
    assert not is_mobile(TARGET+['CB'],GENERATED,TARGET,False,False,False)
    assert is_mobile(TARGET+['CG1'],GENERATED,TARGET,False,False,False)
    assert not is_mobile(TARGET+['CA'],GENERATED,TARGET,False,True,False)


def test_only_target_sidechain_h_is_mobile_and_virtuals_are_protected():
    assert is_mobile(TARGET+['HG11'],GENERATED,TARGET,True,False,False)
    assert not is_mobile(TARGET+['H'],GENERATED,TARGET,True,True,False)
    assert not is_mobile(['A',30,' ','PHE','HD1'],GENERATED,TARGET,True,False,False)
    assert not is_mobile(TARGET+['V1'],GENERATED,TARGET,False,False,True)


def atoms():
    return [dict(key=TARGET+['CB'],xyz=[0.,0.,0.],mobile=False,is_hydrogen=False,is_virtual=False,is_backbone=False),
            dict(key=TARGET+['CG1'],xyz=[1.,0.,0.],mobile=True,is_hydrogen=False,is_virtual=False,is_backbone=False)]


def test_motion_is_reported_not_undone_or_hidden():
    before=atoms()
    after=copy.deepcopy(before)
    after[0]['xyz'][0]+=.01
    after[1]['xyz'][0]+=1.
    result=movement(before,after)
    assert result==dict(protected_max_motion_angstrom=.01,mobile_max_motion_angstrom=1.)
    assert after[0]['xyz'][0]==.01


def test_disappearing_atom_or_changed_mobility_refuses():
    with pytest.raises(ValueError,match='inventory changed'):
        movement(atoms(),atoms()[:1])
    a=atoms()
    b=copy.deepcopy(a)
    b[0]['mobile']=True
    with pytest.raises(ValueError,match='classification changed'):
        movement(a,b)
