"""Independent mobility audit distinguishes protected atoms from generated tips."""
import copy

import pytest
from scripts.audit_generated_chi_relax import mobility_check, snapshot_atoms

CASE=dict(target_key=['A',29,' ','VAL'],generated=[['A',29,' ','VAL',n] for n in ('CB','CG1','CG2')])


def atoms():
    return [dict(key=CASE['target_key']+['CB'],xyz=[0.,0.,0.],mobile=False,is_hydrogen=False,is_virtual=False,is_backbone=False),
            dict(key=CASE['target_key']+['CG1'],xyz=[1.,0.,0.],mobile=True,is_hydrogen=False,is_virtual=False,is_backbone=False)]


def test_generated_cb_motion_fails_even_when_distal_motion_is_allowed():
    a=atoms()
    b=copy.deepcopy(a)
    b[1]['xyz'][0]+=1.
    assert mobility_check(a,b,CASE,1e-9)['mobile_max_motion_angstrom']==1.
    b[0]['xyz'][0]+=.001
    with pytest.raises(ValueError,match='protected atom moved'):
        mobility_check(a,b,CASE,1e-9)


def test_claimed_mobility_cannot_make_cb_or_backbone_h_free():
    a=atoms()
    a[0]['mobile']=True
    with pytest.raises(ValueError,match='unexpected mobility'):
        mobility_check(a,a,CASE,1e-9)


def test_all_atom_binding_reconstructs_hydrogens_without_renaming_residue_identity():
    snapshot=dict(heavy_atoms=[dict(chain='a',resid=29,icode='B',resname='VAL',atom='CA',xyz=[0.,0.,0.])],
        hydrogens=[dict(chain='a',resid=29,icode='B',atom='H',xyz=[1.,0.,0.])],virtual_atoms=[])
    result=snapshot_atoms(snapshot)
    assert result[('a',29,'B','VAL','H')]==[1.,0.,0.]
    assert ('A',29,'B','VAL','H') not in result
