"""Synthetic context inventories distinguish target ownership from scorer state."""
import copy

import pytest
from scripts.audit_mean_comparability import context


def snapshot():
    atoms=[dict(chain='A',resid=1,icode=' ',resname='VAL',atom='CA',xyz=[0.,0.,0.]),
           dict(chain='A',resid=1,icode=' ',resname='VAL',atom='CG1',xyz=[1.,0.,0.])]
    return dict(heavy_atoms=atoms,virtual_atoms=[],hydrogens=[dict(chain='A',resid=1,icode=' ',atom='HG11',xyz=[2.,0.,0.]),
        dict(chain='A',resid=2,icode=' ',atom='H',xyz=[3.,0.,0.])],residues=[],external=[],target_internal=[])


def test_generated_target_coordinates_and_target_h_are_excluded_but_context_h_is_not():
    original=snapshot()
    changed=copy.deepcopy(original)
    generated=[['A',1,' ','VAL','CG1']]
    target=['A',1,' ','VAL']
    changed['heavy_atoms'][1]['xyz'][0]+=1.
    changed['hydrogens'][0]['xyz'][0]+=1.
    assert context(original,generated,target)==context(changed,generated,target)
    changed['hydrogens'][1]['xyz'][0]+=1.
    assert context(original,generated,target)!=context(changed,generated,target)


def test_known_backbone_motion_remains_visible():
    a=snapshot()
    b=copy.deepcopy(a)
    b['heavy_atoms'][0]['xyz'][0]+=.001
    assert context(a,[],['A',1,' ','VAL'])!=context(b,[],['A',1,' ','VAL'])


def test_duplicate_atom_identity_refuses():
    s=snapshot()
    s['heavy_atoms'].append(s['heavy_atoms'][0])
    with pytest.raises(ValueError,match='duplicate physical'):
        context(s,[],['A',1,' ','VAL'])
