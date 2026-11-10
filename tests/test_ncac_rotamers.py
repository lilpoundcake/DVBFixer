"""Synthetic relative rotations, ownership and stereo; no native coordinates."""
import copy
import math

import pytest
from scripts.audit_template_causality import volume
from scripts.research_ncac_rotamers import rotate_branches


def points():
    return dict(N=[-1.,0.,0.],CA=[0.,0.,0.],C=[0.,1.,0.],O=[0.,2.,0.],
                CB=[0.,0.,1.5],OG1=[1.,0.,2.],CG2=[0.,1.,2.])


@pytest.mark.parametrize('angle',[0.,120.,-120.])
def test_relative_rotation_preserves_axis_backbone_internal_distances_and_ca_volume(angle):
    p=points()
    before=copy.deepcopy(p)
    q=rotate_branches(p,angle)
    assert p==before
    for n in ('N','CA','C','O','CB'):
        assert q[n]==p[n]
    for a,b in [('CA','CB'),('CB','OG1'),('CB','CG2'),('OG1','CG2')]:
        assert math.dist(q[a],q[b])==pytest.approx(math.dist(p[a],p[b]),abs=1e-12)
    assert volume(q)==volume(p)
    if angle==0.:
        assert q==p


def test_both_branches_rotate_together_and_inverse_recovers_inputs():
    p=points()
    q=rotate_branches(p,120.)
    assert q['OG1']!=p['OG1'] and q['CG2']!=p['CG2']
    recovered=rotate_branches(q,-120.)
    for n in p:
        assert recovered[n]==pytest.approx(p[n],abs=1e-12)


def test_unsupported_angle_chemistry_and_degenerate_axis_refuse():
    with pytest.raises(ValueError,match='predeclared'):
        rotate_branches(points(),60.)
    p=points()
    p['CD1']=[1.,2.,3.]
    with pytest.raises(ValueError,match='VAL/THR'):
        rotate_branches(p,120.)
    p=points()
    p['CB']=p['CA']
    with pytest.raises(ValueError,match='degenerate'):
        rotate_branches(p,120.)


@pytest.mark.parametrize('angle',[0.,120.,-120.])
def test_independent_angle_readback_accepts_declared_rigid_rotation(angle):
    from scripts.audit_ncac_rotamers import check_rotation

    check_rotation(points(),rotate_branches(points(),angle),angle)


@pytest.mark.parametrize('kind',['axis','angle','one_branch'])
def test_independent_readback_rejects_ownership_or_transport_change(kind):
    from scripts.audit_ncac_rotamers import check_rotation

    p=points()
    q=rotate_branches(p,120.)
    if kind=='axis':
        q['CB'][0]+=.01
    elif kind=='angle':
        q=rotate_branches(p,-120.)
    else:
        q['CG2']=p['CG2']
    with pytest.raises(ValueError):
        check_rotation(p,q,120.)
