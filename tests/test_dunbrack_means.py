"""Synthetic rigid branch transport to arbitrary chi offsets, no library data."""
import copy
import math

import pytest
from scripts.research_dunbrack_means import rotate_distal

from dvbfixer.diagnose.geometry import signed_dihedral_degrees


def points():
    return dict(N=[-1.,0.,0.],CA=[0.,0.,0.],C=[0.,1.,0.],O=[0.,2.,0.],
                CB=[0.,0.,1.5],OG1=[1.,0.,2.],CG2=[0.,1.,2.])


@pytest.mark.parametrize('delta',[-173.1,-61.7,0.,42.8,179.9])
def test_arbitrary_mean_transport_preserves_cb_observations_and_branch_geometry(delta):
    p=points()
    before=copy.deepcopy(p)
    q=rotate_distal(p,delta)
    assert p==before
    for n in ('N','CA','C','O','CB'):
        assert q[n]==p[n]
    for a,b in [('CB','OG1'),('CB','CG2'),('OG1','CG2')]:
        assert math.dist(p[a],p[b])==pytest.approx(math.dist(q[a],q[b]),abs=1e-12)
    old=signed_dihedral_degrees(*(p[n] for n in ('N','CA','CB','OG1')))
    new=signed_dihedral_degrees(*(q[n] for n in ('N','CA','CB','OG1')))
    assert abs((new-old-delta+180.)%360.-180.)<1e-10


def test_nonfinite_rotation_refuses_instead_of_creating_nan_coordinates():
    with pytest.raises(ValueError,match='finite'):
        rotate_distal(points(),float('nan'))
