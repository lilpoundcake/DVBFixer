"""Geometry-only source reconstruction retains chi and immutable anchors."""
import pytest
from scripts.research_start_geometry import regularize

from dvbfixer.protein_completion.branched import named_beta_volume
from dvbfixer.protein_completion.geometry import extend


def example():
    points=dict(N=(0.,1.45,0.),CA=(0.,0.,0.),C=(1.5,-.5,0.),O=(2.,-.5,.1),CB=(.3,-.6,.05))
    points['OG1']=extend(points['N'],points['CA'],points['CB'],1.6,80.,33.)
    points['CG2']=extend(points['N'],points['CA'],points['CB'],1.7,80.,153.)
    return points


@pytest.mark.parametrize('arm',['beta_only','full_generated'])
def test_chi_and_backbone_stay_fixed_while_new_atom_shape_changes(arm):
    points=example()
    rebuilt,chi=regularize(points,'THR',arm)
    assert chi==pytest.approx(33.)
    for name in ('N','CA','C','O'):
        assert rebuilt[name]==points[name]
    if arm=='beta_only':
        assert rebuilt['CB']==points['CB']
    assert named_beta_volume(rebuilt['CA'],rebuilt['CB'],rebuilt['OG1'],rebuilt['CG2'])>0


def test_raw_is_exact_noop_and_val_repair_arms_are_identity_controls():
    points=example()
    assert regularize(points,'THR','raw')[0]==points
    val={n:v for n,v in points.items() if n!='OG1'}
    val['CG1']=points['OG1']
    assert regularize(val,'VAL','full_generated')[0]==val
