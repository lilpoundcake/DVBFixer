"""Legacy torsions use their actual backbone, unlike fixed-source completion."""
import pytest
from scripts.amend_legacy_chi_reporting import actual_chi


def test_actual_chi_is_rigid_translation_invariant():
    points=dict(N=(0.,1.,0.),CA=(0.,0.,0.),CB=(1.,0.,0.),OG1=(1.,1.,1.))
    shifted={n:tuple(x+3 for x in xyz) for n,xyz in points.items()}
    assert actual_chi(points,'THR')==pytest.approx(actual_chi(shifted,'THR'))
    hybrid=dict(points,CB=shifted['CB'],OG1=shifted['OG1'])
    assert actual_chi(hybrid,'THR')!=pytest.approx(actual_chi(shifted,'THR'))
