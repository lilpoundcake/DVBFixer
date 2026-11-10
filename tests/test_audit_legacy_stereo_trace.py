"""Signed-volume readback identifies R/S crossings without reflecting coordinates."""
import pytest
from scripts.audit_legacy_stereo_trace import determinant


def test_determinant_detects_sign_and_preserves_small_nonzero_volume():
    points=dict(CB=(0.,0.,0.),CA=(1.,0.,0.),OG1=(0.,1.,0.),CG2=(0.,0.,.01))
    assert determinant(points,'CB','CA','OG1','CG2')==pytest.approx(.01)
    points['CG2']=(0.,0.,-.01)
    assert determinant(points,'CB','CA','OG1','CG2')==pytest.approx(-.01)


def test_nonfinite_stereo_snapshot_is_not_a_clean_state():
    with pytest.raises(ValueError,match='nonfinite'):
        determinant(dict(CB=(0.,0.,0.),CA=(1.,0.,0.),OG1=(0.,1.,0.),CG2=(0.,0.,float('nan'))),'CB','CA','OG1','CG2')
