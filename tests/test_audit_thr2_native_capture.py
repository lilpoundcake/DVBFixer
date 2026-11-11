"""Saved-only audit helpers: exact inventories and independent volume sign."""
import pytest
from scripts.audit_thr2_native_capture import delta, geometry, verify_files
from scripts.research_dunbrack_prior import sha


def test_geometry_tracks_volume_without_repair():
    pytest.importorskip('rdkit')
    xyz = dict(CA=[0., 0., 0.], N=[1., 0., 0.], C=[0., 1., 0.], CB=[0., 0., 1.],
        O=[0., 2., 0.], OG1=[1., 0., 1.], CG2=[0., 1., 1.])
    rows = [dict(key=['A', '2', '', 'THR', name]) for name in xyz]
    points = list(xyz.values())
    result = geometry(rows, points)
    assert result['target_ca_volume_angstrom3'] == 1000.
    assert result['target_cb_volume_angstrom3'] == -1000.
    assert result['negative_alpha'] == []
    reflected = [[-p[0], *p[1:]] for p in points]
    opposite = geometry(rows, reflected)
    assert opposite['target_ca_volume_angstrom3'] == -1000.
    assert len(opposite['negative_alpha']) == 1
    assert opposite['target_cip']['CB'] != result['target_cip']['CB']
    with pytest.raises(ValueError, match='duplicate'):
        geometry(rows+[rows[0]], points+[points[0]])


def test_delta_units_and_inventory():
    assert delta([[0., 0., 0.]], [[1., 0., 0.]])['maximum_angstrom'] == 10.
    with pytest.raises(ValueError, match='inventory'):
        delta([], [])


def test_changed_capture_is_rejected(tmp_path):
    path = tmp_path/'payload'
    path.write_bytes(b'initial')
    record = dict(artifacts={'payload':sha(path.read_bytes())})
    verify_files(tmp_path, record)
    path.write_bytes(b'changed')
    with pytest.raises(ValueError, match='changed'):
        verify_files(tmp_path, record)
