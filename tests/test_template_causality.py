"""Synthetic pre-execution proper-fit and factorial-enrollment controls."""
import json

import numpy as np
import pytest
from scripts import research_template_causality as module


def test_proper_fit_refuses_reflection_without_modifying_input_anchors():
    q = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    p = q*np.array([-1., 1., 1.])+np.array([2., 3., 4.])
    before = p.copy(), q.copy()
    t, rotation, center = module.proper_overlay(p, q)
    assert np.linalg.det(rotation) == pytest.approx(1.)
    assert np.allclose(rotation.T @ rotation, np.eye(3))
    assert not np.allclose((rotation @ (q+t).T).T+center, p)
    assert np.array_equal(p, before[0]) and np.array_equal(q, before[1])


def test_proper_fit_exact_for_a_proper_rigid_transform_and_signed_volume_preserved():
    q = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    p = (rotation @ q.T).T+np.array([2., 3., 4.])
    t, r, center = module.proper_overlay(p, q)
    assert np.allclose((r @ (q+t).T).T+center, p)
    names = ['CA', 'N', 'C', 'CB']
    assert module.signed_volume(dict(zip(names, q, strict=True))) == pytest.approx(module.signed_volume(dict(zip(names, p, strict=True))))
    with pytest.raises(ValueError, match='at least three'):
        module.proper_overlay(p[:2], q[:2])


def test_enrollment_keeps_all_four_case_roles_and_full_factorial_without_execution(tmp_path, monkeypatch):
    root = tmp_path/'parent'
    root.mkdir()
    (root/'evidence.json').write_text(json.dumps({'environment': {}}))
    target = dict(model=1, chain_id='I', sequence_number='28', insertion_code='', alternate_location='', occurrence=1)
    for name in module.CASES:
        folder = root/name
        folder.mkdir()
        (folder/'input.pdb').write_text('HEADER    SYNTHETIC PLAN ONLY\nEND\n')
        (folder/'source-plan.json').write_text(json.dumps({'request': {'residue': target, 'residue_name': 'VAL'}}))
        (folder/'template-proposal.json').write_text(json.dumps({'snapshots': []}))
    monkeypatch.setattr(module, 'PARENT', root)
    monkeypatch.setattr(module, 'preflight', lambda: {})

    def never(*args, **kwargs):
        raise AssertionError('enrollment must not rebuild or run a control')

    monkeypatch.setattr(module, 'run_arm', never)
    monkeypatch.setattr(module, 'rebuild_missing_atoms_with_retry', never)
    output = tmp_path/'plan.json'
    module.plan(output)
    record = json.loads(output.read_bytes())
    assert len(record['cases']) == len(record['arms']) == 4
    assert {tuple(a.values()) for a in record['arms']} == {('stock', 'normal'), ('stock', 'none'), ('proper', 'normal'), ('proper', 'none')}
    assert not record['public_capability']
    with pytest.raises(FileExistsError):
        module.plan(output)


def test_independent_scalar_geometry_replays_signed_volumes_and_rejects_tampering():
    from scripts.audit_template_causality import audit_points, determinant, volume

    points = {'CA': [0., 0., 0.], 'N': [-1., 0., 0.], 'C': [0., 1., 0.], 'CB': [0., 0., -1.]}
    assert volume(points) == module.signed_volume(points) == 1.
    assert determinant([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]) == 1.
    assert audit_points(points, 1.) == 1.
    with pytest.raises(ValueError, match='signed volume differs'):
        audit_points(points, -1.)


def test_proper_full_template_handedness_does_not_guarantee_hybrid_handedness():
    from scripts.audit_template_causality import determinant, volume

    template = {'CA': [0., 0., 0.], 'N': [-1., 0., 0.], 'C': [0., 1., 0.], 'CB': [0., 0., -1.]}
    proper = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    hybrid = dict(template, N=[1., 0., 0.])
    assert determinant(proper) == 1. and volume(template) > 0 and volume(hybrid) < 0


def test_independent_cip_replay_cannot_accept_a_relabelled_ca():
    from scripts.audit_template_causality import audit_points

    from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE
    from dvbfixer.protein_completion.geometry import extend, l_cb

    points = {'N': (-1.46, 0., 0.), 'CA': (0., 0., 0.), 'C': (.55, 1.43, 0.), 'O': (1.76, 1.66, 0.)}
    points['CB'] = l_cb(points['N'], points['CA'], points['C'])
    branch = BRANCHES['VAL']
    points['CG1'] = extend(points['N'], points['CA'], points['CB'], branch.first_length, TETRAHEDRAL_ANGLE, 60.)
    points['CG2'] = extend(points['N'], points['CA'], points['CB'], branch.second_length, TETRAHEDRAL_ANGLE, 180.)
    value = module.signed_volume(points)
    audit_points(points, value, 'VAL', {'CA': 'S'})
    with pytest.raises(ValueError, match='CIP differs'):
        audit_points(points, value, 'VAL', {'CA': 'R'})
