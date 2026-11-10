"""Geometric invariants of private research controls, not backend qualification."""
import importlib
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def research(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/'scripts'))
    return importlib.import_module('research_chiral_completion')


def test_cb_construction_has_explicit_l_branch(research):
    n, ca, c = np.array([-0.14, 0., 0.]), np.zeros(3), np.array([0.05, 0.14, 0.])
    cb = research.construct_l_cb(n, ca, c)
    assert np.dot(np.cross(n-ca, c-ca), cb-ca) > 0
    assert np.linalg.norm(cb-ca) == pytest.approx(0.153)
    for anchor, expected in ((n, 110.5), (c, 110.1)):
        angle = np.degrees(np.arccos(np.dot(anchor-ca, cb-ca)/np.linalg.norm(anchor-ca)/np.linalg.norm(cb-ca)))
        assert angle == pytest.approx(expected)


def test_cb_is_equivariant_under_proper_motion(research):
    rng = np.random.default_rng(17)
    n, ca, c = np.array([-0.14, 0., 0.]), np.zeros(3), np.array([0.05, 0.14, 0.])
    cb = research.construct_l_cb(n, ca, c)
    for _ in range(100):
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        q[:, 0] *= np.linalg.det(q)
        shift = rng.normal(size=3)
        actual = research.construct_l_cb(n@q+shift, ca@q+shift, c@q+shift)
        np.testing.assert_allclose(actual, cb@q+shift, atol=1e-12)


@pytest.mark.parametrize('n,ca,c', [
    ([0,0,0], [0,0,0], [1,0,0]),
    ([1,0,0], [0,0,0], [-1,0,0]),
    ([float('nan'),0,0], [0,0,0], [1,0,0]),
    ([1,0,0], [0,0,0], [-0.98,0.1,0]),
])
def test_cb_refuses_invalid_or_incompatible_anchors(research, n, ca, c):
    with pytest.raises(ValueError):
        research.construct_l_cb(n, ca, c)


def test_proper_fit_cannot_silently_reflect(research):
    module = importlib.import_module('research_pdbfixer')
    template = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]])
    mirror = template*np.array([-1.,1.,1.])
    shift, rotation, center = module.proper_overlay(mirror, template)
    assert np.linalg.det(rotation) == pytest.approx(1)
    assert np.linalg.norm((template+shift)@rotation.T+center-mirror) > 0.5
