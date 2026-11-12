"""Saved initializer formula and distinct centre/ligand geometry observables."""
import math

import pytest
from scripts.audit_thr2_h_initialization import initializer, target_geometry, tuple_tree


def test_formula_normalization_and_refusal():
    parent = [0., 0., 0.]
    base, perturbed = initializer(parent, [[1., 0., 0.]], [0., 0.5, 0.5])
    assert base == [-0.1, 0., 0.]
    assert math.sqrt(sum(x*x for x in perturbed)) == pytest.approx(0.1)
    assert tuple_tree([1, [2, 3], None]) == (1, (2, 3), None)
    with pytest.raises(ValueError, match='branch'):
        initializer(parent, [], [0., 0., 0.])


def test_hb_alone_can_conflict_with_unchanged_heavy_volume():
    xyz = dict(CB=[0., 0., 0.], CA=[1., 0., 0.], OG1=[0., 1., 0.], CG2=[0., 0., 1.], HB=[-1., -1., -1.])
    rows = [dict(key=['A', '2', '', 'THR', name], xyz=point) for name, point in xyz.items()]
    good = target_geometry(rows)
    rows[-1] = dict(rows[-1], xyz=[1., 1., 1.])
    wrong = target_geometry(rows)
    assert good['heavy_center_volume_angstrom3'] == wrong['heavy_center_volume_angstrom3']
    assert not good['heavy_vs_four_ligand_sign_conflict']
    assert wrong['heavy_vs_four_ligand_sign_conflict']
