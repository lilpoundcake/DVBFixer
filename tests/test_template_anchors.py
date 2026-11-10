"""Synthetic source-bound anchor selection, before sampling the frozen cases."""
import copy

import numpy as np
import pytest
from scripts import research_template_anchors as module
from scripts.research_template_causality import proper_overlay


def anchors():
    q = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [.2, .3, .8]])
    p = q+np.array([1., 2., 3.])
    return ['N', 'CA', 'C', 'O'], p, q


def result_arrays(result):
    return [np.asarray(r) for r in result]


def test_subset_excludes_o_from_fit_but_retains_supplied_anchor_transport_and_inputs():
    names, p, q = anchors()
    before = p.copy(), q.copy()
    changed = p.copy()
    changed[3] += [2., -4., 7.]
    calls = []
    first = module.make_anchor_overlay(names, p, q, 'ncac', proper_overlay, calls)(p, q)
    second = module.make_anchor_overlay(names, changed, q, 'ncac', proper_overlay, calls)(changed, q)
    assert all(np.array_equal(a, b) for a, b in zip(result_arrays(first), result_arrays(second), strict=True))
    assert calls[0]['supplied_names'] == names and calls[0]['used_names'] == ['N', 'CA', 'C']
    assert calls[0]['used_indices'] == [0, 1, 2]
    assert np.array_equal(p, before[0]) and np.array_equal(q, before[1])


def test_all_anchor_control_is_exact_original_fit_and_o_is_not_ignored():
    names, p, q = anchors()
    p[3] += [2., -4., 7.]
    calls = []
    actual = module.make_anchor_overlay(names, p, q, 'all', proper_overlay, calls)(p, q)
    expected = proper_overlay(p, q)
    subset = proper_overlay(p[:3], q[:3])
    assert all(np.array_equal(a, b) for a, b in zip(result_arrays(actual), result_arrays(expected), strict=True))
    assert not np.allclose(actual[1], subset[1])
    assert calls[0]['used_names'] == names


def test_anchor_order_is_identity_bound_not_assumed_to_be_first_three_positions():
    names, p, q = anchors()
    order = [3, 1, 0, 2]
    names, p, q = [names[i] for i in order], p[order], q[order]
    calls = []
    module.make_anchor_overlay(names, p, q, 'ncac', proper_overlay, calls)(p, q)
    assert calls[0]['used_indices'] == [1, 2, 3]
    assert calls[0]['used_names'] == ['CA', 'N', 'C']


def test_unbound_anchors_or_undeclared_policy_fail_instead_of_falling_back():
    names, p, q = anchors()
    calls = []
    hook = module.make_anchor_overlay(names, p, q, 'ncac', proper_overlay, calls)
    bad = copy.deepcopy(p)
    bad[0, 0] += .001
    with pytest.raises(ValueError, match='unexpected overlay'):
        hook(bad, q)
    with pytest.raises(ValueError, match='declared policy'):
        module.make_anchor_overlay(names, p, q, 'native', proper_overlay, calls)
    assert calls == []


def audit_fixture():
    names, p, q = anchors()
    p = p @ np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
    calls = []
    module.make_anchor_overlay(names, p, q, 'ncac', proper_overlay, calls)(p, q)
    call = calls[0]
    overlay = dict(reference_points_nm=p.tolist(), template_points_nm=q.tolist(),
                   rotation=call['rotation'], translate_template_nm=call['translate_template_nm'],
                   translate_reference_nm=call['translate_reference_nm'])
    return names, p, q, call, overlay


def test_independent_auditor_verifies_used_subset_centers_and_proper_fit_optimum():
    from scripts.audit_template_anchors import fit_check

    names, p, q, call, overlay = audit_fixture()
    result = fit_check(call, overlay, names, p, q, 'ncac')
    assert result['determinant'] == pytest.approx(1.)
    assert result['used_anchor_rmsd_angstrom'] < 1e-10


@pytest.mark.parametrize('kind', ['used_names', 'used_coordinates', 'centering', 'nonoptimal_rotation'])
def test_independent_auditor_rejects_hidden_or_incorrect_anchor_fit(kind):
    from scripts.audit_template_anchors import fit_check

    names, p, q, call, overlay = audit_fixture()
    if kind == 'used_names':
        call['used_names'].append('O')
    elif kind == 'used_coordinates':
        call['used_reference_points_nm'][0][0] += .1
    elif kind == 'centering':
        call['translate_reference_nm'][0] += .1
        overlay['translate_reference_nm'] = call['translate_reference_nm']
    else:
        call['rotation'] = np.eye(3).tolist()
        overlay['rotation'] = call['rotation']
    with pytest.raises(ValueError):
        fit_check(call, overlay, names, p, q, 'ncac')
