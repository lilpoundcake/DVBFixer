"""Synthetic controls only; retained scientific states are enrolled separately."""
import random
from dataclasses import replace

import pytest
from scripts.hydrogen_initialization_guard import Atom, GuardRefusal, Input
from scripts.hydrogen_initialization_guard import propose as strict
from scripts.post_hydrogen_local_noop import propose

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def center(radius=0.11, compatible=False, resid='2', offset=0.):
    component = ComponentInstanceRef(1, 'A', resid)
    atoms = tuple(Atom(ExactAtomRef(component, n), 'THR', e) for n, e in
                  [('CB', 'C'), ('CA', 'C'), ('OG1', 'O'), ('CG2', 'C'), ('HB', 'H')])
    points = ((0., 0., 0.), (0.14, 0., 0.01), (-0.07, 0.121, 0.01),
              (-0.07, -0.121, 0.01), (0., 0., -radius if compatible else radius))
    return Input(atoms, tuple((x+offset, y, z) for x, y, z in points),
                 ((0, 1), (0, 2), (0, 3), (0, 4)), frozenset([atoms[4].identity]))


def mixed(a, b):
    n = len(a.atoms)
    return Input(a.atoms+b.atoms, a.positions_nm+b.positions_nm,
                 a.bonds+tuple((x+n, y+n) for x, y in b.bonds), a.new_hydrogens | b.new_hydrogens)


def test_compatible_elongated_is_exact_noop_with_visible_warning_not_strict_pass():
    value = center(0.165, True)
    with pytest.raises(GuardRefusal, match='distance'):
        strict(value)
    result = propose(value)
    assert result.result.positions_nm is value.positions_nm or result.result.positions_nm == value.positions_nm
    assert not result.result.moved and not result.intervention_hydrogens
    assert result.original_new_hydrogens == value.new_hydrogens
    assert len(result.warnings) == 1
    assert result.warnings[0].radius_nm == 0.165
    assert 'not_geometry_certified' in result.warnings[0].reason


@pytest.mark.parametrize('compatible', [False, True])
def test_admitted_states_match_original_exactly(compatible):
    value = center(compatible=compatible)
    assert propose(value).result == strict(value)


@pytest.mark.parametrize('radius', [0.049, 0.151, 0.21])
def test_incompatible_distance_refuses_without_partial_output(radius):
    with pytest.raises(GuardRefusal, match='distance'):
        propose(center(radius))


def test_mixed_keeps_compatible_warning_and_strictly_moves_only_conflict():
    value = mixed(center(), center(0.165, True, '3', 1.))
    proposal = propose(value)
    assert proposal.result.moved == (value.atoms[4].identity,)
    assert proposal.result.positions_nm[5:] == value.positions_nm[5:]
    assert len(proposal.warnings) == 1
    assert proposal.original_new_hydrogens == value.new_hydrogens


@pytest.mark.parametrize('problem', ['flat-heavy', 'flat-four', 'coincident'])
def test_unclassifiable_sign_is_not_compatible_noop(problem):
    value = center(0.165, True)
    p = list(value.positions_nm)
    if problem == 'flat-heavy':
        p[1:4] = [(x, y, 0.) for x, y, z in p[1:4]]
    elif problem == 'flat-four':
        p[4] = (0., 0., 0.01)
    else:
        p[4] = p[1]
    with pytest.raises(GuardRefusal, match='degenerate'):
        propose(replace(value, positions_nm=tuple(p)))


def test_contact_check_keeps_untouched_compatible_center_in_full_inventory():
    # Second center's heavy parent is exactly at the first center's proposal.
    value = mixed(center(), center(0.165, True, '3'))
    p = value.positions_nm[:5]+tuple((x, y, z-0.11) for x, y, z in value.positions_nm[5:])
    with pytest.raises(GuardRefusal, match='overlaps'):
        propose(replace(value, positions_nm=p))


def test_existing_h_is_immutable_even_when_it_is_incompatible():
    value = replace(center(), new_hydrogens=frozenset())
    assert propose(value).result.positions_nm == value.positions_nm
    assert not propose(value).result.moved


def test_alias_supported_and_source_and_rng_unchanged():
    value = center()
    # Canonical SER CA also has the same three-neighbor ownership pattern.
    renames = {'CB': 'CA', 'CA': 'N', 'OG1': 'C', 'CG2': 'CB', 'HB': 'HA'}
    atoms = tuple(replace(a, identity=ExactAtomRef(a.identity.component, renames[a.identity.name]),
                          residue_name='SER') for a in value.atoms)
    value = replace(value, atoms=atoms, new_hydrogens=frozenset([atoms[4].identity]))
    before = repr(value)
    rng = random.getstate()
    assert propose(value) == propose(value)
    assert repr(value) == before and random.getstate() == rng


def test_shared_protonation_alias_is_supported():
    value = center()
    renames = {'CB': 'CA', 'CA': 'N', 'OG1': 'C', 'CG2': 'CB', 'HB': 'HA'}
    atoms = tuple(replace(a, identity=ExactAtomRef(a.identity.component, renames[a.identity.name]),
                          residue_name='GLH') for a in value.atoms)
    aliased = replace(value, atoms=atoms, new_hydrogens=frozenset([atoms[4].identity]))
    assert propose(aliased).result == strict(aliased)


@pytest.mark.parametrize('corruption', ['observation', 'warning', 'mask', 'projection'])
def test_independent_audit_rejects_corrupt_proposals(corruption):
    from scripts.audit_post_h_local_noop import audit

    value = mixed(center(), center(0.165, True, '3', 1.))
    proposal = propose(value)
    assert audit(value, proposal)['passed']
    if corruption == 'warning':
        bad = replace(proposal, warnings=())
    elif corruption == 'mask':
        bad = replace(proposal, intervention_hydrogens=frozenset())
    else:
        p = list(proposal.result.positions_nm)
        i = 0 if corruption == 'observation' else 4
        x, y, z = p[i]
        p[i] = (x+0.001, y, z)
        bad = replace(proposal, result=replace(proposal.result, positions_nm=tuple(p)))
    with pytest.raises(ValueError, match='audit'):
        audit(value, bad)


def test_proper_rotation_translation_preserves_policy_and_projection():
    import numpy as np

    value = mixed(center(), center(0.165, True, '3', 1.))
    rotation = np.array([[0., -1., 0.], [0., 0., 1.], [-1., 0., 0.]])
    assert np.linalg.det(rotation) == 1
    shift = np.array([1., -2., 3.])
    moved = replace(value, positions_nm=tuple(tuple(p) for p in np.asarray(value.positions_nm)@rotation+shift))
    a, b = propose(value), propose(moved)
    assert a.result.moved == b.result.moved
    assert np.asarray(b.result.positions_nm) == pytest.approx(np.asarray(a.result.positions_nm)@rotation+shift, abs=1e-14)
    assert [w.atom for w in a.warnings] == [w.atom for w in b.warnings]
