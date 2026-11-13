"""Ownership/return-transport controls without scientific native execution."""
from dataclasses import replace

import pytest
from scripts.hydrogen_initialization_guard import Atom, GuardRefusal, Input, propose
from scripts.research_post_hydrogen_guard import bind_return

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def source():
    c = ComponentInstanceRef(1, 'A', '2')
    atoms = tuple(Atom(ExactAtomRef(c, name), 'THR', element) for name, element in
                  [('CB', 'C'), ('CA', 'C'), ('OG1', 'O'), ('CG2', 'C'), ('HB', 'H')])
    return Input(atoms, ((0., 0., 0.), (0.14, 0., 0.01), (-0.07, 0.121, 0.01),
                         (-0.07, -0.121, 0.01), (0., 0., 0.1)),
                 ((0, 1), (0, 2), (0, 3), (0, 4)), frozenset([atoms[4].identity]))


def test_post_h_uses_returned_radius_and_retains_original_ownership():
    value = source()
    points = value.positions_nm[:4]+((0., 0., 0.11),)
    returned = bind_return(value, points)
    assert returned.new_hydrogens == value.new_hydrogens
    result = propose(returned)
    assert result.positions_nm[4] == (0., 0., -0.11)
    assert result.positions_nm[:4] == value.positions_nm[:4]
    assert value.positions_nm[4] == (0., 0., 0.1)


@pytest.mark.parametrize('observed_h', [False, True])
def test_moved_observations_in_native_return_refuse(observed_h):
    value = source()
    points = list(value.positions_nm)
    if observed_h:
        value = replace(value, new_hydrogens=frozenset())
        points[4] = (0., 0., 0.11)
    else:
        points[1] = (0.141, 0., 0.01)
    with pytest.raises(GuardRefusal, match='observation'):
        bind_return(value, points)


def test_return_graph_endpoint_and_atom_identity_are_not_ignored():
    value = source()
    with pytest.raises(GuardRefusal, match='graph'):
        bind_return(value, value.positions_nm, returned_bonds=((0, 1), (0, 2), (0, 3), (1, 4)))
    rows = [dict(key=['a', '2', '', a.residue_name, a.identity.name], element=a.element) for a in value.atoms]
    with pytest.raises(GuardRefusal, match='identity'):
        bind_return(value, value.positions_nm, returned_atoms=rows)


def test_bond_ordering_and_exact_noop_return_admit():
    value = source()
    returned = bind_return(value, value.positions_nm, returned_bonds=tuple((b, a) for a, b in reversed(value.bonds)))
    assert returned == value
    observed = replace(value, new_hydrogens=frozenset())
    assert propose(bind_return(observed, observed.positions_nm)).positions_nm == observed.positions_nm


def test_incomplete_or_nonfinite_return_refuses():
    value = source()
    with pytest.raises(GuardRefusal, match='inventory'):
        bind_return(value, value.positions_nm[:-1])
    with pytest.raises(GuardRefusal, match='finite'):
        bind_return(value, value.positions_nm[:4]+((float('nan'), 0., 0.),))
