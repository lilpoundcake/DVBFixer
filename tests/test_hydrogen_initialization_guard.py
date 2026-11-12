"""Private proposal engineering only; no native sampling or minimization."""
from dataclasses import replace

import pytest
from scripts.hydrogen_initialization_guard import Atom, GuardRefusal, Input, propose

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def graph(*, observed=False):
    component = ComponentInstanceRef(1, 'A', '2')
    atoms = tuple(Atom(ExactAtomRef(component, name), 'THR', element)
                  for name, element in [('CB', 'C'), ('CA', 'C'), ('OG1', 'O'), ('CG2', 'C'), ('HB', 'H')])
    # A shallow heavy tetrahedron: the H near +z is on the wrong side.
    positions = ((0., 0., 0.), (0.14, 0., 0.01), (-0.07, 0.121, 0.01),
                 (-0.07, -0.121, 0.01), (0., 0., 0.1))
    return Input(atoms, positions, ((0, 1), (0, 2), (0, 3), (0, 4)),
                 frozenset() if observed else frozenset([atoms[4].identity]))


def test_only_new_h_moves_preserving_radius_and_all_input_arrays():
    value = graph()
    result = propose(value)
    assert result.positions_nm[:4] == value.positions_nm[:4]
    assert result.positions_nm[4] == (0., 0., -0.1)
    assert result.moved == result.supported == (value.atoms[4].identity,)
    assert not result.unsupported
    assert value.positions_nm[4] == (0., 0., 0.1)


def test_valid_is_exact_noop_and_idempotent():
    value = graph()
    good = replace(value, positions_nm=propose(value).positions_nm)
    assert propose(good).positions_nm == good.positions_nm
    assert not propose(good).moved


def test_observed_bad_h_is_never_moved():
    value = graph(observed=True)
    assert propose(value).positions_nm == value.positions_nm
    assert not propose(value).supported


def test_proper_frame_covariance_all_24_rotations_and_translation():
    from scripts.research_thr2_initializer_frames import inverse, rotate, rotations

    value = graph()
    expected = propose(value).positions_nm
    shift = [1.25, -0.75, 0.5]
    for matrix in rotations():
        positions = tuple(tuple(a+b for a, b in zip(rotate(matrix, p), shift, strict=True)) for p in value.positions_nm)
        result = propose(replace(value, positions_nm=positions))
        for p, q in zip(result.positions_nm, expected, strict=True):
            assert inverse(matrix, p, shift) == pytest.approx(q, abs=1e-14)


def test_atom_and_bond_order_independence():
    value = graph()
    order = [4, 3, 1, 0, 2]
    inverse = {old:new for new, old in enumerate(order)}
    shuffled = replace(value, atoms=tuple(value.atoms[i] for i in order),
                       positions_nm=tuple(value.positions_nm[i] for i in order),
                       bonds=tuple((inverse[b], inverse[a]) for a, b in reversed(value.bonds)))
    expected = dict(zip((a.identity for a in value.atoms), propose(value).positions_nm, strict=True))
    assert dict(zip((a.identity for a in shuffled.atoms), propose(shuffled).positions_nm, strict=True)) == expected


@pytest.mark.parametrize('field', ['chain', 'icode', 'occurrence'])
def test_distinct_occurrences_chains_and_insertion_codes_not_merged(field):
    value = graph()
    component = value.atoms[0].identity.component
    other = replace(component, **{'chain_id':'a'} if field == 'chain' else
                    {'insertion_code':'A'} if field == 'icode' else {'occurrence':2})
    atoms = value.atoms+tuple(replace(a, identity=ExactAtomRef(other, a.identity.name)) for a in value.atoms)
    combined = Input(atoms, value.positions_nm+tuple((x+1., y, z) for x, y, z in value.positions_nm),
                     value.bonds+tuple((a+5, b+5) for a, b in value.bonds), value.new_hydrogens)
    result = propose(combined)
    assert result.positions_nm[5:] == combined.positions_nm[5:]
    assert result.moved == (value.atoms[4].identity,)


def test_multi_h_and_nonprotein_are_explicit_unsupported_noops():
    value = graph()
    heterogen = replace(value, atoms=tuple(replace(a, residue_name='LIG') for a in value.atoms))
    assert propose(heterogen).unsupported == tuple(value.new_hydrogens)
    sibling = replace(value.atoms[4], identity=ExactAtomRef(value.atoms[0].identity.component, 'HB2'))
    multi = replace(value, atoms=value.atoms+(sibling,), positions_nm=value.positions_nm+((0.1, 0., 0.),),
                    bonds=value.bonds+((0, 5),))
    assert propose(multi).positions_nm == multi.positions_nm
    assert propose(multi).unsupported == tuple(value.new_hydrogens)


@pytest.mark.parametrize('kind', ['heavy-planar', 'h-planar', 'distance', 'overlap'])
def test_refusals_are_atomic_and_preserve_input(kind):
    value = graph()
    points = list(value.positions_nm)
    if kind == 'heavy-planar':
        points[1:4] = [(x, y, 0.) for x, y, _ in points[1:4]]
    elif kind == 'h-planar':
        points[4] = (0.1, 0., 0.01)
    elif kind == 'distance':
        points[4] = (0., 0., 0.3)
    else:
        obstacle = Atom(ExactAtomRef(ComponentInstanceRef(1, 'B', '1'), 'O'), 'HOH', 'O')
        value = replace(value, atoms=value.atoms+(obstacle,), positions_nm=value.positions_nm+((0., 0., -0.1),))
        points = list(value.positions_nm)
    value = replace(value, positions_nm=tuple(points))
    before = value.positions_nm
    with pytest.raises(GuardRefusal):
        propose(value)
    assert value.positions_nm == before


@pytest.mark.parametrize('kind', ['duplicate', 'nan', 'mask-heavy', 'mask-absent', 'duplicate-bond', 'self-bond'])
def test_invalid_transport_refuses(kind):
    value = graph()
    if kind == 'duplicate':
        value = replace(value, atoms=value.atoms[:4]+(value.atoms[0],))
    elif kind == 'nan':
        value = replace(value, positions_nm=value.positions_nm[:4]+((float('nan'), 0., 0.),))
    elif kind == 'mask-heavy':
        value = replace(value, new_hydrogens=frozenset([value.atoms[0].identity]))
    elif kind == 'mask-absent':
        value = replace(value, new_hydrogens=frozenset([ExactAtomRef(ComponentInstanceRef(1, 'Z', '1'), 'H')]))
    else:
        value = replace(value, bonds=value.bonds+(((1, 0),) if kind == 'duplicate-bond' else ((0, 0),)))
    with pytest.raises(GuardRefusal):
        propose(value)


def test_atomic_refusal_after_an_earlier_valid_proposal():
    value = graph()
    other = ComponentInstanceRef(1, 'B', '2')
    atoms = value.atoms+tuple(replace(a, identity=ExactAtomRef(other, a.identity.name)) for a in value.atoms)
    points = value.positions_nm+tuple((x+1., y, 0. if i in (1, 2, 3) else z)
                                    for i, (x, y, z) in enumerate(value.positions_nm))
    combined = Input(atoms, points, value.bonds+tuple((a+5, b+5) for a, b in value.bonds),
                     value.new_hydrogens | frozenset([atoms[9].identity]))
    with pytest.raises(GuardRefusal, match='degenerate'):
        propose(combined)
    assert combined.positions_nm == points


def test_native_mass_admission_requires_observed_h_pinned():
    openmm = pytest.importorskip('openmm')
    from scripts.qualify_hydrogen_initialization_guard import admit

    value = graph()
    rows = [dict(key=['A', '2', '', a.residue_name, a.identity.name], element=a.element,
                 xyz=[10*x for x in p]) for a, p in zip(value.atoms, value.positions_nm, strict=True)]
    before = dict(atoms=rows[:4], positions_nm=[list(p) for p in value.positions_nm[:4]])
    system = openmm.System()
    for mass in (0., 0., 0., 0., 1.):
        system.addParticle(mass)
    admitted = admit(rows, before, [list(p) for p in value.positions_nm], value.bonds, system)
    assert admitted.new_hydrogens == value.new_hydrogens
    before['atoms'] = rows
    before['positions_nm'] = [list(p) for p in value.positions_nm]
    with pytest.raises(GuardRefusal, match='ownership'):
        admit(rows, before, [list(p) for p in value.positions_nm], value.bonds, system)
    system.setParticleMass(4, 0.)
    assert not admit(rows, before, [list(p) for p in value.positions_nm], value.bonds, system).new_hydrogens


def test_adapter_does_not_guess_legacy_repeated_atom_keys():
    from scripts.qualify_hydrogen_initialization_guard import atoms_from_rows

    rows = [dict(key=['A', '2', '', 'THR', 'CB'], element='C')]*2
    with pytest.raises(GuardRefusal, match='ambiguous'):
        atoms_from_rows(rows)
    rows = [dict(rows[0], chain_index=0, residue_index=0), dict(rows[0], chain_index=1, residue_index=1)]
    atoms = atoms_from_rows(rows)
    assert atoms[0].identity != atoms[1].identity
    assert [a.identity.component.occurrence for a in atoms] == [1, 2]
