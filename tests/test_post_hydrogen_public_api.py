"""Fake public calls test transaction/ownership; no H sampling/minimization."""
import random

import pytest
from scripts.post_hydrogen_public_api import PostHGuardError, apply, snapshot


def model():
    from openmm import Vec3
    from openmm.app import Modeller, Topology, element
    from openmm.unit import nanometer

    topology = Topology()
    residue = topology.addResidue('THR', topology.addChain('A'), '2')
    atoms = [topology.addAtom(n, element.Element.getBySymbol(e), residue)
             for n, e in [('CB', 'C'), ('CA', 'C'), ('OG1', 'O'), ('CG2', 'C')]]
    for neighbor in atoms[1:]:
        topology.addBond(atoms[0], neighbor)
    return Modeller(topology, [Vec3(*p) for p in ((0., 0., 0.), (0.14, 0., 0.01),
                    (-0.07, 0.121, 0.01), (-0.07, -0.121, 0.01))]*nanometer)


def add_fake(working, *args, **kwargs):
    from openmm import Vec3
    from openmm.app import element
    from openmm.unit import nanometer

    atoms = list(working.topology.atoms())
    h = working.topology.addAtom('HB', element.hydrogen, atoms[0].residue)
    working.topology.addBond(atoms[0], h)
    working.positions = list(working.positions.value_in_unit(nanometer))+[Vec3(0., 0., 0.11)]
    working.positions *= nanometer
    return ['unchanged-public-return-value']


def test_guard_after_api_preserves_radius_observations_and_return_value():
    value = model()
    before = snapshot(value)
    events = []
    out = apply(value, add_fake, protected_hydrogens=frozenset(), journal=events.append)
    assert out.returned == ['unchanged-public-return-value']
    assert snapshot(value).positions_nm[:4] == before.positions_nm
    assert snapshot(value).positions_nm[4] == (0., 0., -0.11)
    assert events[0]['raw'].positions_nm[4] == (0., 0., 0.11)


def test_observational_baseline_does_not_apply_proposal():
    value = model()
    out = apply(value, add_fake, protected_hydrogens=frozenset(), enabled=False)
    assert snapshot(value).positions_nm[4] == (0., 0., 0.11)
    assert len(out.result.moved) == 1


def test_api_exception_does_not_modify_original_or_hide_error():
    value = model()
    before = snapshot(value)

    def failed_api(working):
        add_fake(working)
        raise ValueError('native template failure')

    with pytest.raises(ValueError, match='template failure'):
        apply(value, failed_api, protected_hydrogens=frozenset())
    assert snapshot(value) == before


@pytest.mark.parametrize('problem', ['deleted', 'new-heavy', 'moved-observation', 'graph'])
def test_ownership_refusal_is_typed_and_atomic(problem):
    from openmm import Vec3
    from openmm.app import element
    from openmm.unit import nanometer

    value = model()
    before = snapshot(value)

    def bad_api(working):
        add_fake(working)
        atoms = list(working.topology.atoms())
        if problem == 'deleted':
            working.delete([atoms[2]])
        elif problem == 'new-heavy':
            working.topology.addAtom('XX', element.carbon, atoms[0].residue)
            working.positions = list(working.positions.value_in_unit(nanometer))+[Vec3(1., 1., 1.)]
            working.positions *= nanometer
        elif problem == 'moved-observation':
            working.positions[1] = Vec3(1., 1., 1.)*nanometer
        else:
            working.topology.addBond(atoms[1], atoms[2])

    with pytest.raises(PostHGuardError):
        apply(value, bad_api, protected_hydrogens=frozenset())
    assert snapshot(value) == before
    assert not issubclass(PostHGuardError, (ValueError, KeyError))


def test_missing_stripped_original_h_refuses_before_native_call():
    from dvbfixer.domain.structure_identity import ExactAtomRef

    value = model()
    atom = ExactAtomRef(snapshot(value).atoms[0].identity.component, 'HB')
    called = []
    with pytest.raises(PostHGuardError, match='before API'):
        apply(value, lambda w: called.append(w), protected_hydrogens=frozenset([atom]))
    assert not called


def test_existing_h_is_not_new_and_valid_noop_keeps_it_exact():
    value = model()
    add_fake(value)
    before = snapshot(value)
    protected = frozenset([before.atoms[4].identity])
    out = apply(value, lambda w: 42, protected_hydrogens=protected)
    assert out.returned == 42
    assert not out.result.moved
    assert snapshot(value) == before


def test_existing_h_deletion_is_not_relabelled_as_new():
    value = model()
    add_fake(value)
    before = snapshot(value)

    def deletes_h(w):
        w.delete([list(w.topology.atoms())[-1]])

    with pytest.raises(PostHGuardError):
        apply(value, deletes_h, protected_hydrogens=frozenset([before.atoms[-1].identity]))
    assert snapshot(value) == before


def test_adapter_adds_no_random_draws():
    state = random.getstate()
    try:
        random.seed(17)
        expected = random.Random(17)

        def sampling_api(w):
            assert random.random() == expected.random()
            return add_fake(w)

        apply(model(), sampling_api, protected_hydrogens=frozenset())
        assert random.getstate() == expected.getstate()
    finally:
        random.setstate(state)


def test_nonempty_residue_templates_refuse_before_api_and_journal_failure_is_atomic():
    value = model()
    before = snapshot(value)
    with pytest.raises(PostHGuardError, match='residueTemplates'):
        apply(value, lambda w: pytest.fail('must not call'), protected_hydrogens=frozenset(), residueTemplates={object():'ALA'})

    def failed_sink(event):
        raise OSError('diagnostic disk full')

    with pytest.raises(OSError):
        apply(value, add_fake, protected_hydrogens=frozenset(), journal=failed_sink)
    assert snapshot(value) == before


def test_duplicate_return_atom_is_a_geometry_error_not_ff_fallback():
    from openmm import Vec3
    from openmm.app import element
    from openmm.unit import nanometer

    value = model()
    before = snapshot(value)

    def duplicate(w):
        add_fake(w)
        residue = list(w.topology.residues())[0]
        w.topology.addAtom('HB', element.hydrogen, residue)
        w.positions = list(w.positions.value_in_unit(nanometer))+[Vec3(0., 0., -0.11)]
        w.positions *= nanometer

    with pytest.raises(PostHGuardError, match='duplicate'):
        apply(value, duplicate, protected_hydrogens=frozenset())
    assert snapshot(value) == before


def test_canonical_variant_alias_does_not_invent_a_new_observed_atom():
    value = model()
    residue = list(value.topology.residues())[0]
    residue.name = 'CYS'
    # Keep this a pure ownership no-op, not a claim about a THR-shaped CYS graph.
    def aliases(w):
        list(w.topology.residues())[0].name = 'CYX'
        return 'alias'

    out = apply(value, aliases, protected_hydrogens=frozenset())
    assert out.returned == 'alias' and not out.result.moved
    assert list(value.topology.residues())[0].name == 'CYX'


def test_explicit_reference_observation_transport_is_exact():
    from openmm import Platform, Vec3
    from openmm.unit import nanometer

    value = model()
    before = snapshot(value)

    def rounded(w, **kwargs):
        add_fake(w)
        w.positions[1] = Vec3(0.140000001, 0., 0.01)*nanometer

    with pytest.raises(PostHGuardError, match='exact'):
        apply(value, rounded, protected_hydrogens=frozenset(), platform=Platform.getPlatformByName('Reference'))
    assert snapshot(value) == before


def test_changed_existing_bond_order_refuses_even_with_same_graph():
    from openmm.app import Double

    value = model()
    before = snapshot(value)

    def changes_order(w):
        edge = list(w.topology.bonds())[0]
        names = edge[0].name, edge[1].name
        w.delete([edge])
        atoms = {a.name:a for a in w.topology.atoms()}
        w.topology.addBond(atoms[names[0]], atoms[names[1]], Double, 2)

    with pytest.raises(PostHGuardError, match='metadata'):
        apply(value, changes_order, protected_hydrogens=frozenset())
    assert snapshot(value) == before


def test_clone_bonds_point_to_cloned_atoms_not_original_objects():
    from scripts.post_hydrogen_public_api import clone_model

    original = model()
    cloned = clone_model(original)
    original_atoms, cloned_atoms = set(original.topology.atoms()), set(cloned.topology.atoms())
    assert not original_atoms & cloned_atoms
    assert snapshot(original) == snapshot(cloned)
    assert all(a in cloned_atoms and b in cloned_atoms for a, b in cloned.topology.bonds())
    assert all(a not in original_atoms and b not in original_atoms for a, b in cloned.topology.bonds())
