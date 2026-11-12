"""Template-only topology restoration before forcefield-free H fallback."""
import pytest


def glycam_graph(name='UYB', chain='H', resid='232', icode='', extra=()):
    from openmm import Vec3, unit
    from openmm.app import Topology, element

    top = Topology()
    res = top.addResidue(name, top.addChain(chain), resid, icode)
    c1 = top.addAtom('C1', element.carbon, res)
    h1 = top.addAtom('H1', element.hydrogen, res)
    top.addAtom('O5', element.oxygen, res)
    top.addAtom('C2', element.carbon, res)
    coords = [Vec3(0., 0., 0.), Vec3(0.109, 0., 0.), Vec3(0., 0.146, 0.), Vec3(0., 0., 0.156)]
    for atom_name, symbol in extra:
        top.addAtom(atom_name, element.Element.getBySymbol(symbol), res)
        coords.append(Vec3(0.25, 0.01*len(coords), 0.))
    anchor_res = top.addResidue('NLN', top.addChain('a'), '72')
    anchor = top.addAtom('ND2', element.nitrogen, anchor_res)
    top.addBond(c1, anchor)
    positions = (coords+[Vec3(-0.145, 0., 0.)])*unit.nanometer
    return top, positions, c1, h1


def test_unbound_glycam_h_restores_exact_template_links_without_atom_motion():
    pytest.importorskip('openmm')
    from openmm.unit import nanometer

    from dvbfixer.prepare.fallback import restore_glycam_fallback_bonds

    top, pos, c1, h1 = glycam_graph()
    atoms = list(top.atoms())
    before = [list(p) for p in pos.value_in_unit(nanometer)]
    inherited = set(tuple(sorted((a.index, b.index))) for a, b in top.bonds())
    assert restore_glycam_fallback_bonds(top, pos) == 3
    assert list(top.atoms()) == atoms
    assert [list(p) for p in pos.value_in_unit(nanometer)] == before
    bonds = set(tuple(sorted((a.index, b.index))) for a, b in top.bonds())
    assert bonds == inherited | {(c1.index, h1.index), (0, 2), (0, 3)}
    assert restore_glycam_fallback_bonds(top, pos) == 0


def state(top, pos):
    from openmm.unit import nanometer

    return ([(a.residue.chain.index, a.residue.chain.id, a.residue.id,
              a.residue.insertionCode, a.name, a.element) for a in top.atoms()],
            [list(p) for p in pos.value_in_unit(nanometer)],
            [(a.index, b.index, bond.type, bond.order) for bond in top.bonds() for a, b in [bond]])


@pytest.mark.parametrize('defect', ['duplicate', 'missing_parent', 'wrong_element', 'wrong_parent', 'unknown_h'])
def test_invalid_affected_residue_fails_without_live_mutation(defect):
    pytest.importorskip('openmm')
    from openmm.app import Modeller, element

    from dvbfixer.prepare.fallback import GlycamFallbackTopologyError, restore_glycam_fallback_bonds

    extra = {'duplicate': [('H1', 'H')], 'unknown_h': [('HX', 'H')],
             'wrong_parent': [('C6', 'C'), ('H61', 'H')]}.get(defect, [])
    top, pos, c1, h1 = glycam_graph(extra=extra)
    if defect == 'missing_parent':
        model = Modeller(top, pos)
        model.delete([c1])
        top, pos = model.topology, model.positions
    elif defect == 'wrong_element':
        c1.element = element.oxygen
    elif defect == 'wrong_parent':
        top.addBond(h1, list(top.atoms())[2])
        # The other unbound, template-known H selects this residue for repair.
    before = state(top, pos)
    with pytest.raises(GlycamFallbackTopologyError):
        restore_glycam_fallback_bonds(top, pos)
    assert state(top, pos) == before


def test_bonded_input_and_unknown_ligand_do_not_load_templates(monkeypatch):
    pytest.importorskip('openmm')
    import dvbfixer.ffutils
    from dvbfixer.prepare.fallback import restore_glycam_fallback_bonds

    def unexpected(*a, **k):
        raise AssertionError('template loading is outside this control scope')

    monkeypatch.setattr(dvbfixer.ffutils, 'build_glycam_system', unexpected)
    for name, bonded in [('UYB', True), ('LIG', False)]:
        top, pos, c1, h1 = glycam_graph(name)
        if bonded:
            top.addBond(c1, h1)
        before = state(top, pos)
        assert restore_glycam_fallback_bonds(top, pos) == 0
        assert state(top, pos) == before


def test_chain_occurrence_case_icode_and_external_bonds_preserved():
    pytest.importorskip('openmm')
    from openmm.app import Modeller

    from dvbfixer.prepare.fallback import restore_glycam_fallback_bonds

    first, pos, _, _ = glycam_graph()
    model = Modeller(first, pos)
    for chain, icode in [('H', 'A'), ('h', '')]:
        other, coords, _, _ = glycam_graph(chain=chain, icode=icode)
        model.add(other, coords)
    top, pos = model.topology, model.positions
    before = state(top, pos)
    assert restore_glycam_fallback_bonds(top, pos) == 9
    after = state(top, pos)
    assert after[:2] == before[:2]
    assert set(before[2]) <= set(after[2])
    atoms = list(top.atoms())
    assert all(atoms[a].residue is atoms[b].residue for a, b, _, _ in set(after[2])-set(before[2]))


def test_two_residue_validation_is_all_or_nothing():
    pytest.importorskip('openmm')
    from openmm.app import Modeller

    from dvbfixer.prepare.fallback import GlycamFallbackTopologyError, restore_glycam_fallback_bonds

    first, pos, _, _ = glycam_graph()
    model = Modeller(first, pos)
    second, coords, c1, _ = glycam_graph(chain='K')
    other = Modeller(second, coords)
    other.delete([c1])
    model.add(other.topology, other.positions)
    before = state(model.topology, model.positions)
    with pytest.raises(GlycamFallbackTopologyError):
        restore_glycam_fallback_bonds(model.topology, model.positions)
    assert state(model.topology, model.positions) == before


@pytest.mark.parametrize('failure', ['load', 'missing_template', 'wrapper_bonds'])
def test_template_and_wrapper_failures_are_typed_and_atomic(monkeypatch, failure):
    pytest.importorskip('openmm')
    from types import SimpleNamespace

    import dvbfixer.ffutils
    from dvbfixer.prepare.fallback import GlycamFallbackTopologyError, restore_glycam_fallback_bonds

    original = dvbfixer.ffutils.build_glycam_system

    def build(xmls, view, **kwargs):
        if failure == 'load':
            raise ValueError('toy resource failure')
        if failure == 'missing_template':
            return SimpleNamespace(_templates={})
        return original(xmls, view, add_bonds=False, **kwargs)

    monkeypatch.setattr(dvbfixer.ffutils, 'build_glycam_system', build)
    top, pos, _, _ = glycam_graph()
    before = state(top, pos)
    with pytest.raises(GlycamFallbackTopologyError):
        restore_glycam_fallback_bonds(top, pos)
    assert state(top, pos) == before


def test_graph_recovery_cannot_call_native_h_placement_or_create_system(monkeypatch):
    pytest.importorskip('openmm')
    from openmm import LocalEnergyMinimizer
    from openmm.app import ForceField, Modeller

    from dvbfixer.prepare.fallback import restore_glycam_fallback_bonds

    def unexpected(*a, **k):
        raise AssertionError('topology metadata recovery must not sample or create a System')

    monkeypatch.setattr(ForceField, 'createSystem', unexpected)
    monkeypatch.setattr(Modeller, 'addHydrogens', unexpected)
    monkeypatch.setattr(LocalEnergyMinimizer, 'minimize', unexpected)
    top, pos, _, _ = glycam_graph()
    assert restore_glycam_fallback_bonds(top, pos) == 3


@pytest.mark.parametrize('template_failure', [False, True])
def test_prepare_helper_scope_and_typed_error_never_resample(tmp_path, monkeypatch, template_failure):
    pytest.importorskip('openmm')
    pytest.importorskip('pdbfixer')
    from openmm.app import Modeller
    from scripts.research_ser_hg_capture import synthetic_input

    from dvbfixer.prepare import fallback, main

    calls, fallback_calls = [], []

    def h(model, *a, **kw):
        calls.append('H')
        if template_failure:
            raise ValueError('toy template failure')

    def restore(*a, **k):
        fallback_calls.append('restore')
        raise fallback.GlycamFallbackTopologyError('toy ambiguous existing glycan H')

    monkeypatch.setattr(Modeller, 'addHydrogens', h)
    monkeypatch.setattr(fallback, 'restore_glycam_fallback_bonds', restore)
    source, output = tmp_path/'input.pdb', tmp_path/'output.pdb'
    raw = synthetic_input()
    source.write_bytes(raw)
    args = [str(source), '-o', str(output), '--no-propka', '--no-protassign']
    if template_failure:
        with pytest.raises(fallback.GlycamFallbackTopologyError):
            main(args)
        assert fallback_calls == ['restore']
        assert not output.exists() and not output.with_suffix('.dat').exists()
    else:
        main(args)
        assert not fallback_calls
    assert calls == ['H']
    assert source.read_bytes() == raw
