"""Template-only topology restoration before forcefield-free H fallback."""
import pytest


def glycam_graph(name='UYB', chain='H', resid='232', icode=''):
    from openmm import Vec3, unit
    from openmm.app import Topology, element

    top = Topology()
    res = top.addResidue(name, top.addChain(chain), resid, icode)
    c1 = top.addAtom('C1', element.carbon, res)
    h1 = top.addAtom('H1', element.hydrogen, res)
    top.addAtom('O5', element.oxygen, res)
    top.addAtom('C2', element.carbon, res)
    anchor_res = top.addResidue('NLN', top.addChain('a'), '72')
    anchor = top.addAtom('ND2', element.nitrogen, anchor_res)
    top.addBond(c1, anchor)
    positions = [Vec3(0., 0., 0.), Vec3(0.109, 0., 0.), Vec3(0., 0.146, 0.),
        Vec3(0., 0., 0.156), Vec3(-0.145, 0., 0.)]*unit.nanometer
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
