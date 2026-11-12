"""Offline chain-run transport only; never sample the real cohort again."""
import copy

import pytest


def saved_graph():
    from openmm import Vec3, unit
    from openmm.app import Topology, element
    from scripts.research_legacy_start_controls import snapshot_atoms

    top = Topology()
    coords = []
    for index, (chain, resid, icode) in enumerate([('A', '1', ''), ('a', '1', ''), ('A', '1', 'B')]):
        res = top.addResidue('HOH', top.addChain(chain), resid, icode)
        o = top.addAtom('O', element.oxygen, res)
        h = top.addAtom('H1', element.hydrogen, res)
        top.addBond(o, h)
        coords.extend([Vec3(index, 0., 0.), Vec3(index+0.12458133889150464, 0., 0.)])
    pos = coords*unit.nanometer
    return dict(atoms=snapshot_atoms(top, pos), positions_nm=[list(p) for p in coords],
                bonds=[[a.index, b.index] for a, b in top.bonds()])


def test_ordered_chain_runs_preserve_case_icode_graph_and_exact_nm():
    pytest.importorskip('openmm')
    from openmm.unit import nanometer
    from scripts.audit_oh_native_v2 import restore_ordered

    saved = saved_graph()
    top, pos = restore_ordered(saved)
    assert [c.id for c in top.chains()] == ['A', 'a', 'A']
    assert [r.insertionCode for r in top.residues()] == ['', '', 'B']
    assert [list(p) for p in pos.value_in_unit(nanometer)] == saved['positions_nm']
    assert [[a.index, b.index] for a, b in top.bonds()] == saved['bonds']


def test_ambiguous_repeated_residue_is_refused():
    pytest.importorskip('openmm')
    from scripts.audit_oh_native_v2 import restore_ordered

    saved = copy.deepcopy(saved_graph())
    for row in saved['atoms'][-2:]:
        row['key'][2] = ''
    with pytest.raises(ValueError, match='ambiguous'):
        restore_ordered(saved)


def test_ordered_transport_replays_identical_repair():
    pytest.importorskip('openmm')
    from openmm.unit import nanometer
    from scripts.audit_oh_native_v2 import restore_ordered

    from dvbfixer.ffutils.geometry import repair_misplaced_hydrogens

    saved = saved_graph()
    a, p = restore_ordered(saved)
    b, q = restore_ordered(saved)
    assert repair_misplaced_hydrogens(a, p) == repair_misplaced_hydrogens(b, q) == 3
    assert p.value_in_unit(nanometer) == q.value_in_unit(nanometer)
