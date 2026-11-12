"""Occurrence-aware observer transport and topology-only recovery capture."""
import json
import random

import pytest


def test_snapshot_roundtrips_repeated_adjacent_chain_ids_and_typed_bonds(tmp_path):
    pytest.importorskip('openmm')
    from openmm import Vec3, unit
    from openmm.app import Single, Topology, element
    from scripts.audit_oh_native import read
    from scripts.qualify_glycam_fallback import restore, snapshot

    top = Topology()
    for _ in range(2):
        res = top.addResidue('HOH', top.addChain('H'), '1', 'A')
        o = top.addAtom('O', element.oxygen, res)
        h = top.addAtom('H1', element.hydrogen, res)
        top.addBond(o, h, Single, 1)
    pos = [Vec3(0., 0., 0.), Vec3(0.1, 0., 0.), Vec3(1., 0., 0.), Vec3(1.1, 0., 0.)]*unit.nanometer
    snapshot(tmp_path, 'input', top, pos)
    saved = read(tmp_path/'input.json.gz')
    reconstructed, coords = restore(saved)
    assert [c.id for c in reconstructed.chains()] == ['H', 'H']
    assert [r.index for r in reconstructed.residues()] == [0, 1]
    assert [b.type for b in reconstructed.bonds()] == [Single, Single]
    assert [b.order for b in reconstructed.bonds()] == [1, 1]
    assert coords.value_in_unit(unit.nanometer) == pos.value_in_unit(unit.nanometer)
    saved['atoms'][-1]['key'][0] = 'h'
    with pytest.raises(ValueError, match='Inconsistent'):
        restore(saved)


def test_worker_captures_recovery_without_sampling_and_restores_hooks(tmp_path, monkeypatch):
    pytest.importorskip('openmm')
    from scripts import qualify_glycam_fallback as observer

    from dvbfixer.prepare import fallback
    from tests.test_prepare_glycam_fallback import glycam_graph

    monkeypatch.setattr(observer.qualify_oh_native, 'check', lambda _: None)

    def fake_worker(plan, index, output):
        top, pos, _, _ = glycam_graph()
        assert fallback.restore_glycam_fallback_bonds(top, pos) == 3

    monkeypatch.setattr(observer.qualify_oh_native, 'worker', fake_worker)
    original = fallback.restore_glycam_fallback_bonds
    before_rng = random.getstate()
    plan = tmp_path/'plan.json'
    plan.write_text('{}')
    observer.worker(plan, 0, tmp_path)
    assert fallback.restore_glycam_fallback_bonds is original
    assert random.getstate() == before_rng
    capture = json.loads((tmp_path/'recovery.json').read_bytes())
    assert capture['recoveries'] == [dict(number=1, returned=True, added_bonds=3)]
    before = observer.audit_oh_native.read(tmp_path/'recovery-1-before.json.gz')
    after = observer.audit_oh_native.read(tmp_path/'recovery-1-after.json.gz')
    assert before['atoms'] == after['atoms'] and before['positions_nm'] == after['positions_nm']
    assert len(after['bonds']) == len(before['bonds'])+3
