"""Synthetic engineering controls; not original-input scientific enrollment."""
import gzip
import inspect
import json
import random
from functools import wraps

import pytest
from scripts import research_thr2_native_capture as capture


@pytest.mark.parametrize('index', [0, 1])
def test_native_h_and_whole_capture_and_separate_replay(tmp_path, monkeypatch, index):
    openmm = pytest.importorskip('openmm')
    from openmm import unit
    from openmm.app import ForceField, Modeller, Topology, element

    import dvbfixer.cli

    plan = tmp_path/'plan.json'
    record = dict(protocol=dict(h_seeds=[None, 1, 2], labels=['unseeded', 'seed-1', 'seed-2'],
        maximum_native_calls_per_label=8, maximum_callbacks=4096))
    plan.write_text(json.dumps(record))
    folder = tmp_path/'native'
    folder.mkdir()
    monkeypatch.setattr(capture, 'check', lambda record: None)
    monkeypatch.setattr(capture.sys, 'argv', [])
    monkeypatch.setattr(capture, 'HISTORICAL', tmp_path)
    monkeypatch.setattr(dvbfixer.cli, '__file__', str(tmp_path/'src/dvbfixer/cli.py'))
    native_create = ForceField.createSystem
    native_h = Modeller.addHydrogens
    native_min = openmm.LocalEnergyMinimizer.minimize

    def main():
        # Verify wraps preserves ForceField's default arguments for ArgTracker.
        assert inspect.unwrap(ForceField.createSystem) is native_create
        topology = Topology()
        residue = topology.addResidue('HOH', topology.addChain('A'), '1')
        topology.addAtom('O', element.oxygen, residue)
        model = Modeller(topology, [openmm.Vec3(0, 0, 0)]*unit.nanometer)
        ff = ForceField('tip3p.xml')
        model.addHydrogens(ff, platform=openmm.Platform.getPlatformByName('Reference'))
        system = ff.createSystem(model.topology, rigidWater=False)
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        context.setPositions(model.positions)
        openmm.LocalEnergyMinimizer.minimize(context, 10., 5)
        del context

    monkeypatch.setattr(dvbfixer.cli, 'main', main)
    before = random.getstate()
    capture.native_worker(plan, index, folder)
    assert (random.getstate() == before) == (index == 1)
    assert ForceField.createSystem is native_create
    assert Modeller.addHydrogens is native_h
    assert openmm.LocalEnergyMinimizer.minimize is native_min
    result = json.loads((folder/'capture.json').read_text())
    assert [row['stage'] for row in result['calls']] == ['hydrogen', 'whole']
    assert all(row['returned'] for row in result['calls'])
    for mode in ('captured', 'Reference'):
        output = tmp_path/mode
        output.mkdir()
        capture.replay_worker(plan, folder, 2, mode, output)
        path = json.loads(gzip.decompress((output/'path.json.gz').read_bytes()))
        assert path['returned']
        assert len(path['initial_nm']) == 3


def test_launch_retains_failed_exit(tmp_path):
    import sys

    result = capture.launch(dict(python=sys.executable), ['invalid-command'], tmp_path, 30)
    assert result['status'] == 'execution_blocker'
    assert result['returncode'] != 0
    assert 'run.log' in result['artifacts']


def test_native_failure_retains_boundary_and_restores_hooks(tmp_path, monkeypatch):
    openmm = pytest.importorskip('openmm')
    from openmm import unit
    from openmm.app import ForceField, Modeller, Topology, element

    import dvbfixer.cli

    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(protocol=dict(h_seeds=[None, 1, 2], labels=['unseeded', 'seed-1', 'seed-2'],
        maximum_native_calls_per_label=8))))
    folder = tmp_path/'native'
    folder.mkdir()
    monkeypatch.setattr(capture, 'check', lambda record: None)
    monkeypatch.setattr(capture.sys, 'argv', [])
    monkeypatch.setattr(capture, 'HISTORICAL', tmp_path)
    monkeypatch.setattr(dvbfixer.cli, '__file__', str(tmp_path/'src/dvbfixer/cli.py'))
    native_min = openmm.LocalEnergyMinimizer.minimize
    native_create, native_h = ForceField.createSystem, Modeller.addHydrogens

    @wraps(native_min)
    def fail(context, *args, **kwargs):
        raise RuntimeError('synthetic native failure')

    monkeypatch.setattr(openmm.LocalEnergyMinimizer, 'minimize', fail)

    def main():
        topology = Topology()
        residue = topology.addResidue('HOH', topology.addChain('A'), '1')
        topology.addAtom('O', element.oxygen, residue)
        model = Modeller(topology, [openmm.Vec3(0, 0, 0)]*unit.nanometer)
        model.addHydrogens(ForceField('tip3p.xml'), platform=openmm.Platform.getPlatformByName('Reference'))

    monkeypatch.setattr(dvbfixer.cli, 'main', main)
    saved = random.getstate()
    with pytest.raises(RuntimeError, match='synthetic native failure'):
        capture.native_worker(plan, 1, folder)
    assert random.getstate() == saved
    assert ForceField.createSystem is native_create
    assert Modeller.addHydrogens is native_h
    assert openmm.LocalEnergyMinimizer.minimize is fail
    result = json.loads((folder/'capture.json').read_text())
    assert not result['calls'][0]['returned']
    assert not result['h_calls'][0]['returned']
    assert (folder/'call-1-before.chk').exists()
    assert (folder/'call-1-after.chk').exists()
