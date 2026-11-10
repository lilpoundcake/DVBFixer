"""Synthetic H-stage capture and restoration controls, no cohort experiments."""
import json
from pathlib import Path
from unittest.mock import patch

import pytest
from scripts.research_ha_placement_trace import worker


def setup_plan(tmp_path, openmm):
    source = tmp_path/'source'
    source.mkdir()
    for name in ('prepared.pdb', 'prepared.dat', 'input-change.json'):
        (source/name).write_text('{}\n')
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, historical_root=str(Path.cwd()),
        openmm_version=openmm.__version__, original_plan='unused-synthetic.json',
        cases=[dict(id='synthetic-water', seed=1, source_arm='raw', source=str(source))])))
    output = tmp_path/'output'
    output.mkdir()
    return plan, output


def test_capture_forwards_real_native_h_minimizer_and_stops_before_full_phase(tmp_path):
    openmm = pytest.importorskip('openmm')
    from openmm.app import Modeller, Simulation, Topology, element
    from openmm.unit import nanometer

    plan, output = setup_plan(tmp_path, openmm)
    native_h, native_m = Modeller.addHydrogens, openmm.LocalEnergyMinimizer.minimize

    def synthetic_worker(*args):
        topology = Topology()
        chain = topology.addChain('A')
        water = topology.addResidue('HOH', chain, '1')
        topology.addAtom('O', element.oxygen, water)
        model = Modeller(topology, [[0., 0., 0.]]*nanometer)
        model.addHydrogens(platform=openmm.Platform.getPlatformByName('Reference'))
        Simulation.minimizeEnergy(None)
        raise AssertionError('must not proceed past first full-system boundary')

    with patch('scripts.research_start_geometry.worker', synthetic_worker):
        worker(plan, 0, output)
    trace = json.loads((output/'placement-trace.json').read_bytes())
    assert trace['phase1_not_executed']
    assert trace['h_calls'][0]['platform'] == 'Reference'
    assert trace['h_calls'][0]['maximum_iterations'] == 50
    assert (output/'h-internal-before.json.gz').exists()
    assert (output/'h-internal-after.json.gz').exists()
    assert (output/'h-internal-system.xml').exists()
    assert Modeller.addHydrogens is native_h
    assert openmm.LocalEnergyMinimizer.minimize is native_m


def test_unexpected_h_failure_restores_hooks_and_does_not_emit_success(tmp_path):
    openmm = pytest.importorskip('openmm')
    from openmm.app import Modeller

    plan, output = setup_plan(tmp_path, openmm)
    native_h, native_m = Modeller.addHydrogens, openmm.LocalEnergyMinimizer.minimize

    def broken_worker(*args):
        raise RuntimeError('retained synthetic failure')

    with patch('scripts.research_start_geometry.worker', broken_worker), pytest.raises(RuntimeError, match='synthetic failure'):
        worker(plan, 0, output)
    assert Modeller.addHydrogens is native_h
    assert openmm.LocalEnergyMinimizer.minimize is native_m
    assert not (output/'placement-trace.json').exists()
