"""Toy observer and deterministic guard-gap controls, not SER capture samples."""
import inspect
import json
import random
from pathlib import Path

import pytest
from scripts import research_ser_hg_capture as capture


def test_original_literal_and_historical_detector_gap():
    pytest.importorskip('openmm')
    raw = capture.synthetic_input()
    assert b'OG  SER B 126' in raw
    assert b'HG  SER B 126' not in raw
    # This research assertion describes the pre-fix, source-bound observation,
    # not a requirement that current production preserve the old detector gap.
    # Keep the recorded expected values; executable candidate behavior lives in
    # test_hydrogen_postflight.py and the matched old/new qualification.
    record = json.loads((Path(__file__).resolve().parents[1]/'docs/research/ser-hg-capture-evidence.json').read_bytes())
    sweep = record['static_detector']
    observed_length = next(r for r in sweep['rows'] if r['input_og_hg_angstrom'] == 1.2458133889150464)
    assert not observed_length['test_distance_predicate']
    assert not observed_length['findings']
    assert observed_length['repairs'] == 0
    assert observed_length['coordinates_unchanged']
    assert sweep['distance_only_cutoff_angstrom'] == pytest.approx(1.552)
    far = sweep['rows'][-1]
    assert far['repairs'] == 1
    assert far['final_og_hg_angstrom'] == pytest.approx(0.97)


def test_toy_observer_preserves_h_rng_and_scientific_hooks(tmp_path, monkeypatch):
    openmm = pytest.importorskip('openmm')
    from openmm import unit
    from openmm.app import ForceField, Modeller, Topology, element

    import dvbfixer.prepare
    from dvbfixer.ffutils import geometry

    monkeypatch.setattr(capture, 'check', lambda record: None)
    originals = ForceField.createSystem, Modeller.addHydrogens, openmm.LocalEnergyMinimizer.minimize, geometry.repair_misplaced_hydrogens

    def main(args):
        assert inspect.unwrap(ForceField.createSystem) is originals[0]
        top = Topology()
        residue = top.addResidue('HOH', top.addChain('A'), '1')
        top.addAtom('O', element.oxygen, residue)
        model = Modeller(top, [openmm.Vec3(0, 0, 0)]*unit.nanometer)
        model.addHydrogens(ForceField('tip3p.xml'), platform=openmm.Platform.getPlatformByName('Reference'))
        geometry.repair_misplaced_hydrogens(model.topology, model.positions)

    monkeypatch.setattr(dvbfixer.prepare, 'main', main)
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(cases=[dict(label='seed-1', seed=1)],
        protocol=dict(maximum_h_calls=8, maximum_native_h_calls=8))))
    output = tmp_path/'capture'
    output.mkdir()
    before = random.getstate()
    capture.worker(plan, 0, output)
    assert random.getstate() == before
    assert (ForceField.createSystem, Modeller.addHydrogens, openmm.LocalEnergyMinimizer.minimize, geometry.repair_misplaced_hydrogens) == originals
    report = json.loads((output/'capture.json').read_text())
    assert len(report['native']) == len(report['h_calls']) == len(report['repairs']) == 1
    assert report['native'][0]['returned']
    assert report['repairs'][0]['returned']
