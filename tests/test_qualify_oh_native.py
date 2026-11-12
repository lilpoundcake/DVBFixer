"""Observer engineering controls only; no real cohort input or native sampling."""
import json
import random

import pytest


def test_exact_nm_transport_and_independent_measurement(tmp_path):
    pytest.importorskip('openmm')
    from openmm.unit import nanometer
    from scripts.audit_oh_native import measure, read, replay, restore
    from scripts.qualify_hydrogen_postflight import load
    from scripts.qualify_oh_native import BASELINE, snapshot

    from dvbfixer.ffutils import geometry
    from tests.test_hydrogen_postflight import oh_graph

    top, positions = oh_graph()
    snapshot(tmp_path, 'toy', top, positions)
    saved = read(tmp_path/'toy.json.gz')
    _, restored = restore(saved)
    assert [list(p) for p in restored.value_in_unit(nanometer)] == saved['positions_nm']
    assert not measure(saved)['passed']
    baseline = load('oh_native_engineering_baseline', BASELINE)
    old, _ = replay(saved, baseline)
    new, _ = replay(saved, geometry)
    assert not old['measurement']['passed']
    assert new['measurement']['passed']
    assert new['heavy_preserved'] and new['valid_h_preserved']


@pytest.mark.parametrize('mode', ['normal', 'fallback', 'geometry_error', 'budget', 'native_budget', 'native'])
def test_observer_restores_hooks_rng_and_does_not_retry_bounds(tmp_path, monkeypatch, mode):
    pytest.importorskip('openmm')
    pytest.importorskip('pdbfixer')
    from types import SimpleNamespace

    import openmm
    from openmm.app import Modeller
    from scripts import qualify_oh_native as observer

    import dvbfixer.ffutils
    import dvbfixer.prepare
    from dvbfixer.ffutils import geometry
    from tests.test_hydrogen_postflight import oh_graph

    calls = []
    ff = object()
    monkeypatch.setattr(observer, 'check', lambda record: None)
    monkeypatch.setattr(dvbfixer.ffutils, 'build_glycam_system', lambda xmls, *a, **k: ff)
    native_calls = []

    def fake_min(context, tolerance=10., maxIterations=0, reporter=None):
        native_calls.append(context)

    monkeypatch.setattr(openmm.LocalEnergyMinimizer, 'minimize', fake_min)
    monkeypatch.setattr(observer, 'capture_context', lambda *a: dict(platform={'name': 'toy-not-native'}))

    def fake_h(model, forcefield=None, pH=7., variants=None):
        calls.append(forcefield)
        random.random()
        if mode == 'native':
            newTopology = model.topology  # noqa: N806 -- exact upstream local name used by observer.
            system = SimpleNamespace(this=123, getNumParticles=lambda: newTopology.getNumAtoms())
            openmm.LocalEnergyMinimizer.minimize(SimpleNamespace(getSystem=lambda: system))
        if mode == 'native_budget':
            openmm.LocalEnergyMinimizer.minimize(None)
        if mode == 'fallback' and forcefield is not None:
            raise ValueError('toy missing template')

    monkeypatch.setattr(Modeller, 'addHydrogens', fake_h)
    if mode == 'geometry_error':
        def refused(*a, **k):
            raise geometry.HydrogenGeometryError([])
        monkeypatch.setattr(geometry, 'repair_misplaced_hydrogens', refused)

    def fake_main(args):
        top, positions = oh_graph(distance=0.97)
        model = Modeller(top, positions)
        built = dvbfixer.ffutils.build_glycam_system(['toy.xml'], top)
        try:
            model.addHydrogens(built)
        except ValueError:
            model.addHydrogens(None)
        geometry.repair_misplaced_hydrogens(model.topology, model.positions)

    monkeypatch.setattr(dvbfixer.prepare, 'main', fake_main)
    record = dict(cases=[dict(seed=1, requested_ff='amber', ph=7.)],
        protocol=dict(maximum_h_calls=0 if mode == 'budget' else 8,
                      maximum_native_h_calls=0 if mode == 'native_budget' else 8))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(record))
    output = tmp_path/'run'
    output.mkdir()
    before = random.getstate()
    hooks = Modeller.addHydrogens, geometry.repair_misplaced_hydrogens, dvbfixer.ffutils.build_glycam_system
    if mode in ('budget', 'native_budget'):
        with pytest.raises(SystemExit, match='bound exceeded'):
            observer.worker(plan, 0, output)
        assert calls == ([ff] if mode == 'native_budget' else [])
        assert not native_calls
    elif mode == 'geometry_error':
        with pytest.raises(geometry.HydrogenGeometryError):
            observer.worker(plan, 0, output)
        assert calls == [ff]
    else:
        observer.worker(plan, 0, output)
        capture = json.loads((output/'capture.json').read_bytes())
        assert capture['h_calls'][0]['forcefield_xmls'] == ['toy.xml']
        assert calls == ([ff, None] if mode == 'fallback' else [ff])
        if mode == 'fallback':
            assert capture['h_calls'][1]['forcefield_provided'] is False
        if mode == 'native':
            assert len(native_calls) == 1
            assert capture['native'][0]['returned']
            assert capture['native'][0]['before']['platform']['name'] == 'toy-not-native'
    assert random.getstate() == before
    assert hooks == (Modeller.addHydrogens, geometry.repair_misplaced_hydrogens, dvbfixer.ffutils.build_glycam_system)
