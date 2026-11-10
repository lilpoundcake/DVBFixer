"""Native whole baselines precede all fixed-arm continuations, no safety repair."""
import gzip
import json

import pytest
from scripts.research_carbon_H_whole import preflight, worker


def fixture(tmp_path):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    rows = [dict(key=['d', '82', 'A', 'ALA', n], element=e) for n, e in [('CA', 'C'), ('HA', 'H')]]
    points, bonds = [[0., 0., 0.], [.109, 0., 0.]], [[0, 1]]
    def save(name, data):
        path = tmp_path/name
        raw = json.dumps(data).encode()
        path.write_bytes(gzip.compress(raw) if path.suffix == '.gz' else raw)
        return str(path)
    topology = save('topology.json', dict(atoms=rows, bonds=bonds))
    case = dict(id='synthetic', seed=1, conditioning='prepared-view', status='admitted',
        context_kind='prepared-protein-view', atoms=topology, topology=topology,
        pre_H=save('pre.json.gz', dict(atoms=[dict(rows[0], xyz=[0., 0., 0.])])),
        ownership=save('ownership.json', dict(rows=rows, bonds=bonds, new_h=[1], mobile=[1])),
        baseline_path=save('H-baseline.json.gz', dict(initial_nm=points, final_nm=points)),
        system=str(tmp_path/'H-system.xml'))
    system = openmm.System()
    for mass in (0., 1.):
        system.addParticle(mass)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, .109, 1000.)
    system.addForce(bond)
    (tmp_path/'H-system.xml').write_text(openmm.XmlSerializer.serialize(system))
    system.setParticleMass(0, 12.)
    positional = openmm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    for name in ('k', 'x0', 'y0', 'z0'):
        positional.addPerParticleParameter(name)
    positional.addParticle(0, [100., 0., 0., 0.])
    system.addForce(positional)
    pull = openmm.CustomExternalForce('-10*x')
    pull.addParticle(1, [])
    system.addForce(pull)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(points*nanometer)
        state = context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
        for name, obj in [('system', system), ('integrator', integrator), ('state', state)]:
            (tmp_path/f'whole-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
    finally:
        del context
        del integrator
    spec = dict(rows=rows, bonds=bonds, mapping=[0, 1], system=str(tmp_path/'whole-system.xml'),
        state=str(tmp_path/'whole-state.xml'), integrator=str(tmp_path/'whole-integrator.xml'))
    case['whole'] = spec
    h_root = tmp_path/'H-run'
    h_folder = h_root/'synthetic/seed-1/prepared-view/unchanged'
    h_folder.mkdir(parents=True)
    (h_folder/'path.json.gz').write_bytes(gzip.compress(json.dumps(dict(final_nm=points)).encode()))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, cases=[case], h_root=str(h_root),
        openmm_version=openmm.__version__, policy=dict(arms=['unchanged']))))
    return plan, case, dict(rows=rows, bonds=bonds, final_nm=points)


def test_native_baseline_and_direct_continuation_agree_at_both_endpoints(tmp_path):
    plan, case, data = fixture(tmp_path)
    assert preflight(case, case['whole'], data)['exact_original_nm_start']
    baseline_root = tmp_path/'baseline'
    baseline = baseline_root/'synthetic/seed-1/prepared-view'
    baseline.mkdir(parents=True)
    worker(plan, 0, 'native_baseline', None, baseline)
    output = tmp_path/'continuation'
    output.mkdir()
    worker(plan, 0, 'unchanged', baseline_root, output)
    trace = json.loads((output/'trace.json').read_bytes())
    assert all(f['exact_baseline_nm'] for f in trace['frames'])
    assert not trace['emergency_reflection'] and not trace['observations_restored']
    assert not list(output.glob('*.pdb'))


def test_original_State_start_drift_refuses_before_whole_sampling(tmp_path):
    _, case, data = fixture(tmp_path)
    data['final_nm'][1][0] += .001
    with pytest.raises(ValueError, match='original native nm State'):
        preflight(case, case['whole'], data)


def test_legacy_label_alias_is_not_guessed_at_whole_bridge(tmp_path):
    _, case, data = fixture(tmp_path)
    case['context_kind'] = 'geometry-reconstruction-stress'
    data['rows'] = [dict(r, key=['D', *r['key'][1:]]) for r in data['rows']]
    with pytest.raises(ValueError, match='legacy full labels'):
        preflight(case, case['whole'], data)


def test_historical_Angstrom_comparison_is_not_relabelled_as_nm(tmp_path):
    plan, case, _ = fixture(tmp_path)
    historical = tmp_path/'historical.json.gz'
    # Deliberately incorrect initial Å positions; the worker must retain its
    # measured nm evidence and refuse rather than weakening the comparison.
    historical.write_bytes(gzip.compress(json.dumps(dict(atoms=[dict(r, xyz=[0., 0., 0.])
        for r in case['whole']['rows']])).encode()))
    record = json.loads(plan.read_bytes())
    record['cases'][0]['whole']['historical_angstrom'] = dict(initial=str(historical))
    plan.write_text(json.dumps(record))
    output = tmp_path/'baseline'
    output.mkdir()
    with pytest.raises(ValueError, match='historical Angstrom'):
        worker(plan, 0, 'native_baseline', None, output)
    assert (output/'initial.json.gz').exists() and not (output/'phase1.json.gz').exists()
