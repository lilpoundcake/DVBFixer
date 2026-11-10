"""Native reporters preserve both endpoints and retain interrupted callback data."""
import gzip
import json

import pytest
from scripts.research_carbon_H_whole_path import worker


def fixture(tmp_path):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    source = tmp_path/'source'
    source.mkdir()
    system = openmm.System()
    system.addParticle(12.)
    system.addParticle(1.)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, .109, 1000.)
    system.addForce(bond)
    positional = openmm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    for name in ('k', 'x0', 'y0', 'z0'):
        positional.addPerParticleParameter(name)
    positional.addParticle(0, [100., 0., 0., 0.])
    system.addForce(positional)
    pull = openmm.CustomExternalForce('-100*x')
    pull.addParticle(1, [])
    system.addForce(pull)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions([[0., 0., 0.], [.109, 0., 0.]]*nanometer)
        state = context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
        (source/'state.xml').write_text(openmm.XmlSerializer.serialize(state))
        (source/'integrator.xml').write_text(openmm.XmlSerializer.serialize(integrator))
        for stage in ('initial', 'phase1', 'phase2'):
            if stage == 'phase2':
                particle, parameters = positional.getParticleParameters(0)
                positional.setParticleParameters(0, particle, [parameters[0]/10, *parameters[1:]])
                positional.updateParametersInContext(context)
            if stage != 'initial':
                openmm.LocalEnergyMinimizer.minimize(context, 10., 1000)
            points = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
            (source/f'{stage}.json.gz').write_bytes(gzip.compress(json.dumps(dict(positions_nm=points)).encode()))
            (source/f'{stage}-system.xml').write_text(openmm.XmlSerializer.serialize(system))
    finally:
        del context
        del integrator
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, openmm_version=openmm.__version__,
        policy=dict(maximum_callbacks_per_phase=2048, tolerance_kj_mol_nm=10.),
        cases=[dict(id='synthetic', seed=1, conditioning='raw', arm='unchanged', source=str(source),
            state=str(source/'state.xml'), integrator=str(source/'integrator.xml'), centres=[],
            rows=[dict(key=['d', '82', 'A', 'ALA', n], element=e) for n, e in [('CA', 'C'), ('HA', 'H')]])])))
    output = tmp_path/'run/synthetic/seed-1/raw/unchanged'
    output.mkdir(parents=True)
    return openmm, plan, source, output


def test_reporter_preserves_endpoints_and_saved_only_gradient_audit(tmp_path, monkeypatch):
    openmm, plan, _, output = fixture(tmp_path)
    worker(plan, 0, output)
    trace = json.loads((output/'trace.json').read_bytes())
    assert all(p['exact_original_nm_endpoint'] for p in trace['phases'])
    assert not trace['reporter_requests_early_exit'] and not trace['new_coordinates_or_physics']
    for phase in ('phase1', 'phase2'):
        data = json.loads(gzip.decompress((output/f'{phase}-path.json.gz').read_bytes()))
        assert data['native_call_returned']
        assert all(cb['ordinal'] == i for i, cb in enumerate(data['callbacks']))
    from scripts.audit_carbon_H_whole_path import audit
    from scripts.research_dunbrack_prior import sha

    # Two-atom IO/gradient contract fixture; four-atom focal/guard mathematics
    # and exact case/insertion-code resolution are tested separately.
    root = tmp_path/'run'
    (root/'executed-plan.json').write_bytes(plan.read_bytes())
    artifacts = {p.name:sha(p.read_bytes()) for p in output.iterdir()}
    (root/'evidence.json').write_text(json.dumps(dict(plan_sha256=sha(plan.read_bytes()), cases=[
        dict(id='synthetic', seed=1, conditioning='raw', arm='unchanged', status='completed', returncode=0, artifacts=artifacts)])))
    def forbidden(*args, **kwargs):
        raise AssertionError('native callback auditor must not minimize')
    monkeypatch.setattr(openmm.LocalEnergyMinimizer, 'minimize', forbidden)
    report = tmp_path/'audit.json'
    audit(plan, root, report)
    result = json.loads(report.read_bytes())
    assert result['summary']['statuses'] == {'completed':1}
    assert result['summary']['phases'] == 2
    assert result['summary']['maximum_gradient_readback_error_kj_mol_nm'] == 0.


def test_interrupted_minimizer_retains_collected_callbacks_without_success_trace(tmp_path, monkeypatch):
    openmm, plan, _, output = fixture(tmp_path)
    native = openmm.LocalEnergyMinimizer.minimize
    def interrupted(*args, **kwargs):
        native(*args, **kwargs)
        raise RuntimeError('controlled observation interruption')
    monkeypatch.setattr(openmm.LocalEnergyMinimizer, 'minimize', interrupted)
    with pytest.raises(RuntimeError, match='observation interruption'):
        worker(plan, 0, output)
    data = json.loads(gzip.decompress((output/'phase1-path.json.gz').read_bytes()))
    assert not data['native_call_returned'] and data['callbacks']
    assert not (output/'trace.json').exists() and not (output/'phase2-path.json.gz').exists()


def test_changed_endpoint_is_retained_as_blocker_not_a_successful_replay(tmp_path):
    _, plan, source, output = fixture(tmp_path)
    path = source/'phase1.json.gz'
    data = json.loads(gzip.decompress(path.read_bytes()))
    data['positions_nm'][0][0] += .001
    path.write_bytes(gzip.compress(json.dumps(data).encode()))
    with pytest.raises(ValueError, match='endpoint'):
        worker(plan, 0, output)
    assert (output/'phase1-path.json.gz').exists() and not (output/'trace.json').exists()
