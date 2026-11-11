"""Toy saved-force decomposition; no THR2 scientific sample or optimization."""
import json

import pytest
from scripts import audit_thr2_force_terms as audit
from scripts.research_legacy_start_controls import compressed


def test_saved_only_force_sum_and_no_optimizer(tmp_path, monkeypatch):
    openmm = pytest.importorskip('openmm')
    from openmm import unit

    root = tmp_path/'run'
    source = root/'toy'
    folder = root/'replays/toy/call-1/Reference'
    source.mkdir(parents=True)
    folder.mkdir(parents=True)
    monkeypatch.setattr(audit, 'check', lambda record: None)

    def forbidden(*args, **kwargs):
        pytest.fail('saved-force audit called an optimizer')

    monkeypatch.setattr(openmm.LocalEnergyMinimizer, 'minimize', forbidden)
    points = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [0., 2., 0.], [1., 0., 1.], [0., 1., 1.]]
    names = ['CA', 'N', 'C', 'CB', 'O', 'OG1', 'CG2']
    atoms = [dict(key=['A', '2', '', 'THR', name]) for name in names]
    system = openmm.System()
    for _ in names:
        system.addParticle(12.)
    for strength in (2., 3.):
        force = openmm.CustomExternalForce(f'{strength}*(x*x+y*y+z*z)')
        for index in range(len(names)):
            force.addParticle(index, [])
        system.addForce(force)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    context.setPositions(points*unit.nanometer)
    state = context.getState(positions=True, parameters=True)
    for label, obj in (('system', system), ('integrator', integrator), ('state', state)):
        (source/f'before-{label}.xml').write_text(openmm.XmlSerializer.serialize(obj))
    del context
    compressed(source/'before-atoms.json.gz', dict(atoms=atoms))
    compressed(folder/'path.json.gz', dict(initial_nm=points, final_nm=points, callbacks=[]))
    (source/'capture.json').write_text(json.dumps(dict(calls=[dict(number=1, stage='whole',
        before=dict(prefix='before', platform=dict(name='Reference', constructor_properties={})))])))
    (root/'evidence.json').write_text(json.dumps(dict(replays=[dict(label='toy', number=1, mode='Reference', status='completed')])))
    parent = tmp_path/'parent.json'
    parent.write_text('{}')
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, capture_plan=str(parent), root=str(root), limits='toy')))
    private, public = tmp_path/'private.json', tmp_path/'public.json'
    audit.audit(plan, private, public)
    result = json.loads(public.read_text())
    assert result['summary']['paths'] == 1
    assert result['summary']['nodes'] == 2
    assert max(result['summary']['maximum_sum_errors'].values()) < 1e-10
    assert not result['extra_minimization']
