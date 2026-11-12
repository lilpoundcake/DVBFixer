"""Toy transplant inventory and native phase transition tests, not THR2 sampling."""
import json

import pytest
from scripts import research_thr2_hswap as swap
from scripts.research_legacy_start_controls import compressed


def test_partition_closure_and_immutable_heavy():
    a = [[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]
    b = [[0., 0., 0.], [3., 0., 0.], [4., 0., 0.]]
    keys = [['A', '2', '', 'THR', name] for name in ('CB', 'HB', 'HA')]
    elements = ['C', 'H', 'H']
    first, selected = swap.swap_points(a, b, elements, keys, 'HB')
    assert selected == [1]
    assert first == [a[0], b[1], a[2]]
    combined, selected = swap.swap_points(first, b, elements, keys, 'other-H')
    assert selected == [2]
    assert combined == b
    assert a[1] == [1., 0., 0.]
    with pytest.raises(ValueError, match='uniquely'):
        swap.swap_points(a, b, elements, [keys[0], keys[1], keys[1]], 'HB')


@pytest.mark.parametrize('mask', ['unchanged', 'HB', 'other-H'])
def test_synthetic_worker_transports_and_chains_exactly(tmp_path, monkeypatch, mask):
    openmm = pytest.importorskip('openmm')
    from openmm import unit

    monkeypatch.setattr(swap, 'check_plan', lambda record: None)
    starts = []
    for label, h in (('recipient', 0.1), ('donor', 0.2)):
        source = tmp_path/label
        source.mkdir()
        system = openmm.System()
        force = openmm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        for name in ('k', 'x0', 'y0', 'z0'):
            force.addPerParticleParameter(name)
        for i in range(3):
            system.addParticle(1.)
            force.addParticle(i, [100., 0., 0., 0.])
        system.addForce(force)
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        context.setPositions([[0., 0.1, 0.], [h, 0., 0.], [0., h, 0.]]*unit.nanometer)
        state = context.getState(positions=True, parameters=True)
        for name, obj in (('system', system), ('integrator', integrator), ('state', state)):
            (source/f'p1-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
        del context
        for i in range(3):
            force.setParticleParameters(i, i, [10., 0., 0., 0.])
        (source/'p2-system.xml').write_text(openmm.XmlSerializer.serialize(system))
        compressed(source/'p1-atoms.json.gz', dict(atoms=[dict(key=['A', '2', '', 'THR', name], element=element)
            for name, element in (('CB', 'C'), ('HB', 'H'), ('HA', 'H'))]))
        starts.append(dict(label=label, source=str(source), calls=[dict(before=dict(prefix='p1', platform={})), dict(before=dict(prefix='p2'))]))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(cases=[dict(recipient='recipient', donor='donor', mask=mask, mode='Reference')],
        starts=starts, policy=dict(maximum_callbacks_per_phase=4096))))
    output = tmp_path/'output'
    output.mkdir()
    swap.worker(plan, 0, output)
    first, second = [swap.read(output/f'phase-{i}/path.json.gz') for i in (1, 2)]
    assert first['returned'] and second['returned']
    assert first['final_nm'] == second['initial_nm']
    assert first['initial_nm'][0] == [0., 0.1, 0.]
    assert first['initial_nm'][1][0] == (0.2 if mask == 'HB' else 0.1)
    assert first['initial_nm'][2][1] == (0.2 if mask == 'other-H' else 0.1)
