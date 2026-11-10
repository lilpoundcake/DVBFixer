"""Native CH2 worker reproduces its bound baseline without touching other atoms."""
import gzip
import json

import pytest
from scripts.research_ASP49_methylene_control import ARMS, worker


@pytest.mark.parametrize('arm', ARMS)
def test_native_CH2_worker_has_exact_unchanged_control_and_fixed_heavy_atoms(tmp_path, arm):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    points = [[0., 0., 0.], [.15, 0., 0.], [-.05, .1414, 0.], [.001, .001, -.1], [.002, .001, -.11]]
    group = dict(parent=0, parent_key=['A', '48', '', 'LEU', 'CB'], hydrogens=[3, 4], heavy=[1, 2])
    system = openmm.System()
    for i in range(5):
        system.addParticle(1. if i in (3, 4) else 0.)
    bond = openmm.HarmonicBondForce()
    for i in (3, 4):
        bond.addBond(0, i, .109, 1000.)
    system.addForce(bond)
    source = tmp_path/'source'
    source.mkdir()
    (source/'H-system.xml').write_text(openmm.XmlSerializer.serialize(system))
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(points*nanometer)
        openmm.LocalEnergyMinimizer.minimize(context, 1., 50)
        final = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
    finally:
        del context
        del integrator
    path = source/'source-path.json.gz'
    path.write_bytes(gzip.compress(json.dumps(dict(initial_nm=points, final_nm=final)).encode()))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, openmm_version=openmm.__version__, policy=dict(maximum_callbacks=512),
        cases=[dict(status='admitted', source=str(source), path=str(path), local=[group], groups=[group])])))
    output = tmp_path/'result'
    output.mkdir()
    worker(plan, 0, arm, output)
    trace = json.loads((output/'trace.json').read_bytes())
    if arm == 'unchanged':
        assert trace['exact_original_nm_endpoint']
    data = json.loads(gzip.decompress((output/'path.json.gz').read_bytes()))
    assert all(data['initial_nm'][i] == points[i] and data['final_nm'][i] == points[i] for i in (0, 1, 2))
    assert not trace['whole_minimization_executed'] and not trace['production_changed']
