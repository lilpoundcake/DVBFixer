"""Native worker preserves exact baseline and uses separate scoring contexts."""
import gzip
import json

import pytest
from scripts.graph_methyl_geometry import discover
from scripts.research_graph_methyl_control import ARMS, worker


@pytest.mark.parametrize('arm', ARMS)
def test_graph_worker_uses_fixed_system_and_exact_unchanged_baseline(tmp_path, arm):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    rows = [dict(key=['a', '82', 'A', 'THR', n], element=e) for n, e in
        [('CG2', 'C'), ('CB', 'C'), ('HG21', 'H'), ('HG22', 'H'), ('HG23', 'H'), ('CA', 'C')]]
    bonds = [[0, 1], [0, 2], [0, 3], [0, 4], [1, 5]]
    points = [[0., 0., 0.], [0., 0., .15], [.01, .001, -.1], [.002, .01, -.11], [.01, .01, -.09], [.1, 0., .2]]
    system = openmm.System()
    for row in rows:
        system.addParticle(1. if row['element'] == 'H' else 0.)
    force = openmm.HarmonicBondForce()
    for i in (2, 3, 4):
        force.addBond(0, i, .109, 1000.)
    system.addForce(force)
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
    (source/'H-path.json.gz').write_bytes(gzip.compress(json.dumps(dict(initial_nm=points, final_nm=final)).encode()))
    graph = discover(rows, bonds, {2, 3, 4}, {2, 3, 4})
    (source/'H-ownership.json').write_text(json.dumps(dict(graph=graph)))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, openmm_version=openmm.__version__,
        cases=[dict(status='captured', groups=1, source=str(source))], policy=dict(maximum_callbacks=512))))
    output = tmp_path/'output'
    output.mkdir()
    worker(plan, 0, arm, output)
    trace = json.loads((output/'trace.json').read_bytes())
    assert trace['groups'] == 1
    assert trace['exact_original_endpoint'] if arm == 'unchanged' else True
    assert len(trace['grids']) == (1 if arm == 'graph_tetrahedral_energy_selected' else 0)
    actual = json.loads(gzip.decompress((output/'path.json.gz').read_bytes()))
    assert all(actual['initial_nm'][i] == points[i] for i in (0, 1, 5))
    assert all(actual['final_nm'][i] == points[i] for i in (0, 1, 5))
