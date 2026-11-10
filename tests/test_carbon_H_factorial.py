"""Fixed factorial has disjoint ownership, no phase search and exact controls."""
import copy
import gzip
import json

import pytest
from scripts.enroll_carbon_H_factorial import ARMS
from scripts.research_carbon_H_factorial import initialize, worker


def fixture():
    # CH3: indices0..5; CH2: indices6..10. Separate heavy frames.
    points = [[0., 0., 0.], [.15, 0., 0.], [.15, .15, 0.],
        [0., 0., .1], [0., .01, .11], [0., .02, .09],
        [1., 0., 0.], [1.15, 0., 0.], [.95, .14, 0.],
        [1., 0., .1], [1., .01, .11], [9., 8., 7.]]
    case = dict(methyl=dict(groups=[dict(parent=0, parent_key=['A', '1', '', 'ALA', 'CB'],
        axis=1, frame=2, hydrogens=[3, 4, 5])]),
        methylene=dict(groups=[dict(parent=6, parent_key=['A', '2', '', 'LEU', 'CB'],
            heavy=[7, 8], hydrogens=[9, 10])]))
    return points, case


@pytest.mark.parametrize('arm', ARMS)
def test_fixed_factorial_only_changes_declared_initial_H(arm):
    points, case = fixture()
    before = copy.deepcopy(points)
    result = initialize(points, case, arm)
    owned = set()
    if arm in {'all_methyl_zero', 'combined_zero'}:
        owned.update([3, 4, 5])
    if arm in {'all_methylene', 'combined_zero'}:
        owned.update([9, 10])
    assert all(result[i] == points[i] for i in range(len(points)) if i not in owned)
    assert all(result[i] != points[i] for i in owned)
    assert points == before


def test_combined_is_exact_disjoint_union_not_sequential_phase_optimization():
    points, case = fixture()
    methyl = initialize(points, case, 'all_methyl_zero')
    methylene = initialize(points, case, 'all_methylene')
    combined = initialize(points, case, 'combined_zero')
    assert combined == [methyl[i] if i in {3, 4, 5} else methylene[i] for i in range(len(points))]
    case['methylene']['groups'][0]['hydrogens'] = [3, 10]
    with pytest.raises(ValueError, match='not disjoint'):
        initialize(points, case, 'combined_zero')


def test_native_worker_exact_control_and_immutable_observation(tmp_path):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    rows = [dict(key=['d', '82', 'A', 'ALA', n], element=e) for n, e in [('CA', 'C'), ('HA', 'H')]]
    points = [[0., 0., 0.], [.1, 0., 0.]]
    system = openmm.System()
    system.addParticle(0.)
    system.addParticle(1.)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, .109, 1000.)
    system.addForce(bond)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(points*nanometer)
        openmm.LocalEnergyMinimizer.minimize(context, 1., 50)
        final = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
    finally:
        del context
        del integrator
    def save(name, data):
        path = tmp_path/name
        raw = json.dumps(data).encode()
        path.write_bytes(gzip.compress(raw) if path.suffix == '.gz' else raw)
        return str(path)
    topology = save('topology.json', dict(atoms=rows, bonds=[[0, 1]]))
    case = dict(id='synthetic', seed=1, conditioning='prepared-view', status='admitted',
        context_kind='prepared-protein-view', atoms=topology, topology=topology,
        pre_H=save('pre.json.gz', dict(atoms=[dict(rows[0], xyz=[0., 0., 0.])])),
        ownership=save('ownership.json', dict(rows=rows, bonds=[[0, 1]], new_h=[1], mobile=[1])),
        baseline_path=save('baseline.json.gz', dict(initial_nm=points, final_nm=final)),
        system=str(tmp_path/'system.xml'), methyl=dict(groups=[]), methylene=dict(groups=[]))
    (tmp_path/'system.xml').write_text(openmm.XmlSerializer.serialize(system))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, cases=[case],
        openmm_version=openmm.__version__, policy=dict(maximum_callbacks=512))))
    output = tmp_path/'run'
    output.mkdir()
    worker(plan, 0, 'unchanged', output)
    trace = json.loads((output/'trace.json').read_bytes())
    assert trace['exact_original_nm_endpoint'] and not trace['production_changed']
    path = json.loads(gzip.decompress((output/'path.json.gz').read_bytes()))
    assert path['final_nm'] == final and path['final_nm'][0] == points[0]
