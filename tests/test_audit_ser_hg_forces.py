"""Bond graph distance and toy saved-force readback, without optimization."""
import json

import pytest
from scripts import audit_ser_hg_forces as audit
from scripts.audit_ser_hg_forces import graph_distance
from scripts.research_legacy_start_controls import compressed


def test_graph_distance_uses_only_declared_bonds():
    assert graph_distance([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]], 0, 5) == 5
    assert graph_distance([[0, 1], [1, 2]], 0, 3) is None
    assert graph_distance([], 2, 2) == 0


def test_saved_force_inventory_and_no_optimizer(tmp_path, monkeypatch):
    openmm = pytest.importorskip('openmm')
    from openmm import unit

    monkeypatch.setattr(audit, 'check', lambda record: None)

    def forbidden(*args, **kwargs):
        pytest.fail('saved-only audit invoked a minimizer')

    monkeypatch.setattr(openmm.LocalEnergyMinimizer, 'minimize', forbidden)
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(cases=[dict(label='toy', seed=1)])))
    root = tmp_path/'run'
    source = root/'toy'
    source.mkdir(parents=True)
    (root/'evidence.json').write_text(json.dumps(dict(plan_sha256=audit.sha(plan.read_bytes()), cases=[dict(label='toy', seed=1, status='completed')])))
    atoms = [dict(key=['B', '126', '', 'SER', name], element=element)
        for name, element in (('HG', 'H'), ('OG', 'O'), ('CB', 'C'), ('CA', 'C'), ('C', 'C'), ('OXT', 'O'))]
    system = openmm.System()
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, 0.096, 400000.)
    system.addForce(bond)
    nonbonded = openmm.NonbondedForce()
    for i in range(6):
        system.addParticle(1. if i == 0 else 0.)
        nonbonded.addParticle(0.4 if i == 0 else -0.8 if i == 5 else 0., 0.2, 0.)
    system.addForce(nonbonded)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    positions = [[0.096, 0., 0.], [0., 0., 0.], [0.1, 0.1, 0.], [0.2, 0.2, 0.], [0.3, 0.3, 0.], [0.4, 0., 0.]]
    for side in ('before', 'after'):
        if side == 'after':
            positions[0] = [0.399, 0., 0.]
        context.setPositions(positions*unit.nanometer)
        state = context.getState(positions=True, parameters=True)
        (source/f'{side}-state.xml').write_text(openmm.XmlSerializer.serialize(state))
    del context
    (source/'before-system.xml').write_text(openmm.XmlSerializer.serialize(system))
    compressed(source/'before-atoms.json.gz', dict(atoms=atoms))
    (source/'native-1-topology.json').write_text(json.dumps(dict(bonds=[[i, i+1] for i in range(5)])))
    (source/'capture.json').write_text(json.dumps(dict(native=[dict(number=1, tolerance=1.,
        before=dict(prefix='before', platform=dict(name='Reference', constructor_properties={})), after=dict(prefix='after'))])))
    output = tmp_path/'audit.json'
    audit.audit(plan, root, output)
    result = json.loads(output.read_text())
    assert result['summary']['readback_states'] == 4
    assert result['summary']['pair_limit_cases'] == 1
    assert result['cases'][0]['native'][0]['hg_oxt_graph_distance'] == 5
    assert not result['additional_native_minimization']
