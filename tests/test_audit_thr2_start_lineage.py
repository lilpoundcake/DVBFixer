"""Matched-start reports must expose, not conceal, heavy-coordinate changes."""
import json

import pytest
from scripts import audit_thr2_start_lineage as audit
from scripts.research_legacy_start_controls import compressed


@pytest.mark.parametrize('move_heavy', [False, True])
def test_lineage_reports_exact_system_and_actual_motion(tmp_path, monkeypatch, move_heavy):
    openmm = pytest.importorskip('openmm')
    from openmm import unit

    monkeypatch.setattr(audit, 'check', lambda record: None)
    plan = tmp_path/'plan.json'
    labels = ['unseeded', 'seed-1', 'seed-2']
    plan.write_text(json.dumps(dict(protocol=dict(labels=labels))))
    root = tmp_path/'run'
    root.mkdir()
    system = openmm.System()
    system.addParticle(12.)
    system.addParticle(1.)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    atoms = [dict(key=['A', '2', '', 'THR', name], element=element)
        for name, element in (('CB', 'C'), ('HB', 'H'))]
    for index, label in enumerate(labels):
        folder = root/label
        folder.mkdir()
        context.setPositions([[float(move_heavy and index == 2), 0., 0.], [0., float(index), 0.]]*unit.nanometer)
        state = context.getState(positions=True, parameters=True)
        for name, obj in (('state', state), ('system', system), ('integrator', integrator)):
            (folder/f'before-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
        compressed(folder/'before-atoms.json.gz', dict(atoms=atoms))
        (folder/'call-2-topology.json').write_text('{"bonds": [[0, 1]]}')
        (folder/'capture.json').write_text(json.dumps(dict(calls=[dict(number=2, stage='whole', before=dict(prefix='before'))])))
    del context
    (root/'evidence.json').write_text(json.dumps(dict(plan_sha256=audit.sha(plan.read_bytes()),
        native=[dict(label=label) for label in labels])))
    output = tmp_path/'audit.json'
    audit.audit(plan, root, output)
    result = json.loads(output.read_text())
    last = result['comparisons'][-1]
    assert last['exact_system_xml']
    assert last['exact_bond_inventory']
    assert last['heavy']['exact'] == (not move_heavy)
    assert last['heavy']['maximum_angstrom'] == (10. if move_heavy else 0.)
    assert last['hydrogen']['maximum_angstrom'] == 20.
    assert not result['extra_minimization']
