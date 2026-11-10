"""Saved-State native baseline and repaired unchanged continuation agree exactly."""
import gzip
import json

import pytest
from scripts.research_graph_methyl_whole import worker


def fixture(tmp_path):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    rows = [dict(key=['a', '82', 'A', 'ALA', n], element=e) for n, e in [('CA', 'C'), ('HA', 'H')]]
    bonds, points = [[0, 1]], [[0., 0., 0.], [.109, 0., 0.]]
    source = tmp_path/'capture'
    source.mkdir()
    (source/'whole-topology.json').write_text(json.dumps(dict(atoms=rows, bonds=bonds, H_to_whole_indices=[0, 1])))
    (source/'H-topology.json').write_text(json.dumps(dict(atoms=rows, bonds=bonds)))
    system = openmm.System()
    for mass in (12., 1.):
        system.addParticle(mass)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, .109, 1000.)
    system.addForce(bond)
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
            (source/f'whole-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
    finally:
        del context
        del integrator
    h_root = tmp_path/'H-run'
    h_folder = h_root/'synthetic'/'seed-1'/'unchanged'
    h_folder.mkdir(parents=True)
    (h_folder/'path.json.gz').write_bytes(gzip.compress(json.dumps(dict(final_nm=points)).encode()))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, openmm_version=openmm.__version__, h_root=str(h_root),
        cases=[dict(id='synthetic', seed=1, status='captured', source=str(source))], policy=dict(arms=['unchanged']))))
    baseline_root = tmp_path/'baseline'
    baseline = baseline_root/'synthetic'/'seed-1'
    baseline.mkdir(parents=True)
    worker(plan, 0, 'native_baseline', None, baseline)
    return plan, baseline_root


def test_unchanged_matches_native_baseline_start_and_both_endpoints(tmp_path):
    plan, baseline = fixture(tmp_path)
    output = tmp_path/'continuation'
    output.mkdir()
    worker(plan, 0, 'unchanged', baseline, output)
    trace = json.loads((output/'trace.json').read_bytes())
    assert all(f['exact_baseline_nm_state'] for f in trace['frames'])
    assert trace['post_H_repairs'] == 0
    assert not trace['emergency_reflection'] and not trace['observations_restored']
    assert not list(output.glob('*.pdb'))


def test_baseline_parameter_drift_refuses_and_retains_evidence(tmp_path):
    plan, baseline = fixture(tmp_path)
    path = baseline/'synthetic'/'seed-1'/'initial-system.xml'
    path.write_text(path.read_text().replace('mass="12"', 'mass="13"'))
    output = tmp_path/'continuation'
    output.mkdir()
    with pytest.raises(ValueError, match='force/anchor/phase'):
        worker(plan, 0, 'unchanged', baseline, output)
    assert (output/'initial.json.gz').exists()
    assert not (output/'phase1.json.gz').exists()


def test_saved_only_auditor_reads_worker_baseline_and_continuation_end_to_end(tmp_path):
    from scripts.audit_graph_methyl_whole import audit
    from scripts.research_dunbrack_prior import sha

    from dvbfixer.ffutils.dat import DatRecord

    plan, baseline = fixture(tmp_path)
    record = json.loads(plan.read_bytes())
    source = tmp_path/'capture'
    DatRecord().save(source/'prepared.dat', verbose=False)
    root = tmp_path/'run'
    folder = root/'synthetic'/'seed-1'/'unchanged'
    folder.mkdir(parents=True)
    worker(plan, 0, 'unchanged', baseline, folder)
    for directory, mode in [(baseline, 'native_baseline'), (root, 'continuations')]:
        (directory/'executed-plan.json').write_bytes(plan.read_bytes())
        cell_folder = directory/'synthetic'/'seed-1'
        if mode == 'continuations':
            cell_folder = cell_folder/'unchanged'
        artifacts = {p.name: sha(p.read_bytes()) for p in cell_folder.iterdir() if p.is_file()}
        evidence = dict(mode=mode, plan_sha256=sha(plan.read_bytes()),
            cases=[dict(id='synthetic', seed=1, arm='native_baseline' if mode == 'native_baseline' else 'unchanged',
                status='completed', returncode=0, artifacts=artifacts)])
        if mode == 'continuations':
            evidence['baseline_run_sha256'] = sha((baseline/'evidence.json').read_bytes())
        (directory/'evidence.json').write_text(json.dumps(evidence))
    output = tmp_path/'audit.json'
    audit(plan, root, baseline, output)
    result = json.loads(output.read_bytes())
    assert result['summary']['arms']['unchanged']['exact_native_baseline_cells'] == 1
    assert result['summary']['arms']['unchanged']['final_negative_CA'] == 0
    assert result['policy'] == record['policy']
