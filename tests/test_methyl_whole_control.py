"""Native two-phase worker matches an independently constructed tiny baseline."""
import gzip
import json

import pytest
from scripts.research_methyl_whole_control import worker


def fixture(tmp_path):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import angstrom, nanometer

    rows = [dict(key=['A', '16', '', 'GLN', n], element=e, xyz=p)
            for n, e, p in [('CA', 'C', [0., 0., 0.]), ('HA', 'H', [1.09, 0., 0.])]]
    points = [[0., 0., 0.], [.109, 0., 0.]]
    source, original = tmp_path/'H-source', tmp_path/'whole-source'
    source.mkdir()
    original.mkdir()
    (source/'h-internal-before.json.gz').write_bytes(gzip.compress(json.dumps(dict(atoms=rows)).encode()))
    (source/'topology.json').write_text(json.dumps(dict(bonds=[[0, 1]])))
    system = openmm.System()
    system.addParticle(12)
    system.addParticle(1)
    bonds = openmm.HarmonicBondForce()
    bonds.addBond(0, 1, .109, 1000.)
    system.addForce(bonds)
    restraint = openmm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    for name in ('k', 'x0', 'y0', 'z0'):
        restraint.addPerParticleParameter(name)
    restraint.addParticle(0, [100., 0., 0., 0.])
    system.addForce(restraint)
    pull = openmm.CustomExternalForce('-10*x')
    pull.addParticle(1, [])
    system.addForce(pull)
    (original/'phase-1-before-system.xml').write_text(openmm.XmlSerializer.serialize(system))
    integrator = openmm.VerletIntegrator(0.)
    (original/'initial-integrator.xml').write_text(openmm.XmlSerializer.serialize(integrator))
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(points*nanometer)
        state = context.getState(positions=True, velocities=True, parameters=True, energy=True)
        (original/'initial-state.xml').write_text(openmm.XmlSerializer.serialize(state))
        for phase in (0, 1, 2):
            if phase == 2:
                restraint.setParticleParameters(0, 0, [10., 0., 0., 0.])
                restraint.updateParametersInContext(context)
                (original/'phase-2-before-system.xml').write_text(openmm.XmlSerializer.serialize(system))
            if phase:
                openmm.LocalEnergyMinimizer.minimize(context, 10., 1000)
            xyz = context.getState(positions=True).getPositions().value_in_unit(angstrom)
            material = [dict(r, xyz=list(p)) for r, p in zip(rows, xyz, strict=True)]
            name = 'phase-1-before.json.gz' if phase == 0 else f'phase-{phase}-after.json.gz'
            (original/name).write_bytes(gzip.compress(json.dumps(dict(atoms=material)).encode()))
    finally:
        del context
        del integrator
    h_run = tmp_path/'H-run'
    endpoint = h_run/'synthetic'/'seed-1'/'raw'/'unchanged'
    endpoint.mkdir(parents=True)
    (endpoint/'path.json.gz').write_bytes(gzip.compress(json.dumps(dict(final_nm=points)).encode()))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, openmm_version=openmm.__version__, h_run=str(h_run),
        cases=[dict(id='synthetic', seed=1, source_arm='raw', source=str(source), original_whole=str(original))])))
    output = tmp_path/'output'
    output.mkdir()
    return plan, output, original


def test_worker_reproduces_both_native_baseline_phases(tmp_path):
    plan, output, _ = fixture(tmp_path)
    worker(plan, 0, 'unchanged', output)
    trace = json.loads((output/'trace.json').read_bytes())
    assert trace['original_parameters_exact']
    assert trace['post_H_repairs'] == 0
    assert [f['stage'] for f in trace['frames']] == ['initial', 'phase1', 'phase2']
    assert all(f['baseline_difference'][kind]['exact'] for f in trace['frames'] for kind in ('heavy', 'hydrogen'))
    assert not trace['emergency_reflection'] and not trace['restored_observations']


def test_phase_parameter_mismatch_retains_partial_evidence_and_fails(tmp_path):
    openmm = pytest.importorskip('openmm')
    plan, output, original = fixture(tmp_path)
    path = original/'phase-2-before-system.xml'
    bad_system = openmm.XmlSerializer.deserialize(path.read_text())
    bad_system.setParticleMass(1, 2.)
    path.write_text(openmm.XmlSerializer.serialize(bad_system))
    with pytest.raises(ValueError, match='parameters/anchors/reference'):
        worker(plan, 0, 'unchanged', output)
    assert (output/'initial.json.gz').exists()
    assert (output/'phase1.json.gz').exists()
    assert not (output/'phase2.json.gz').exists()
