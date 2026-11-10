"""Downstream CH2 worker matches independently called native baseline phases."""
import gzip
import json

import pytest
from scripts.research_graph_methyl_whole import worker as baseline_worker
from scripts.research_methylene_whole_control import worker


def fixture(tmp_path):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    source = tmp_path/'source'
    source.mkdir()
    rows = [dict(key=['A','1','','ALA',n],element=e) for n,e in [('CA','C'),('HA','H')]]
    bonds,points = [[0,1]],[[0.,0.,0.],[.109,0.,0.]]
    (source/'whole-topology.json').write_text(json.dumps(dict(atoms=rows,bonds=bonds,H_to_whole_indices=[0,1])))
    (source/'H-topology.json').write_text(json.dumps(dict(atoms=rows,bonds=bonds)))
    system = openmm.System()
    system.addParticle(12.)
    system.addParticle(1.)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0,1,.109,1000.)
    system.addForce(bond)
    restraint = openmm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    for name in ('k','x0','y0','z0'):
        restraint.addPerParticleParameter(name)
    restraint.addParticle(0,[100.,0.,0.,0.])
    system.addForce(restraint)
    pull = openmm.CustomExternalForce('-10*x')
    pull.addParticle(1,[])
    system.addForce(pull)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system,integrator,openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(points*nanometer)
        state = context.getState(positions=True,velocities=True,parameters=True,integratorParameters=True)
        for name,obj in [('system',system),('integrator',integrator),('state',state)]:
            (source/f'whole-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
    finally:
        del context
        del integrator
    h_root = tmp_path/'H-run'
    h_folder = h_root/'synthetic/seed-1/source/unchanged'
    h_folder.mkdir(parents=True)
    (h_folder/'path.json.gz').write_bytes(gzip.compress(json.dumps(dict(final_nm=points)).encode()))
    baseline = tmp_path/'baseline'
    baseline.mkdir()
    plan = tmp_path/'plan.json'
    case = dict(id='synthetic',seed=1,status='captured',source_arm='source',source=str(source),baseline=str(baseline))
    plan.write_text(json.dumps(dict(inputs={},source_hashes={},openmm_version=openmm.__version__,h_root=str(h_root),cases=[case],policy=dict(arms=['unchanged']))))
    baseline_worker(plan,0,'native_baseline',None,baseline)
    output = tmp_path/'continuation'
    output.mkdir()
    return plan,baseline,output


def test_six_control_contract_reduces_to_exact_tiny_native_phases(tmp_path):
    plan,_,output = fixture(tmp_path)
    worker(plan,0,'unchanged',output)
    trace = json.loads((output/'trace.json').read_bytes())
    assert all(f['exact_baseline_nm'] for f in trace['frames'])
    assert trace['post_H_repairs'] == 0 and not trace['emergency_reflection']
    assert not trace['observations_restored'] and not list(output.glob('*.pdb'))


def test_phase_parameter_change_is_not_hidden_by_reproduced_start(tmp_path):
    plan,baseline,output = fixture(tmp_path)
    path = baseline/'phase2-system.xml'
    path.write_text(path.read_text().replace('mass="12"','mass="13"'))
    with pytest.raises(ValueError,match='forces/anchors/phase'):
        worker(plan,0,'unchanged',output)
    assert (output/'initial.json.gz').exists() and (output/'phase1.json.gz').exists()
    assert not (output/'phase2.json.gz').exists()


def test_saved_only_audit_reads_worker_payloads_without_reminimizing(tmp_path,monkeypatch):
    from scripts import audit_methylene_whole_control as auditor
    from scripts.research_dunbrack_prior import sha

    from dvbfixer.ffutils.dat import DatRecord

    plan,_,_ = fixture(tmp_path)
    record = json.loads(plan.read_bytes())
    DatRecord().save(tmp_path/'source/prepared.dat',verbose=False)
    root = tmp_path/'audited-run'
    folder = root/'synthetic/seed-1/source/unchanged'
    folder.mkdir(parents=True)
    worker(plan,0,'unchanged',folder)
    (root/'executed-plan.json').write_bytes(plan.read_bytes())
    artifacts = {p.name:sha(p.read_bytes()) for p in folder.iterdir() if p.is_file()}
    (root/'evidence.json').write_text(json.dumps(dict(plan_sha256=sha(plan.read_bytes()),cases=[
        dict(id='synthetic',seed=1,source_arm='source',arm='unchanged',status='completed',returncode=0,artifacts=artifacts)])))
    # This two-atom fixture exercises IO/force/bridge contracts. Actual focal
    # geometry and exact case/insertion-code matching have separate unit tests.
    sites = []
    def site(rows,points,target):
        sites.append(target)
        return dict(residue=target,synthetic_contract_only=True)
    monkeypatch.setattr(auditor,'site_geometry',site)
    import openmm
    def forbidden(*args,**kwargs):
        raise AssertionError('saved-only whole audit must not minimize')
    monkeypatch.setattr(openmm.LocalEnergyMinimizer,'minimize',forbidden)
    output = tmp_path/'whole-audit.json'
    auditor.audit(plan,root,output)
    result = json.loads(output.read_bytes())
    assert result['summary']['arms']['unchanged']['exact_native_baseline_cells'] == 1
    assert result['policy'] == record['policy'] and len(sites) == 6
    assert sites == [['A','49','','ASP'],['A','151','','ALA']]*3
