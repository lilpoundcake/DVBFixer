"""All18 CH2 outputs through the same two pre-emergency native whole phases.

Six unchanged controls bind to the previous methyl-source-arm whole trajectories.
No PDB, solvent, emergency repair, new seed or force-field modification.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import time
from pathlib import Path

from scripts.preflight_graph_methyl_transport import mapped_start
from scripts.research_ASP49_methylene_control import ARMS
from scripts.research_beta_wall import restraint
from scripts.research_dunbrack_prior import sha, write
from scripts.research_h_stage_transport import topology_from_rows
from scripts.research_ha_minimizer_path import verify


def plan(h_plan, h_root, h_audit, baseline_root, baseline_audit, output):
    raw = h_plan.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    run_raw = (h_root/'evidence.json').read_bytes()
    h_run = json.loads(run_raw)
    audit_raw = h_audit.read_bytes()
    audit = json.loads(audit_raw)
    baseline_raw = (baseline_root/'evidence.json').read_bytes()
    baseline = json.loads(baseline_raw)
    baseline_audit_raw = baseline_audit.read_bytes()
    prior = json.loads(baseline_audit_raw)
    if audit['plan_sha256'] != sha(raw) or audit['run_sha256'] != sha(run_raw) or audit['summary']['statuses'] != {'completed':18} or audit['summary']['arms']['unchanged']['exact_original_nm_endpoints'] != 6:
        raise ValueError('complete exact CH2 H-stage cohort not accepted')
    if prior['run_sha256'] != sha(baseline_raw) or prior['summary']['statuses'] != {'completed':30, 'refused_source':18}:
        raise ValueError('previous whole native baseline cohort differs')
    inputs = {**parent['inputs'], str(h_plan):sha(raw), str(h_root/'evidence.json'):sha(run_raw), str(h_audit):sha(audit_raw),
        str(baseline_root/'evidence.json'):sha(baseline_raw), str(baseline_audit):sha(baseline_audit_raw)}
    for row in h_run['cases']:
        folder = h_root/row['id']/f"seed-{row['seed']}"/row['source_arm']/row['arm']
        for name, digest in row['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('CH2 H output changed')
            inputs[str(folder/name)] = digest
    cases = []
    for case in parent['cases']:
        if case['status'] != 'admitted':
            raise ValueError('whole continuation requires complete declared six-context H cohort')
        original = next(c for c in baseline['cases'] if (c['id'],c['seed'],c['arm']) == (case['id'],case['seed'],case['source_arm']))
        if original['status'] != 'completed':
            raise ValueError('corresponding original whole baseline has blocker')
        folder = baseline_root/case['id']/f"seed-{case['seed']}"/case['source_arm']
        for name,digest in original['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('original whole baseline changed')
            inputs[str(folder/name)] = digest
        source = Path(case['source'])
        topology = json.loads((source/'whole-topology.json').read_bytes())
        h_topology = json.loads((source/'H-topology.json').read_bytes())
        h_path = json.loads(gzip.decompress((h_root/case['id']/f"seed-{case['seed']}"/case['source_arm']/'unchanged/path.json.gz').read_bytes()))
        points,_ = mapped_start(h_topology['atoms'],h_topology['bonds'],h_path['final_nm'],topology['H_to_whole_indices'])
        expected = json.loads(gzip.decompress((folder/'initial.json.gz').read_bytes()))['positions_nm']
        if points != expected:
            raise ValueError('original H-label repair/lineage nm start differs from corresponding whole baseline')
        cases.append(dict(id=case['id'],seed=case['seed'],source_arm=case['source_arm'],source=case['source'],baseline=str(folder),exact_start_preflight=True))
    sources = dict(parent['source_hashes'])
    for name in ('scripts/research_methylene_whole_control.py','scripts/preflight_graph_methyl_transport.py',
                 'scripts/research_h_stage_transport.py','scripts/research_beta_wall.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output,dict(schema='dvbfixer.methylene-whole-control-plan.v1',cases=cases,inputs=inputs,source_hashes=sources,
        h_root=str(h_root),python=parent['python'],openmm_version=parent['openmm_version'],
        policy=dict(arms=ARMS,per_cell_seconds=1800,total_seconds=10800,platform='Reference',no_solvent=True,
            tolerance_kj_mol_nm=10.,maximum_iterations_per_phase=1000,
            enrollment='all18 CH2 outputs, six original matched source-arm whole baselines, no extra source/seed/retry/phase choice',
            bridge='unchanged original H-label production repair once, object-index mapping to full labels; original starts all exact in nm',
            physics='original captured full System/State/integrator, all anchors/parameters unchanged; only native tenfold positional-weight reduction in phase2',
            audits='exact six unchanged starts/both endpoints; all raw/rounded CA/THR-beta volumes including ASP49 and ALA151, H geometry and original-heavy movement',
            limits='two native phases BEFORE production emergency safety net; no reflection/observation restore/public candidate/PDB or qualification'),production_changed=False))


def worker(plan_path,index,arm,output):
    import openmm
    from openmm.app import Simulation
    from openmm.unit import kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if arm not in ARMS or openmm.__version__ != record['openmm_version']:
        raise ValueError('undeclared whole CH2 arm/runtime')
    case = record['cases'][index]
    source,baseline = Path(case['source']),Path(case['baseline'])
    whole = json.loads((source/'whole-topology.json').read_bytes())
    h_topology = json.loads((source/'H-topology.json').read_bytes())
    h_folder = Path(record['h_root'])/case['id']/f"seed-{case['seed']}"/case['source_arm']/arm
    h_path = json.loads(gzip.decompress((h_folder/'path.json.gz').read_bytes()))
    points,repairs = mapped_start(h_topology['atoms'],h_topology['bonds'],h_path['final_nm'],whole['H_to_whole_indices'])
    system = openmm.XmlSerializer.deserialize((source/'whole-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((source/'whole-integrator.xml').read_text())
    state = openmm.XmlSerializer.deserialize((source/'whole-state.xml').read_text())
    topology = topology_from_rows(whole['atoms'],whole['bonds'])
    simulation = Simulation(topology,system,integrator,openmm.Platform.getPlatformByName('Reference'))
    frames = []
    try:
        simulation.context.setState(state)
        simulation.context.setPositions(points*nanometer)
        positional = restraint(system)
        for stage in ('initial','phase1','phase2'):
            if stage == 'phase2':
                for i in range(positional.getNumParticles()):
                    particle,parameters = positional.getParticleParameters(i)
                    positional.setParticleParameters(i,particle,[parameters[0]/10,*parameters[1:]])
                positional.updateParametersInContext(simulation.context)
            xml = openmm.XmlSerializer.serialize(system)
            if xml != (baseline/f'{stage}-system.xml').read_text():
                raise ValueError('whole CH2 forces/anchors/phase parameters differ from baseline')
            if stage != 'initial':
                openmm.LocalEnergyMinimizer.minimize(simulation.context,10*kilojoules_per_mole/nanometer,1000)
            result = simulation.context.getState(positions=True,energy=True)
            current = [list(p) for p in result.getPositions().value_in_unit(nanometer)]
            expected = json.loads(gzip.decompress((baseline/f'{stage}.json.gz').read_bytes()))['positions_nm']
            exact = current == expected
            data = dict(positions_nm=current,energy_kj_mol=result.getPotentialEnergy().value_in_unit(kilojoules_per_mole))
            payload = gzip.compress(json.dumps(data,sort_keys=True,allow_nan=False).encode(),mtime=0)
            with (output/f'{stage}.json.gz').open('xb') as stream:
                stream.write(payload)
            with (output/f'{stage}-system.xml').open('x') as stream:
                stream.write(xml)
            frames.append(dict(stage=stage,coordinate_sha256=sha(payload),system_sha256=sha(xml.encode()),exact_baseline_nm=exact))
            if arm == 'unchanged' and not exact:
                raise ValueError('unchanged whole CH2 baseline endpoint differs; retain evidence')
        write(output/'trace.json',dict(case=case,arm=arm,plan_sha256=sha(raw),frames=frames,post_H_repairs=repairs,
            emergency_reflection=False,observations_restored=False,production_changed=False))
    finally:
        del simulation
        del integrator


def run(plan_path,output):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    env = dict(os.environ,PYTHONPATH=os.pathsep.join((str(Path('src').resolve()),str(Path('.').resolve()))))
    results,started = [],time.monotonic()
    for index,case in enumerate(record['cases']):
        for arm in ARMS:
            folder = (output/case['id']/f"seed-{case['seed']}"/case['source_arm']/arm).resolve()
            folder.mkdir(parents=True,exist_ok=False)
            row = dict(id=case['id'],seed=case['seed'],source_arm=case['source_arm'],arm=arm,status='budget_not_executed',returncode=None)
            remaining = record['policy']['total_seconds']-(time.monotonic()-started)
            if remaining > 0:
                command = [record['python'],str(Path(__file__).resolve()),'worker','--plan',str(plan_path.resolve()),
                    '--index',str(index),'--arm',arm,'--output',str(folder)]
                with (folder/'run.log').open('xb') as log:
                    try:
                        process = subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,
                            timeout=min(remaining,record['policy']['per_cell_seconds']),check=False)
                        row.update(status='completed' if process.returncode == 0 else 'execution_blocker',returncode=process.returncode)
                    except subprocess.TimeoutExpired:
                        row['status'] = 'timeout'
            row['artifacts'] = {p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(row)
            write(folder/'result.json',row)
            print(json.dumps({k:v for k,v in row.items() if k != 'artifacts'}),flush=True)
    verify(record)
    write(output/'evidence.json',dict(schema='dvbfixer.methylene-whole-control-run.v1',cases=results,plan_sha256=sha(raw),production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    p = sub.add_parser('plan')
    p.add_argument('--h-plan',type=Path,required=True)
    p.add_argument('--h-root',type=Path,required=True)
    p.add_argument('--h-audit',type=Path,required=True)
    p.add_argument('--baseline-root',type=Path,required=True)
    p.add_argument('--baseline-audit',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    for name in ('run','worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan',type=Path,required=True)
        p.add_argument('--output',type=Path,required=True)
        if name == 'worker':
            p.add_argument('--index',type=int,required=True)
            p.add_argument('--arm',choices=ARMS,required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.h_plan,args.h_root,args.h_audit,args.baseline_root,args.baseline_audit,args.output)
    elif args.command == 'run':
        run(args.plan,args.output)
    else:
        worker(args.plan,args.index,args.arm,args.output)
