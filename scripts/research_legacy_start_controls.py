"""Private H-start controls and complete native-state capture; no production patch.

Two prospectively fixed Python H-initialization seeds on all six original inputs.
Native minimizers/settings remain unchanged. Same-platform checkpoint repeats and
Reference XML-start controls are a separate, bounded saved-state stage.
"""
from __future__ import annotations

import argparse
import gzip
import inspect
import json
import os
import random
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify

PARENT=Path('tmp/pdbfixer-research/legacy-stereo-trace-plan-01.json')
FOCAL=('1CRN-THR2-CB-OG1-CG2','1ENH-THR27-CB-OG1-CG2')
SEEDS=(1,2)
WRITABLE={'Reference':(), 'CPU':('Threads','DeterministicForces'),
    'OpenCL':('DeviceIndex','OpenCLPlatformIndex','Precision','UseCpuPme','DisablePmeStream')}


def platform_properties(platform,context):
    name=platform.getName()
    if name not in WRITABLE:
        raise ValueError('undeclared native replay platform')
    all_values={k:platform.getPropertyValue(context,k) for k in platform.getPropertyNames()}
    return dict(name=name,all_properties=all_values,
                constructor_properties={k:all_values[k] for k in WRITABLE[name]})


def plan(output):
    import openmm
    import openmm.app.modeller as module
    from openmm.app import Modeller, Simulation

    parent=json.loads(PARENT.read_bytes())
    verify(parent)
    record=dict(parent)
    record['schema']='dvbfixer.legacy-start-controls-plan.v1'
    record['inputs']=dict(parent['inputs'])
    for p in (PARENT,Path(module.__file__)):
        record['inputs'][str(p)]=sha(p.read_bytes())
    record['source_hashes']=dict(parent['source_hashes'])
    record['source_hashes']['scripts/research_legacy_start_controls.py']=sha(Path(__file__).read_bytes())
    record['environment']=dict(openmm_version=openmm.__version__,
        add_hydrogens_source_sha256=sha(inspect.getsource(Modeller.addHydrogens).encode()),
        simulation_minimize_source_sha256=sha(inspect.getsource(Simulation.minimizeEnergy).encode()))
    record['protocol']=dict(seeds=SEEDS,enrollment='ALL six retained prepared inputs/dat, each seed once; no native-error or success filter',
        intervention='save/restore Python random state and set seed ONLY around native Modeller.addHydrogens; all other calls/forces/settings unchanged',
        capture='full atom map/bonds, System/Integrator/State XML and same-platform checkpoint at every native phase boundary, H-before/after coordinates and platform properties',
        cli='minimize --no-solvent --atom-naming standard; same input/dat/default physics and automatic platform',
        replay=dict(focal_cases=FOCAL,phase=1,seeds=SEEDS,repeats=2,platforms=('captured','Reference'),
            captured_start='native checkpoint on same declared platform/properties; record initialization agreement',
            reference_start='identical serialized System/Integrator/State; no FF/restraint edits; explicit research numerical control only',
            max_native_calls=16,limits='native phase-1 only, not full pipeline output or cross-platform bit identity'),
        no_adaptation='no additional seeds/repeats/cases/budgets to recover a chosen stereo outcome',
        timeout_seconds_per_case=1800,native_accuracy=False,production_integration=False)
    write(output,record)


def snapshot_atoms(topology,positions):
    from openmm.unit import angstrom

    xyz=positions.value_in_unit(angstrom)
    return [dict(key=[a.residue.chain.id,str(a.residue.id),a.residue.insertionCode.strip(),a.residue.name,a.name],
        element=a.element.symbol if a.element else None,xyz=[float(v) for v in p])
        for a,p in zip(topology.atoms(),xyz,strict=True)]


def compressed(path,record):
    path.write_bytes(gzip.compress(json.dumps(record,sort_keys=True,allow_nan=False).encode(),mtime=0))


def worker(plan_path,case_id,seed,output):
    import openmm
    from openmm.app import Modeller, Simulation
    from openmm.unit import kilojoules_per_mole, nanometer

    from dvbfixer.cli import main

    record=json.loads(plan_path.read_bytes())
    verify(record)
    if seed not in record['protocol']['seeds']:
        raise ValueError('undeclared H-start seed')
    case=next(c for c in record['cases'] if c['id']==case_id)
    original_h=Modeller.addHydrogens
    original_m=Simulation.minimizeEnergy
    h_calls=[]
    phases=[]

    def hydrogen(model,*args,**kwargs):
        i=len(h_calls)+1
        compressed(output/f'h-{i}-before.json.gz',dict(atoms=snapshot_atoms(model.topology,model.positions)))
        saved=random.getstate()
        random.seed(seed)
        try:
            result=original_h(model,*args,**kwargs)
        finally:
            random.setstate(saved)
        compressed(output/f'h-{i}-after.json.gz',dict(atoms=snapshot_atoms(model.topology,model.positions)))
        h_calls.append(dict(call=i,seed=seed,python_random_state_restored=random.getstate()==saved))
        return result

    def capture(simulation,phase,stage):
        prefix=f'phase-{phase}-{stage}'
        state=simulation.context.getState(positions=True,velocities=True,parameters=True,integratorParameters=True)
        for name,obj in (('system',simulation.system),('integrator',simulation.integrator),('state',state)):
            (output/f'{prefix}-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
        (output/f'{prefix}.chk').write_bytes(simulation.context.createCheckpoint())
        compressed(output/f'{prefix}-atoms.json.gz',dict(atoms=snapshot_atoms(simulation.topology,state.getPositions())))
        bonds=[[a.index,b.index] for a,b in simulation.topology.bonds()]
        write(output/f'{prefix}-topology.json',dict(bonds=bonds))
        return dict(prefix=prefix,platform=platform_properties(simulation.context.getPlatform(),simulation.context))

    def minimize(simulation,*args,**kwargs):
        phase=len(phases)+1
        bound=inspect.signature(original_m).bind(simulation,*args,**kwargs)
        bound.apply_defaults()
        if bound.arguments['reporter'] is not None:
            raise ValueError('unexpected native reporter')
        settings=dict(tolerance_kj_mol_nm=bound.arguments['tolerance'].value_in_unit(kilojoules_per_mole/nanometer),
                      max_iterations=bound.arguments['maxIterations'])
        before=capture(simulation,phase,'before')
        result=original_m(simulation,*args,**kwargs)
        after=capture(simulation,phase,'after')
        phases.append(dict(phase=phase,before=before,after=after,settings=settings))
        return result

    sys.argv=['dvbfixer','minimize',str(output/'prepared.pdb'),'--dat',str(output/'prepared.dat'),
              '--no-solvent','--atom-naming','standard','-o',str(output/'minimized.pdb')]
    with patch.object(Modeller,'addHydrogens',hydrogen),patch.object(Simulation,'minimizeEnergy',minimize):
        main()
    if Modeller.addHydrogens is not original_h or Simulation.minimizeEnergy is not original_m:
        raise ValueError('native start-control hooks not restored')
    verify(record)
    write(output/'capture.json',dict(id=case_id,seed=seed,h_calls=h_calls,phases=phases,hooks_restored=True,
        plan_sha256=sha(plan_path.read_bytes()),native_coordinates_read=False,prepared_source=case['prepared']))


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    verify(record)
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    env=dict(os.environ,PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH',''),
        PYTHONPATH=os.pathsep.join((str(Path('src').resolve()),str(Path('.').resolve()))))
    results=[]
    for case in record['cases']:
        for seed in record['protocol']['seeds']:
            folder=output/case['id']/f'seed-{seed}'
            folder.mkdir(parents=True,exist_ok=False)
            for name,source in (('prepared.pdb',case['prepared']),('prepared.dat',case['dat'])):
                (folder/name).write_bytes(Path(source).read_bytes())
            argv=[record['python'],str(Path(__file__).resolve()),'worker','--plan',str(plan_path.resolve()),
                  '--case',case['id'],'--seed',str(seed),'--output',str(folder.resolve())]
            with (folder/'minimize.log').open('xb') as log:
                try:
                    proc=subprocess.run(argv,env=env,stdout=log,stderr=subprocess.STDOUT,
                        timeout=record['protocol']['timeout_seconds_per_case'],check=False)
                    result=dict(id=case['id'],seed=seed,status='completed' if proc.returncode==0 else 'execution_blocker',returncode=proc.returncode)
                except subprocess.TimeoutExpired:
                    result=dict(id=case['id'],seed=seed,status='timeout',returncode=None)
            result['artifacts']={p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(result)
            write(folder/'case-result.json',result)
    verify(record)
    write(output/'evidence.json',dict(schema='dvbfixer.legacy-start-capture-run.v1',cases=results,
        plan_sha256=sha(plan_path.read_bytes()),native_coordinates_read=False,production_changed=False))
    print(json.dumps(dict(enrolled=len(results),completed=sum(c['status']=='completed' for c in results)),indent=2))


def freeze_replay(plan_path,capture_root,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    verify(record)
    run=json.loads((capture_root/'evidence.json').read_bytes())
    if run['plan_sha256']!=sha(raw):
        raise ValueError('complete-state capture plan differs')
    bound={str(plan_path):sha(raw),str(capture_root/'evidence.json'):sha((capture_root/'evidence.json').read_bytes())}
    for result in run['cases']:
        folder=capture_root/result['id']/f"seed-{result['seed']}"
        for name,digest in result['artifacts'].items():
            if sha((folder/name).read_bytes())!=digest:
                raise ValueError('captured native artifact changed')
            bound[str(folder/name)]=digest
    write(output,dict(schema='dvbfixer.legacy-start-replay-plan.v1',parent_plan=str(plan_path),capture_root=str(capture_root),
        protocol=record['protocol']['replay'],inputs=bound,source_hashes=record['source_hashes'],native_coordinates_read=False))


def replay(plan_path,output):
    import openmm
    from openmm.unit import angstrom, kilojoules_per_mole, nanometer

    record=json.loads(plan_path.read_bytes())
    verify(record)
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    outcomes=[]
    for case_id in record['protocol']['focal_cases']:
        for seed in record['protocol']['seeds']:
            source=Path(record['capture_root'])/case_id/f'seed-{seed}'
            capture=json.loads((source/'capture.json').read_bytes())
            phase=next(p for p in capture['phases'] if p['phase']==record['protocol']['phase'])
            prefix=phase['before']['prefix']
            for mode in record['protocol']['platforms']:
                for repeat in range(1,record['protocol']['repeats']+1):
                    folder=output/case_id/f'seed-{seed}'/mode/f'repeat-{repeat}'
                    folder.mkdir(parents=True,exist_ok=False)
                    row=dict(id=case_id,seed=seed,mode=mode,repeat=repeat)
                    try:
                        system=openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
                        integrator=openmm.XmlSerializer.deserialize((source/f'{prefix}-integrator.xml').read_text())
                        if mode=='captured':
                            declared=phase['before']['platform']
                            platform=openmm.Platform.getPlatformByName(declared['name'])
                            context=openmm.Context(system,integrator,platform,declared['constructor_properties'])
                            context.loadCheckpoint((source/f'{prefix}.chk').read_bytes())
                            if platform_properties(platform,context)!=declared:
                                raise ValueError('checkpoint replay platform properties differ')
                        elif mode=='Reference':
                            platform=openmm.Platform.getPlatformByName('Reference')
                            context=openmm.Context(system,integrator,platform)
                            context.setState(openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()))
                        else:
                            raise ValueError('undeclared replay mode')
                        for stage in ('before','after'):
                            state=context.getState(positions=True,energy=True)
                            compressed(folder/f'{stage}.json.gz',dict(
                                positions_angstrom=state.getPositions(asNumpy=True).value_in_unit(angstrom).tolist(),
                                energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)))
                            if stage=='before':
                                openmm.LocalEnergyMinimizer.minimize(context,
                                    phase['settings']['tolerance_kj_mol_nm']*kilojoules_per_mole/nanometer,
                                    phase['settings']['max_iterations'])
                        row.update(status='completed',platform=platform.getName(),settings=phase['settings'])
                        del context,integrator,system
                    except Exception as error:
                        row.update(status='execution_blocker',error_type=type(error).__name__,error=str(error))
                    row['artifacts']={p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
                    outcomes.append(row)
                    write(folder/'result.json',row)
    verify(record)
    write(output/'evidence.json',dict(schema='dvbfixer.legacy-start-replay-run.v1',cases=outcomes,
        plan_sha256=sha(plan_path.read_bytes()),native_coordinates_read=False,production_changed=False))
    print(json.dumps(dict(enrolled=len(outcomes),completed=sum(r['status']=='completed' for r in outcomes)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--output',type=Path,required=True)
    for name in ('run','worker','freeze-replay','replay'):
        p=sub.add_parser(name)
        p.add_argument('--plan',type=Path,required=True)
        p.add_argument('--output',type=Path,required=True)
        if name=='worker':
            p.add_argument('--case',required=True)
            p.add_argument('--seed',type=int,required=True)
        if name=='freeze-replay':
            p.add_argument('--capture',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.output)
    elif args.command=='run':
        run(args.plan,args.output)
    elif args.command=='worker':
        worker(args.plan,args.case,args.seed,args.output)
    elif args.command=='freeze-replay':
        freeze_replay(args.plan,args.capture,args.output)
    else:
        replay(args.plan,args.output)
