"""Logging-only phase boundaries in unchanged production legacy minimize.

Replay ALL six retained prepared inputs once, no settings/geometry repair/native
metrics. Supported Python Simulation.minimizeEnergy observer forwards args/results
unchanged and snapshots actual Context coordinates, not a C++ virtual callback.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify

PARENT=Path('tmp/pdbfixer-research/legacy-comparator-plan-01.json')
RUN=Path('tmp/pdbfixer-research/legacy-comparator-run-01')


def worker_argv(python,plan_path,case,folder):
    return [python,str(Path(__file__).resolve()),'worker','--plan',str(plan_path.resolve()),
            '--case',case,'--output',str(folder.resolve())]


def plan(output):
    parent=json.loads(PARENT.read_bytes())
    verify(parent)
    inputs={str(PARENT):sha(PARENT.read_bytes())}
    cases=[]
    for case in parent['cases']:
        folder=RUN/case['id']
        for name in ('prepared.pdb','prepared.dat','minimized.pdb','case-result.json'):
            inputs[str(folder/name)]=sha((folder/name).read_bytes())
        cases.append(dict(id=case['id'],target=case['target'],prepared=str(folder/'prepared.pdb'),
                          dat=str(folder/'prepared.dat'),baseline=str(folder/'minimized.pdb')))
    hashes=dict(parent['source_hashes'])
    hashes['scripts/research_legacy_stereo_trace.py']=sha(Path(__file__).read_bytes())
    write(output,dict(schema='dvbfixer.legacy-stereo-trace-plan.v1',cases=cases,inputs=inputs,
        source_hashes=hashes,python=parent['python'],native_coordinates_read=False,
        protocol=dict(enrollment='ALL six original prepared cases, not only final stereo failures',
            intervention='logging-only Simulation.minimizeEnergy wrapper; forward original self/args/kwargs/return untouched',
            snapshots='all actual topology atom positions before/after every native minimizeEnergy call, plus raw output',
            settings='same production minimize --no-solvent --atom-naming standard with copied dat; default physics/platform',
            repetitions=1,timeout_seconds_per_case=1800,geometry_repair=False,
            limits='previously native-evaluated inputs; stage mechanism only, no native metrics or optimizer profile choice; automatic-platform endpoint differences retained, no repeat-until-match')))


def worker(plan_path,case_id,output):
    record=json.loads(plan_path.read_bytes())
    verify(record)
    case=next(c for c in record['cases'] if c['id']==case_id)
    from openmm.app import Simulation
    from openmm.unit import angstrom

    from dvbfixer.cli import main

    original=Simulation.minimizeEnergy
    events=[]
    count=0

    def snapshot(simulation,call,stage):
        coordinates=simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(angstrom)
        atoms=[]
        for atom,xyz in zip(simulation.topology.atoms(),coordinates,strict=True):
            residue=atom.residue
            atoms.append(dict(key=[residue.chain.id,str(residue.id),residue.insertionCode.strip(),residue.name,atom.name],
                element=atom.element.symbol if atom.element else None,xyz=[float(x) for x in xyz]))
        path=output/f'call-{call}-{stage}.json.gz'
        raw=json.dumps(dict(atoms=atoms),sort_keys=True,allow_nan=False).encode()
        path.write_bytes(gzip.compress(raw,mtime=0))
        events.append(dict(call=call,stage=stage,snapshot=path.name,sha256=sha(path.read_bytes())))

    def observe(simulation,*args,**kwargs):
        nonlocal count
        count+=1
        call=count
        snapshot(simulation,call,'before')
        result=original(simulation,*args,**kwargs)
        snapshot(simulation,call,'after')
        return result

    sys.argv=['dvbfixer','minimize',str(output/'prepared.pdb'),'--dat',str(output/'prepared.dat'),
              '--no-solvent','--atom-naming','standard','-o',str(output/'minimized.pdb')]
    with patch.object(Simulation,'minimizeEnergy',observe):
        main()
    if Simulation.minimizeEnergy is not original:
        raise ValueError('native phase observer not restored')
    verify(record)
    write(output/'trace.json',dict(id=case_id,events=events,calls=count,hook_restored=True,
        plan_sha256=sha(plan_path.read_bytes()),native_coordinates_read=False,
        args_forwarded_unchanged=True,coordinates_restored=False,prepared_source=case['prepared']))


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    verify(record)
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    env=dict(os.environ,PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH',''),
             PYTHONPATH=os.pathsep.join((str(Path('src').resolve()),str(Path('.').resolve()))))
    results=[]
    for case in record['cases']:
        folder=output/case['id']
        folder.mkdir(exist_ok=False)
        for name,path in (('prepared.pdb',case['prepared']),('prepared.dat',case['dat'])):
            (folder/name).write_bytes(Path(path).read_bytes())
        argv=worker_argv(record['python'],plan_path,case['id'],folder)
        with (folder/'minimize.log').open('xb') as log:
            try:
                proc=subprocess.run(argv,env=env,stdout=log,stderr=subprocess.STDOUT,
                    timeout=record['protocol']['timeout_seconds_per_case'],check=False)
                result=dict(id=case['id'],returncode=proc.returncode,
                    status='completed' if proc.returncode==0 else 'execution_blocker')
            except subprocess.TimeoutExpired:
                result=dict(id=case['id'],returncode=None,status='timeout')
        result['artifacts']={p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
        results.append(result)
        write(folder/'case-result.json',result)
    verify(record)
    write(output/'evidence.json',dict(schema='dvbfixer.legacy-stereo-trace-run.v1',cases=results,
        plan_sha256=sha(plan_path.read_bytes()),native_coordinates_read=False,production_changed=False))
    print(json.dumps(dict(enrolled=len(results),completed=sum(r['status']=='completed' for r in results)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--output',type=Path,required=True)
    for name in ('run','worker'):
        p=sub.add_parser(name)
        p.add_argument('--plan',type=Path,required=True)
        p.add_argument('--output',type=Path,required=True)
        if name=='worker':
            p.add_argument('--case',required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.output)
    elif args.command=='run':
        run(args.plan,args.output)
    else:
        worker(args.plan,args.case,args.output)
