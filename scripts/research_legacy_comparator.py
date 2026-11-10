"""Frozen private full legacy prepare/minimize comparator; real observation motion.

Run production CLI unchanged on fresh owned copies of ALL shared-admitted masks.
No reference coordinates, restoration, scientific thresholds or production edits.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write

MASKS=Path('tmp/pdbfixer-research/prospective-chi-guard-masks-01')
LEDGER=Path('docs/research/prospective-chi-guard-evidence-selection.json')
ENTRY='from dvbfixer.cli import main; main()'


def commands(python,folder):
    common=[python,'-c',ENTRY]
    return [common+['prepare',str(folder/'input.pdb'),'--backend','legacy','--atom-naming','standard',
                    '-o',str(folder/'prepared.pdb'),'--dat',str(folder/'prepared.dat')],
            common+['minimize',str(folder/'prepared.pdb'),'--no-solvent','--atom-naming','standard',
                    '--dat',str(folder/'prepared.dat'),'-o',str(folder/'minimized.pdb')]]


def plan(output):
    masks=json.loads((MASKS/'plan.json').read_bytes())
    ledger=json.loads(LEDGER.read_bytes())
    inputs={str(p):sha(p.read_bytes()) for p in (MASKS/'plan.json',LEDGER,
        Path('docs/reference/prepare.md'),Path('docs/reference/minimize.md'))}
    cases=[]
    for case in masks['cases']:
        if ledger['cases'][case['id']]['status']!='audited':
            continue
        source=MASKS/case['source']
        inputs[str(source)]=sha(source.read_bytes())
        cases.append(dict(id=case['id'],source=str(source),source_sha256=sha(source.read_bytes()),
                          target=case['target'],mask=case['mask'],reference=case['reference']))
    hashes={str(p):sha(p.read_bytes()) for p in sorted(Path('src/dvbfixer').rglob('*.py'))
            if 'protein_completion' not in p.parts}
    hashes['scripts/research_legacy_comparator.py']=sha(Path(__file__).read_bytes())
    write(output,dict(schema='dvbfixer.legacy-comparator-plan.v1',cases=cases,inputs=inputs,source_hashes=hashes,
        python=sys.executable,mask_root=str(MASKS),
        protocol=dict(enrollment='ALL six shared-admitted new-cohort masks; no native-error/stock-success filtering',
            commands='commands() in bound script; prepare legacy then minimize --no-solvent, .dat handoff',
            naming='explicit standard at both stages; no alternate naming/coordinate postprocessing',
            physics='unchanged pH7, default FF/platform, original100/weak5 restraints, max1000 per phase; no solvent only explicit dev convention',
            repetitions=1,timeout_seconds_per_stage=1800,continue_other_cases_on_failure=True,
            child_path='prepend selected interpreter bin only in child; absolute project PYTHONPATH',
            retention='fresh exclusive case directories, raw outputs/logs/sidecars and exit statuses retained',
            audit='actual relaxed coordinates; original-heavy identity/motion and target original-frame RMSD, geometry in real relaxed frame; no coordinate restore/alignment',
            limits='retrospective declared comparator, not an independent new-native accuracy cohort or qualification'),
        native_coordinates_read=False,production_changed=False))


def verify(record):
    for name,digest in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(name).read_bytes())!=digest:
            raise ValueError('frozen legacy comparator source changed')


def run(plan_path,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    verify(record)
    output=output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    env=dict(os.environ,PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH',''),
             PYTHONPATH=str(Path('src').resolve()))
    results=[]
    for case in record['cases']:
        folder=output/case['id']
        folder.mkdir(exist_ok=False)
        source=Path(case['source']).read_bytes()
        (folder/'input.pdb').write_bytes(source)
        stages=[]
        for stage,argv in zip(('prepare','minimize'),commands(record['python'],folder),strict=True):
            row=dict(stage=stage,argv=argv)
            with (folder/(stage+'.log')).open('xb') as log:
                try:
                    proc=subprocess.run(argv,env=env,stdout=log,stderr=subprocess.STDOUT,
                        timeout=record['protocol']['timeout_seconds_per_stage'],check=False)
                    row.update(status='completed' if proc.returncode==0 else 'command_failed',returncode=proc.returncode)
                except subprocess.TimeoutExpired:
                    row.update(status='timeout',returncode=None)
            stages.append(row)
            write(folder/(stage+'-execution.json'),row)
            if row['status']!='completed':
                break
        if Path(case['source']).read_bytes()!=source or (folder/'input.pdb').read_bytes()!=source:
            raise ValueError('legacy execution mutated source input')
        artifacts={str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        result=dict(id=case['id'],stages=stages,artifacts=artifacts,
            status='completed' if len(stages)==2 and all(s['status']=='completed' for s in stages) else 'execution_blocker')
        results.append(result)
        write(folder/'case-result.json',result)
    verify(record)
    write(output/'evidence.json',dict(schema='dvbfixer.legacy-comparator-run.v1',plan_sha256=sha(raw),
        cases=results,native_coordinates_read=False,observations_restored=False,production_changed=False))
    print(json.dumps(dict(enrolled=len(results),completed=sum(r['status']=='completed' for r in results),
        blockers=sum(r['status']!='completed' for r in results)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('run')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.output)
    else:
        run(args.plan,args.output)
