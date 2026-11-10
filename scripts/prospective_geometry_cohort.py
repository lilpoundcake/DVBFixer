"""Pre-download new-input qualification of frozen source-geometry interventions.

No native-based case/arm selection, source cleanup, replacement structures or new
geometry/seed/force settings. Keep H/altloc/links and every refusal/shortfall.
"""
from __future__ import annotations

import argparse
import inspect
import json
import os
import random
import subprocess
import sys
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from unittest.mock import patch
from urllib.request import urlopen

from scripts.branched_template_baseline import preflight
from scripts.prospective_independent_cohort import project_first, selected_targets
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.branched_source import plan_branched_source

SYSTEMS={'1FNA':'fibronectin type-III domain, crystal input',
         '1SHG':'alpha-spectrin SH3 domain, crystal input',
         '1VII':'villin headpiece, deposited H/ensemble retention control'}


def enroll(output):
    output.mkdir(parents=True,exist_ok=False)
    hashes={str(p):sha(p.read_bytes()) for p in Path('src/dvbfixer').rglob('*.py')}
    for filename in ('scripts/prospective_geometry_cohort.py','scripts/research_start_geometry.py',
        'scripts/prospective_independent_cohort.py','scripts/research_legacy_start_controls.py',
        'scripts/audit_legacy_start_controls.py','scripts/audit_start_geometry.py'):
        hashes[filename]=sha(Path(filename).read_bytes())
    write(output/'enrollment.json',dict(schema='dvbfixer.prospective-geometry-enrollment.v1',
        systems=[dict(id=i,description=d,url=f'https://files.rcsb.org/download/{i}.pdb') for i,d in SYSTEMS.items()],
        source_hashes=hashes,inputs={'tmp/pdbfixer-research/start-geometry-plan-01.json':
            sha(Path('tmp/pdbfixer-research/start-geometry-plan-01.json').read_bytes())},environment=preflight(),python=sys.executable,
        policy=dict(projection='unchanged project_first: first coordinate model/first ATOM chain, preserve H/altloc/occupancy/touching metadata endpoints',
            masks='first two THR and first two VAL identities in file order, whole CB+both-branch masks only',
            admission='unchanged plan_branched_source; no H removal/altloc selection/link scrub to obtain admission',
            prepare='unchanged legacy prepare/standard names/default pH/PROPKA/Reduce; PDBFixer platform Reference, H seed1/Reference in isolated research child',
            refinement='unchanged research_start_geometry materialize/worker, all raw/beta_only/full_generated arms, seeds1,2, forced H rebuild/Reference inside and outside',
            geometry='unchanged raw/rounded final target geometry plus global CA/THR-beta/CIP; retain initial failures, measured input chi and original observation motion',
            native='only AFTER all geometry outcomes commit; all completed output target RMSDs including geometry failures, plus predeclared common-passing matched arm pairs; never choose arm/case by accuracy',
            original_frame=True,replacement_sources=False,source_normalization=False,
            scientific_independence='prospective source/kernel enrollment, deliberate small proteins not blinded/random/training-independent; H control may refuse, do not count as executed success',
            max_masks=12,max_refinement_cells=72,download_limit_bytes=20*1024*1024,timeout_seconds_per_cell=1800),
        structures_downloaded=False,production_changed=False,native_coordinates_read=False))
    print('Enrolled 1FNA/1SHG/1VII before downloads; fixed kernels and all three arms.')


def fetch(enrollment_path,output):
    raw=enrollment_path.read_bytes()
    record=json.loads(raw)
    verify(record)
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for system in record['systems']:
        row=dict(id=system['id'],url=system['url'],retrieved_utc=datetime.now(UTC).isoformat())
        try:
            limit=record['policy']['download_limit_bytes']
            with urlopen(system['url'],timeout=60) as response:
                data=response.read(limit+1)
                row.update(http_status=response.status,resolved_url=response.geturl())
            if len(data)>limit:
                raise ValueError('download exceeds predeclared bound')
            (output/(system['id']+'.pdb')).write_bytes(data)
            row.update(sha256=sha(data),bytes=len(data))
            if row['http_status']!=200 or not any(line.startswith(b'ATOM  ') for line in data.splitlines()):
                raise ValueError('official response is not a coordinate PDB')
            row['status']='downloaded'
        except Exception as error:
            row.update(status='download_blocker',error_type=type(error).__name__,error=str(error))
        rows.append(row)
        write(output/(system['id']+'-provenance.json'),row)
    write(output/'provenance.json',dict(systems=rows,enrollment_sha256=sha(raw),replacements=False))
    if any(r['status']!='downloaded' for r in rows):
        raise RuntimeError('source download shortfall retained, no replacement')


def freeze(enrollment_path,inputs,output):
    raw=enrollment_path.read_bytes()
    enrollment=json.loads(raw)
    verify(enrollment)
    provenance=json.loads((inputs/'provenance.json').read_bytes())
    if provenance['enrollment_sha256']!=sha(raw):
        raise ValueError('new-source provenance enrollment differs')
    output.mkdir(parents=True,exist_ok=False)
    artifacts={}
    cases=[]
    systems=[]
    for system,download in zip(enrollment['systems'],provenance['systems'],strict=True):
        source=(inputs/(system['id']+'.pdb')).read_bytes()
        if download['id']!=system['id'] or sha(source)!=download['sha256']:
            raise ValueError('downloaded original bytes differ')
        try:
            reference,chain=project_first(source)
            reference_name=system['id']+'-reference.pdb'
            (output/reference_name).write_bytes(reference)
            artifacts[reference_name]=sha(reference)
            targets,counts=selected_targets(reference)
            systems.append(dict(id=system['id'],status='projected',chain=chain,counts=counts,
                count_shortfall=counts!=dict(VAL=2,THR=2)))
            for residue,target in targets:
                branch=BRANCHES[residue]
                mask=('CB',branch.first,branch.second)
                name=f'{system["id"]}-{residue}{target.sequence_number}{target.insertion_code}'
                masked=b''.join(line for line in reference.splitlines(keepends=True) if not(
                    line.startswith(b'ATOM  ') and (line[21:22].decode(),line[22:26].strip().decode(),line[26:27].strip().decode())
                    ==(target.chain_id,target.sequence_number,target.insertion_code) and line[12:16].strip().decode() in mask))
                filename=name+'.pdb'
                (output/filename).write_bytes(masked)
                artifacts[filename]=sha(masked)
                row=dict(id=name,target=asdict(target),residue=residue,mask=mask,source=filename,reference=reference_name)
                try:
                    request=plan_branched_source(masked,target).request
                    if set(request.missing_atoms)!=set(mask):
                        raise ValueError('selected mask does not exactly match missing atoms')
                    row['status']='admitted'
                except ValueError as error:
                    row.update(status='shared_admission_refused',reason=str(error))
                cases.append(row)
        except ValueError as error:
            systems.append(dict(id=system['id'],status='projection_refused',reason=str(error)))
    write(output/'plan.json',dict(schema='dvbfixer.prospective-geometry-mask-plan.v1',cases=cases,systems=systems,
        artifacts=artifacts,source_hashes=enrollment['source_hashes'],inputs={**enrollment['inputs'],str(enrollment_path):sha(raw),
        str(inputs/'provenance.json'):sha((inputs/'provenance.json').read_bytes())},enrollment=enrollment,
        mask_root=str(output),native_coordinates_read=False,production_changed=False))
    print(json.dumps(dict(masks=len(cases),admitted=sum(c['status']=='admitted' for c in cases)),indent=2))


def prepare_worker(plan_path,case_id,folder):
    from openmm import Platform
    from openmm.app import Modeller
    from pdbfixer import PDBFixer

    from dvbfixer.cli import main

    record=json.loads(plan_path.read_bytes())
    verify(record)
    old_init,old_h=PDBFixer.__init__,Modeller.addHydrogens
    fixes=[]
    hydrogens=[]
    reference=Platform.getPlatformByName('Reference')

    def init(fixer,*args,**kwargs):
        old_init(fixer,*args,**kwargs)
        fixer.platform=reference
        fixes.append(dict(platform='Reference'))

    def hydrogen(model,*args,**kwargs):
        bound=inspect.signature(old_h).bind(model,*args,**kwargs)
        bound.apply_defaults()
        if bound.arguments['platform'] is not None:
            raise ValueError('unexpected preparation H platform')
        saved=random.getstate()
        random.seed(1)
        try:
            result=old_h(model,*args,**dict(kwargs,platform=reference))
        finally:
            random.setstate(saved)
        hydrogens.append(dict(seed=1,platform='Reference',state_restored=random.getstate()==saved))
        return result

    sys.argv=['dvbfixer','prepare',str(folder/'input.pdb'),'--backend','legacy','--atom-naming','standard',
        '-o',str(folder/'prepared.pdb'),'--dat',str(folder/'prepared.dat')]
    with patch.object(PDBFixer,'__init__',init),patch.object(Modeller,'addHydrogens',hydrogen):
        main()
    if PDBFixer.__init__ is not old_init or Modeller.addHydrogens is not old_h:
        raise ValueError('research preparation hooks not restored')
    write(folder/'prepare-controls.json',dict(id=case_id,fixer_instances=fixes,H_calls=hydrogens,hooks_restored=True,
        native_coordinates_read=False,plan_sha256=sha(plan_path.read_bytes())))


def prepare(plan_path,output):
    raw=plan_path.read_bytes()
    record=json.loads(raw)
    verify(record)
    if preflight()!=record['enrollment']['environment']:
        raise ValueError('enrolled preparation environment changed')
    output.mkdir(parents=True,exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    env=dict(os.environ,PATH=str(Path(record['enrollment']['python']).parent)+os.pathsep+os.environ.get('PATH',''),
        PYTHONPATH=os.pathsep.join((str(Path('src').resolve()),str(Path('.').resolve()))))
    results=[]
    for case in record['cases']:
        folder=output/case['id']
        folder.mkdir()
        source=Path(record['mask_root'])/case['source']
        if sha(source.read_bytes())!=record['artifacts'][case['source']]:
            raise ValueError('masked source changed before preparation')
        (folder/'input.pdb').write_bytes(source.read_bytes())
        row=dict(id=case['id'],status=case['status'])
        if case['status']=='admitted':
            argv=[record['enrollment']['python'],str(Path(__file__).resolve()),'prepare-worker','--plan',str(plan_path.resolve()),
                  '--case',case['id'],'--output',str(folder.resolve())]
            with (folder/'prepare.log').open('xb') as log:
                try:
                    process=subprocess.run(argv,env=env,stdout=log,stderr=subprocess.STDOUT,check=False,
                        timeout=record['enrollment']['policy']['timeout_seconds_per_cell'])
                    row.update(status='prepared' if process.returncode==0 else 'preparation_blocker',returncode=process.returncode)
                except subprocess.TimeoutExpired:
                    row.update(status='preparation_timeout',returncode=None)
        else:
            row['reason']=case.get('reason')
        if source.read_bytes()!=(folder/'input.pdb').read_bytes():
            raise ValueError('preparation mutated input copy')
        row['artifacts']={p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
        write(folder/'result.json',row)
        results.append(row)
    verify(record)
    write(output/'evidence.json',dict(schema='dvbfixer.prospective-geometry-preparation.v1',cases=results,
        plan_sha256=sha(raw),native_coordinates_read=False,production_changed=False))


def refinement_plan(mask_plan,prepared_root,output):
    raw=mask_plan.read_bytes()
    record=json.loads(raw)
    verify(record)
    preparation_raw=(prepared_root/'evidence.json').read_bytes()
    preparation=json.loads(preparation_raw)
    if preparation['plan_sha256']!=sha(raw):
        raise ValueError('prospective preparation execution differs')
    inputs=dict(record['inputs'])
    inputs[str(mask_plan)]=sha(raw)
    inputs[str(prepared_root/'evidence.json')]=sha(preparation_raw)
    by_id={c['id']:c for c in record['cases']}
    cases=[]
    inherited=[]
    for status in preparation['cases']:
        case=by_id[status['id']]
        folder=prepared_root/case['id']
        for name,digest in status['artifacts'].items():
            if sha((folder/name).read_bytes())!=digest:
                raise ValueError('new preparation artifact changed')
            inputs[str(folder/name)]=digest
        if status['status']!='prepared':
            inherited.append(dict(id=case['id'],status=status['status'],reason=status.get('reason')))
            continue
        owned=DatRecord.load(folder/'prepared.dat').added_keys()
        target=ComponentInstanceRef(**case['target'])
        if any((target.chain_id,target.sequence_number,target.insertion_code,n) not in owned for n in case['mask']):
            inherited.append(dict(id=case['id'],status='ownership_refused'))
            continue
        cases.append(dict(case,prepared=str(folder/'prepared.pdb'),dat=str(folder/'prepared.dat'),
            original_masked_source=str(Path(record['mask_root'])/case['source'])))
    protocol=json.loads(Path('tmp/pdbfixer-research/start-geometry-plan-01.json').read_bytes())['protocol']
    protocol=dict(protocol,source='ALL successfully prepared, provenance-admitted cases of the frozen new cohort; preserve all inherited nonexecutions')
    write(output,dict(schema='dvbfixer.private-start-geometry-plan.v1',cases=cases,inputs=inputs,
        source_hashes=record['source_hashes'],python=record['enrollment']['python'],protocol=protocol,
        inherited_nonexecutions=inherited,cohort_mask_plan=str(mask_plan),native_coordinates_read=False,production_changed=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ('enroll','fetch','freeze','prepare','prepare-worker','refinement-plan'):
        p=sub.add_parser(name)
        p.add_argument('--output',type=Path,required=True)
        if name in ('fetch','freeze'):
            p.add_argument('--enrollment',type=Path,required=True)
        if name=='freeze':
            p.add_argument('--inputs',type=Path,required=True)
        if name in ('prepare','prepare-worker','refinement-plan'):
            p.add_argument('--plan',type=Path,required=True)
        if name=='prepare-worker':
            p.add_argument('--case',required=True)
        if name=='refinement-plan':
            p.add_argument('--prepared',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='enroll':
        enroll(args.output)
    elif args.command=='fetch':
        fetch(args.enrollment,args.output)
    elif args.command=='freeze':
        freeze(args.enrollment,args.inputs,args.output)
    elif args.command=='prepare':
        prepare(args.plan,args.output)
    elif args.command=='prepare-worker':
        prepare_worker(args.plan,args.case,args.output)
    else:
        refinement_plan(args.plan,args.prepared,args.output)
