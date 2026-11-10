"""Source-bound transport for pre-enrolled new-input query/mean/control views.

Scientific kernels and complete-pool choice policy were frozen before download.
This adapter transports masked inputs/retained states only; no native references,
new angles/weights/numerical grid or production integration.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.prospective_chi_guard_cohort import verify
from scripts.research_dunbrack_prior import FLAGS, REVISION, sha, write

GROUP='prospective_new_three'
C_POLICY=Path('tmp/pdbfixer-research/chi-stopping-plans-01/absolute-1e-6-gmax1e-4.json')
GUARD_POLICY=Path('tmp/pdbfixer-research/chi-endpoint-guard-plan-01.json')
LIBRARY_ENV=Path('tmp/pdbfixer-research/dunbrack-prior-run-02/environment.json')


def scored_pool_status(cells,scored):
    if [c['id'] for c in cells]!=[c['id'] for c in scored]:
        raise ValueError('complete scored starting pool inventory differs')
    if not cells:
        return 'empty_geometry_pool'
    if any(c['score_status']!='scored' or not c['import_audit']['accepted'] for c in scored):
        return 'score_import_blocker'
    return 'ready'


def bind(inputs,path):
    inputs[str(path)]=sha(path.read_bytes())


def query_plan(mask_root,preparation_root,output):
    from dvbfixer.diagnose.geometry import signed_dihedral_degrees
    from dvbfixer.domain.structure_identity import ComponentInstanceRef
    from dvbfixer.protein_completion.backbone import backbone_angles
    from dvbfixer.protein_completion.branched_source import plan_branched_source

    raw=(mask_root/'plan.json').read_bytes()
    masks=json.loads(raw)
    verify(masks['enrollment'])
    prepared=json.loads((preparation_root/'evidence.json').read_bytes())
    if prepared['plan_sha256']!=sha(raw) or (preparation_root/'executed-plan.json').read_bytes()!=raw or prepared['native_coordinates_read']:
        raise ValueError('new-input preparation binding differs')
    inputs={}
    for path in (mask_root/'plan.json',preparation_root/'evidence.json',preparation_root/'executed-plan.json',LIBRARY_ENV):
        bind(inputs,path)
    cases=[]
    for original in masks['cases']:
        name=original['id']
        row=prepared['results'][name]
        if row['status']!='admitted' or row.get('ncac_status')!='evaluated':
            cases.append(dict(id=name,group=GROUP,status='preparation_not_executed',parent_status=row['status'],
                reason=row.get('reason') or row.get('error') or row.get('ncac_status')))
            continue
        source_path=preparation_root/name/'input.pdb'
        proposal_path=preparation_root/name/'ncac-proposals.json'
        for path in (source_path,proposal_path):
            bind(inputs,path)
        source=source_path.read_bytes()
        if sha(source)!=masks['artifacts'][original['source']]:
            raise ValueError('masked new-input source differs')
        request=plan_branched_source(source,ComponentInstanceRef(**original['target'])).request
        try:
            phi,psi=backbone_angles(request)
        except ValueError as error:
            cases.append(dict(id=name,group=GROUP,status='backbone_refused',reason=str(error)))
            continue
        candidates=json.loads(proposal_path.read_bytes())['candidates']
        if [c['index'] for c in candidates]!=[0,1,2]:
            raise ValueError('new-input parent pool must retain three declared alternatives')
        measured=[]
        for c in candidates:
            points=c['points']
            tip='CG1' if request.residue_name=='VAL' else 'OG1'
            chi=signed_dihedral_degrees(*(points[n] for n in ('N','CA','CB',tip)))
            if chi is None:
                raise ValueError('new-input parent chi degenerate')
            measured.append(dict(index=c['index'],chi_degrees=chi,common_accepted=c['common_audit']['accepted']))
        cases.append(dict(id=name,group=GROUP,status='query_ready',residue_name=request.residue_name,
            phi_degrees=phi,psi_degrees=psi,candidates=measured,source_path=str(source_path),target=original['target'],
            proposal_path=str(proposal_path)))
    hashes=dict(masks['source_hashes'])
    hashes['scripts/research_prospective_guard_views.py']=sha(Path(__file__).read_bytes())
    write(output,dict(schema='dvbfixer.pyrosetta-prior-raw-weight-plan.v2',cases=cases,inputs=inputs,
        source_hashes=hashes,init_flags=FLAGS,expected_revision=REVISION,native_coordinates_read=False,
        selection='unchanged raw three-family-weight selector; samples for frozen mean generation, no native metrics',
        prospective_enrollment_sha256=masks['enrollment_sha256'],production_enabled=False))


def query_run(plan_path,output):
    from scripts.research_dunbrack_prior_weights import run

    run(plan_path,output)
    environment=json.loads((output/'environment.json').read_bytes())
    expected=json.loads(LIBRARY_ENV.read_bytes())
    if environment['binary_sha256']!=expected['binary_sha256']:
        raise ValueError('new-input query library differs from enrolled installed binary')
    write(output/'prospective-library-binding.json',dict(binary_sha256=environment['binary_sha256'],
        historical_environment_sha256=sha(LIBRARY_ENV.read_bytes()),native_coordinates_read=False))


def mean_plan(query_path,query_root,output):
    query_raw=query_path.read_bytes()
    query=json.loads(query_raw)
    for p,d in {**query['inputs'],**query['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen prospective query state changed')
    rows=json.loads((query_root/'selection-before-native.json').read_bytes())
    if rows['plan_sha256']!=sha(query_raw) or rows['native_coordinates_read']:
        raise ValueError('new-input query execution binding differs')
    binding=json.loads((query_root/'prospective-library-binding.json').read_bytes())
    if binding['binary_sha256']!=json.loads(LIBRARY_ENV.read_bytes())['binary_sha256']:
        raise ValueError('new-input library binding differs')
    inputs=dict(query['inputs'])
    for path in (query_path,query_root/'selection-before-native.json',query_root/'environment.json',
                 query_root/'weight-semantics.json',query_root/'prospective-library-binding.json'):
        bind(inputs,path)
    cases=[]
    for case in query['cases']:
        if case['status']!='query_ready':
            cases.append(dict(case,status='parent_not_executed'))
            continue
        row=rows['results'][case['id']]
        if row['status']!='queried':
            raise ValueError('ready new-input library query omitted')
        base=json.loads(Path(case['proposal_path']).read_bytes())['candidates'][0]
        cases.append(dict(case,status='generation_ready',parent_points=base['points'],samples=row['samples'],
            parent_chi_degrees=case['candidates'][0]['chi_degrees']))
    write(output,dict(schema='dvbfixer.dunbrack-means-plan.v1',cases=cases,inputs=inputs,
        source_hashes=query['source_hashes'],native_coordinates_read=False,production_enabled=False,
        protocol='unchanged frozen three API mean-chi transports; observations/CB fixed; no sigma perturbation',
        prospective_enrollment_sha256=query['prospective_enrollment_sha256']))


def control_plan(mean_root,score_root,output):
    score_plan_root=mean_root/'score-plan'
    score_raw=(score_plan_root/'plan.json').read_bytes()
    scores=json.loads(score_raw)
    executed=json.loads((score_root/'evidence.json').read_bytes())
    if executed['plan_sha256']!=sha(score_raw) or (score_root/'executed-plan.json').read_bytes()!=score_raw:
        raise ValueError('new-input mean scoring binding differs')
    record=json.loads(C_POLICY.read_bytes())
    record['inputs']={}
    for path in (C_POLICY,score_plan_root/'plan.json',score_root/'evidence.json',score_root/'executed-plan.json',
                 score_root/'environment.json',mean_root/'geometry-before-scoring.json',LIBRARY_ENV):
        bind(record['inputs'],path)
    geometry=json.loads((mean_root/'geometry-before-scoring.json').read_bytes())
    cases,inherited=[],[]
    views={c['id']:c for c in scores['cases']}
    for name,g in geometry['results'].items():
        if g['status']!='generated':
            inherited.append(dict(id=name,group=GROUP,status=g['status']))
            continue
        case=views[name+'-means']
        row=executed['results'][case['id']]
        status=scored_pool_status(case['cells'],row['cells'])
        if status!='ready':
            inherited.append(dict(id=name,group=GROUP,status=status))
            continue
        source=score_plan_root/case['id']/'source.pdb'
        bind(record['inputs'],source)
        cells=[]
        for declaration,saved in zip(case['cells'],row['cells'],strict=True):
            snapshot=score_root/case['id']/declaration['id']/'pose-snapshot.json.gz'
            candidate=score_plan_root/declaration['pdb']
            if sha(candidate.read_bytes())!=scores['artifacts'][declaration['pdb']] or sha(snapshot.read_bytes())!=saved['snapshot_sha256']:
                raise ValueError('new-input scored candidate/snapshot changed')
            for path in (candidate,snapshot):
                bind(record['inputs'],path)
            cells.append(dict(declaration=declaration,saved=saved,snapshot=str(snapshot)))
        cases.append(dict(id=name,group=GROUP,case=case,source=str(source),cells=cells))
    record.update(cases=cases,inherited_nonexecutions=inherited,score_plan_root=str(score_plan_root),
        source_hashes=scores['source_hashes'],prospective_new_inputs=True)
    write(output,record)


def guard_plan(control_path,control_root,output):
    raw=control_path.read_bytes()
    record=json.loads(raw)
    if (control_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('prospective C execution binding differs')
    run=json.loads((control_root/'evidence.json').read_bytes())
    if run['plan_sha256']!=sha(raw):
        raise ValueError('prospective C result binding differs')
    original_policy=json.loads(GUARD_POLICY.read_bytes())
    record['schema']='dvbfixer.chi-endpoint-guard-plan.v1'
    record['guard_protocol']=dict(original_policy['guard_protocol'],parent_run=str(control_root))
    for path in (control_path,control_root/'evidence.json',control_root/'executed-plan.json',GUARD_POLICY):
        bind(record['inputs'],path)
    cases=[]
    for case in record['cases']:
        rows=run['results'][case['id']]['cells']
        if [r['index'] for r in rows]!=[c['declaration']['index'] for c in case['cells']]:
            raise ValueError('prospective C starting inventory differs')
        if any(r['status']!='refined' for r in rows):
            record['inherited_nonexecutions'].append(dict(id=case['id'],group=GROUP,status='incomplete_C_pool'))
            continue
        for cell in case['cells']:
            folder=control_root/case['id']/cell['declaration']['id']
            for name in ('before-atoms.json','after-atoms.json','result.json'):
                bind(record['inputs'],folder/name)
        cases.append(case)
    record['cases']=cases
    write(output,record)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('query-plan')
    p.add_argument('--masks',type=Path,required=True)
    p.add_argument('--preparation',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('query')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('mean-plan')
    p.add_argument('--query-plan',type=Path,required=True)
    p.add_argument('--query-run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('control-plan')
    p.add_argument('--mean-root',type=Path,required=True)
    p.add_argument('--score-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('guard-plan')
    p.add_argument('--control-plan',type=Path,required=True)
    p.add_argument('--control-run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='query-plan':
        query_plan(args.masks,args.preparation,args.output)
    elif args.command=='query':
        query_run(args.plan,args.output)
    elif args.command=='mean-plan':
        mean_plan(args.query_plan,args.query_run,args.output)
    elif args.command=='control-plan':
        control_plan(args.mean_root,args.score_root,args.output)
    else:
        guard_plan(args.control_plan,args.control_run,args.output)
