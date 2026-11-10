"""Prospective geometry ledger, followed by separate all-output native reporting.

Reuse frozen geometry and rewrite auditors. No new coordinates, minimization,
source normalization, selected-profile search or native-dependent output choice.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import statistics
from pathlib import Path

from scripts.audit_legacy_comparator import atoms, motion
from scripts.audit_legacy_start_controls import load_atoms
from scripts.audit_start_geometry import geometry, heavy_map, input_audit
from scripts.compare_branched_completion import measure
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_comparator import verify

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def check_run(plan_path,root):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    verify(plan)
    run_raw=(root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(raw) or (root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('prospective geometry executed plan differs')
    for row in run['cases']:
        folder=root/row['id']
        if 'seed' in row:
            folder=folder/f"seed-{row['seed']}"/row['arm']
        for name,digest in row['artifacts'].items():
            if sha((folder/name).read_bytes())!=digest:
                raise ValueError('prospective retained artifact changed')
    return plan,run,run_raw


def select(mask_plan,prepared_root,refinement_plan,refined_root,output):
    masks,prepared,prepared_raw=check_run(mask_plan,prepared_root)
    plan,run,run_raw=check_run(refinement_plan,refined_root)
    if plan['cohort_mask_plan']!=str(mask_plan):
        raise ValueError('prospective refinement refers to another cohort')
    if [r['id'] for r in prepared['cases']]!=[c['id'] for c in masks['cases']]:
        raise ValueError('complete source preparation enrollment differs')
    maskroot=Path(masks['mask_root'])
    source_reports=[]
    for case,prep in zip(masks['cases'],prepared['cases'],strict=True):
        source=(maskroot/case['source']).read_bytes()
        if sha(source)!=masks['artifacts'][case['source']] or source!=(prepared_root/case['id']/'input.pdb').read_bytes():
            raise ValueError('prospective original masked input changed')
        try:
            request=plan_branched_source(source,ComponentInstanceRef(**case['target'])).request
            if set(request.missing_atoms)!=set(case['mask']):
                raise ValueError('prospective exact missing mask changed')
        except ValueError as error:
            if case['status']!='shared_admission_refused' or str(error)!=case['reason']:
                raise ValueError('prospective source refusal no longer reproduced') from error
        else:
            if case['status']!='admitted':
                raise ValueError('prospective source unexpectedly admitted')
        row=dict(id=case['id'],admission=case['status'],preparation_status=prep['status'],reason=case.get('reason'),
            retained_source_H_records=sum(line.startswith((b'ATOM  ',b'HETATM')) and line[76:78].strip() in (b'H',b'D') for line in source.splitlines()))
        if prep['status']=='prepared':
            folder=prepared_root/case['id']
            controls=json.loads((folder/'prepare-controls.json').read_bytes())
            if not controls['hooks_restored'] or any(c['platform']!='Reference' for c in controls['fixer_instances']):
                raise ValueError('prospective preparation heavy-atom platform differs')
            if any(c['platform']!='Reference' or c['seed']!=1 or not c['state_restored'] for c in controls['H_calls']):
                raise ValueError('prospective preparation H controls differ')
            row['preparation_original_observation_motion']=motion(source,(folder/'prepared.pdb').read_bytes())
        source_reports.append(row)
    by_id={c['id']:c for c in plan['cases']}
    masks_by_id={c['id']:c for c in masks['cases']}
    expected=[(c['id'],s,a) for c in plan['cases'] for s in plan['protocol']['seeds'] for a in plan['protocol']['arms']]
    if [(r['id'],r['seed'],r['arm']) for r in run['cases']]!=expected:
        raise ValueError('all prospective refinement cells not retained')
    reports=[]
    for row in run['cases']:
        result=dict(id=row['id'],seed=row['seed'],arm=row['arm'],status=row['status'],residue=by_id[row['id']]['residue'])
        if row['status']!='completed':
            result['error']=row.get('error')
            reports.append(result)
            continue
        case=by_id[row['id']]
        folder=refined_root/row['id']/f"seed-{row['seed']}"/row['arm']
        change=input_audit(case,folder,row['arm'])
        trace=json.loads((folder/'trace.json').read_bytes())
        if trace['plan_sha256']!=sha(refinement_plan.read_bytes()) or not trace['hooks_restored'] or len(trace['h_calls'])!=1:
            raise ValueError('prospective refinement hook contract differs')
        if any(h['platform']!='Reference' or h['seed']!=row['seed'] or not h['python_random_state_restored'] for h in trace['h_calls']):
            raise ValueError('prospective controlled H placement differs')
        bonds=json.loads((folder/'topology.json').read_bytes())['bonds']
        phases=[]
        for phase in trace['phases']:
            if any(phase[s]['platform']['name']!='Reference' for s in ('before','after')):
                raise ValueError('prospective native phase platform differs')
            after=load_atoms(folder/f"phase-{phase['phase']}-after.json.gz")
            phases.append(dict(phase=phase['phase'],raw=geometry(after,bonds,case['target'],case['residue']),
                rounded=geometry(after,bonds,case['target'],case['residue'],rounded=True)))
        initial=load_atoms(folder/'phase-1-before.json.gz')
        final_raw=(folder/'minimized.pdb').read_bytes()
        target_atoms=atoms(final_raw)
        actual_heavy={(a.chain_id,a.sequence_number,a.insertion_code,n):atom for (a,n),atom in target_atoms.items()}
        if actual_heavy.keys()!=heavy_map(initial).keys():
            raise ValueError('prospective final heavy identities differ')
        actual=[dict(a,xyz=actual_heavy[tuple(a['key'][:3]+[a['key'][-1]])].position) if a['element'] not in {'H','D',None} else a for a in after]
        final=geometry(actual,bonds,case['target'],case['residue'],rounded=True)
        source=(maskroot/masks_by_id[row['id']]['source']).read_bytes()
        result.update(source_rewrite=change,initial_geometry=geometry(initial,bonds,case['target'],case['residue']),
            phases=phases,final_geometry=final,final_geometry_accepted=final['accepted'],
            original_observation_motion=motion(source,final_raw),final_pdb_sha256=sha(final_raw))
        reports.append(result)
    summary={}
    for arm in plan['protocol']['arms']:
        rows=[r for r in reports if r['arm']==arm]
        finished=[r for r in rows if r['status']=='completed']
        summary[arm]=dict(enrolled_masks=len(masks['cases']),executed_cells=len(finished),
            geometry_passes=sum(r['final_geometry_accepted'] for r in finished),
            THR_executed=sum(r['residue']=='THR' for r in finished),
            THR_geometry_passes=sum(r['residue']=='THR' and r['final_geometry_accepted'] for r in finished),
            execution_blockers=len(rows)-len(finished))
    ledger=dict(schema='dvbfixer.prospective-geometry-ledger.v1',sources=source_reports,cases=reports,summary=summary,
        systems=masks['systems'],inherited_nonexecutions=plan['inherited_nonexecutions'],
        mask_plan_sha256=sha(mask_plan.read_bytes()),refinement_plan_sha256=sha(refinement_plan.read_bytes()),
        preparation_run_sha256=sha(prepared_raw),refinement_run_sha256=sha(run_raw),
        audit_script_sha256=sha(Path(__file__).read_bytes()),native_reference_coordinates_parsed=False,
        native_accuracy_reported=False,profile_selected=False,production_changed=False)
    write(output,ledger)
    print(json.dumps(summary,indent=2))


def target_errors(points,truth,names,residue):
    if any(n not in truth or n not in points for n in names):
        raise ValueError('native reference lacks an exact held-out atom')
    errors=[math.dist(points[n],truth[n]) for n in names]
    exact=math.sqrt(sum(x*x for x in errors)/len(errors))
    symmetry=exact
    if residue=='VAL':
        swapped=dict(truth,CG1=truth['CG2'],CG2=truth['CG1'])
        symmetry=min(exact,math.sqrt(sum(math.dist(points[n],swapped[n])**2 for n in names)/len(names)))
    return dict(exact_name_rmsd_angstrom=exact,generated_VAL_symmetry_rmsd_angstrom=symmetry,
        maximum_error_angstrom=max(errors),aligned=False,used_for_selection=False)


def report(mask_plan,refinement_plan,refined_root,ledger_path,output):
    masks=json.loads(mask_plan.read_bytes())
    plan=json.loads(refinement_plan.read_bytes())
    verify(masks)
    verify(plan)
    ledger_raw=ledger_path.read_bytes()
    ledger=json.loads(ledger_raw)
    if ledger['audit_script_sha256']!=sha(Path(__file__).read_bytes()) or ledger['refinement_plan_sha256']!=sha(refinement_plan.read_bytes()):
        raise ValueError('geometry ledger changed before native reporting')
    by_id={c['id']:c for c in masks['cases']}
    results=[]
    native_controls={}
    for cell in ledger['cases']:
        if cell['status']!='completed':
            continue
        case=by_id[cell['id']]
        target=ComponentInstanceRef(**case['target'])
        folder=refined_root/cell['id']/f"seed-{cell['seed']}"/cell['arm']
        actual=(folder/'minimized.pdb').read_bytes()
        if sha(actual)!=cell['final_pdb_sha256']:
            raise ValueError('committed output changed before native accuracy')
        reference=Path(masks['mask_root'])/case['reference']
        native=reference.read_bytes()
        if sha(native)!=masks['artifacts'][case['reference']]:
            raise ValueError('held-out reference bytes changed')
        truth={n:a.position for (r,n),a in atoms(native).items() if r==target}
        points={n:a.position for (r,n),a in atoms(actual).items() if r==target}
        row=dict(id=cell['id'],seed=cell['seed'],arm=cell['arm'],residue=case['residue'],geometry_accepted=cell['final_geometry_accepted'])
        try:
            row['accuracy']=target_errors(points,truth,case['mask'],case['residue'])
            if cell['id'] not in native_controls:
                source=(Path(masks['mask_root'])/case['source']).read_bytes()
                request=plan_branched_source(source,target).request
                additions=tuple(Addition(ExactAtomRef(target,n),n[0],truth[n]) for n in case['mask'])
                measured,_=measure(source,request,additions,anchored=False)
                native_controls[cell['id']]=dict(accepted=measured.accepted,raw_findings=measured.raw['findings'])
        except ValueError as error:
            row['accuracy_unavailable']=str(error)
        results.append(row)
    by_cell={(r['id'],r['seed'],r['arm']):r for r in results if 'accuracy' in r}
    pairs=[]
    for first,second in itertools.combinations(plan['protocol']['arms'],2):
        for stratum in ('all','THR','VAL'):
            for mode in ('all_completed','both_geometry_pass'):
                matched=[]
                for cell in results:
                    key=(cell['id'],cell['seed'],second)
                    if cell['arm']!=first or 'accuracy' not in cell or key not in by_cell:
                        continue
                    other=by_cell[key]
                    if stratum!='all' and cell['residue']!=stratum:
                        continue
                    if mode=='both_geometry_pass' and not (cell['geometry_accepted'] and other['geometry_accepted']):
                        continue
                    matched.append((cell,other))
                medians={}
                for index,arm in enumerate((first,second)):
                    values=[p[index]['accuracy']['exact_name_rmsd_angstrom'] for p in matched]
                    medians[arm]=statistics.median(values) if values else None
                pairs.append(dict(first=first,second=second,stratum=stratum,mode=mode,matched_cells=len(matched),
                    median_exact_rmsd_angstrom=medians,ids=[dict(id=a['id'],seed=a['seed']) for a,_ in matched]))
    write(output,dict(schema='dvbfixer.prospective-geometry-native-report.v1',cases=results,matched_pairs=pairs,
        native_geometry_controls=native_controls,ledger_sha256=sha(ledger_raw),
        native_reporting_after_geometry_commit=True,geometry_failures_in_all_completed=True,
        profile_selected=False,scientific_qualification=False,production_changed=False))
    print(json.dumps(dict(completed_outputs_reported=len(results),native_controls=len(native_controls),matched_rows=len(pairs)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ('select','report'):
        p=sub.add_parser(name)
        p.add_argument('--masks',type=Path,required=True)
        p.add_argument('--refinement-plan',type=Path,required=True)
        p.add_argument('--refined',type=Path,required=True)
        p.add_argument('--output',type=Path,required=True)
        if name=='select':
            p.add_argument('--prepared',type=Path,required=True)
        else:
            p.add_argument('--ledger',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='select':
        select(args.masks,args.prepared,args.refinement_plan,args.refined,args.output)
    else:
        report(args.masks,args.refinement_plan,args.refined,args.ledger,args.output)
