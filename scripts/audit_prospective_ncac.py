"""Independent new-input stage/geometry/score readback, then native reporting.

No coordinates generated or scored. Keep original denominators, including shared
refusals and NCAC blockers; never condition enrollment on a favorable outcome.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import statistics
from dataclasses import asdict
from pathlib import Path

import pdbfixer.pdbfixer as pf
from scripts.audit_ncac_expansion import accuracy
from scripts.audit_ncac_rotamers import check_rotation
from scripts.audit_template_anchors import fit_check
from scripts.audit_template_causality import audit_points, read_pdb
from scripts.compare_branched_completion import measure
from scripts.research_rosetta_joint import sha, verify
from scripts.research_rosetta_repack import identity, observed_atoms, write_json
from scripts.rosetta_fixed_contract import choose_scored_own, import_audit
from scripts.summarize_rosetta_fixed_scoring import source_records_preserved

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def canonical(value):
    return json.loads(json.dumps(value))


def snapshot_check(plan_root,run_root,case,row):
    """Read-only physical/topology import and score-term replay."""
    folder=run_root/case['id']/row['id']
    if json.loads((folder/'score.json').read_bytes())!=row:
        raise ValueError('score record differs')
    raw=(folder/'pose-snapshot.json.gz').read_bytes()
    if sha(raw)!=row['snapshot_sha256']:
        raise ValueError('snapshot digest differs')
    s=json.loads(gzip.decompress(raw))
    expected=observed_atoms((plan_root/row['pdb']).read_bytes())
    ids={identity(a) for a in expected}
    outgoing={tuple(a if a[-1]=='C' else b)[:4] for a,b in case['expected_external'] if {a[-1],b[-1]}=={'C','N'}}
    terminal={(r['chain'],r['resid'],r['icode']) for r in s['residues'] if 'CtermProteinFull' in r['name']}
    allowed=[identity(a) for a in s['heavy_atoms'] if a['atom']=='OXT' and identity(a)[:3] in terminal
             and identity(a)[:4] not in outgoing and identity(a) not in ids]
    audit=import_audit(expected,s['heavy_atoms'],case['expected_external'],s['external'],case['target_internal'],s['target_internal'],allowed)
    if canonical(audit)!=row['import_audit']:
        raise ValueError('independent scorer import audit differs')
    expected_status='scored' if audit['accepted'] else 'import_contract_blocker'
    if row['score_status'] not in (expected_status,'inconsistent_type_topology_blocker'):
        raise ValueError('scorer import status differs')
    if any(not v['virtual_by_rosetta_type'] for v in s['virtual_atoms']) or {identity(v) for v in s['virtual_atoms']} & {identity(a) for a in s['heavy_atoms']}:
        raise ValueError('virtual/physical inventory differs')
    if not row['virtual_coordinates_unchanged_by_scoring']:
        raise ValueError('scoring moved virtual coordinates')
    if row['score_status']=='scored':
        if not math.isfinite(row['total_score']) or not math.isclose(sum(t['raw']*t['weight'] for t in row['score_terms'].values()),row['total_score'],rel_tol=1e-10,abs_tol=1e-8):
            raise ValueError('score-term total differs')
    elif row['total_score'] is not None:
        raise ValueError('blocked score must not be usable')
    return dict(extras=sorted((identity(a),a['xyz']) for a in s['heavy_atoms'] if identity(a) in allowed),
        virtual=sorted((identity(a),a['xyz']) for a in s['virtual_atoms']),
        context_h=sorted(((h['chain'],h['resid'],h['icode'],h['atom']),h['xyz']) for h in s['hydrogens'] if (h['chain'],h['resid'],h['icode'])!=tuple(case['target_key'][:3])))


def additions(request,points):
    return tuple(Addition(ExactAtomRef(request.residue,n),n[0],tuple(points[n])) for n in request.missing_atoms)


def audit(mask_root,science_root,score_run_root,output):
    masks_raw=(mask_root/'plan.json').read_bytes()
    masks=json.loads(masks_raw)
    # Input-only verification here; native/source reference checks are deferred.
    for p,d in masks['source_hashes'].items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('enrolled implementation changed')
    if sha((mask_root/'enrollment.json').read_bytes())!=masks['enrollment_sha256']:
        raise ValueError('enrollment binding changed')
    science_raw=(science_root/'evidence.json').read_bytes()
    science=json.loads(science_raw)
    if science['plan_sha256']!=sha(masks_raw) or (science_root/'executed-plan.json').read_bytes()!=masks_raw or science['native_coordinates_read']:
        raise ValueError('scientific plan/run binding differs')
    score_root=science_root/'score-plan'
    plan_raw=(score_root/'plan.json').read_bytes()
    plan=json.loads(plan_raw)
    verify(plan,score_root)
    scored_raw=(score_run_root/'evidence.json').read_bytes()
    scored=json.loads(scored_raw)
    if scored['plan_sha256']!=sha(plan_raw) or (score_run_root/'executed-plan.json').read_bytes()!=plan_raw:
        raise ValueError('score plan/run binding differs')
    if set(science['results'])!={c['id'] for c in masks['cases']}:
        raise ValueError('scientific case omitted')
    score_cases={c['id']:c for c in plan['cases']}
    if set(scored['results'])!=set(score_cases):
        raise ValueError('scored view omitted')
    ledger,replay={},{}
    snapshots,stage_count,seed_count,fit_count=0,0,0,0
    score_statuses={}
    for case in masks['cases']:
        name=case['id']
        result=science['results'][name]
        folder=science_root/name
        if json.loads((folder/'result.json').read_bytes())!=result:
            raise ValueError('scientific result record changed')
        source=(mask_root/case['source']).read_bytes()
        if sha(source)!=masks['artifacts'][case['source']] or source!=(folder/'input.pdb').read_bytes():
            raise ValueError('source input changed')
        if result['status']=='shared_admission_refused':
            try:
                plan_branched_source(source,ComponentInstanceRef(**case['target']))
            except ValueError as error:
                if str(error)!=result['reason']:
                    raise ValueError('shared refusal differs') from error
            else:
                raise ValueError('shared refusal no longer reproduced')
            if any(name+'-'+a in score_cases for a in ('baseline','ncac')):
                raise ValueError('refused case scored')
            ledger[name]=dict(result,retained_files={str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()})
            continue
        request=plan_branched_source(source,ComponentInstanceRef(**case['target'])).request
        if set(request.missing_atoms)!=set(case['mask']):
            raise ValueError('missing mask changed')
        target=request.residue
        key=(target.chain_id,target.sequence_number,target.insertion_code)
        own=json.loads((folder/'own-attempts.json').read_bytes())
        proposals={}
        own_points=[{a['atom']['name']:a['xyz'] for a in x['candidate']['additions']} for x in own['attempts']]
        proposals['baseline']=own_points
        cell_ledger=dict(status='admitted',ncac_status=result['ncac_status'],methods={})
        cell_ledger['retained_files']={str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        if result['stock']['status']=='evaluated':
            stock=json.loads((folder/'stock-proposal.json').read_bytes())
            stock_added=tuple(Addition(ExactAtomRef(target,a['atom']['name']),a['element'],tuple(a['xyz'])) for a in stock['additions'])
            measured,pdb=measure(source,request,stock_added,anchored=False)
            if canonical(asdict(measured))!=result['stock']['common_audit'] or pdb!=(folder/'stock-proposal.pdb').read_bytes():
                raise ValueError('stock common audit differs')
            cell_ledger['stock']=dict(common_accepted=measured.accepted,common_audit=asdict(measured))
        else:
            cell_ledger['stock']=result['stock']
        if result['ncac_status']=='evaluated':
            trace_raw=(folder/'ncac-trace/evidence.json').read_bytes()
            if sha(trace_raw)!=result['trace_sha256']:
                raise ValueError('NCAC trace digest differs')
            trace=json.loads(trace_raw)
            transport=json.loads((folder/'anchor-policy.json').read_bytes())
            if not trace['hooks_restored'] or not transport['research_hook_restored']:
                raise ValueError('NCAC hooks not restored')
            observed={tuple(a['key']):a['xyz'] for a in read_pdb(source)}
            for event in trace['events']:
                payload=(folder/'ncac-trace'/event['coordinates']).read_bytes()
                if sha(payload)!=event['coordinate_sha256']:
                    raise ValueError('stage digest differs')
                atoms=json.loads(gzip.decompress(payload))['atoms']
                actual={tuple(a['key']):a['xyz'] for a in atoms}
                if len(actual)!=len(atoms) or not observed.keys()<=actual.keys() or any(math.dist(actual[k],p)>1e-6 for k,p in observed.items()):
                    raise ValueError('NCAC stage observation/identity changed')
                points={k[3]:p for k,p in actual.items() if k[:3]==key}
                if points!=event['target_points_angstrom']:
                    raise ValueError('target/full-stage coordinates differ')
                audit_points(points,event['ca_volume_angstrom3'],case['residue_name'],event['target_cip'])
                stage_count+=1
                seed_count+=event['stage']=='returned-seed'
            data=transport['anchor_data']
            template_raw=(Path(pf.__file__).parent/data['resource']).read_bytes()
            if sha(template_raw)!=data['resource_sha256']:
                raise ValueError('canonical fit template changed')
            resource=read_pdb(template_raw)
            target_observed={k[3]:p for k,p in observed.items() if k[:3]==key}
            names=[a['key'][3] for a in resource if a['key'][3] in target_observed]
            expected_p=[[x/10. for x in target_observed[n]] for n in names]
            expected_q=[[x/10. for x in a['xyz']] for a in resource if a['key'][3] in target_observed]
            if data['names']!=names or data['reference_points_nm']!=expected_p or data['template_points_nm']!=expected_q:
                raise ValueError('fit anchors not bound to observed source/template')
            if len(transport['calls'])!=len(trace['overlays']):
                raise ValueError('used fit inventory differs')
            for used,supplied in zip(transport['calls'],trace['overlays'],strict=True):
                fit_check(used,supplied,data['names'],data['reference_points_nm'],data['template_points_nm'],'ncac')
                fit_count+=1
            for step in trace['minimizer_iterations']:
                audit_points(step['target_points_angstrom'],step['ca_volume_angstrom3'])
                if step['observed_max_movement_angstrom']>1e-6:
                    raise ValueError('NCAC minimizer moved observations')
            parent=[e for e in trace['events'] if e['stage']=='returned-seed'][-1]['target_points_angstrom']
            ncac=json.loads((folder/'ncac-proposals.json').read_bytes())['candidates']
            if [c['relative_angle_degrees'] for c in ncac]!=[0.,120.,-120.]:
                raise ValueError('relative-angle pool changed')
            for p in ncac:
                check_rotation(parent,p['points'],p['relative_angle_degrees'])
            proposals['ncac']=[p['points'] for p in ncac]
        else:
            cell_ledger['ncac_blocker']=result
        method_replay={}
        for arm,points_pool in proposals.items():
            identifier=name+'-'+arm
            declaration=score_cases[identifier]
            row=scored['results'][identifier]
            feasible,audits,materialized=[],[],[]
            for i,points in enumerate(points_pool):
                added=additions(request,points)
                common,pdb=measure(source,request,added,anchored=False)
                filename=f'{identifier}/own-{i}.pdb'
                if pdb!=(score_root/filename).read_bytes() or not source_records_preserved(source,pdb,declaration['generated']):
                    raise ValueError('candidate/source materialization differs')
                if arm=='ncac' and canonical(asdict(common))!=ncac[i]['common_audit']:
                    raise ValueError('NCAC common gate differs')
                audits.append(dict(index=i,common_audit=asdict(common)))
                materialized.append((added,pdb))
                if common.accepted:
                    feasible.append(i)
            if [c['index'] for c in declaration['cells']]!=feasible or [c['index'] for c in row['cells']]!=feasible:
                raise ValueError('complete feasible score inventory differs')
            contexts=[]
            for c,d in zip(row['cells'],declaration['cells'],strict=True):
                if any(c[k]!=v for k,v in d.items()):
                    raise ValueError('scored declaration changed')
                contexts.append(snapshot_check(score_root,score_run_root,declaration,c))
                score_statuses[c['score_status']]=score_statuses.get(c['score_status'],0)+1
                snapshots+=1
            selection=choose_scored_own(row['cells'])
            if selection!=row['selection'] or selection!=json.loads((score_run_root/identifier/'selection-before-native.json').read_bytes()):
                raise ValueError('independent selection differs')
            cell_ledger['methods'][arm]=dict(candidates=audits,selection=selection,
                context_coordinates_identical=not contexts or all(c==contexts[0] for c in contexts))
            method_replay[arm]=materialized
        ledger[name]=cell_ledger
        replay[name]=(request,source,method_replay)
    selection_path=output.with_name(output.stem+'-selection.json')
    if output.exists() or selection_path.exists():
        raise FileExistsError('fresh publication destinations required')
    write_json(selection_path,dict(schema='dvbfixer.prospective-ncac-selection.v1',cases=ledger,
        native_coordinates_read=False,validated_snapshots=snapshots,stage_snapshots=stage_count,
        ncac_seed_states=seed_count,used_ncac_fits=fit_count,mask_plan_sha256=sha(masks_raw),
        score_statuses=score_statuses,
        science_run_sha256=sha(science_raw),score_run_sha256=sha(scored_raw),audit_script_sha256=sha(Path(__file__).read_bytes())))
    # All selections now committed. Only reporting below reads reference coordinates.
    reports={}
    for case in masks['cases']:
        name=case['id']
        if name not in replay:
            reports[name]=ledger[name]
            continue
        request,source,methods=replay[name]
        raw=(mask_root/case['reference']).read_bytes()
        if sha(raw)!=masks['artifacts'][case['reference']]:
            raise ValueError('native reference changed')
        truth={a.name:a.position for a in _parse_atom_lines(raw.decode().splitlines(keepends=True)) if _component_ref(a)==request.residue}
        native_audit,_=measure(source,request,additions(request,truth),anchored=False)
        rows={}
        for arm,pool in methods.items():
            selected=ledger[name]['methods'][arm]['selection']['selected_index']
            rows[arm]=dict(selected_index=selected,context_coordinates_identical=ledger[name]['methods'][arm]['context_coordinates_identical'],
                selected_accuracy=None if selected is None else accuracy(request,*pool[selected],truth),
                zero_accuracy=accuracy(request,*pool[0],truth))
        reports[name]=dict(status='evaluated',methods=rows,native_common_audit=asdict(native_audit))
    mutual=[r for r in reports.values() if r['status']=='evaluated' and all(
        r['methods'].get(a,{}).get('selected_index') is not None and r['methods'][a]['context_coordinates_identical'] for a in ('baseline','ncac'))]
    ncac_selected=[r for r in reports.values() if r['status']=='evaluated' and r['methods'].get('ncac',{}).get('selected_index') is not None]
    aggregate=dict(planned=len(masks['cases']),shared_refusals=sum(r['status']=='shared_admission_refused' for r in reports.values()),
        ncac_blockers=sum(r.get('ncac_status')=='execution_or_protocol_blocker' for r in ledger.values()),
        validated_snapshots=snapshots,mutual_selected_cases=len(mutual),
        score_statuses=score_statuses,
        common_selected={a:sum(r['status']=='evaluated' and r['methods'].get(a,{}).get('selected_index') is not None for r in reports.values()) for a in ('baseline','ncac')},
        native_common_passes=sum(r.get('native_common_audit',{}).get('accepted',False) for r in reports.values()),
        ncac_changed_choices=sum(r['methods']['ncac']['selected_index']!=0 for r in ncac_selected),
        ncac_zero_rounded_exact_median=statistics.median(r['methods']['ncac']['zero_accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in ncac_selected) if ncac_selected else None,
        mutual_selected_rounded_exact_medians={a:statistics.median(r['methods'][a]['selected_accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in mutual) if mutual else None for a in ('baseline','ncac')})
    write_json(output,dict(schema='dvbfixer.prospective-ncac-audit.v1',cases=reports,aggregate=aggregate,
        selection_ledger_name=selection_path.name,selection_ledger_sha256=sha(selection_path.read_bytes()),
        mask_plan_sha256=sha(masks_raw),science_run_sha256=sha(science_raw),score_run_sha256=sha(scored_raw),
        audit_script_sha256=sha(Path(__file__).read_bytes()),production_changed=False,scientific_qualification=False))
    print(json.dumps(aggregate,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--masks',type=Path,required=True)
    parser.add_argument('--science-run',type=Path,required=True)
    parser.add_argument('--score-run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.masks,args.science_run,args.score_run,args.output)
