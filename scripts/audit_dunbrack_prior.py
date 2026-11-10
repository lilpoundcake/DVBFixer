"""Read-only source/selection/API-repeat audit, then native reporting.

Public evidence excludes query samples, means, sigmas and per-family weights.
Those bounded research responses remain private; no bundled table is published.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path

from scripts.audit_ncac_expansion import accuracy
from scripts.compare_branched_completion import measure
from scripts.research_dunbrack_prior import sha, write

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.backbone import backbone_angles
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def ranking(candidates,samples):
    weights={s['well']:s['probability'] for s in samples}
    if len(samples)!=3 or set(weights)!={1,2,3} or any(not math.isfinite(v) or v<0 for v in weights.values()) or sum(weights.values())<=0:
        raise ValueError('invalid raw family weights')
    if [c['index'] for c in candidates]!=[0,1,2] or {c['well'] for c in candidates}!={1,2,3}:
        raise ValueError('candidate/well inventory differs')
    feasible=sorted((c for c in candidates if c['common_accepted']),key=lambda c:(-weights[c['well']],c['index']))
    return feasible[0]['index'] if feasible else None


def audit(plan_path,run_root,repeat_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen inputs/implementation changed')
    for name in ('executed-plan.json','environment.json','selection-before-native.json','weight-semantics.json'):
        if (run_root/name).read_bytes()!=(repeat_root/name).read_bytes():
            raise ValueError('complete independent API repeat not byte-identical')
    if (run_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('executed plan differs')
    env_raw=(run_root/'environment.json').read_bytes()
    env=json.loads(env_raw)
    result_raw=(run_root/'selection-before-native.json').read_bytes()
    result=json.loads(result_raw)
    semantics=json.loads((run_root/'weight-semantics.json').read_bytes())
    if result['plan_sha256']!=sha(raw) or result['environment_sha256']!=sha(env_raw) or result['native_coordinates_read']:
        raise ValueError('selection input binding differs')
    if semantics['normalized_weights_used_for_selection'] or not semantics['selection_adapter_restored'] or semantics['native_coordinates_read']:
        raise ValueError('weight semantics differ')
    if set(result['results'])!={c['id'] for c in plan['cases']}:
        raise ValueError('case inventory differs')
    old_plan=json.loads(Path('tmp/pdbfixer-research/ncac-expansion-plan-01.json').read_bytes())
    new_plan=json.loads(Path('tmp/pdbfixer-research/ncac-new-system-masks-01/plan.json').read_bytes())
    old_candidates=json.loads(Path('tmp/pdbfixer-research/ncac-rotamers-score-plan-01/geometry.json').read_bytes())
    original_cases={c['id']:c for c in old_plan['cases']+new_plan['cases']}
    ledger,replay,totals={}, {}, []
    for declaration in plan['cases']:
        name=declaration['id']
        row=result['results'][name]
        if declaration['status']!='query_ready':
            if row['status']!=declaration['status']:
                raise ValueError('parent refusal status differs')
            ledger[name]=dict(group=declaration['group'],status=row['status'],parent_status=row.get('parent_status'),reason=row.get('reason'))
            continue
        source=Path(declaration['source_path']).read_bytes()
        request=plan_branched_source(source,ComponentInstanceRef(**declaration['target'])).request
        phi,psi=backbone_angles(request)
        if phi!=declaration['phi_degrees'] or psi!=declaration['psi_degrees'] or row['phi_degrees']!=phi or row['psi_degrees']!=psi:
            raise ValueError('API backbone input not source-bound')
        expected_chi=['N','CA','CB','CG1' if request.residue_name=='VAL' else 'OG1']
        if row['chi_atoms']!=expected_chi or len(row['periodic_controls'])!=2:
            raise ValueError('chi convention or periodic control count differs')
        for control in row['periodic_controls']:
            if len(control)!=3:
                raise ValueError('incomplete periodic library view')
            for a,b in zip(row['samples'],control,strict=True):
                if a['well']!=b['well'] or any(not math.isclose(a[k],b[k],abs_tol=1e-6,rel_tol=0.) for k in ('probability','mean_degrees','sigma_degrees')):
                    raise ValueError('periodic API repeat differs')
        selected=ranking(row['candidates'],row['samples'])
        if selected!=row['selected_index']:
            raise ValueError('independent feasible raw-weight choice differs')
        total=sum(s['probability'] for s in row['samples'])
        if not math.isclose(semantics['cases'][name]['raw_total'],total,rel_tol=0.,abs_tol=1e-14):
            raise ValueError('raw weight total differs')
        for sample in row['samples']:
            displayed=semantics['cases'][name]['normalized_display_weights'][str(sample['well'])]
            if not math.isclose(displayed,sample['probability']/total,rel_tol=0.,abs_tol=1e-14):
                raise ValueError('diagnostic normalized display differs')
        totals.append(total)
        if declaration['group']=='old_eleven':
            proposals=old_candidates['cases'][name]['candidates']
        else:
            proposals=json.loads((Path(declaration['source_path']).parent/'ncac-proposals.json').read_bytes())['candidates']
        materialized=[]
        for i,c in enumerate(proposals):
            added=tuple(Addition(ExactAtomRef(request.residue,n),n[0],tuple(c['points'][n])) for n in request.missing_atoms)
            common,pdb=measure(source,request,added,anchored=False)
            if common.accepted!=row['candidates'][i]['common_accepted']:
                raise ValueError('common geometry feasibility differs')
            materialized.append((added,pdb))
        ledger[name]=dict(group=declaration['group'],status='queried',selected_index=selected,
                         common_feasible_indices=[c['index'] for c in row['candidates'] if c['common_accepted']],
                         chi_convention_verified=True,periodic_controls_pass=True)
        replay[name]=(request,source,materialized)
    ledger_path=output.with_name(output.stem+'-selection.json')
    if output.exists() or ledger_path.exists():
        raise FileExistsError('fresh audit/selection destinations required')
    write(ledger_path,dict(schema='dvbfixer.installed-dunbrack-prior-selection.v1',cases=ledger,
        plan_sha256=sha(raw),private_query_sha256=sha(result_raw),environment=env,native_coordinates_read=False,
        complete_api_repeat_byte_identical=True,table_data_included=False,audit_script_sha256=sha(Path(__file__).read_bytes())))
    # All independently verified selections commit before opening native reporting.
    old_reporting=json.loads(Path('docs/research/ncac-relative-rotamers-evidence.json').read_bytes())
    new_reporting=json.loads(Path('docs/research/ncac-new-system-evidence.json').read_bytes())
    reports={}
    for name,(request,source,pool) in replay.items():
        original=original_cases[name]
        group=ledger[name]['group']
        mask_root=Path(old_plan['mask_root']) if group=='old_eleven' else Path('tmp/pdbfixer-research/ncac-new-system-masks-01')
        mask_plan=json.loads((mask_root/'plan.json').read_bytes())
        reference=(mask_root/original['reference']).read_bytes()
        if sha(reference)!=mask_plan['artifacts'][original['reference']]:
            raise ValueError('native reporting reference digest differs')
        truth={a.name:a.position for a in _parse_atom_lines(reference.decode().splitlines(keepends=True)) if _component_ref(a)==request.residue}
        selected=ledger[name]['selected_index']
        prior_error=None if selected is None else accuracy(request,*pool[selected],truth)
        zero=accuracy(request,*pool[0],truth)
        if group=='old_eleven':
            parent_row=old_reporting['cases'][name]
            scored_index=parent_row['selected_index']
            scored_rmsd=parent_row['selected_rmsd']
        else:
            parent_row=new_reporting['cases'][name]['methods']['ncac']
            scored_index=parent_row['selected_index']
            scored_rmsd=parent_row['selected_accuracy']['rounded']['exact_name_rmsd_angstrom']
        scored_check=accuracy(request,*pool[scored_index],truth)['rounded']['exact_name_rmsd_angstrom']
        if scored_check!=scored_rmsd:
            raise ValueError('parent scored accuracy does not repeat saved candidate')
        reports[name]=dict(group=group,selected_index=selected,prior_accuracy=prior_error,zero_accuracy=zero,
            parent_scored_index=scored_index,parent_scored_rmsd=scored_rmsd)
    aggregate={}
    for group in ('old_eleven','new_inputs_eleven'):
        subset=[r for r in reports.values() if r['group']==group and r['selected_index'] is not None]
        aggregate[group]=dict(planned=sum(c['group']==group for c in plan['cases']),queried=len(subset),
            inherited_nonexecutions=sum(c['group']==group and c['status']!='query_ready' for c in plan['cases']),
            changed_from_zero=sum(r['selected_index']!=0 for r in subset),
            changed_from_parent_scored=sum(r['selected_index']!=r['parent_scored_index'] for r in subset),
            prior_rounded_exact_median=statistics.median(r['prior_accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in subset) if subset else None,
            zero_rounded_exact_median=statistics.median(r['zero_accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in subset) if subset else None,
            parent_scored_rounded_exact_median=statistics.median(r['parent_scored_rmsd'] for r in subset) if subset else None)
    write(output,dict(schema='dvbfixer.installed-dunbrack-prior-audit.v1',aggregate=aggregate,cases=reports,
        selection_ledger_name=ledger_path.name,selection_ledger_sha256=sha(ledger_path.read_bytes()),
        api_raw_weight_total_range=[min(totals),max(totals)],normalized_probability_claim=False,
        table_data_included=False,production_changed=False,scientific_qualification=False,
        audit_script_sha256=sha(Path(__file__).read_bytes())))
    print(json.dumps(aggregate,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--repeat',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.repeat,args.output)
