"""Independent saved bounded-guard policy, ownership and every-cycle geometry audit.

No optimization, PyRosetta, native reference or proposal/profile selection.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict
from pathlib import Path

from scripts.audit_generated_chi_relax import mobility_check
from scripts.compare_branched_completion import measure
from scripts.research_dunbrack_prior import sha, write

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def expected_decision(maximum,cycle,limit,cutoff):
    if not math.isfinite(maximum) or maximum<0 or cycle<1 or cycle>limit:
        raise ValueError('invalid guarded endpoint gradient/cycle')
    if maximum<=cutoff:
        return 'endpoint_gradient_met'
    if cycle==limit:
        return 'budget_exhausted'
    return 'continue'


def audit(plan_path,run_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen endpoint guard state changed')
    if (run_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('endpoint guard executed plan differs')
    index_raw=(run_root/'guard-evidence.json').read_bytes()
    index=json.loads(index_raw)
    run_raw=(run_root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    for record in (index,run):
        if record['plan_sha256']!=sha(raw) or record['native_coordinates_read'] or record['weights_modified']:
            raise ValueError('endpoint guard execution binding differs')
    if not index['hook_restored']:
        raise ValueError('endpoint guard adapter not restored')
    inventory=[(case,cell) for case in plan['cases'] for cell in case['cells']]
    if set(run['results'])!={case['id'] for case in plan['cases']}:
        raise ValueError('guarded run pool inventory differs')
    for case in plan['cases']:
        if [r['index'] for r in run['results'][case['id']]['cells']]!=[c['declaration']['index'] for c in case['cells']]:
            raise ValueError('guarded run complete start inventory differs')
    if index['enrolled_starts']!=len(inventory) or index['created_minimizers']!=len(inventory):
        raise ValueError('endpoint guard starting inventory differs')
    traces={(t['id'],t['index']):t for t in index['traces']}
    if len(traces)!=len(index['traces']):
        raise ValueError('duplicate endpoint guard trace')
    if not traces.keys()<={(case['id'],cell['declaration']['index']) for case,cell in inventory}:
        raise ValueError('unexpected endpoint guard trace identity')
    prior=json.loads(Path('docs/research/chi-checkpoint-components-evidence.json').read_bytes())
    dimensions={(r['id'],r['index']):len(r['last_analytic_components_per_degree']) for r in prior['cases']}
    settings=plan['guard_protocol']
    parent=Path(settings['parent_run'])
    reports=[]
    for case,cell in inventory:
        key=(case['id'],cell['declaration']['index'])
        folder=run_root/case['id']/cell['declaration']['id']
        baseline=parent/case['id']/cell['declaration']['id']
        result=json.loads((folder/'result.json').read_bytes())
        declared=next(c for c in run['results'][case['id']]['cells'] if c['index']==key[1])
        if result!=declared:
            raise ValueError('guarded retained final result differs')
        if result['status']!='refined' or key not in traces:
            reports.append(dict(id=key[0],index=key[1],status='execution_blocker',error=result.get('error')))
            continue
        trace=traces[key]
        if trace['actual_settings']!=plan['minimizer'] or not trace['first_cycle_exact']:
            raise ValueError('guarded inner settings or initial-cycle replay differ')
        for filename in ('before-atoms.json','pose-snapshot.json.gz'):
            if (folder/filename).read_bytes()!=(baseline/filename).read_bytes():
                raise ValueError('guarded original start not exact')
        original=json.loads((folder/'before-atoms.json').read_bytes())['atoms']
        if sha((folder/'before-atoms.json').read_bytes())!=result['before_sha256'] or sha((folder/'after-atoms.json').read_bytes())!=result['after_sha256']:
            raise ValueError('guarded final atom trace digest differs')
        source=Path(case['source']).read_bytes()
        target=case['case']['target_key']
        request=plan_branched_source(source,ComponentInstanceRef(1,target[0],str(target[1]),target[2].strip())).request
        cycles=[]
        phases=trace['phases']
        if not phases or len(phases)>settings['total_cycle_limit']:
            raise ValueError('guarded cycle budget violated')
        for cycle,descriptor in enumerate(phases,1):
            phase_raw=(folder/f'cycle-{cycle}.json').read_bytes()
            phase=json.loads(phase_raw)
            if sha(phase_raw)!=descriptor['phase_sha256'] or phase['cycle']!=cycle or descriptor['cycle']!=cycle:
                raise ValueError('guarded phase inventory/digest differs')
            atoms=phase['atoms']
            if cycle==1:
                if atoms!=json.loads((baseline/'after-atoms.json').read_bytes())['atoms'] or phase['raw_score']!=json.loads((baseline/'result.json').read_bytes())['final_raw_score']:
                    raise ValueError('first guarded inner cycle not exact mode C')
            motion=mobility_check(original,atoms,case['case'],plan['protected_numerical_bound_angstrom'])
            gradient=phase['gradient']
            values=gradient['components']
            if sorted(v['chi'] for v in values)!=list(range(1,dimensions[key]+1)) or not all(math.isfinite(v['analytic_per_degree']) for v in values):
                raise ValueError('guarded endpoint gradient inventory differs')
            maximum=max(abs(v['analytic_per_degree']) for v in values)
            if maximum!=gradient['maximum_per_degree'] or not gradient['original_pose_unchanged'] or max(gradient['mapped_motion'].values())>settings['mapped_endpoint_numerical_bound_angstrom']:
                raise ValueError('actual endpoint derivative scope differs')
            if not math.isfinite(phase['raw_score']) or abs(gradient['raw_score']-phase['raw_score'])>settings['center_score_absolute_bound'] or abs(gradient['minimizer_center_score']-gradient['raw_score'])>settings['center_score_absolute_bound']:
                raise ValueError('actual endpoint raw/gradient objective differs')
            decision=expected_decision(maximum,cycle,settings['total_cycle_limit'],settings['cutoff_per_degree'])
            if phase['decision']!=decision or descriptor['decision']!=decision or (cycle<len(phases) and decision!='continue') or (cycle==len(phases) and decision=='continue'):
                raise ValueError('bounded actual-endpoint continuation policy differs')
            points={a['key'][-1]:a['xyz'] for a in atoms if a['key'][:4]==target and not a['is_hydrogen'] and not a['is_virtual']}
            additions=tuple(Addition(ExactAtomRef(request.residue,name),name[0],tuple(points[name])) for name in request.missing_atoms)
            common,_=measure(source,request,additions,anchored=False)
            cycles.append(dict(cycle=cycle,raw_score=phase['raw_score'],gradient_maximum_per_degree=maximum,
                decision=decision,movement=motion,common_audit=asdict(common),phase_sha256=sha(phase_raw)))
        if trace['final_decision']!=cycles[-1]['decision']:
            raise ValueError('guarded final stop decision differs')
        final_phase=json.loads((folder/f'cycle-{len(phases)}.json').read_bytes())
        if final_phase['atoms']!=json.loads((folder/'after-atoms.json').read_bytes())['atoms'] or final_phase['raw_score']!=result['final_raw_score']:
            raise ValueError('guarded final result not bound to last audited phase')
        if not result['original_pose_exact'] or result['convergence_claim'] or result['coordinate_restore_performed']:
            raise ValueError('guarded scope/convergence/restore claim differs')
        if not result['graph_unchanged'] or not result['only_target_chi_enabled']:
            raise ValueError('guarded graph or target-only scope differs')
        weights=cell['saved']['score_terms']
        if set(weights)!=set(result['score_terms']) or any(v['weight']!=weights[k]['weight'] for k,v in result['score_terms'].items()):
            raise ValueError('guarded score weights changed')
        if not math.isclose(sum(t['raw']*t['weight'] for t in result['score_terms'].values()),result['final_raw_score'],rel_tol=1e-10,abs_tol=1e-8):
            raise ValueError('guarded final score terms do not sum to total')
        reports.append(dict(id=key[0],index=key[1],status='audited',first_cycle_exact=True,cycles=cycles,
            final_decision=cycles[-1]['decision'],final_common_accepted=cycles[-1]['common_audit']['accepted'],
            final_raw_score_delta_vs_mode_c=cycles[-1]['raw_score']-cycles[0]['raw_score']))
    audited=[r for r in reports if r['status']=='audited']
    all_cycles=[c for r in audited for c in r['cycles']]
    summary=dict(starts=len(reports),audited=len(audited),execution_blockers=len(reports)-len(audited),
        inner_cycles=len(all_cycles),extra_cycles=sum(len(r['cycles'])-1 for r in audited),
        starts_continued=sum(len(r['cycles'])>1 for r in audited),
        actual_endpoint_gradient_met=sum(r['final_decision']=='endpoint_gradient_met' for r in audited),
        budget_exhausted=sum(r['final_decision']=='budget_exhausted' for r in audited),
        final_common_passes=sum(r['final_common_accepted'] for r in audited),
        final_gradient_and_common_passes=sum(r['final_common_accepted'] and r['final_decision']=='endpoint_gradient_met' for r in audited),
        protected_max_motion_angstrom=max((c['movement']['protected_max_motion_angstrom'] for c in all_cycles),default=None),
        cycles_with_global_stereo_findings=sum(bool(c['common_audit']['global_stereo_findings']) for c in all_cycles))
    write(output,dict(schema='dvbfixer.chi-endpoint-guard-audit.v1',summary=summary,cases=reports,
        inherited_nonexecutions=plan['inherited_nonexecutions'],plan_sha256=sha(raw),guard_index_sha256=sha(index_raw),
        run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),native_coordinates_read=False,
        native_accuracy_reported=False,profile_selected=False,convergence_certified=False,production_changed=False,
        interpretation='bounded actual-endpoint-gradient criterion, separate every-cycle raw/rounded geometry, no ranking/native qualification'))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
