"""Native-free stopping-policy readback: no profile selection or reweighting.

Read all start/finish traces, enforce protected ownership and raw/rounded geometry,
then report every predeclared profile. No native coordinates or accuracy files.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from dataclasses import asdict
from pathlib import Path

from scripts.audit_generated_chi_relax import mobility_check, snapshot_atoms
from scripts.compare_branched_completion import measure
from scripts.research_dunbrack_prior import sha, write

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def select_complete(rows):
    if any(r['status']!='audited' for r in rows):
        return None
    feasible=[r for r in rows if r['common_audit']['accepted']]
    return min(feasible,key=lambda r:(r['final_raw_score'],r['index']))['index'] if feasible else None


def audit(plans_root,runs_prefix,output):
    suite_raw=(plans_root/'suite.json').read_bytes()
    suite=json.loads(suite_raw)
    if suite['native_accuracy_used'] or suite['cross_profile_selection']:
        raise ValueError('native-free suite contract differs')
    baseline_path=Path('docs/research/generated-chi-relax-evidence-selection.json')
    baseline_raw=baseline_path.read_bytes()
    baseline=json.loads(baseline_raw)
    if baseline['native_coordinates_read'] or baseline['plan_sha256']!=suite['baseline_plan_sha256']:
        raise ValueError('baseline ledger contains native reporting')
    profiles={}
    for declaration in suite['profiles']:
        plan_path=plans_root/declaration['filename']
        raw=plan_path.read_bytes()
        if sha(raw)!=declaration['sha256']:
            raise ValueError('suite profile binding differs')
        plan=json.loads(raw)
        for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
            if sha(Path(p).read_bytes())!=d:
                raise ValueError('frozen stopping input changed')
        run_root=Path(str(runs_prefix)+declaration['id']+'-run-01')
        run_raw=(run_root/'evidence.json').read_bytes()
        run=json.loads(run_raw)
        if run['plan_sha256']!=sha(raw) or (run_root/'executed-plan.json').read_bytes()!=raw or run['native_coordinates_read'] or run['weights_modified']:
            raise ValueError('stopping profile execution contract differs')
        settings=json.loads((run_root/'stopping-control-settings.json').read_bytes())
        if not settings['hook_restored'] or settings['native_coordinates_read'] or settings['weights_modified'] or settings['cross_profile_selection']:
            raise ValueError('control adapter not restored or scope differs')
        expected_starts=sum(len(c['cells']) for c in plan['cases'])
        if len(settings['actual_minimizers'])!=expected_starts:
            raise ValueError('runtime minimizer inventory differs')
        for actual in settings['actual_minimizers']:
            if any(actual[k]!=plan['minimizer'][k] for k in actual):
                raise ValueError('actual stopping policy differs from frozen plan')
        if set(run['results'])!={c['id'] for c in plan['cases']}:
            raise ValueError('profile pool omitted')
        pools={}
        energy_deltas=[]
        for case in plan['cases']:
            name=case['id']
            results=run['results'][name]['cells']
            if [c['index'] for c in results]!=[c['declaration']['index'] for c in case['cells']]:
                raise ValueError('profile start inventory differs')
            source=Path(case['source']).read_bytes()
            key=case['case']['target_key']
            request=plan_branched_source(source,ComponentInstanceRef(1,key[0],str(key[1]),key[2].strip())).request
            reports=[]
            for original,result in zip(case['cells'],results,strict=True):
                folder=run_root/name/original['declaration']['id']
                if json.loads((folder/'result.json').read_bytes())!=result or result['coordinate_restore_performed']:
                    raise ValueError('retained profile result differs/restoration used')
                if result['status']!='refined':
                    reports.append(dict(index=result['index'],status=result['status'],error=result['error']))
                    continue
                old=Path(original['snapshot']).read_bytes()
                if old!=(folder/'pose-snapshot.json.gz').read_bytes() or not result['original_pose_exact']:
                    raise ValueError('complete original start not exactly reproduced')
                before_raw=(folder/'before-atoms.json').read_bytes()
                after_raw=(folder/'after-atoms.json').read_bytes()
                if sha(before_raw)!=result['before_sha256'] or sha(after_raw)!=result['after_sha256']:
                    raise ValueError('profile atom trace digest differs')
                before=json.loads(before_raw)['atoms']
                after=json.loads(after_raw)['atoms']
                if {tuple(a['key']):a['xyz'] for a in before}!=snapshot_atoms(json.loads(gzip.decompress(old))):
                    raise ValueError('profile initial atom trace not original-state bound')
                moves=mobility_check(before,after,case['case'],plan['protected_numerical_bound_angstrom'])
                if result['initial_score']!=original['saved']['total_score'] or not result['graph_unchanged'] or not result['only_target_chi_enabled'] or result['convergence_claim']:
                    raise ValueError('profile score/scope/replay claim differs')
                if set(result['score_terms'])!=set(original['saved']['score_terms']) or any(t['weight']!=original['saved']['score_terms'][n]['weight'] for n,t in result['score_terms'].items()):
                    raise ValueError('profile changed score weights/terms')
                if not math.isfinite(result['final_raw_score']) or not math.isclose(sum(t['raw']*t['weight'] for t in result['score_terms'].values()),result['final_raw_score'],rel_tol=1e-10,abs_tol=1e-8):
                    raise ValueError('profile final score accounting differs')
                points={a['key'][-1]:a['xyz'] for a in after if a['key'][:4]==key and not a['is_hydrogen'] and not a['is_virtual']}
                added=tuple(Addition(ExactAtomRef(request.residue,n),n[0],tuple(points[n])) for n in request.missing_atoms)
                common,pdb=measure(source,request,added,anchored=False)
                before_map={tuple(a['key']):a for a in before}
                heavy_motion=max((math.dist(a['xyz'],before_map[tuple(a['key'])]['xyz']) for a in after if a['mobile'] and not a['is_hydrogen']),default=0.)
                prior_row=next(c for c in baseline['cases'][name]['candidates'] if c['index']==result['index'])
                delta=result['final_raw_score']-prior_row['final_raw_score']
                energy_deltas.append(delta)
                reports.append(dict(index=result['index'],status='audited',common_audit=asdict(common),movement=moves,
                    heavy_max_motion_angstrom=heavy_motion,final_raw_score=result['final_raw_score'],score_delta_vs_baseline=delta,
                    rounded_heavy_output_sha256=sha(pdb)))
            selected=select_complete(reports)
            baseline_index=baseline['cases'][name]['selected_index']
            baseline_start=next(c for c in case['cells'] if c['declaration']['index']==baseline_index)
            baseline_candidate=next(c for c in baseline['cases'][name]['candidates'] if c['index']==baseline_index)
            if not baseline_candidate['rounded_heavy_output_equals_start']:
                raise ValueError('baseline selected-output hash requires saved unchanged-start control')
            baseline_hash=sha((Path(plan['score_plan_root'])/baseline_start['declaration']['pdb']).read_bytes())
            selected_hash=next((c['rounded_heavy_output_sha256'] for c in reports if c['index']==selected),None)
            pools[name]=dict(group=case['group'],candidates=reports,selected_index=selected,
                baseline_selected_index=baseline_index,selected_rounded_output_sha256=selected_hash,
                baseline_rounded_output_sha256=baseline_hash)
        candidates=[c for p in pools.values() for c in p['candidates']]
        audited=[c for c in candidates if c['status']=='audited']
        summary=dict(starts=len(candidates),audited=len(audited),execution_blockers=len(candidates)-len(audited),
            common_passes=sum(c['common_audit']['accepted'] for c in audited),
            protected_max_motion_angstrom=max((c['movement']['protected_max_motion_angstrom'] for c in audited),default=None),
            selected_pools=sum(p['selected_index'] is not None for p in pools.values()),
            changed_selected_indices=sum(p['selected_index']!=p['baseline_selected_index'] for p in pools.values()),
            changed_selected_rounded_outputs=sum(p['selected_rounded_output_sha256']!=p['baseline_rounded_output_sha256'] for p in pools.values()),
            final_raw_score_delta_vs_baseline_range=[min(energy_deltas),max(energy_deltas)] if energy_deltas else None,
            maximum_heavy_motion_angstrom=max((c['heavy_max_motion_angstrom'] for c in audited),default=None))
        profiles[declaration['id']]=dict(summary=summary,cases=pools,minimizer_settings=plan['minimizer'],
            plan_sha256=sha(raw),run_sha256=sha(run_raw),settings_sha256=sha((run_root/'stopping-control-settings.json').read_bytes()),
            inherited_nonexecutions=plan['inherited_nonexecutions'])
    sequential=[]
    order=[p['id'] for p in suite['profiles']]
    for previous,current in zip(order,order[1:]):
        a,b=profiles[previous]['cases'],profiles[current]['cases']
        score_changes=[]
        for name,pool in a.items():
            other={c['index']:c for c in b[name]['candidates']}
            for c in pool['candidates']:
                d=other[c['index']]
                if c['status']=='audited' and d['status']=='audited':
                    score_changes.append(abs(c['final_raw_score']-d['final_raw_score']))
        sequential.append(dict(previous=previous,current=current,compared_states=len(score_changes),
            maximum_absolute_score_change=max(score_changes,default=None),
            changed_selected_indices=sum(p['selected_index']!=b[n]['selected_index'] for n,p in a.items()),
            changed_selected_rounded_outputs=sum(p['selected_rounded_output_sha256']!=b[n]['selected_rounded_output_sha256'] for n,p in a.items())))
    write(output,dict(schema='dvbfixer.chi-stopping-controls-audit.v1',profiles=profiles,sequential_comparisons=sequential,
        suite_sha256=sha(suite_raw),baseline_selection_sha256=sha(baseline_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,native_accuracy_reported=False,cross_profile_selection=False,
        convergence_certified=False,weights_modified=False,production_changed=False,table_data_included=False))
    print(json.dumps({k:v['summary'] for k,v in profiles.items()},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plans',type=Path,required=True)
    parser.add_argument('--runs-prefix',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plans,args.runs_prefix,args.output)
