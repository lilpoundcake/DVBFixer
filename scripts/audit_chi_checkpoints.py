"""Saved-only native derivative-check checkpoint audit; no inferred zero rows.

Report every exact/changed replay. Native L2 checkpoints are not a full function,
iteration or return-branch trace and do not certify final endpoint convergence.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_generated_chi_relax import mobility_check
from scripts.research_dunbrack_prior import sha, write


def checkpoint_bounds(norm,dimensions,cutoff):
    if not math.isfinite(norm) or norm<0 or dimensions<1:
        raise ValueError('invalid native analytic norm/dimension')
    lower=norm/math.sqrt(dimensions)
    return dict(component_max_lower_bound=lower,component_max_upper_bound=norm,
        certainly_at_or_below_gmax=norm<=cutoff,certainly_above_gmax=lower>cutoff)


def audit(plan_path,run_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen native checkpoint input changed')
    if (run_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('native checkpoint executed plan differs')
    index_raw=(run_root/'checkpoint-evidence.json').read_bytes()
    index=json.loads(index_raw)
    run=json.loads((run_root/'evidence.json').read_bytes())
    if index['plan_sha256']!=sha(raw) or not index['hook_restored'] or index['native_coordinates_read'] or index['weights_modified']:
        raise ValueError('native checkpoint run contract differs')
    if run['plan_sha256']!=sha(raw) or run['native_coordinates_read'] or run['weights_modified']:
        raise ValueError('native checkpoint underlying run contract differs')
    inventory=[(case,cell) for case in plan['cases'] for cell in case['cells']]
    if len(index['comparisons'])!=len(inventory):
        raise ValueError('checkpoint case inventory incomplete')
    endpoint_record=Path('docs/research/chi-derivatives-evidence.json').read_bytes()
    endpoints={(r['id'],r['index']):r for r in json.loads(endpoint_record)['endpoints']}
    root=Path(plan['checkpoint_protocol']['parent_run'])
    reports=[]
    protected=[]
    for (case,cell),comparison in zip(inventory,index['comparisons'],strict=True):
        key=(case['id'],cell['declaration']['index'])
        if key!=(comparison['id'],comparison['index']):
            raise ValueError('checkpoint endpoint identity differs')
        folder=run_root/case['id']/cell['declaration']['id']
        parent=root/case['id']/cell['declaration']['id']
        exact={name:(parent/name).exists() and (folder/name).exists() and (parent/name).read_bytes()==(folder/name).read_bytes()
               for name in ('before-atoms.json','after-atoms.json','result.json','pose-snapshot.json.gz')}
        if comparison['exact']!=exact or comparison['endpoint_exact']!=all(exact.values()):
            raise ValueError('independent exact replay comparison differs')
        retained=json.loads((folder/'result.json').read_bytes())
        declared=next(r for r in run['results'][case['id']]['cells'] if r['index']==key[1])
        if retained!=declared:
            raise ValueError('retained native checkpoint result differs')
        if retained['status']!='refined' or comparison['checkpoint_sha256'] is None:
            reports.append(dict(id=key[0],index=key[1],status='execution_blocker',exact=exact,error=retained.get('error')))
            continue
        before=json.loads((folder/'before-atoms.json').read_bytes())['atoms']
        after=json.loads((folder/'after-atoms.json').read_bytes())['atoms']
        motion=mobility_check(before,after,case['case'],plan['protected_numerical_bound_angstrom'])
        protected.append(motion['protected_max_motion_angstrom'])
        trace_raw=(folder/'checkpoints.json').read_bytes()
        if sha(trace_raw)!=comparison['checkpoint_sha256']:
            raise ValueError('native checkpoint trace digest differs')
        trace=json.loads(trace_raw)
        if not trace['deriv_check'] or trace['deriv_check_to_stdout'] or trace['python_virtual_callback']:
            raise ValueError('native checker diagnostic flags differ')
        if trace['settings']!=plan['minimizer']:
            raise ValueError('native checker numerical settings changed')
        rows=trace['checkpoints']
        if not rows or [r['checkpoint'] for r in rows]!=list(range(1,len(rows)+1)):
            raise ValueError('native derivative checkpoint sequence missing')
        dimensions=len(endpoints[key]['chis'])
        bounds=[]
        for row in rows:
            if row['dimensions']!=dimensions or row['steps']!=5:
                raise ValueError('native derivative checkpoint vector inventory differs')
            norm=row['analytic_l2_norm_per_degree']
            bound=checkpoint_bounds(norm,dimensions,plan['minimizer']['gradient_cutoff'])
            bounds.append(bound)
            sparse={}
            for point in row['sparse_populated_points']:
                pair=(point['vector_index'],point['step_index'])
                if pair in sparse or not 1<=pair[0]<=dimensions or not 1<=pair[1]<=5:
                    raise ValueError('duplicate/invalid sparse native diagnostic row')
                values=[point[k] for k in ('analytic_per_degree','numeric_per_degree','dof_value','f00','f11','f22')]
                if not all(math.isfinite(v) for v in values):
                    raise ValueError('nonfinite populated native diagnostic row')
                if not(abs(point['analytic_per_degree'])>.001 or abs(point['numeric_per_degree'])>.001):
                    raise ValueError('native sparse-row population policy differs')
                step=.0005*2**pair[1]
                numeric=(point['f11']-point['f22'])/(2*step)
                if not math.isclose(numeric,point['numeric_per_degree'],rel_tol=1e-9,abs_tol=1e-9):
                    raise ValueError('native sparse finite-difference arithmetic differs')
                sparse[pair]=point
            # Sparse defaults are absent, not numerical zeros. Check only present entries.
            analytic_by_dof={}
            for point in sparse.values():
                dof=point['vector_index']
                analytic=point['analytic_per_degree']
                if dof in analytic_by_dof and analytic_by_dof[dof]!=analytic:
                    raise ValueError('native per-step analytic derivative differs')
                analytic_by_dof[dof]=analytic
            if math.sqrt(sum(v*v for v in analytic_by_dof.values()))>norm+1e-10:
                raise ValueError('populated analytic subset exceeds complete L2 norm')
        last=rows[-1]['analytic_l2_norm_per_degree']
        endpoint=endpoints[key]
        reports.append(dict(id=key[0],index=key[1],status='audited',endpoint_exact=all(exact.values()),exact=exact,
            checkpoint_count=len(rows),first_analytic_l2_norm_per_degree=rows[0]['analytic_l2_norm_per_degree'],
            last_analytic_l2_norm_per_degree=last,last_component_bounds=bounds[-1],
            endpoint_analytic_maximum_per_degree=endpoint['maximum_analytic_derivative_per_degree'],
            endpoint_above_mode_c_gmax=endpoint['above_mode_c_gmax'],
            explicit_max200_warning=endpoint['explicit_max200_warning'],movement=motion,
            native_trace_sha256=sha(trace_raw)))
    audited=[r for r in reports if r['status']=='audited']
    summary=dict(starts=len(reports),audited=len(audited),execution_blockers=len(reports)-len(audited),
        exact_endpoint_repeats=sum(r['endpoint_exact'] for r in audited),
        changed_endpoint_replays=sum(not r['endpoint_exact'] for r in audited),
        derivative_checkpoints=sum(r['checkpoint_count'] for r in audited),
        exact_above_endpoint_but_last_norm_below_gmax=sum(r['endpoint_exact'] and r['endpoint_above_mode_c_gmax'] and r['last_component_bounds']['certainly_at_or_below_gmax'] for r in audited),
        protected_max_motion_angstrom=max(protected,default=None))
    write(output,dict(schema='dvbfixer.chi-checkpoints-audit.v1',summary=summary,cases=reports,
        plan_sha256=sha(raw),checkpoint_index_sha256=sha(index_raw),
        endpoint_derivative_evidence_sha256=sha(endpoint_record),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,native_accuracy_reported=False,profile_selected=False,
        full_function_trace=False,return_branch_recorded=False,convergence_certified=False,production_changed=False,
        interpretation='native derivative-check norm sequence, exact-versus-changed replay retained, sparse rows absent not zero'))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
