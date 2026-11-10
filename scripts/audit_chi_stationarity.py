"""Independent saved finite-chi sensitivity audit; no native data or new scores."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

from scripts.audit_generated_chi_relax import mobility_check
from scripts.research_dunbrack_prior import sha, write


def finite_summary(center,probes,chi_count,steps,bound):
    expected={(chi,h,sign) for chi in range(1,chi_count+1) for h in steps for sign in (-1,1)}
    lookup={(p['chi'],p['step_degrees'],p['sign']):p for p in probes}
    if len(lookup)!=len(probes) or lookup.keys()!=expected or chi_count<1:
        raise ValueError('complete signed probe inventory differs')
    if not math.isfinite(center) or any(not math.isfinite(p['score']) for p in probes):
        raise ValueError('nonfinite sensitivity score')
    rows=[]
    for chi in range(1,chi_count+1):
        slopes=[]
        for h in steps:
            minus=lookup[(chi,h,-1)]['score']
            plus=lookup[(chi,h,1)]['score']
            slopes.append(dict(step_degrees=h,slope_per_radian=(plus-minus)/(2*math.radians(h)),
                curvature_per_radian_squared=((plus-center)+(minus-center))/math.radians(h)**2,
                best_signed_energy_change=min(minus-center,plus-center),
                resolved_descent=min(minus-center,plus-center)<-bound))
        rows.append(dict(chi=chi,steps=slopes,slope_step_discrepancy=abs(slopes[0]['slope_per_radian']-slopes[1]['slope_per_radian'])))
    return dict(chis=rows,resolved_descent=any(s['resolved_descent'] for r in rows for s in r['steps']),
        maximum_absolute_slope_per_radian=max(abs(s['slope_per_radian']) for r in rows for s in r['steps']),
        best_signed_energy_change=min(s['best_signed_energy_change'] for r in rows for s in r['steps']))


def audit(plan_path,run_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen stationarity input changed')
    run_raw=(run_root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if (run_root/'executed-plan.json').read_bytes()!=raw or run['plan_sha256']!=sha(raw) or run['native_coordinates_read'] or run['weights_modified']:
        raise ValueError('stationarity execution binding differs')
    if [(r['id'],r['index']) for r in run['results']]!=[(r['id'],r['result']['index']) for r in plan['endpoints']]:
        raise ValueError('all endpoint enrollment differs')
    reports=[]
    inputs={str(plan_path):sha(raw),str(run_root/'evidence.json'):sha(run_raw)}
    max_protected=0.
    for declaration,result in zip(plan['endpoints'],run['results'],strict=True):
        folder=run_root/declaration['id']/declaration['cell']['declaration']['id']
        if json.loads((folder/'result.json').read_bytes())!=result:
            raise ValueError('retained endpoint result differs')
        if result['status']!='probed':
            reports.append(dict(id=result['id'],index=result['index'],status=result['status'],error=result['error']))
            continue
        if not result['endpoint_unchanged'] or not result['initial_exact_replay'] or not result['endpoint_atoms_exact']:
            raise ValueError('endpoint import/materialization scope differs')
        if (folder/'pose-snapshot.json.gz').read_bytes()!=Path(declaration['cell']['snapshot']).read_bytes():
            raise ValueError('original pose snapshot not exact')
        if result['saved_endpoint_score']!=declaration['result']['final_raw_score'] or abs(result['endpoint_score']-result['saved_endpoint_score'])>plan['materialized_score_absolute_bound']:
            raise ValueError('endpoint materialized score differs')
        endpoint=json.loads((Path(declaration['folder'])/'after-atoms.json').read_bytes())['atoms']
        endpoint_map={tuple(a['key']):a for a in endpoint}
        heavy_motion_by_chi={}
        for probe in result['probes']:
            path=folder/probe['atoms_file']
            probe_raw=path.read_bytes()
            if sha(probe_raw)!=probe['atoms_sha256']:
                raise ValueError('probe atom evidence changed')
            inputs[str(path)]=sha(probe_raw)
            atoms=json.loads(gzip.decompress(probe_raw))
            motion=mobility_check(endpoint,atoms,declaration['case'],plan['protected_bound_angstrom'])
            if any(not math.isclose(v,probe['movement'][k],rel_tol=0,abs_tol=1e-12) for k,v in motion.items()):
                raise ValueError('probe movement summary differs')
            max_protected=max(max_protected,motion['protected_max_motion_angstrom'])
            heavy_motion=max((math.dist(a['xyz'],endpoint_map[tuple(a['key'])]['xyz']) for a in atoms
                              if a['mobile'] and not a['is_hydrogen']),default=0.)
            heavy_motion_by_chi[probe['chi']]=max(heavy_motion_by_chi.get(probe['chi'],0.),heavy_motion)
        sensitivity=finite_summary(result['endpoint_score'],result['probes'],result['target_chi_count'],plan['steps_degrees'],plan['descent_reporting_bound'])
        for chi in sensitivity['chis']:
            chi['maximum_generated_heavy_probe_motion_angstrom']=heavy_motion_by_chi[chi['chi']]
            chi['heavy_moving_above_1e-9']=heavy_motion_by_chi[chi['chi']]>1e-9
        reports.append(dict(id=result['id'],index=result['index'],status='audited',sensitivity=sensitivity,
            materialized_score_difference=result['endpoint_score']-result['saved_endpoint_score'],
            target_chi_count=result['target_chi_count'],probes=len(result['probes'])))
    audited=[r for r in reports if r['status']=='audited']
    summary=dict(enrolled=len(reports),audited=len(audited),blockers=len(reports)-len(audited),
        probes=sum(r['probes'] for r in audited),target_chis=sum(r['target_chi_count'] for r in audited),
        endpoints_with_resolved_descent=sum(r['sensitivity']['resolved_descent'] for r in audited),
        endpoints_with_heavy_moving_chi_descent=sum(any(c['heavy_moving_above_1e-9'] and any(s['resolved_descent'] for s in c['steps'])
                                                     for c in r['sensitivity']['chis']) for r in audited),
        protected_max_motion_angstrom=max_protected,
        maximum_absolute_slope_per_radian=max((r['sensitivity']['maximum_absolute_slope_per_radian'] for r in audited),default=None),
        minimum_signed_energy_change=min((r['sensitivity']['best_signed_energy_change'] for r in audited),default=None),
        maximum_absolute_materialized_score_difference=max((abs(r['materialized_score_difference']) for r in audited),default=None))
    write(output,dict(schema='dvbfixer.chi-stationarity-audit.v1',summary=summary,endpoints=reports,
        plan_sha256=sha(raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        retained_probe_binding_sha256=sha(json.dumps(inputs,sort_keys=True).encode()),
        native_coordinates_read=False,native_accuracy_reported=False,profile_selected=False,
        convergence_certified=False,analytic_derivatives_validated=False,production_changed=False,table_data_included=False))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
