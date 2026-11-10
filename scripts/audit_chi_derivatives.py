"""Independent saved-only analytic/finite derivative and unit audit.

No PyRosetta, optimization, native accuracy, profile selection or reweighting.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

from scripts.audit_generated_chi_relax import mobility_check
from scripts.research_dunbrack_prior import sha, write


def compare(analytic,numeric,absolute,relative):
    if not all(math.isfinite(v) for v in (analytic,numeric)):
        raise ValueError('nonfinite derivative readback')
    error=abs(analytic-numeric)
    bound=absolute+relative*max(abs(analytic),abs(numeric))
    return dict(absolute_error_per_degree=error,bound_per_degree=bound,within_reporting_bound=error<=bound)


def audit(plan_path,run_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen derivative input changed')
    if (run_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('derivative executed plan differs')
    run_raw=(run_root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(raw) or any(run[k] for k in ('native_coordinates_read','weights_modified','minimizer_executed')):
        raise ValueError('derivative execution scope differs')
    rows=run['results']
    if [(r['id'],r['index']) for r in rows]!=[(r['id'],r['result']['index']) for r in plan['endpoints']]:
        raise ValueError('complete derivative enrollment differs')
    finite=json.loads(Path('docs/research/chi-stationarity-evidence.json').read_bytes())
    original={(r['id'],r['index']):r for r in finite['endpoints']}
    termination=json.loads(Path('docs/research/chi-termination-evidence.json').read_bytes())
    logs={(r['id'],r['index']):r for r in termination['endpoints']}
    settings=plan['derivative_protocol']
    reports=[]
    protected=[]
    for declaration,row in zip(plan['endpoints'],rows,strict=True):
        folder=run_root/row['id']/declaration['cell']['declaration']['id']
        if json.loads((folder/'result.json').read_bytes())!=row:
            raise ValueError('retained derivative result differs')
        if row['status']!='evaluated':
            reports.append(dict(id=row['id'],index=row['index'],status=row['status'],error=row['error']))
            continue
        if not all(row[k] for k in ('initial_exact_replay','endpoint_atoms_exact','endpoint_unchanged')):
            raise ValueError('derivative endpoint replay differs')
        if (folder/'pose-snapshot.json.gz').read_bytes()!=Path(declaration['cell']['snapshot']).read_bytes():
            raise ValueError('initial derivative pose snapshot differs')
        endpoint=json.loads((Path(declaration['folder'])/'after-atoms.json').read_bytes())['atoms']
        mapped_raw=(folder/'mapped-endpoint-atoms.json.gz').read_bytes()
        if sha(mapped_raw)!=row['mapped_endpoint_atoms_sha256']:
            raise ValueError('mapped endpoint digest differs')
        mapped_atoms=json.loads(gzip.decompress(mapped_raw))
        motion=mobility_check(endpoint,mapped_atoms,declaration['case'],plan['protected_bound_angstrom'])
        if max(motion.values())>settings['endpoint_atoms_bound_angstrom']:
            raise ValueError('base DOF mapping moved endpoint beyond declared bound')
        protected.append(motion['protected_max_motion_angstrom'])
        if abs(row['minimizer_center_score']-row['endpoint_score'])>plan['materialized_score_absolute_bound'] or abs(row['endpoint_score']-declaration['result']['final_raw_score'])>plan['materialized_score_absolute_bound']:
            raise ValueError('derivative center objective differs')
        count=len(row['mapping'])
        if count!=len(row['analytic_derivatives_per_degree']) or count!=len(row['variables']):
            raise ValueError('analytic-vector mapping inventory differs')
        if sorted(m['chi'] for m in row['mapping'])!=list(range(1,count+1)) or [m['vector_index'] for m in row['mapping']]!=list(range(1,count+1)):
            raise ValueError('chi/vector mapping incomplete or duplicated')
        probes={(p['vector_index'],p['step_degrees'],p['sign']):p for p in row['probes']}
        expected={(i,h,sign) for i in range(1,count+1) for h in settings['numeric_steps_degrees'] for sign in (-1,1)}
        if len(probes)!=len(row['probes']) or probes.keys()!=expected:
            raise ValueError('complete analytic finite probe inventory differs')
        for probe in row['probes']:
            probe_raw=(folder/probe['atoms_file']).read_bytes()
            if sha(probe_raw)!=probe['atoms_sha256']:
                raise ValueError('analytic finite probe digest differs')
            atoms=json.loads(gzip.decompress(probe_raw))
            moves=mobility_check(endpoint,atoms,declaration['case'],plan['protected_bound_angstrom'])
            if any(not math.isclose(v,probe['movement'][k],rel_tol=0.,abs_tol=1e-12) for k,v in moves.items()):
                raise ValueError('reported derivative probe motion differs')
            protected.append(moves['protected_max_motion_angstrom'])
        key=(row['id'],row['index'])
        prior={r['chi']:r for r in original[key]['sensitivity']['chis']}
        chis=[]
        for mapped,analytic in zip(row['mapping'],row['analytic_derivatives_per_degree'],strict=True):
            if not math.isclose(mapped['scale_factor'],180/math.pi,rel_tol=0.,abs_tol=1e-12):
                raise ValueError('analytic torsion scaling differs')
            step_reports=[]
            for h in settings['numeric_steps_degrees']:
                minus=probes[(mapped['vector_index'],h,-1)]['score']
                plus=probes[(mapped['vector_index'],h,1)]['score']
                numeric=(plus-minus)/(2*h)
                external=next(s['slope_per_radian'] for s in prior[mapped['chi']]['steps'] if s['step_degrees']==h)*math.pi/180
                step_reports.append(dict(step_degrees=h,numeric_per_degree=numeric,
                    analytic_comparison=compare(analytic,numeric,settings['comparison_absolute'],settings['comparison_relative']),
                    external_chi_numeric_per_degree=external,
                    numeric_path_comparison=compare(external,numeric,settings['comparison_absolute'],settings['comparison_relative'])))
            chis.append(dict(chi=mapped['chi'],analytic_per_degree=analytic,analytic_per_radian=analytic*180/math.pi,steps=step_reports))
        maximum=max(abs(c['analytic_per_degree']) for c in chis)
        reports.append(dict(id=row['id'],index=row['index'],status='audited',chis=chis,
            maximum_analytic_derivative_per_degree=maximum,above_mode_c_gmax=maximum>settings['cutoff_per_degree'],
            explicit_max200_warning=logs[key]['explicit_max200_warning'],
            finite_resolved_descent=original[key]['sensitivity']['resolved_descent']))
    audited=[r for r in reports if r['status']=='audited']
    comparisons=[s for r in audited for c in r['chis'] for s in c['steps']]
    summary=dict(enrolled=len(reports),audited=len(audited),blockers=len(reports)-len(audited),
        analytic_dofs=sum(len(r['chis']) for r in audited),finite_comparisons=len(comparisons),
        analytic_numeric_mismatches=sum(not s['analytic_comparison']['within_reporting_bound'] for s in comparisons),
        external_vs_minmap_numeric_mismatches=sum(not s['numeric_path_comparison']['within_reporting_bound'] for s in comparisons),
        maximum_analytic_numeric_error_per_degree=max((s['analytic_comparison']['absolute_error_per_degree'] for s in comparisons),default=None),
        endpoints_above_mode_c_gmax=sum(r['above_mode_c_gmax'] for r in audited),
        above_gmax_with_explicit_max200=sum(r['above_mode_c_gmax'] and r['explicit_max200_warning'] for r in audited),
        protected_max_motion_angstrom=max(protected,default=None))
    write(output,dict(schema='dvbfixer.chi-derivatives-audit.v1',summary=summary,endpoints=reports,
        plan_sha256=sha(raw),run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),
        native_coordinates_read=False,native_accuracy_reported=False,profile_selected=False,
        convergence_certified=False,production_changed=False,table_data_included=False))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
