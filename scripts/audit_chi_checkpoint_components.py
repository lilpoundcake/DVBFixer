"""Saved-only full native component readback with print precision bounds.

Native verbose traces are not an optimizer, native reference, full function trace
or recorded return-branch identifier.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write


def parse_components(text,dimensions):
    if dimensions<1:
        raise ValueError('positive native component dimension required')
    rows=[]
    for line in text.splitlines():
        if 'deriv_dev:' not in line:
            continue
        fields=line.split('deriv_dev:',1)[1].split()
        if len(fields)!=9:
            raise ValueError('native component text format differs')
        i,n=int(fields[0]),int(fields[1])
        values=[float(v) for v in fields[2:]]
        if n!=dimensions or i!=len(rows)%dimensions+1 or not all(math.isfinite(v) for v in values):
            raise ValueError('native component sequence, dimension or finite value differs')
        rows.append(dict(vector_index=i,dimensions=n,score=values[0],dof_type=int(fields[3]),
            atom=int(fields[4]),residue=int(fields[5]),analytic_per_degree=values[4],minimum_absolute_deviation=values[5]))
    if not rows or len(rows)%dimensions:
        raise ValueError('incomplete native component checkpoints')
    return [rows[i:i+dimensions] for i in range(0,len(rows),dimensions)]


def component_cutoff(components,cutoff):
    maximum=max(abs(v) for v in components)
    # Text precision accounting only, not a derivative/scientific acceptance gate.
    uncertainty=1e-11+1e-8*maximum
    return dict(maximum_per_degree=maximum,text_uncertainty=uncertainty,
        certainly_at_or_below_gmax=maximum+uncertainty<=cutoff,
        certainly_above_gmax=maximum-uncertainty>cutoff)


def audit(plan_path,run_root,norm_report_path,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    if 'verbose_checkpoint_amendment' not in plan:
        raise ValueError('separate verbose checkpoint plan required')
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen verbose checkpoint state changed')
    norm_raw=norm_report_path.read_bytes()
    norms=json.loads(norm_raw)
    if norms['plan_sha256']!=sha(raw) or norms['native_coordinates_read'] or norms['audit_script_sha256']!=sha(Path('scripts/audit_chi_checkpoints.py').read_bytes()):
        raise ValueError('independent verbose norm-audit binding differs')
    index_raw=(run_root/'checkpoint-evidence.json').read_bytes()
    if norms['checkpoint_index_sha256']!=sha(index_raw) or (run_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('norm/native executed checkpoint bindings differ')
    adapter_raw=(run_root/'verbose-adapter.json').read_bytes()
    adapter=json.loads(adapter_raw)
    if not adapter['hook_restored'] or any(adapter[k] for k in ('numerical_settings_changed','native_coordinates_read','python_virtual_callback')):
        raise ValueError('verbose native adapter scope differs')
    inventory=[(c,cell) for c in plan['cases'] for cell in c['cells']]
    if len(inventory)!=len(adapter['captures']) or len(inventory)!=len(norms['cases']):
        raise ValueError('complete verbose inventory differs')
    quiet_root=Path('tmp/pdbfixer-research/chi-checkpoints-run-01')
    reports=[]
    total_components=0
    for (case,cell),capture,norm in zip(inventory,adapter['captures'],norms['cases'],strict=True):
        key=(case['id'],cell['declaration']['index'])
        if key!=(capture['id'],capture['index']) or key!=(norm['id'],norm['index']) or norm['status']!='audited' or not norm['endpoint_exact']:
            raise ValueError('exact verbose native replay/identity required')
        if capture['numerical_options_before']!=capture['numerical_options_after'] or set(capture['numerical_options_before'])!=set(plan['verbose_checkpoint_amendment']['numerical_getters']):
            raise ValueError('verbose native numeric option fingerprint differs')
        if capture['actual_native_flags']!=dict(deriv_check=True,deriv_check_verbose=True,deriv_check_to_stdout=True):
            raise ValueError('actual native verbose flags differ')
        folder=run_root/case['id']/cell['declaration']['id']
        log_raw=(folder/'native-components.log').read_bytes()
        if sha(log_raw)!=capture['log_sha256']:
            raise ValueError('native component log digest differs')
        text=log_raw.decode()
        if text.count('DVBFIXER_SYNTHETIC_COMPONENT_CAPTURE_BEGIN')!=1 or text.count('DVBFIXER_SYNTHETIC_COMPONENT_CAPTURE_END')!=1:
            raise ValueError('native component capture controls missing/duplicated')
        trace_raw=(folder/'checkpoints.json').read_bytes()
        if trace_raw!=(quiet_root/case['id']/cell['declaration']['id']/'checkpoints.json').read_bytes():
            raise ValueError('verbose diagnostics changed native checkpoint values')
        trace=json.loads(trace_raw)['checkpoints']
        dimensions=trace[0]['dimensions']
        groups=parse_components(text,dimensions)
        if len(groups)!=len(trace) or len(groups)!=norm['checkpoint_count']:
            raise ValueError('native component checkpoint count differs')
        identities=[(r['dof_type'],r['atom'],r['residue']) for r in groups[0]]
        for group,native in zip(groups,trace,strict=True):
            if [(r['dof_type'],r['atom'],r['residue']) for r in group]!=identities:
                raise ValueError('native component identity/order changed')
            vector=[r['analytic_per_degree'] for r in group]
            parsed_norm=math.sqrt(sum(v*v for v in vector))
            if not math.isclose(parsed_norm,native['analytic_l2_norm_per_degree'],rel_tol=1e-8,abs_tol=2e-12):
                raise ValueError('printed component vector disagrees with full native norm')
            for point in native['sparse_populated_points']:
                parsed=group[point['vector_index']-1]
                if parsed['atom']!=point['dof_id']['atom'] or parsed['residue']!=point['dof_id']['residue'] or not math.isclose(parsed['analytic_per_degree'],point['analytic_per_degree'],rel_tol=1e-8,abs_tol=2e-12):
                    raise ValueError('printed component differs from populated native point')
        total_components+=sum(len(g) for g in groups)
        last=[r['analytic_per_degree'] for r in groups[-1]]
        cutoff=component_cutoff(last,plan['minimizer']['gradient_cutoff'])
        reports.append(dict(id=key[0],index=key[1],endpoint_exact=True,
            native_checkpoint_values_unchanged=True,checkpoint_count=len(groups),
            last_analytic_components_per_degree=last,last_component_cutoff=cutoff,
            endpoint_analytic_maximum_per_degree=norm['endpoint_analytic_maximum_per_degree'],
            endpoint_above_mode_c_gmax=norm['endpoint_above_mode_c_gmax'],
            explicit_max200_warning=norm['explicit_max200_warning'],native_log_sha256=sha(log_raw)))
    summary=dict(starts=len(reports),exact_endpoint_repeats=len(reports),
        native_checkpoint_records=sum(r['checkpoint_count'] for r in reports),
        parsed_component_values=total_components,
        endpoint_above_but_last_components_below=sum(r['endpoint_above_mode_c_gmax'] and r['last_component_cutoff']['certainly_at_or_below_gmax'] for r in reports),
        endpoint_above_and_last_components_above=sum(r['endpoint_above_mode_c_gmax'] and r['last_component_cutoff']['certainly_above_gmax'] for r in reports))
    write(output,dict(schema='dvbfixer.chi-checkpoint-components-audit.v1',summary=summary,cases=reports,
        plan_sha256=sha(raw),norm_audit_sha256=sha(norm_raw),verbose_adapter_sha256=sha(adapter_raw),
        audit_script_sha256=sha(Path(__file__).read_bytes()),native_coordinates_read=False,native_accuracy_reported=False,
        profile_selected=False,full_function_trace=False,return_branch_recorded=False,convergence_certified=False,
        production_changed=False,table_data_included=False,
        interpretation='native derivative component checkpoints and pinned pre-step-gradient check source, not a direct return-branch record'))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--norm-audit',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.norm_audit,args.output)
