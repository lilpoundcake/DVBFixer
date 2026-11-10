"""Saved-only independent diagnostic readback and exact logging-replay audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write


def audit(plan_path,run_root,output):
    raw=plan_path.read_bytes()
    plan=json.loads(raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen termination input changed')
    if (run_root/'executed-plan.json').read_bytes()!=raw:
        raise ValueError('logging replay plan differs')
    parent=Path(plan['termination_replay']['parent_run'])
    parent_plan=json.loads((parent/'executed-plan.json').read_bytes())
    if plan['minimizer']!=parent_plan['minimizer'] or plan['init_flags']!=parent_plan['init_flags']+' -unmute core.optimization.Minimizer':
        raise ValueError('nonlogging settings changed')
    evidence_raw=(run_root/'termination-evidence.json').read_bytes()
    evidence=json.loads(evidence_raw)
    if 'native_flush_amendment' not in plan:
        raise ValueError('original fd capture lacks native-buffer boundary controls; no endpoint attribution')
    adapter=json.loads((run_root/'native-flush-adapter.json').read_bytes())
    if not adapter['hook_restored'] or adapter['numerical_settings_changed'] or adapter['synthetic_optimizer_warnings']:
        raise ValueError('native flush amendment contract differs')
    if evidence['plan_sha256']!=sha(raw) or not evidence['hook_restored'] or evidence['native_coordinates_read'] or evidence['weights_modified']:
        raise ValueError('termination scope differs')
    inventory=[(c,cell) for c in plan['cases'] for cell in c['cells']]
    if len(inventory)!=len(evidence['traces']) or len(inventory)!=len(evidence['comparisons']):
        raise ValueError('incomplete logging replay inventory')
    if len(adapter['captures'])!=len(inventory):
        raise ValueError('incomplete controlled capture inventory')
    rows=[]
    for (case,cell),trace,comparison,capture in zip(inventory,evidence['traces'],evidence['comparisons'],adapter['captures'],strict=True):
        key=(case['id'],cell['declaration']['index'])
        if key!=(trace['id'],trace['index']) or key!=(comparison['id'],comparison['index']):
            raise ValueError('logging endpoint identity differs')
        old=parent/case['id']/cell['declaration']['id']
        new=run_root/case['id']/cell['declaration']['id']
        for filename in ('before-atoms.json','after-atoms.json','result.json','pose-snapshot.json.gz'):
            if (old/filename).read_bytes()!=(new/filename).read_bytes() or not comparison['exact'][filename]:
                raise ValueError('logging replay altered exact endpoint state/result')
        log=(new/'minimizer.log').read_bytes()
        if sha(log)!=trace['sha256']:
            raise ValueError('native diagnostic log digest differs')
        if sha(log)!=capture['sha256'] or capture['path']!=trace['path'] or not capture['synthetic_markers_verified']:
            raise ValueError('native boundary control binding differs')
        text=log.decode()
        if text.count('DVBFIXER_SYNTHETIC_CAPTURE_CONTROL_BEGIN')!=1 or text.count('DVBFIXER_SYNTHETIC_CAPTURE_CONTROL_END')!=1:
            raise ValueError('native tracer markers missing or duplicated')
        maximum='LBFGS MAX CYCLES 200 EXCEEDED, BUT FUNC NOT CONVERGED!' in text
        failed='Line search failed even after resetting Hessian; aborting at iter#' in text
        if maximum!=trace['explicit_max200_warning'] or failed!=trace['explicit_line_search_abort'] or trace['termination_certified']:
            raise ValueError('independent diagnostic interpretation differs')
        if any(trace['actual_settings'][k]!=v for k,v in plan['minimizer'].items()):
            raise ValueError('logging replay numerical settings differ')
        rows.append(dict(id=key[0],index=key[1],endpoint_exact=True,log_sha256=sha(log),log_bytes=len(log),
                         explicit_max200_warning=maximum,explicit_line_search_abort=failed))
    # Saved-only association across ALL endpoints, not an intervention/selection.
    finite_path=Path('docs/research/chi-stationarity-evidence.json')
    finite_raw=finite_path.read_bytes()
    finite=json.loads(finite_raw)
    if finite['native_coordinates_read'] or finite['native_accuracy_reported']:
        raise ValueError('finite probe record is not native-free')
    lookup={(r['id'],r['index']):r for r in finite['endpoints']}
    if set(lookup)!={(r['id'],r['index']) for r in rows}:
        raise ValueError('finite-probe association inventory differs')
    for row in rows:
        probe=lookup[(row['id'],row['index'])]
        row['finite_resolved_descent']=probe['sensitivity']['resolved_descent']
    summary=dict(endpoints=len(rows),exact_endpoint_repeats=sum(r['endpoint_exact'] for r in rows),
        explicit_max200_warnings=sum(r['explicit_max200_warning'] for r in rows),
        explicit_line_search_aborts=sum(r['explicit_line_search_abort'] for r in rows),
        no_explicit_termination_message=sum(not r['explicit_max200_warning'] and not r['explicit_line_search_abort'] for r in rows),
        finite_descent_and_max200=sum(r['finite_resolved_descent'] and r['explicit_max200_warning'] for r in rows),
        finite_descent_and_line_search_abort=sum(r['finite_resolved_descent'] and r['explicit_line_search_abort'] for r in rows))
    write(output,dict(schema='dvbfixer.chi-termination-audit.v1',summary=summary,endpoints=rows,
        plan_sha256=sha(raw),run_evidence_sha256=sha(evidence_raw),finite_probe_evidence_sha256=sha(finite_raw),
        native_flush_adapter_sha256=sha((run_root/'native-flush-adapter.json').read_bytes()),
        audit_script_sha256=sha(Path(__file__).read_bytes()),native_coordinates_read=False,native_accuracy_reported=False,
        numerical_settings_changed=False,profile_selected=False,convergence_certified=False,production_changed=False,
        interpretation='explicit diagnostics and saved-only finite-descent association, not causal attribution or convergence proof'))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
