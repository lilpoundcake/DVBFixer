"""Independent saved residue-pair/onebody accounting and exact-state readback.

Never import PyRosetta, generate coordinates, evaluate scores, read native data or
change choices/weights. Interpret model arithmetic only where accounting closes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write


def pair_map(rows,terms):
    result={}
    for row in rows:
        first,second=tuple(row['first']),tuple(row['second'])
        if first==second:
            raise ValueError('self edge in pair accounting')
        key=tuple(sorted((first,second)))
        if key in result or set(row['weighted_terms'])!=set(terms):
            raise ValueError('duplicate pair or term inventory differs')
        if any(not math.isfinite(v) for v in row['weighted_terms'].values()):
            raise ValueError('nonfinite pair energy')
        result[key]=row['weighted_terms']
    return result


def difference(prior,selected,terms,target):
    a,b=pair_map(prior,terms),pair_map(selected,terms)
    reports=[]
    for key in sorted(a.keys()|b.keys()):
        delta={t:b.get(key,{}).get(t,0.)-a.get(key,{}).get(t,0.) for t in terms}
        if any(v!=0. for v in delta.values()):
            reports.append(dict(first=key[0],second=key[1],contains_target=tuple(target) in key,weighted_deltas=delta))
    return reports


def accounting(values,expected,reported):
    total=math.fsum(values)
    bound=1e-5+1e-7*math.fsum(abs(v) for v in values)
    closed=abs(total-expected)<=bound
    if total!=reported['pair_sum'] or expected!=reported['recorded_total'] or bound!=reported['numerical_bound'] or closed!=reported['closed']:
        raise ValueError('independent accounting closure differs')
    return dict(residual=total-expected,bound=bound,closed=closed)


def audit(plan_path,run_root,output):
    plan_raw=plan_path.read_bytes()
    plan=json.loads(plan_raw)
    for p,d in {**plan['inputs'],**plan['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('frozen pair input changed')
    if (run_root/'executed-plan.json').read_bytes()!=plan_raw:
        raise ValueError('executed pair plan differs')
    if plan['schema']=='dvbfixer.mean-pair-accounting-plan.v2-json-transport':
        adapter=json.loads((run_root/'transport-adapter.json').read_bytes())
        if not adapter['hook_restored'] or adapter['numeric_tolerance_added'] or not adapter['complete_snapshot_byte_equality_required']:
            raise ValueError('JSON transport amendment weakened replay')
    run_raw=(run_root/'evidence.json').read_bytes()
    run=json.loads(run_raw)
    if run['plan_sha256']!=sha(plan_raw) or run['native_coordinates_read'] or run['weights_modified']:
        raise ValueError('pair-run contract differs')
    if set(run['results'])!={r['id'] for r in plan['cases']}:
        raise ValueError('disagreement case omitted')
    results={}
    for case in plan['cases']:
        cells={}
        checks={}
        for declaration in case['cells']:
            role=declaration['role']
            index=run['results'][case['id']][role]
            path=Path(index['path'])
            raw=path.read_bytes()
            if sha(raw)!=index['sha256']:
                raise ValueError('pair-state evidence digest differs')
            cell=json.loads(raw)
            old=Path(declaration['original_snapshot']).read_bytes()
            repeated=(path.parent/'pose-snapshot.json.gz').read_bytes()
            if old!=repeated or sha(old)!=cell['original_snapshot_sha256'] or sha(repeated)!=cell['repeated_snapshot_sha256']:
                raise ValueError('complete repeated state not byte-identical')
            if not cell['exact_score_repeat'] or not cell['atom_coordinates_unchanged']:
                raise ValueError('exact repeat/coordinate invariant differs')
            checks[role]={}
            for term in plan['pair_terms']:
                saved=declaration['saved']['score_terms']['ScoreType.'+term]
                expected=saved['raw']*saved['weight']
                checks[role][term]=accounting([p['weighted_terms'][term] for p in cell['residue_pairs']],expected,cell['accounting'][term])
            saved=declaration['saved']['score_terms']['ScoreType.fa_dun']
            checks[role]['fa_dun']=accounting([p['weighted_fa_dun'] for p in cell['fa_dun_onebody']],saved['raw']*saved['weight'],cell['accounting']['fa_dun'])
            cells[role]=cell
        target=case['case']['target_key']
        deltas=difference(cells['prior']['residue_pairs'],cells['selected']['residue_pairs'],plan['pair_terms'],target)
        one_a={tuple(r['residue']):r['weighted_fa_dun'] for r in cells['prior']['fa_dun_onebody']}
        one_b={tuple(r['residue']):r['weighted_fa_dun'] for r in cells['selected']['fa_dun_onebody']}
        if one_a.keys()!=one_b.keys():
            raise ValueError('onebody residue inventory differs')
        one_delta=[dict(residue=k,is_target=k==tuple(target),weighted_delta=one_b[k]-one_a[k]) for k in sorted(one_a) if one_a[k]!=one_b[k]]
        totals={}
        for term in plan['pair_terms']:
            a=next(c for c in case['cells'] if c['role']=='prior')['saved']['score_terms']['ScoreType.'+term]
            b=next(c for c in case['cells'] if c['role']=='selected')['saved']['score_terms']['ScoreType.'+term]
            if a['weight']!=b['weight']:
                raise ValueError('original score weights differ')
            totals[term]=dict(recorded_total_delta=(b['raw']-a['raw'])*a['weight'],
                target_pair_delta=math.fsum(p['weighted_deltas'][term] for p in deltas if p['contains_target']),
                other_pair_delta=math.fsum(p['weighted_deltas'][term] for p in deltas if not p['contains_target']))
        results[case['id']]=dict(target=target,exact_states_repeated=2,accounting=checks,term_totals=totals,
            changed_pairs=deltas,fa_dun_onebody_deltas=one_delta,
            non_target_changed_pairs=sum(not p['contains_target'] for p in deltas),
            complete_accounting_closed=all(t['closed'] for c in checks.values() for t in c.values()))
    write(output,dict(schema='dvbfixer.mean-pair-accounting-audit.v1',cases=results,plan_sha256=sha(plan_raw),
        run_sha256=sha(run_raw),audit_script_sha256=sha(Path(__file__).read_bytes()),native_coordinates_read=False,
        new_candidates=False,weights_modified=False,production_changed=False,table_data_included=False,
        interpretation='recorded model residue-pair attribution, not physical/native causality or a fitted selector'))
    print(json.dumps({n:dict(complete_accounting_closed=r['complete_accounting_closed'],
        non_target_changed_pairs=r['non_target_changed_pairs'],term_totals=r['term_totals']) for n,r in results.items()},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    audit(args.plan,args.run,args.output)
