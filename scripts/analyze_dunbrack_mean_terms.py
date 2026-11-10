"""Finite saved-score decomposition for ALL prior/score disagreements, no tuning.

No native files, score evaluation, new weights, generation or context restoration.
Energy terms explain the recorded ranking only, not physical/native causality.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write


def deltas(prior,selected):
    if set(prior['score_terms'])!=set(selected['score_terms']):
        raise ValueError('score-term inventory differs')
    rows=[]
    for term,a in prior['score_terms'].items():
        b=selected['score_terms'][term]
        if a['weight']!=b['weight'] or not all(math.isfinite(x) for x in (a['raw'],b['raw'],a['weight'])):
            raise ValueError('nonfinite values or changed score weights')
        rows.append(dict(term=term,weight=a['weight'],weighted_delta=(b['raw']-a['raw'])*a['weight']))
    difference=selected['total_score']-prior['total_score']
    if not math.isclose(sum(r['weighted_delta'] for r in rows),difference,rel_tol=1e-10,abs_tol=1e-8):
        raise ValueError('term delta sum does not match recorded total change')
    return dict(total_delta=difference,terms=sorted(rows,key=lambda r:(r['weighted_delta'],r['term'])))


def analyze(generation_root,score_root,output):
    geometry_raw=(generation_root/'geometry-before-scoring.json').read_bytes()
    geometry=json.loads(geometry_raw)
    score_raw=(score_root/'evidence.json').read_bytes()
    scores=json.loads(score_raw)
    results={}
    for name,row in geometry['results'].items():
        if row['status']!='generated':
            continue
        prior=row['prior_selected_index']
        score=scores['results'][name+'-means']
        chosen=score['selection']['selected_index']
        if prior is None or chosen is None or prior==chosen:
            continue
        by_index={c['index']:c for c in score['cells']}
        results[name]=dict(prior_index=prior,scored_index=chosen,**deltas(by_index[prior],by_index[chosen]))
    write(output,dict(schema='dvbfixer.dunbrack-means-term-deltas.v1',cases=results,
        geometry_sha256=sha(geometry_raw),score_run_sha256=sha(score_raw),script_sha256=sha(Path(__file__).read_bytes()),
        enrollment='ALL nonempty prior-versus-scored mean-pool disagreements, no native-error filter',
        native_coordinates_read=False,new_scores_evaluated=False,weights_modified=False))
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generation',type=Path,required=True)
    parser.add_argument('--score-run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    analyze(args.generation,args.score_run,args.output)
