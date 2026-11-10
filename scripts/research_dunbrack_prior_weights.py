"""Explicit v2 amendment for raw installed-library interpolation weights.

Original v1 source/failed run stay immutable. API probability() values are retained
verbatim; ranking uses finite nonnegative raw weights, not a claimed normalized
distribution. No native-dependent change, geometry change or production route.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from unittest.mock import patch

from scripts import research_dunbrack_prior as parent


def select_weights(candidates,weights):
    if set(weights)!={1,2,3} or any(not math.isfinite(v) or v<0 for v in weights.values()) or sum(weights.values())<=0:
        raise ValueError('three nonnegative finite family weights with positive total required')
    if [c['index'] for c in candidates]!=[0,1,2] or {c['well'] for c in candidates}!={1,2,3}:
        raise ValueError('exactly one saved candidate per family required')
    feasible=[c for c in candidates if c['common_accepted']]
    return max(feasible,key=lambda c:(weights[c['well']],-c['index']))['index'] if feasible else None


def plan(parent_path,output):
    raw=parent_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if parent.sha(Path(p).read_bytes())!=d:
            raise ValueError('original source-bound plan changed')
    record['schema']='dvbfixer.pyrosetta-prior-raw-weight-plan.v2'
    record['inputs'][str(parent_path)]=parent.sha(raw)
    record['source_hashes']['scripts/research_dunbrack_prior_weights.py']=parent.sha(Path(__file__).read_bytes())
    record['amendment']='v1 stopped before any completed selection ledger because API probability() weights were not normalized; raw inter-bin values may exceed 1. Preserve original values; no calibrated probability or cause claim.'
    record['selection']='maximum finite nonnegative API raw family weight among common-feasible candidates, stable original index tie; no normalized-probability requirement; no coordinate/geometry/angle/score-weight change'
    record['weight_semantics']='samples[].probability is the literal API accessor output, NOT a normalized empirical distribution; separate normalized display weights are diagnostic only'
    parent.write(output,record)


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    if record['schema']!='dvbfixer.pyrosetta-prior-raw-weight-plan.v2':
        raise ValueError('separately frozen raw-weight amendment required')
    original=parent.select_prior
    with patch.object(parent,'select_prior',select_weights):
        parent.run(plan_path,output)
    if parent.select_prior is not original:
        raise ValueError('process-local ranking adapter not restored')
    rows=json.loads((output/'selection-before-native.json').read_bytes())['results']
    view={}
    for n,row in rows.items():
        if row['status']!='queried':
            continue
        total=sum(s['probability'] for s in row['samples'])
        view[n]=dict(raw_total=total,normalized_display_weights={str(s['well']):s['probability']/total for s in row['samples']})
    parent.write(output/'weight-semantics.json',dict(schema='dvbfixer.pyrosetta-prior-weight-semantics.v2',
        cases=view,normalized_weights_used_for_selection=False,selection_adapter_restored=True,
        raw_accessor='DunbrackRotamerSampleData.probability()',native_coordinates_read=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('run')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.parent,args.output)
    else:
        run(args.plan,args.output)
