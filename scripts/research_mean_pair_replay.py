"""Explicit JSON-transport amendment for exact pair-accounting replay.

Original run compared live tuples in import_audit against JSON lists and stopped
on the documented terminal OXT inventory. Normalize JSON types only; all numbers,
score outputs and compressed complete-pose snapshot bytes must still match exactly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

from scripts import research_mean_pair_energies as parent
from scripts.research_dunbrack_prior import sha, write


def transport(value):
    return json.loads(json.dumps(value,allow_nan=False))


def plan(parent_path,output):
    raw=parent_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('original frozen state changed')
    record['schema']='dvbfixer.mean-pair-accounting-plan.v2-json-transport'
    record['inputs'][str(parent_path)]=sha(raw)
    for name in ('tmp/pdbfixer-research/mean-pair-energies-run-01/execution-error.json',
                 'tmp/pdbfixer-research/mean-pair-energies-run-02/execution-error.json',
                 'tmp/pdbfixer-research/mean-pair-instrumented-diagnostic-01.log'):
        record['inputs'][name]=sha(Path(name).read_bytes())
    record['source_hashes']['scripts/research_mean_pair_replay.py']=sha(Path(__file__).read_bytes())
    record['amendment']='canonical JSON tuple/list transport for repeated score metadata only; no numeric rounding/tolerance, no snapshot/score/coordinate/weight change; original failed runs retained'
    write(output,record)


def run(plan_path,output):
    from scripts import research_rosetta_fixed_scoring as scoring

    record=json.loads(plan_path.read_bytes())
    if record['schema']!='dvbfixer.mean-pair-accounting-plan.v2-json-transport':
        raise ValueError('separate JSON-transport plan required')
    original=scoring.score_cell

    def normalized(*args,**kwargs):
        return transport(original(*args,**kwargs))

    with patch.object(scoring,'score_cell',normalized):
        parent.run(plan_path,output)
    if scoring.score_cell is not original:
        raise ValueError('research JSON transport hook not restored')
    write(output/'transport-adapter.json',dict(schema='dvbfixer.mean-pair-json-transport.v1',
        hook_restored=True,numeric_tolerance_added=False,complete_snapshot_byte_equality_required=True))


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
