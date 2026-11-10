"""Explicit dependency-isolation amendment for the prospective library query.

Original view adapter imported standard-science enrollment/OpenMM at module load.
Keep its source/failed attempt immutable; this executable imports only the existing
stdlib/PyRosetta query driver. No sample/angle/weight/selection changes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write

LIBRARY_ENV=Path('tmp/pdbfixer-research/dunbrack-prior-run-02/environment.json')


def plan(parent_path,error_log,output):
    raw=parent_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('original prospective query state changed')
    record['inputs'][str(parent_path)]=sha(raw)
    record['inputs'][str(error_log)]=sha(error_log.read_bytes())
    record['source_hashes']['scripts/research_prospective_guard_query.py']=sha(Path(__file__).read_bytes())
    record['transport_amendment']='stdlib/PyRosetta-only entry point; avoid importing OpenMM-dependent enrollment module; original failed import retained, no library/coordinate query executed by that attempt; scientific query/weights/selection unchanged'
    write(output,record)


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    if 'transport_amendment' not in record:
        raise ValueError('separate prospective dependency-isolation plan required')
    from scripts.research_dunbrack_prior_weights import run as run_weights

    run_weights(plan_path,output)
    environment=json.loads((output/'environment.json').read_bytes())
    expected=json.loads(LIBRARY_ENV.read_bytes())
    if environment['binary_sha256']!=expected['binary_sha256']:
        raise ValueError('prospective library differs from enrolled installed binary')
    write(output/'prospective-library-binding.json',dict(binary_sha256=environment['binary_sha256'],
        historical_environment_sha256=sha(LIBRARY_ENV.read_bytes()),native_coordinates_read=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--error-log',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('run')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.parent,args.error_log,args.output)
    else:
        run(args.plan,args.output)
