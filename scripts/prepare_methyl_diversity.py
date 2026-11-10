"""One unchanged legacy preparation per admitted graph-methyl source view."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from scripts.branched_template_baseline import preflight
from scripts.prospective_geometry_cohort import prepare_worker
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify

from dvbfixer.ffutils.dat import DatRecord


def worker(plan_path, label, output):
    record = json.loads(plan_path.read_bytes())
    verify(record)
    if preflight() != record['runtime']:
        raise ValueError('canonical preparation runtime/resources changed since admission')
    case = next(c for c in record['cases'] if c['id'] == label)
    if case['status'] != 'admitted':
        raise ValueError('refused source must not enter preparation')
    source = Path(case['projection']).read_bytes()
    with (output/'input.pdb').open('xb') as stream:
        stream.write(source)
    prepare_worker(plan_path, label, output)
    if (output/'input.pdb').read_bytes() != source:
        raise ValueError('legacy preparation changed/deleted its borrowed input')
    DatRecord.load(output/'prepared.dat')
    if not (output/'prepared.pdb').is_file():
        raise ValueError('legacy preparation output missing')


def run(plan_path, output):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if preflight() != record['runtime']:
        raise ValueError('canonical source runtime/resources changed before preparation')
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-source-plan.json').write_bytes(raw)
    contract = dict(source_plan_sha256=sha(raw), orchestrator_sha256=sha(Path(__file__).read_bytes()),
        python=sys.executable, per_case_seconds=1800, total_seconds=5400,
        policy='one legacy prepare per admitted view, Reference PDBFixer/H seed1 through unchanged historical prepare_worker; no retries/source cleanup/replacement')
    write(output/'execution-contract.json', contract)
    env = dict(os.environ, PATH=str(Path(sys.executable).parent)+os.pathsep+os.environ.get('PATH', ''),
        PYTHONPATH=os.pathsep.join((str(Path('src').resolve()), str(Path('.').resolve()))))
    started, results = time.monotonic(), []
    for case in record['cases']:
        folder = (output/case['id']).resolve()
        folder.mkdir(parents=True, exist_ok=False)
        row = dict(id=case['id'], source_status=case['status'], status='refused_source', returncode=None)
        if case['status'] != 'admitted':
            row['refusal'] = case.get('refusal')
        else:
            remaining = contract['total_seconds']-(time.monotonic()-started)
            if remaining <= 0:
                row['status'] = 'budget_not_executed'
            else:
                command = [sys.executable, str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                    '--id', case['id'], '--output', str(folder)]
                with (folder/'run.log').open('xb') as log:
                    try:
                        child = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                            timeout=min(remaining, contract['per_case_seconds']), check=False)
                        row.update(status='prepared' if child.returncode == 0 else 'execution_blocker', returncode=child.returncode)
                    except subprocess.TimeoutExpired:
                        row['status'] = 'timeout'
        row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
        results.append(row)
        write(folder/'result.json', row)
        print(json.dumps({k: v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.methyl-diversity-preparation.v1', cases=results,
        source_plan_sha256=sha(raw), execution_contract_sha256=sha((output/'execution-contract.json').read_bytes()),
        production_changed=False, native_coordinates_read=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for command in ('run', 'worker'):
        p = sub.add_parser(command)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if command == 'worker':
            p.add_argument('--id', required=True)
    args = parser.parse_args()
    if args.command == 'worker':
        worker(args.plan, args.id, args.output)
    else:
        run(args.plan, args.output)
