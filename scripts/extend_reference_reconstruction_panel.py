#!/usr/bin/env python3
"""Prospective bounded completion of the two reference-case development panels.

Thirty serial calls: five remaining approved seeds x two existing cases x three
backends. Original seed-7 results and every raw failure remain unchanged. This is
not the reserved validation cohort, an anchor sweep or a capability decision.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import traceback
from dataclasses import replace
from pathlib import Path

from scripts.audit_mac_reconstruction_smoke import audit
from scripts.remap_reference_model_output import run as remap

from dvbfixer.model.diffusion.contract import DiffusionRequest

SEEDS = (42, 101, 102, 103, 104)
PROJECT = Path(__file__).parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def declare(root, output):
    original_path = root / 'reference-pilot/plan.json'
    original = json.loads(original_path.read_text())
    output.mkdir(exist_ok=False)
    jobs = []
    for seed in SEEDS:
        for old in original['jobs']:
            case = Path(old['case'])
            work = output / f'{case.name}-{old["arm"]}-seed{seed}'
            work.mkdir()
            command = list(old['command'])
            command[0] = sys.executable
            command = [part.replace(old['work'], str(work)) for part in command]
            boundary = command.index('--')
            args = command[boundary+1:]
            if old['arm'] == 'rfd3':
                args = [f'seed={seed}' if part == 'seed=7' else part for part in args]
            elif old['arm'] == 'boltz':
                args[args.index('--seed')+1] = str(seed)
            else:
                request = DiffusionRequest.from_json((case / 'evaluation-mask.json').read_text())
                (work / 'request.json').write_text(replace(request, seeds=(seed,)).to_json())
                (work / request.normalized_pdb.path).write_bytes((case / 'source.pdb').read_bytes())
                args[args.index('--request')+1] = str(work / 'request.json')
                args += ['--seed', str(seed)]
            jobs.append({**old, 'work': str(work), 'seed': seed,
                         'command': command[:boundary+1]+args})
    files = [Path(__file__), PROJECT / 'scripts/audit_mac_reconstruction_smoke.py',
             PROJECT / 'scripts/remap_reference_model_output.py', PROJECT / 'scripts/run_mac_reconstruction_guarded.py',
             PROJECT / 'scripts/rfd3_mps_seeded_research.py', PROJECT / 'scripts/trace_rfd3_seed_transport.py',
             PROJECT / 'deploy/protenix-v1/mps_research_smoke.py', PROJECT / 'deploy/boltz-2/atom_mapping.py']
    for case in sorted({Path(j['case']) for j in jobs}):
        files.extend(p for p in case.rglob('*') if p.is_file())
    files.extend(Path(j['work']) / p for j in jobs if j['arm'] == 'protenix' for p in ('request.json', 'source.pdb'))
    engines = PROJECT / '.artifacts/mac-reconstruction-2026-11-11/engines'
    for engine in ('boltz', 'rfd3', 'protenix-v1'):
        diff = subprocess.check_output(['git', 'diff'], cwd=engines / engine)
        (output / f'{engine}-source.diff').write_bytes(diff)
    plan = {'source_pilot_plan_sha256': digest(original_path), 'seeds': list(SEEDS),
        'seed7': 'reuse original retained pilots; not resampled', 'reserved_seeds': '211/223/227/229/233 unused',
        'calls': len(jobs), 'serial_ml': True, 'automatic_retries': False,
        'same_source_masks_and_per_backend_sampling_settings': True,
        'boltz_source_amendment': 'owned fork now has optional post-update callback patch; stock CLI invokes no callback; separate synthetic disabled/identity control passed; original pilot remains unchanged',
        'rfd3_seed_contract': 'private-mps-seed-at-model-v1',
        'qualification': False, 'anchor_sweep': False, 'reference_used_only_for_scoring': True,
        'files': {str(p): digest(p) for p in files}, 'jobs': jobs}
    (output / 'plan.json').write_text(json.dumps(plan, indent=2)+'\n')
    print(f'Declared {len(jobs)} serial reference-panel extension calls')


def run(root, output):
    plan = json.loads((output / 'plan.json').read_text())
    for name, expected in plan['files'].items():
        if digest(Path(name)) != expected:
            raise ValueError(f'declared reference-panel input changed: {name}')
    rows = []
    with (output / 'orchestration.log').open('xb') as log:
        for job in plan['jobs']:
            work, case = Path(job['work']), Path(job['case'])
            paths = [PROJECT, PROJECT / 'src', PROJECT / 'scripts']
            env = os.environ.copy()
            if job['arm'] == 'protenix':
                paths += [PROJECT / 'deploy/protenix-v1', PROJECT / '.artifacts/mac-reconstruction-2026-11-11/engines/protenix-v1']
            env.update(PYTHONPATH=os.pathsep.join(str(p) for p in paths),
                PROTENIX_ROOT_DIR=str(PROJECT / '.artifacts/mac-reconstruction-2026-11-11/assets/protenix'), LAYERNORM_TYPE='torch')
            if job['arm'] == 'rfd3':
                env.update(DVBFIXER_RFD3_MODEL_SEED=str(job['seed']),
                    DVBFIXER_RFD3_TRANSPORT_OUTPUT=str(work / 'run/transport'),
                    DVBFIXER_RFD3_EXPECTED_FIXED=str(job['expected_fixed']),
                    DVBFIXER_RFD3_EXPECTED_GENERATED=str(job['expected_generated']))
            code = subprocess.run(job['command'], env=env, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT).returncode
            row = {'case': case.name, 'arm': job['arm'], 'seed': job['seed'], 'returncode': code}
            if code == 0:
                try:
                    pattern = {'rfd3': '*.cif.gz', 'boltz': '*_model_0.pdb', 'protenix': 'candidate.pdb'}[job['arm']]
                    predictions = list((work / 'run/predictions').rglob(pattern))
                    if len(predictions) != 1:
                        raise ValueError('native output missing/ambiguous')
                    prediction = predictions[0]
                    if job['arm'] != 'protenix':
                        prediction = remap(case, prediction, work / 'native-mapping', model=job['arm'])
                    row['audit'] = audit(case, prediction, work / 'audit', seed=job['seed'])
                except Exception:
                    row['audit_failure'] = traceback.format_exc()
            execution_path = work / 'run/execution.json'
            execution = json.loads(execution_path.read_text()) if execution_path.exists() else {}
            row['execution'] = execution
            rows.append(row)
            (output / 'results.json').write_text(json.dumps(rows, indent=2)+'\n')
            print(json.dumps({k: row[k] for k in ('case', 'arm', 'seed', 'returncode')}), flush=True)
            if not execution or execution.get('guard_stop_reason'):
                raise RuntimeError('reference extension stopped on guard/orchestration failure; no implicit retry')
    (output / 'completed.json').write_text(json.dumps({'executed': len(rows), 'qualified': False})+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('declare', 'run'))
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    (declare if args.mode == 'declare' else run)(args.root.absolute(), args.output.absolute())
