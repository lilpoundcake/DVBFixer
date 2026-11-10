"""Private faithful replay of the original H-placement boundary, ALL24 cells.

Observe before/after OpenMM's existing 50-iteration H-only minimization. Abort
before whole-system phase1. No alternative H placement, parameter search or fix.
Workers run against the preserved historical source, never modifying that tree.
"""
from __future__ import annotations

import argparse
import functools
import gzip
import inspect
import json
import os
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from scripts.research_dunbrack_prior import sha, write


class BoundaryReached(Exception):
    """Expected stop before any full-system minimization."""


def verify(record):
    for name, digest in {**record['inputs'], **record['source_hashes']}.items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError(f'placement trace binding changed: {name}')


def plan(force_path, output):
    import openmm
    import openmm.app.modeller as module
    from scripts.audit_stereo_force_projection import historical_verify

    force = json.loads(force_path.read_bytes())
    verify(force)
    baseline = Path(force['historical_root'])
    parent = baseline/force['parent_plan']
    historical_verify(json.loads(parent.read_bytes()), baseline)
    inputs = {str(force_path.resolve()): sha(force_path.read_bytes()), str(parent): sha(parent.read_bytes())}
    for p in (baseline/'src/dvbfixer').rglob('*.py'):
        inputs[str(p)] = sha(p.read_bytes())
    for p in (baseline/'scripts').glob('*.py'):
        inputs[str(p)] = sha(p.read_bytes())
    modeller_path = Path(inspect.getsourcefile(module.Modeller))
    inputs[str(modeller_path)] = sha(modeller_path.read_bytes())
    cases = []
    for c in force['executions']:
        source = Path(force['run_root'])/c['id']/f"seed-{c['seed']}"/c['arm']
        for name in ('prepared.pdb', 'prepared.dat', 'input-change.json', 'h-1-after.json.gz',
                     'phase-1-before.json.gz', 'phase-1-before-system.xml'):
            path = (source/name).resolve()
            inputs[str(path)] = sha(path.read_bytes())
        cases.append(dict(id=c['id'], seed=c['seed'], source_arm=c['arm'], source=str(source.resolve())))
    if len(cases) != 24:
        raise ValueError('complete saved cohort not enrolled')
    write(output, dict(schema='dvbfixer.ha-placement-trace-plan.v1', cases=cases, inputs=inputs,
        source_hashes={str(Path(__file__).resolve()): sha(Path(__file__).read_bytes())},
        historical_root=str(baseline), original_plan=str(parent), python=sys.executable,
        openmm_version=openmm.__version__,
        policy=dict(cells=24, per_cell_seconds=900, total_seconds=3600,
            enrollment='ALL24 original starts and original H seeds, no new seeds or alternative placement',
            observation='capture existing H-only Context immediately before/after original 1.0-tolerance/50-iteration minimize; forward exact original call',
            stop='expected boundary stop before full-system phase1; retain exceptions/nonexecutions distinctly',
            invariants='original H endpoint, first whole-system start coordinates and System must reproduce saved originals exactly; existing atoms immobile in H System',
            interpretation='provenance of known HA only; not a repair, native-error score, force isolation, original THR2 replay or general sampling qualification'),
        production_changed=False, native_coordinates_read=False))


def worker(plan_path, index, folder):
    import openmm
    from openmm.app import Modeller, Simulation
    from openmm.unit import angstrom
    from scripts.research_legacy_start_controls import snapshot_atoms
    from scripts.research_start_geometry import worker as historical_worker

    record = json.loads(plan_path.read_bytes())
    verify(record)
    if Path.cwd() != Path(record['historical_root']) or openmm.__version__ != record['openmm_version']:
        raise ValueError('historical cwd or OpenMM version differs')
    case = record['cases'][index]
    source = Path(case['source'])
    for name in ('prepared.pdb', 'prepared.dat', 'input-change.json'):
        with (folder/name).open('xb') as stream:
            stream.write((source/name).read_bytes())
    native_h = Modeller.addHydrogens
    native_min = openmm.LocalEnergyMinimizer.minimize
    captures = []
    h_calls = []
    active = False

    def observe(context, *args, **kwargs):
        if not active or args != (1.0, 50) or kwargs:
            raise ValueError('unexpected H-placement minimizer call')
        before = context.getState(positions=True).getPositions().value_in_unit(angstrom)
        system_xml = openmm.XmlSerializer.serialize(context.getSystem())
        result = native_min(context, *args, **kwargs)
        after = context.getState(positions=True).getPositions().value_in_unit(angstrom)
        if openmm.XmlSerializer.serialize(context.getSystem()) != system_xml:
            raise ValueError('H-stage observation changed System')
        captures.append(dict(before=[list(p) for p in before], after=[list(p) for p in after],
            system_xml=system_xml, platform=context.getPlatform().getName()))
        return result

    @functools.wraps(native_h)
    def hydrogen(model, *args, **kwargs):
        nonlocal active
        if h_calls or active:
            raise ValueError('unexpected second/reentrant H placement')
        original = snapshot_atoms(model.topology, model.positions)
        active = True
        try:
            with patch.object(openmm.LocalEnergyMinimizer, 'minimize', observe):
                result = native_h(model, *args, **kwargs)
        finally:
            active = False
        final = snapshot_atoms(model.topology, model.positions)
        if len(captures) != 1 or len(final) != len(captures[0]['after']):
            raise ValueError('H-stage capture inventory differs')
        event = captures[0]
        for stage in ('before', 'after'):
            material = [dict(row, xyz=point) for row, point in zip(final, event[stage], strict=True)]
            payload = json.dumps(dict(atoms=material), sort_keys=True, allow_nan=False).encode()
            with (folder/f'h-internal-{stage}.json.gz').open('xb') as stream:
                stream.write(gzip.compress(payload, mtime=0))
        with (folder/'h-internal-system.xml').open('x') as stream:
            stream.write(event['system_xml'])
        payload = json.dumps(dict(atoms=original), sort_keys=True, allow_nan=False).encode()
        with (folder/'h-input.json.gz').open('xb') as stream:
            stream.write(gzip.compress(payload, mtime=0))
        h_calls.append(dict(platform=event['platform'], existing_atoms=len(original), final_atoms=len(final),
            original_call_forwarded=True, tolerance=1., maximum_iterations=50))
        return result

    def stop_before_phase(*args, **kwargs):
        raise BoundaryReached()

    stopped = False
    with patch.object(Modeller, 'addHydrogens', hydrogen), patch.object(Simulation, 'minimizeEnergy', stop_before_phase):
        try:
            historical_worker(Path(record['original_plan']), case['id'], case['seed'], case['source_arm'], folder)
        except BoundaryReached:
            stopped = True
    if not stopped or len(h_calls) != 1 or Modeller.addHydrogens is not native_h or openmm.LocalEnergyMinimizer.minimize is not native_min:
        raise ValueError('H trace boundary or hook restoration differs')
    write(folder/'placement-trace.json', dict(case=case, h_calls=h_calls, phase1_not_executed=True,
        plan_sha256=sha(plan_path.read_bytes()), production_changed=False, native_coordinates_read=False))


def run(plan_path, output):
    import time

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    baseline = Path(record['historical_root'])
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(baseline/'src'), str(baseline))))
    started = time.monotonic()
    results = []
    for i, case in enumerate(record['cases']):
        folder = (output/case['id']/f"seed-{case['seed']}"/case['source_arm']).resolve()
        folder.mkdir(parents=True, exist_ok=False)
        row = {k: case[k] for k in ('id', 'seed', 'source_arm')}
        remaining = record['policy']['total_seconds']-(time.monotonic()-started)
        if remaining <= 0:
            row.update(status='budget_not_executed', returncode=None)
        else:
            command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                '--index', str(i), '--output', str(folder)]
            with (folder/'run.log').open('xb') as log:
                try:
                    process = subprocess.run(command, cwd=baseline, env=env, stdout=log, stderr=subprocess.STDOUT,
                        timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                    row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                except subprocess.TimeoutExpired:
                    row.update(status='timeout', returncode=None)
        row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
        results.append(row)
        write(folder/'result.json', row)
        print(json.dumps({k: v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.ha-placement-trace-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, native_coordinates_read=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--force-plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.force_plan, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.output)
