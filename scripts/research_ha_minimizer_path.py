"""Private ALL24 H-only replay with a read-only native iteration reporter.

Same serialized System/start/Reference/1.0 tolerance/50 iterations. No new H,
seeds, restraints, integration, whole-system minimization or repair output.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.research_dunbrack_prior import sha, write


def verify(record):
    for name, digest in {**record['inputs'], **record['source_hashes']}.items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError(f'H-minimizer path binding changed: {name}')


def reporter_type(openmm):
    class Recorder(openmm.MinimizationReporter):
        def __init__(self, particles, maximum=512):
            super().__init__()
            self.particles = particles
            self.maximum = maximum
            self.records = []

        def report(self, iteration, x, grad, args):
            if len(self.records) >= self.maximum:
                raise ValueError('declared callback storage bound exceeded; retain blocker, no early success')
            if len(x) != 3*self.particles or len(grad) != len(x):
                raise ValueError('native reporter particle/gradient inventory differs')
            points, gradients = list(x), list(grad)
            stats = dict(args)
            if not all(math.isfinite(v) for v in points+gradients+list(stats.values())):
                raise ValueError('native reporter returned nonfinite values')
            self.records.append(dict(ordinal=len(self.records), iteration=iteration,
                positions_nm=points, gradient_kj_mol_nm=gradients, statistics=stats))
            return False  # Never request an early exit or mutate native arrays.
    return Recorder


def plan(parent_path, source_root, output):
    import openmm

    parent_raw = parent_path.read_bytes()
    parent = json.loads(parent_raw)
    verify(parent)
    run_raw = (source_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(parent_raw):
        raise ValueError('original H-placement enrollment differs')
    expected = [(c['id'], c['seed'], c['source_arm']) for c in parent['cases']]
    if [(c['id'], c['seed'], c['source_arm']) for c in run['cases']] != expected or len(expected) != 24:
        raise ValueError('complete H-placement cohort not bound')
    inputs = {**parent['inputs'], str(parent_path.resolve()): sha(parent_raw), str((source_root/'evidence.json').resolve()): sha(run_raw)}
    driver = Path('scripts/research_ha_placement_trace.py').resolve()
    cases = []
    for c in run['cases']:
        if c['status'] != 'completed':
            raise ValueError('original H trace has execution blockers')
        source = source_root/c['id']/f"seed-{c['seed']}"/c['source_arm']
        for name, digest in c['artifacts'].items():
            path = source/name
            if sha(path.read_bytes()) != digest:
                raise ValueError('original H trace artifact changed')
            inputs[str(path.resolve())] = digest
        system = openmm.XmlSerializer.deserialize((source/'h-internal-system.xml').read_text())
        if system.getNumConstraints() or system.usesPeriodicBoundaryConditions() or any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
            raise ValueError('constraints/PBC/virtual sites outside the declared replay')
        cases.append(dict(id=c['id'], seed=c['seed'], source_arm=c['source_arm'], source=str(source.resolve())))
    sources = {**parent['source_hashes']}
    for name in ('scripts/research_ha_minimizer_path.py', 'scripts/research_ha_placement_trace.py',
                 'scripts/research_dunbrack_prior.py', 'scripts/audit_legacy_start_controls.py'):
        sources[str(Path(name).resolve())] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.ha-minimizer-path-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, python=sys.executable, openmm_version=openmm.__version__,
        historical_root=parent['historical_root'], placement_driver=str(driver), original_plan=str(parent_path.resolve()),
        policy=dict(cells=24, tolerance_kj_mol_nm=1., maximum_iterations=50, platform='Reference',
            per_cell_seconds=900, total_seconds=3600, maximum_callbacks_per_cell=512,
            enrollment='faithful native H-stage reproduction for ALL24 original H seeds, with reporter added only at the original call boundary; no alternative H placement/seed/strength/iteration search',
            reporter='read-only copies of native callback coordinates/objective gradients/stats, always return False; preserve ordinal and native iteration separately',
            invariants='zero constraints/PBC/virtual sites; identical original serialized System, fixed-existing coordinates and final H endpoint required',
            retention='all callbacks privately, failures/timeouts/budget nonexecutions separate, no retry/adoption',
            audit='initial/every callback/final geometry; static force group readback along fixed heavy outward direction and its HA-sphere tangent; callback gradient vs physical HA force',
            limits='callbacks are optimizer iteration observations, not all line-search evaluations, physical dynamics, integrated work, convergence or a production repair'),
        production_changed=False, native_coordinates_read=False))


def worker(plan_path, index, output):
    from unittest.mock import patch

    import openmm
    from openmm.unit import angstrom, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if openmm.__version__ != record['openmm_version'] or Path.cwd() != Path(record['historical_root']):
        raise ValueError('OpenMM version or historical cwd differs')
    case = record['cases'][index]
    source = Path(case['source'])
    initial = load_atoms(source/'h-internal-before.json.gz')
    spec = importlib.util.spec_from_file_location('frozen_placement_driver', record['placement_driver'])
    if spec is None or spec.loader is None:
        raise ValueError('cannot load frozen placement driver')
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    native = openmm.LocalEnergyMinimizer.minimize
    recorder = reporter_type(openmm)(len(initial), record['policy']['maximum_callbacks_per_cell'])
    captures = []

    def observe(context, *args, **kwargs):
        if captures or args != (1., 50) or kwargs:
            raise ValueError('unexpected native H-minimizer call')
        system = context.getSystem()
        if system.getNumParticles() != len(initial) or system.getNumConstraints():
            raise ValueError('particle or zero-constraint contract differs')
        serialized = openmm.XmlSerializer.serialize(system)
        first_state = context.getState(positions=True)
        first_nm = [list(p) for p in first_state.getPositions().value_in_unit(nanometer)]
        before = [dict(r, xyz=list(p)) for r, p in zip(initial, first_state.getPositions().value_in_unit(angstrom), strict=True)]
        result = native(context, 1., 50, recorder)
        last_state = context.getState(positions=True)
        final_nm = [list(p) for p in last_state.getPositions().value_in_unit(nanometer)]
        final = [dict(r, xyz=list(p)) for r, p in zip(initial, last_state.getPositions().value_in_unit(angstrom), strict=True)]
        if openmm.XmlSerializer.serialize(system) != serialized:
            raise ValueError('native iteration recording changed System')
        captures.append(dict(initial_atoms=before, final_atoms=final,
            initial_nm=first_nm, final_nm=final_nm, callbacks=recorder.records))
        return result

    with patch.object(openmm.LocalEnergyMinimizer, 'minimize', observe):
        driver.worker(Path(record['original_plan']), index, output)
    if len(captures) != 1:
        raise ValueError('native H minimizer was not observed exactly once')
    payload = json.dumps(captures[0], sort_keys=True, allow_nan=False).encode()
    with (output/'path.json.gz').open('xb') as stream:
        stream.write(gzip.compress(payload, mtime=0))
    delta = coordinate_delta(load_atoms(source/'h-internal-after.json.gz'), captures[0]['final_atoms'])
    start = coordinate_delta(initial, captures[0]['initial_atoms'])
    write(output/'path-trace.json', dict(case=case, plan_sha256=sha(raw),
        callbacks=len(recorder.records), saved_endpoint_difference=delta,
        saved_start_difference=start, reporter_early_stop=False,
        integration=False, production_changed=False, native_coordinates_read=False))
    if any(not d[kind]['exact'] for d in (delta, start) for kind in ('heavy', 'hydrogen')):
        raise ValueError('reported replay differs from original start/endpoint; retain artifacts, no causal path claim')


def run(plan_path, output):
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
    write(output/'evidence.json', dict(schema='dvbfixer.ha-minimizer-path-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, native_coordinates_read=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--parent', type=Path, required=True)
    p.add_argument('--source-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.parent, args.source_root, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.output)
