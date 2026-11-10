"""Private four-arm CH2/CH3 factorial on ALL40 enrolled capture entries.

Original native physics, no phase/energy search, no heavy/kept-H rewriting.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from scripts.carbon_H_capture_adapter import adapt
from scripts.enroll_carbon_H_factorial import ARMS
from scripts.graph_methyl_geometry import initialize as methyl_initialize
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_methylene_geometry import initialize as methylene_initialize


def initialize(points, case, arm):
    if arm not in ARMS:
        raise ValueError('undeclared factorial arm')
    methyl = case['methyl']['groups'] if arm in {'all_methyl_zero', 'combined_zero'} else []
    methylene = case['methylene']['groups'] if arm in {'all_methylene', 'combined_zero'} else []
    first = {i for g in methyl for i in g['hydrogens']}
    second = {i for g in methylene for i in g['hydrogens']}
    if first & second:
        raise ValueError('combined factorial ownership is not disjoint')
    # Both constructions read the original radii/frame; never one another's H.
    methyl_points, grids = methyl_initialize(points, methyl)
    methylene_points = methylene_initialize(points, methylene)
    if grids:
        raise ValueError('fixed-zero arm must not evaluate any phase grid')
    result = [list(p) for p in points]
    for i in first:
        result[i] = methyl_points[i]
    for i in second:
        result[i] = methylene_points[i]
    return result


def plan(admission_path, output):
    import openmm

    raw = admission_path.read_bytes()
    admission = json.loads(raw)
    verify(admission)
    if (admission['summary']['statuses'] != {'admitted': 34, 'refused_source': 6}
            or admission['arms_evaluated'] or admission['openmm_version'] != openmm.__version__):
        raise ValueError('complete pre-execution factorial admission differs')
    sources = dict(admission['source_hashes'])
    for name in ('scripts/research_carbon_H_factorial.py', 'scripts/research_ha_minimizer_path.py',
                 'scripts/enroll_carbon_H_factorial.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.carbon-H-factorial-plan.v1', cases=admission['cases'],
        inputs={**admission['inputs'], str(admission_path): sha(raw)}, source_hashes=sources,
        python=sys.executable, openmm_version=openmm.__version__,
        policy=dict(admission['policy'], maximum_callbacks=512,
            execution='four arms in enrolled context order, each fresh Reference Context, original System masses/forces; no retries or extra seeds',
            controls='all34 unchanged native nm endpoints must reproduce exactly; known historical intersections remain labelled replays'),
        production_changed=False, whole_minimization_executed=False))


def worker(plan_path, index, arm, output):
    import openmm
    from openmm.unit import nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    case = record['cases'][index]
    if case['status'] != 'admitted' or openmm.__version__ != record['openmm_version']:
        raise ValueError('refused factorial capture or runtime must not execute')
    data = adapt(case)
    initial = initialize(data['initial_nm'], case, arm)
    system = data['system']
    xml = openmm.XmlSerializer.serialize(system)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    reporter = reporter_type(openmm)(len(initial), record['policy']['maximum_callbacks'])
    try:
        context.setPositions(initial*nanometer)
        openmm.LocalEnergyMinimizer.minimize(context, 1., 50, reporter)
        final = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
        if openmm.XmlSerializer.serialize(system) != xml:
            raise ValueError('factorial changed original System')
        exact = final == data['final_nm']
        with (output/'path.json.gz').open('xb') as stream:
            stream.write(gzip.compress(json.dumps(dict(initial_nm=initial, final_nm=final, callbacks=reporter.records),
                sort_keys=True, allow_nan=False).encode(), mtime=0))
        write(output/'trace.json', dict(case_id=[case['id'], case['seed'], case['conditioning']],
            arm=arm, plan_sha256=sha(raw), exact_original_nm_endpoint=exact,
            system_sha256=sha(xml.encode()), callbacks=len(reporter.records),
            production_changed=False, whole_minimization_executed=False))
        if arm == 'unchanged' and not exact:
            raise ValueError('factorial unchanged native endpoint differs; retain evidence')
    finally:
        del context
        del integrator


def run(plan_path, output):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(Path('src').resolve()), str(Path('.').resolve()))))
    results, started = [], time.monotonic()
    for index, case in enumerate(record['cases']):
        for arm in ARMS:
            folder = (output/case['id']/f"seed-{case['seed']}"/case['conditioning']/arm).resolve()
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=case['seed'], conditioning=case['conditioning'], arm=arm,
                status=case['status'], returncode=None)
            if case['status'] == 'admitted':
                remaining = record['policy']['total_seconds']-(time.monotonic()-started)
                if remaining <= 0:
                    row['status'] = 'budget_not_executed'
                else:
                    command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                        '--index', str(index), '--arm', arm, '--output', str(folder)]
                    with (folder/'run.log').open('xb') as log:
                        try:
                            process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                            row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                        except subprocess.TimeoutExpired:
                            row['status'] = 'timeout'
            row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(row)
            write(folder/'result.json', row)
            print(json.dumps({k: v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.carbon-H-factorial-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, whole_minimization_executed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--admission', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
            p.add_argument('--arm', choices=ARMS, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.admission, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.arm, args.output)
