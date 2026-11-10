"""Private source-diverse three-arm H initializer study; fixed native physics."""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import time
from pathlib import Path

from scripts.graph_methyl_geometry import canonical, initialize
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify

ARMS = ('unchanged', 'graph_tetrahedral_zero', 'graph_tetrahedral_energy_selected')


def plan(capture_plan, capture_root, capture_audit, output):
    raw = capture_plan.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    run_raw = (capture_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    audit_raw = capture_audit.read_bytes()
    audit = json.loads(audit_raw)
    if run['plan_sha256'] != sha(raw) or audit['plan_sha256'] != sha(raw) or audit['run_sha256'] != sha(run_raw):
        raise ValueError('capture/audit binding differs')
    if [(c['id'], c['seed'], c['status']) for c in run['cases']] != [(c['id'], c['seed'], c['status']) for c in audit['cases']]:
        raise ValueError('complete capture/audit enrollment differs')
    inputs = {**parent['inputs'], str(capture_plan): sha(raw), str(capture_root/'evidence.json'): sha(run_raw), str(capture_audit): sha(audit_raw)}
    cases = []
    for cell in run['cases']:
        source = capture_root/cell['id']/f"seed-{cell['seed']}"
        for name, digest in cell['artifacts'].items():
            path = source/name
            if sha(path.read_bytes()) != digest:
                raise ValueError('native capture artifact changed')
            inputs[str(path)] = digest
        case = dict(id=cell['id'], seed=cell['seed'], source=str(source), status=cell['status'])
        if cell['status'] == 'captured':
            ownership = json.loads((source/'H-ownership.json').read_bytes())
            points = json.loads(gzip.decompress((source/'H-path.json.gz').read_bytes()))['initial_nm']
            try:
                for group in ownership['graph']['groups']:
                    for phase in range(0, 360, 15):
                        canonical(points, group, phase)
                case['groups'] = len(ownership['graph']['groups'])
            except ValueError as error:
                case.update(status='refused_geometry', refusal=str(error))
        cases.append(case)
    sources = dict(parent['source_hashes'])
    for name in ('scripts/research_graph_methyl_control.py', 'scripts/research_ha_minimizer_path.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.graph-methyl-control-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, python=parent['python'], openmm_version=parent['openmm_version'],
        policy=dict(arms=ARMS, per_cell_seconds=1800, total_seconds=10800, maximum_callbacks=512,
            platform='Reference', tolerance_kj_mol_nm=1., maximum_iterations=50,
            protocol='immutable eight-source enrollment; all source/capture/geometry refusals retained in all arms; no replacements/retries/new seeds',
            selection='all methyls zero first; single lexicographic parent pass over phases0..345/15, original initial System energy, 1e-8 ties then shortest signed phase; no revisit/final/native score',
            context='separate scoring/minimizing Contexts; original forces/parameters/masses; exact unchanged nm endpoint required',
            limits='H-only research; no whole minimization, accepted candidate/PDB, production repair or general qualification'),
        production_changed=False, whole_minimization_executed=False))


def worker(plan_path, index, arm, output):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if arm not in ARMS or openmm.__version__ != record['openmm_version']:
        raise ValueError('undeclared arm/runtime')
    case = record['cases'][index]
    if case['status'] != 'captured':
        raise ValueError('refused case must not execute')
    source = Path(case['source'])
    path = json.loads(gzip.decompress((source/'H-path.json.gz').read_bytes()))
    ownership = json.loads((source/'H-ownership.json').read_bytes())
    points = path['initial_nm']
    system = openmm.XmlSerializer.deserialize((source/'H-system.xml').read_text())
    xml = openmm.XmlSerializer.serialize(system)
    platform = openmm.Platform.getPlatformByName('Reference')
    selected = arm == 'graph_tetrahedral_energy_selected'
    scoring, scoring_integrator = None, None
    grids = []
    try:
        if selected:
            scoring_integrator = openmm.VerletIntegrator(0.)
            scoring = openmm.Context(system, scoring_integrator, platform)
        def energy(candidate):
            scoring.setPositions(candidate*nanometer)
            return scoring.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        if arm == 'unchanged':
            initial = [p[:] for p in points]
        else:
            initial, grids = initialize(points, ownership['graph']['groups'], selected=selected, energy=energy if selected else None)
    finally:
        if scoring is not None:
            del scoring
        if scoring_integrator is not None:
            del scoring_integrator
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, platform)
    reporter = reporter_type(openmm)(len(points), record['policy']['maximum_callbacks'])
    try:
        context.setPositions(initial*nanometer)
        openmm.LocalEnergyMinimizer.minimize(context, 1., 50, reporter)
        final = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
        exact = final == path['final_nm']
        if openmm.XmlSerializer.serialize(system) != xml:
            raise ValueError('graph initializer changed System')
        data = dict(initial_nm=initial, final_nm=final, callbacks=reporter.records)
        with (output/'path.json.gz').open('xb') as stream:
            stream.write(gzip.compress(json.dumps(data, sort_keys=True, allow_nan=False).encode(), mtime=0))
        write(output/'trace.json', dict(case=case, arm=arm, plan_sha256=sha(raw),
            groups=case['groups'], grids=grids, exact_original_endpoint=exact,
            system_sha256=sha(xml.encode()), callbacks=len(reporter.records), production_changed=False))
        if arm == 'unchanged' and not exact:
            raise ValueError('unchanged original H endpoint differs in nm')
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
    started, results = time.monotonic(), []
    for index, case in enumerate(record['cases']):
        for arm in ARMS:
            folder = (output/case['id']/f"seed-{case['seed']}"/arm).resolve()
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=case['seed'], arm=arm, status=case['status'], returncode=None)
            if case['status'] == 'captured':
                remaining = record['policy']['total_seconds']-(time.monotonic()-started)
                if remaining <= 0:
                    row['status'] = 'budget_not_executed'
                else:
                    command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                        '--index', str(index), '--arm', arm, '--output', str(folder)]
                    with (folder/'run.log').open('xb') as log:
                        try:
                            child = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                            row.update(status='completed' if child.returncode == 0 else 'execution_blocker', returncode=child.returncode)
                        except subprocess.TimeoutExpired:
                            row['status'] = 'timeout'
            row['artifacts'] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(row)
            write(folder/'result.json', row)
            print(json.dumps({k: v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.graph-methyl-control-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, whole_minimization_executed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--capture-plan', type=Path, required=True)
    p.add_argument('--capture-root', type=Path, required=True)
    p.add_argument('--capture-audit', type=Path, required=True)
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
        plan(args.capture_plan, args.capture_root, args.capture_audit, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.arm, args.output)
