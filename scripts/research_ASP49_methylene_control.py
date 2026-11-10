"""Prospective local/all-CH2 initial-condition control on six bound 1TIM starts.

Only generated CH2 H coordinates change. Original System, all other starting
coordinates, Reference/1.0/50 physics and seeds remain fixed; no phase selection.
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

from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_methylene_geometry import canonical, discover, initialize

ARMS = ('unchanged', 'LEU48_CH2', 'all_carbon_CH2')
LOCAL = ['A', '48', '', 'LEU', 'CB']


def plan(attribution_plan, attribution_root, attribution_audit, output):
    import openmm

    raw = attribution_plan.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    prior_raw = attribution_audit.read_bytes()
    prior = json.loads(prior_raw)
    run_raw = (attribution_root/'evidence.json').read_bytes()
    if prior['plan_sha256'] != sha(raw) or prior['run_sha256'] != sha(run_raw) or prior['summary']['contexts'] != 6:
        raise ValueError('complete six-context attribution not bound')
    inputs = {**parent['inputs'], str(attribution_plan): sha(raw), str(attribution_audit): sha(prior_raw),
        str(attribution_root/'evidence.json'): sha(run_raw)}
    for result in prior['cases']:
        path = attribution_root/result['raw_artifact']
        if sha(path.read_bytes()) != result['raw_sha256']:
            raise ValueError('retained attribution pairs changed')
        inputs[str(path)] = result['raw_sha256']
    cases = []
    for entry in parent['cases']:
        source = Path(entry['source'])
        ownership = json.loads((source/'H-ownership.json').read_bytes())
        rows, bonds = ownership['rows'], ownership['bonds']
        path = json.loads(gzip.decompress(Path(entry['path']).read_bytes()))
        case = dict(id=entry['id'], seed=entry['seed'], source_arm=entry['arm'], source=entry['source'],
            path=entry['path'], status='admitted')
        try:
            result = discover(rows, bonds, set(ownership['new_h']), set(ownership['mobile']))
            groups = result['groups']
            local = [g for g in groups if g['parent_key'] == LOCAL]
            if len(local) != 1 or {rows[i]['key'][-1] for i in local[0]['hydrogens']} != {'HB2', 'HB3'}:
                raise ValueError('exact generated mobile LEU48 CB H pair unavailable')
            for group in groups:
                canonical(path['initial_nm'], group)
            case.update(groups=groups, local=local, exclusions=result['exclusions'])
        except ValueError as error:
            case.update(status='refused_geometry', refusal=str(error))
        system = openmm.XmlSerializer.deserialize((source/'H-system.xml').read_text())
        if system.getNumConstraints() or system.usesPeriodicBoundaryConditions() or any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
            raise ValueError('undeclared H-System physics')
        cases.append(case)
    sources = dict(parent['source_hashes'])
    for name in ('scripts/research_ASP49_methylene_control.py', 'scripts/research_methylene_geometry.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.ASP49-methylene-control-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, openmm_version=openmm.__version__, python=sys.executable,
        policy=dict(arms=ARMS, per_cell_seconds=1800, total_seconds=7200, maximum_callbacks=512,
            platform='Reference', tolerance_kj_mol_nm=1., maximum_iterations=50,
            timing='prospective control after saved-only LEU48 HB2 attribution; freeze before any CH2 arm evaluation',
            enrollment='all three original 1TIM source arms and existing seeds1/2, three matched CH2 arms each; no new seeds/retries/geometry cleanup',
            ownership='carbon with exactly2 newly mobile H and2 heavy neighbors in canonical same-residue context; measured/kept H never rewritten',
            construction='lexicographic heavy neighbors u,v define bisector b and normal n; H directions a*b +/-sqrt(1-a^2)*n with a=-1/(3*cos(theta/2)); lexicographic first H takes +normal, original C-H radius retained',
            geometry='both H-heavy cosines -1/3; fixed heavy-heavy angle retained and H-H angle implied; refuse zero/collinear/infeasible frames, no arbitrary fallback or label/energy search',
            invariants='only local LEU48 HB2/HB3 or all generated carbon-CH2 H initially move; ASP49 HA/all heavy/other H exact; all original masses/forces unchanged',
            acceptance='every unchanged source-arm nm endpoint must reproduce; all native callbacks and adverse geometry retained',
            limits='H-stage mechanism only, not new proteins or whole-system/production qualification; no force ablation, integration or extra relaxation'),
        production_changed=False))


def worker(plan_path, index, arm, output):
    import openmm
    from openmm.unit import nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if arm not in ARMS or openmm.__version__ != record['openmm_version']:
        raise ValueError('undeclared arm/runtime')
    case = record['cases'][index]
    if case['status'] != 'admitted':
        raise ValueError('refused geometry must not execute')
    path = json.loads(gzip.decompress(Path(case['path']).read_bytes()))
    groups = [] if arm == 'unchanged' else case['local'] if arm == 'LEU48_CH2' else case['groups']
    initial = initialize(path['initial_nm'], groups)
    system = openmm.XmlSerializer.deserialize((Path(case['source'])/'H-system.xml').read_text())
    xml = openmm.XmlSerializer.serialize(system)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    reporter = reporter_type(openmm)(len(initial), record['policy']['maximum_callbacks'])
    try:
        context.setPositions(initial*nanometer)
        openmm.LocalEnergyMinimizer.minimize(context, 1., 50, reporter)
        final = [list(p) for p in context.getState(positions=True).getPositions().value_in_unit(nanometer)]
        if openmm.XmlSerializer.serialize(system) != xml:
            raise ValueError('CH2 construction changed original System')
        exact = final == path['final_nm']
        with (output/'path.json.gz').open('xb') as stream:
            stream.write(gzip.compress(json.dumps(dict(initial_nm=initial, final_nm=final, callbacks=reporter.records),
                sort_keys=True, allow_nan=False).encode(), mtime=0))
        write(output/'trace.json', dict(case=case, arm=arm, plan_sha256=sha(raw), groups=groups,
            exact_original_nm_endpoint=exact, system_sha256=sha(xml.encode()), callbacks=len(reporter.records),
            production_changed=False, whole_minimization_executed=False))
        if arm == 'unchanged' and not exact:
            raise ValueError('unchanged CH2 control endpoint differs; retain blocker')
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
            folder = (output/case['id']/f"seed-{case['seed']}"/case['source_arm']/arm).resolve()
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=case['seed'], source_arm=case['source_arm'], arm=arm, status=case['status'], returncode=None)
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
    write(output/'evidence.json', dict(schema='dvbfixer.ASP49-methylene-control-run.v1', cases=results,
        plan_sha256=sha(raw), production_changed=False, whole_minimization_executed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--attribution-plan', type=Path, required=True)
    p.add_argument('--attribution-root', type=Path, required=True)
    p.add_argument('--attribution-audit', type=Path, required=True)
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
        plan(args.attribution_plan, args.attribution_root, args.attribution_audit, args.output)
    elif args.command == 'run':
        run(args.plan, args.output)
    else:
        worker(args.plan, args.index, args.arm, args.output)
