"""Private all-context whole continuation after the fixed carbon-H factorial.

Independent native baselines precede every arm. Original no-solvent Systems,
States, integrators, anchors and two native phases only, before emergency repair.
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

from scripts.carbon_H_capture_adapter import adapt, labels, read_json
from scripts.preflight_graph_methyl_transport import mapped_start
from scripts.research_beta_wall import restraint
from scripts.research_dunbrack_prior import sha, write
from scripts.research_h_stage_transport import topology_from_rows
from scripts.research_ha_minimizer_path import verify


def specification(case, legacy):
    source = Path(case['system']).parent
    if case['context_kind'] == 'geometry-reconstruction-stress':
        prior = next(c for c in legacy['cases'] if (c['id'], c['seed'], c['source_arm'])
            == (case['id'], case['seed'], case['conditioning']))
        original = Path(prior['original_whole'])
        rows = labels(read_json(original/'phase-1-before.json.gz')['atoms'])
        spec = dict(rows=rows, bonds=read_json(source/'topology.json')['bonds'], mapping=list(range(len(rows))),
            system=str(original/'phase-1-before-system.xml'), state=str(original/'initial-state.xml'),
            integrator=str(original/'initial-integrator.xml'), dat=str(source/'prepared.dat'),
            historical_angstrom=dict(initial=str(original/'phase-1-before.json.gz'),
                phase1=str(original/'phase-1-after.json.gz'), phase2=str(original/'phase-2-after.json.gz')),
            target=prior['target'], residue=prior['residue'], observed_source=prior['observed_source'])
    else:
        whole = read_json(source/'whole-topology.json')
        spec = dict(rows=whole['atoms'], bonds=whole['bonds'], mapping=whole['H_to_whole_indices'],
            system=str(source/'whole-system.xml'), state=str(source/'whole-state.xml'),
            integrator=str(source/'whole-integrator.xml'), dat=str(source/'prepared.dat'))
    return spec


def preflight(case, spec, data):
    import openmm
    from openmm.unit import nanometer

    mapping = spec['mapping']
    # The V2 map is object-derived. Legacy is accepted only under exact labels.
    if case['context_kind'] == 'geometry-reconstruction-stress' and spec['rows'] != data['rows']:
        raise ValueError('legacy full labels differ; no guessed whole aliases')
    points, count = mapped_start(data['rows'], data['bonds'], data['final_nm'], mapping)
    graph = {tuple(sorted((mapping[a], mapping[b]))) for a, b in data['bonds']}
    if graph != {tuple(sorted(b)) for b in spec['bonds']}:
        raise ValueError('H-to-whole native bond map differs')
    system = openmm.XmlSerializer.deserialize(Path(spec['system']).read_text())
    if system.getNumParticles() != len(spec['rows']) or system.usesPeriodicBoundaryConditions():
        raise ValueError('whole particle/no-solvent contract differs')
    state = openmm.XmlSerializer.deserialize(Path(spec['state']).read_text())
    if points != [list(p) for p in state.getPositions().value_in_unit(nanometer)]:
        raise ValueError('unchanged post-H whole start differs from original native nm State')
    restraint(system)
    topology_from_rows(spec['rows'], spec['bonds'])
    return dict(exact_original_nm_start=True, native_bond_mapping_exact=True, post_H_repairs=count)


def plan(h_plan, h_root, h_audit, legacy_plan, output):
    import openmm

    raw = h_plan.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    audit_raw = h_audit.read_bytes()
    audit = json.loads(audit_raw)
    run_raw = (h_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if (audit['plan_sha256'] != sha(raw) or audit['run_sha256'] != sha(run_raw)
            or audit['summary']['statuses'] != {'completed':136, 'refused_source':24}
            or audit['summary']['complete_matched_contexts'] != 34
            or audit['summary']['maximum_mobile_gradient_error_kj_mol_nm'] > 1e-4):
        raise ValueError('complete matched factorial H audit required')
    legacy_raw = legacy_plan.read_bytes()
    legacy = json.loads(legacy_raw)
    verify(legacy)
    inputs = {**parent['inputs'], **legacy['inputs'], str(h_plan):sha(raw), str(h_audit):sha(audit_raw),
        str(h_root/'evidence.json'):sha(run_raw), str(legacy_plan):sha(legacy_raw)}
    for cell in run['cases']:
        folder = h_root/cell['id']/f"seed-{cell['seed']}"/cell['conditioning']/cell['arm']
        for name, digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('factorial H artifacts changed')
            inputs[str(folder/name)] = digest
    cases = []
    for case in parent['cases']:
        row = dict(case)
        if case['status'] != 'admitted':
            cases.append(row)
            continue
        data = adapt(case)
        spec = specification(case, legacy)
        row.update(whole=spec, whole_preflight=preflight(case, spec, data))
        for name in ('system', 'state', 'integrator', 'dat', 'observed_source'):
            if name in spec:
                inputs[spec[name]] = sha(Path(spec[name]).read_bytes())
        for path in spec.get('historical_angstrom', {}).values():
            inputs[path] = sha(Path(path).read_bytes())
        cases.append(row)
    sources = dict(parent['source_hashes'])
    for name in ('scripts/research_carbon_H_whole.py', 'scripts/preflight_graph_methyl_transport.py',
                 'scripts/research_h_stage_transport.py', 'scripts/research_beta_wall.py',
                 'src/dvbfixer/ffutils/geometry.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.carbon-H-whole-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, h_root=str(h_root), python=sys.executable, openmm_version=openmm.__version__,
        policy=dict(arms=parent['policy']['arms'], baseline_cells=34, continuation_cells=136,
            inherited_refused_arm_cells=24, platform='Reference', no_solvent=True,
            tolerance_kj_mol_nm=10., maximum_iterations_per_phase=1000,
            per_cell_seconds=1800, baseline_total_seconds=7200, continuation_total_seconds=21600,
            enrollment='ALL40 original contexts, all four factorial arms, complete source refusals; no new seed/retry/success-selected subset',
            baseline='independent native Simulation.minimizeEnergy from original saved State/System/integrator, both nm endpoints retained before any arm continuation',
            bridge='native H-label repair once then captured object-index map; original34 nm starts all exact; no alias guessing/observation restore',
            physics='original full no-solvent Systems, forces and anchors; only original tenfold positional-weight reduction in phase2',
            audit='exact unchanged starts and both independent native nm endpoints, honest legacy historical Angstrom readback; all raw/rounded global/focal shapes, existing masked target gate, H strain and original-heavy movement',
            limits='two native phases BEFORE emergency recovery; no reflection, native-accuracy score, new force/phase/H/seed/budget tuning or public candidate/PDB/production qualification'),
        production_changed=False))


def worker(plan_path, index, arm, baseline_root, output):
    import openmm
    from openmm.app import Simulation
    from openmm.unit import kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    case = record['cases'][index]
    if (case['status'] != 'admitted' or openmm.__version__ != record['openmm_version']
            or arm not in ('native_baseline', *record['policy']['arms'])):
        raise ValueError('undeclared whole case/arm/runtime')
    spec = case['whole']
    system = openmm.XmlSerializer.deserialize(Path(spec['system']).read_text())
    integrator = openmm.XmlSerializer.deserialize(Path(spec['integrator']).read_text())
    state = openmm.XmlSerializer.deserialize(Path(spec['state']).read_text())
    initial = [list(p) for p in state.getPositions().value_in_unit(nanometer)]
    repairs = None
    if arm != 'native_baseline':
        data = adapt(case)
        path = read_json(Path(record['h_root'])/case['id']/f"seed-{case['seed']}"/case['conditioning']/arm/'path.json.gz')
        initial, repairs = mapped_start(data['rows'], data['bonds'], path['final_nm'], spec['mapping'])
    simulation = Simulation(topology_from_rows(spec['rows'], spec['bonds']), system, integrator,
        openmm.Platform.getPlatformByName('Reference'))
    frames = []
    try:
        simulation.context.setState(state)
        simulation.context.setPositions(initial*nanometer)
        positional = restraint(system)
        for stage in ('initial', 'phase1', 'phase2'):
            if stage == 'phase2':
                for i in range(positional.getNumParticles()):
                    particle, parameters = positional.getParticleParameters(i)
                    positional.setParticleParameters(i, particle, [parameters[0]/10, *parameters[1:]])
                positional.updateParametersInContext(simulation.context)
            if stage != 'initial':
                if arm == 'native_baseline':
                    simulation.minimizeEnergy(maxIterations=1000)
                else:
                    openmm.LocalEnergyMinimizer.minimize(simulation.context, 10*kilojoules_per_mole/nanometer, 1000)
            current = simulation.context.getState(positions=True, energy=True)
            points = [list(p) for p in current.getPositions().value_in_unit(nanometer)]
            xml = openmm.XmlSerializer.serialize(system)
            payload = gzip.compress(json.dumps(dict(positions_nm=points,
                energy_kj_mol=current.getPotentialEnergy().value_in_unit(kilojoules_per_mole)),
                sort_keys=True, allow_nan=False).encode(), mtime=0)
            with (output/f'{stage}.json.gz').open('xb') as stream:
                stream.write(payload)
            with (output/f'{stage}-system.xml').open('x') as stream:
                stream.write(xml)
            frame = dict(stage=stage, coordinate_sha256=sha(payload), system_sha256=sha(xml.encode()))
            if arm != 'native_baseline':
                baseline = baseline_root/case['id']/f"seed-{case['seed']}"/case['conditioning']
                if xml != (baseline/f'{stage}-system.xml').read_text():
                    raise ValueError('whole native baseline forces/anchors differ')
                exact = points == read_json(baseline/f'{stage}.json.gz')['positions_nm']
                frame['exact_baseline_nm'] = exact
                if arm == 'unchanged' and not exact:
                    raise ValueError('unchanged whole native nm baseline differs; retain evidence')
            if 'historical_angstrom' in spec:
                prior = read_json(spec['historical_angstrom'][stage])['atoms']
                exact_A = labels(prior) == spec['rows'] and [r['xyz'] for r in prior] == [[10*x for x in p] for p in points]
                frame['exact_historical_angstrom'] = exact_A
                if arm == 'native_baseline' and not exact_A:
                    raise ValueError('legacy native baseline differs from historical Angstrom endpoint')
            frames.append(frame)
        write(output/'trace.json', dict(case_id=[case['id'], case['seed'], case['conditioning']],
            arm=arm, plan_sha256=sha(raw), frames=frames, post_H_repairs=repairs,
            emergency_reflection=False, observations_restored=False, production_changed=False))
    finally:
        del simulation
        del integrator


def run(plan_path, output, baseline_root=None):
    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    baseline_mode = baseline_root is None
    baseline_binding = None
    if not baseline_mode:
        base_raw = (baseline_root/'evidence.json').read_bytes()
        base = json.loads(base_raw)
        if (base['plan_sha256'] != sha(raw) or base['mode'] != 'native_baseline'
                or [(c['id'], c['seed'], c['conditioning']) for c in base['cases']]
                    != [(c['id'], c['seed'], c['conditioning']) for c in record['cases']]):
            raise ValueError('complete independent whole baseline differs')
        for c in base['cases']:
            if c['status'] not in {'completed', 'refused_source'}:
                raise ValueError('whole baseline has blocker; no continuations')
            folder = baseline_root/c['id']/f"seed-{c['seed']}"/c['conditioning']
            for name, digest in c['artifacts'].items():
                if sha((folder/name).read_bytes()) != digest:
                    raise ValueError('whole baseline artifact changed')
        baseline_binding = sha(base_raw)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(Path('src').resolve()), str(Path('.').resolve()))))
    results, started = [], time.monotonic()
    total = record['policy']['baseline_total_seconds' if baseline_mode else 'continuation_total_seconds']
    for index, case in enumerate(record['cases']):
        for arm in ('native_baseline',) if baseline_mode else record['policy']['arms']:
            folder = output/case['id']/f"seed-{case['seed']}"/case['conditioning']
            if not baseline_mode:
                folder = folder/arm
            folder = folder.resolve()
            folder.mkdir(parents=True, exist_ok=False)
            row = dict(id=case['id'], seed=case['seed'], conditioning=case['conditioning'], arm=arm,
                status=case['status'], returncode=None)
            if case['status'] == 'admitted':
                remaining = total-(time.monotonic()-started)
                if remaining <= 0:
                    row['status'] = 'budget_not_executed'
                else:
                    command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path.resolve()),
                        '--index', str(index), '--arm', arm, '--output', str(folder)]
                    if not baseline_mode:
                        command += ['--baseline-root', str(baseline_root.resolve())]
                    with (folder/'run.log').open('xb') as log:
                        try:
                            process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                timeout=min(remaining, record['policy']['per_cell_seconds']), check=False)
                            row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                        except subprocess.TimeoutExpired:
                            row['status'] = 'timeout'
            row['artifacts'] = {p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
            results.append(row)
            write(folder/'result.json', row)
            print(json.dumps({k:v for k, v in row.items() if k != 'artifacts'}), flush=True)
    verify(record)
    write(output/'evidence.json', dict(schema='dvbfixer.carbon-H-whole-run.v1', cases=results,
        mode='native_baseline' if baseline_mode else 'continuations', plan_sha256=sha(raw),
        baseline_run_sha256=baseline_binding, production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    for name in ('h-plan', 'h-root', 'h-audit', 'legacy-plan', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        p.add_argument('--baseline-root', type=Path)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
            p.add_argument('--arm', required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.h_plan, args.h_root, args.h_audit, args.legacy_plan, args.output)
    elif args.command == 'run':
        run(args.plan, args.output, args.baseline_root)
    else:
        worker(args.plan, args.index, args.arm, args.baseline_root, args.output)
