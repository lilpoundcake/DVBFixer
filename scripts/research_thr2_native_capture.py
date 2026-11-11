"""Bounded historical THR2 capture and separate saved-state observational replays.

Private coordinates/checkpoints stay below tmp/. Never modifies production source.
"""
from __future__ import annotations

import argparse
import functools
import inspect
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_legacy_start_controls import compressed, platform_properties, snapshot_atoms

ROOT = Path(__file__).resolve().parents[1]
PRIVATE = ROOT/'tmp/pdbfixer-research'
HISTORICAL = PRIVATE/'pre-stereo-baseline-01'
CASE = '1CRN-THR2-CB-OG1-CG2'


def update_owned(path, record):
    """Atomically refresh a progress file inside this run's exclusive directory."""
    staged = path.with_name(path.name+'.next')
    write(staged, record)
    staged.replace(path)


def freeze(output):
    import openmm
    import openmm.version
    from openmm.app import ForceField, Modeller, Simulation

    parent_path = PRIVATE/'legacy-comparator-plan-01.json'
    parent = json.loads(parent_path.read_bytes())
    inputs = {str(parent_path): sha(parent_path.read_bytes())}
    for name, digest in parent['source_hashes'].items():
        path = HISTORICAL/name
        if sha(path.read_bytes()) != digest:
            raise ValueError(f'historical source mismatch: {path}')
        inputs[str(path)] = digest
    for name in ('prepared.pdb', 'prepared.dat', 'minimize-execution.json', 'minimize.log', 'minimized.pdb'):
        path = PRIVATE/'legacy-comparator-run-01'/CASE/name
        inputs[str(path)] = sha(path.read_bytes())
    # Bind transitive research helpers too, without importing production parsers.
    sources = {str(path):sha(path.read_bytes()) for path in sorted((ROOT/'scripts').glob('*.py'))}
    for obj in (Modeller.addHydrogens, ForceField.createSystem, Simulation.minimizeEnergy,
                openmm.LocalEnergyMinimizer.minimize):
        path = Path(inspect.getfile(obj)).resolve()
        inputs[str(path)] = sha(path.read_bytes())
    write(output, dict(schema='dvbfixer.thr2-native-capture-plan.v1', inputs=inputs,
        source_hashes=sources, python=parent['python'], historical_root=str(HISTORICAL),
        original_case=next(c for c in parent['cases'] if c['id'] == CASE),
        runtime=dict(version=openmm.__version__, full_version=openmm.version.version,
            revision=openmm.version.git_revision),
        protocol=dict(labels=['unseeded', 'seed-1', 'seed-2'], h_seeds=[None, 1, 2],
            native='each label ONCE, historical CLI/default automatic platform; --no-solvent --atom-naming standard',
            observer='no native reporter; serialize before/after each LocalEnergyMinimizer call, including H; no energy/force readback in native Context',
            replay='each captured native call once per captured-checkpoint and Reference-XML mode; independent starts, not chained phases',
            maximum_native_calls_per_label=8, maximum_callbacks=4096,
            per_job_seconds=1800, total_seconds=14400,
            no_adaptation='retain failures/nonreproduction; no retries, new seeds, phase extensions or outcome-dependent enrollment',
            numerical='OpenCL is not bit reproducible; report initial/endpoint deltas, not Reference equivalence',
            limits='fresh bounded captures do not recreate uncaptured original failing trajectory; no production qualification'),
        production_changed=False))


def check(record):
    import openmm
    import openmm.version

    verify(record)
    if (openmm.version.version != record['runtime']['full_version']
            or openmm.version.git_revision != record['runtime']['revision']):
        raise ValueError('runtime provenance changed')


def capture_context(context, output, prefix, topology):
    import openmm

    state = context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
    for name, obj in (('system', context.getSystem()), ('integrator', context.getIntegrator()), ('state', state)):
        (output/f'{prefix}-{name}.xml').write_text(openmm.XmlSerializer.serialize(obj))
    (output/f'{prefix}.chk').write_bytes(context.createCheckpoint())
    compressed(output/f'{prefix}-atoms.json.gz', dict(atoms=snapshot_atoms(topology, state.getPositions())))
    return dict(prefix=prefix, platform=platform_properties(context.getPlatform(), context))


def native_worker(plan_path, index, output):
    import openmm
    from openmm.app import ForceField, Modeller
    from openmm.unit import kilojoules_per_mole, nanometer

    record = json.loads(plan_path.read_bytes())
    check(record)
    import dvbfixer.cli
    if Path(dvbfixer.cli.__file__).resolve() != HISTORICAL/'src/dvbfixer/cli.py':
        raise ValueError('not executing frozen historical scientific source')
    seed = record['protocol']['h_seeds'][index]
    original_h, original_system, original_min = Modeller.addHydrogens, ForceField.createSystem, openmm.LocalEnergyMinimizer.minimize
    registry, calls, h_calls = {}, [], []
    scope = {'hydrogen': False}

    def persist():
        update_owned(output/'capture.json', dict(plan_sha256=sha(plan_path.read_bytes()),
            label=record['protocol']['labels'][index], calls=calls, h_calls=h_calls,
            scientific_module=str(Path(dvbfixer.cli.__file__).resolve())))

    @functools.wraps(original_system)
    def create(ff, topology, *args, **kwargs):
        system = original_system(ff, topology, *args, **kwargs)
        registry[int(system.this)] = (system, topology)
        return system

    @functools.wraps(original_h)
    def hydrogen(model, *args, **kwargs):
        ordinal = len(h_calls)+1
        before = random.getstate()
        compressed(output/f'h-{ordinal}-before.json.gz', dict(atoms=snapshot_atoms(model.topology, model.positions)))
        if seed is not None:
            random.seed(seed)
        row = dict(seed=seed, before_rng=before, initialized_rng=random.getstate(), returned=False)
        h_calls.append(row)
        scope['hydrogen'] = True
        persist()
        try:
            result = original_h(model, *args, **kwargs)
            row['returned'] = True
            return result
        finally:
            row['after_rng'] = random.getstate()
            if seed is not None:
                random.setstate(before)
            row['restored_rng'] = random.getstate()
            scope['hydrogen'] = False
            compressed(output/f'h-{ordinal}-after.json.gz', dict(atoms=snapshot_atoms(model.topology, model.positions)))
            persist()

    @functools.wraps(original_min)
    def minimize(context, *args, **kwargs):
        if len(calls) >= record['protocol']['maximum_native_calls_per_label']:
            raise ValueError('prospective native-call bound exceeded')
        bound = inspect.signature(original_min).bind(context, *args, **kwargs)
        bound.apply_defaults()
        if bound.arguments['reporter'] is not None:
            raise ValueError('unexpected existing native reporter')
        tolerance = bound.arguments['tolerance']
        if hasattr(tolerance, 'value_in_unit'):
            tolerance = tolerance.value_in_unit(kilojoules_per_mole/nanometer)
        system = context.getSystem()
        topology = registry[int(system.this)][1]
        number = len(calls)+1
        write(output/f'call-{number}-topology.json', dict(
            bonds=[[a.index, b.index] for a, b in topology.bonds()]))
        row = dict(number=number, stage='hydrogen' if scope['hydrogen'] else 'whole',
            tolerance=float(tolerance), max_iterations=bound.arguments['maxIterations'], returned=False)
        calls.append(row)
        row['before'] = capture_context(context, output, f'call-{number}-before', topology)
        persist()
        try:
            result = original_min(context, *args, **kwargs)
            row['returned'] = True
            return result
        finally:
            row['after'] = capture_context(context, output, f'call-{number}-after', topology)
            persist()

    sys.argv = ['dvbfixer', 'minimize', str(output/'prepared.pdb'), '--dat', str(output/'prepared.dat'),
        '--no-solvent', '--atom-naming', 'standard', '-o', str(output/'minimized.pdb')]
    with patch.object(ForceField, 'createSystem', create), patch.object(Modeller, 'addHydrogens', hydrogen), \
            patch.object(openmm.LocalEnergyMinimizer, 'minimize', minimize):
        dvbfixer.cli.main()
    check(record)
    persist()


def replay_worker(plan_path, source, number, mode, output):
    import openmm
    from openmm.unit import nanometer

    record = json.loads(plan_path.read_bytes())
    check(record)
    capture = json.loads((source/'capture.json').read_bytes())
    row = next(c for c in capture['calls'] if c['number'] == number)
    prefix = row['before']['prefix']
    system = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((source/f'{prefix}-integrator.xml').read_text())
    if mode == 'captured':
        platform = row['before']['platform']
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(platform['name']), platform['constructor_properties'])
        context.loadCheckpoint((source/f'{prefix}.chk').read_bytes())
    else:
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        context.setState(openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()))
    recorder = reporter_type(openmm)(system.getNumParticles(), record['protocol']['maximum_callbacks'])
    initial = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    returned = False
    try:
        openmm.LocalEnergyMinimizer.minimize(context, row['tolerance'], row['max_iterations'], recorder)
        returned = True
    finally:
        final = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        compressed(output/'path.json.gz', dict(initial_nm=initial, final_nm=final,
            callbacks=recorder.records, returned=returned))
        write(output/'replay.json', dict(mode=mode, number=number, stage=row['stage'], returned=returned,
            callbacks=len(recorder.records), platform=platform_properties(context.getPlatform(), context)))
        del context


def launch(record, command, output, timeout):
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(HISTORICAL/'src'), str(ROOT))),
        PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH', ''))
    with (output/'run.log').open('xb') as log:
        try:
            process = subprocess.run([record['python'], str(Path(__file__).resolve()), *command],
                env=env, cwd=HISTORICAL, stdout=log, stderr=subprocess.STDOUT, timeout=timeout, check=False)
            result = dict(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
        except subprocess.TimeoutExpired:
            result = dict(status='timeout', returncode=None)
    result['artifacts'] = {p.name:sha(p.read_bytes()) for p in sorted(output.iterdir()) if p.is_file()}
    write(output/'result.json', result)
    return result


def run(plan_path, output):
    record = json.loads(plan_path.read_bytes())
    check(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    started = time.monotonic()
    native, replays = [], []

    def execute(command, folder):
        remaining = record['protocol']['total_seconds']-(time.monotonic()-started)
        if remaining <= 0:
            result = dict(status='budget_not_executed', returncode=None)
            write(folder/'result.json', result)
            return result
        return launch(record, command, folder, min(remaining, record['protocol']['per_job_seconds']))

    for index, label in enumerate(record['protocol']['labels']):
        folder = output/label
        folder.mkdir()
        for name in ('prepared.pdb', 'prepared.dat'):
            (folder/name).write_bytes((PRIVATE/'legacy-comparator-run-01'/CASE/name).read_bytes())
        command = ['native-worker', '--plan', str(plan_path), '--index', str(index), '--output', str(folder)]
        native.append(dict(label=label, **execute(command, folder)))
        update_owned(output/'native-evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=native))
    # Native cohort is finished before any replay. Bind all captured payloads first.
    bindings = {str(p):sha(p.read_bytes()) for label in record['protocol']['labels']
        for p in sorted((output/label).iterdir()) if p.is_file()}
    write(output/'replay-inputs.json', dict(plan_sha256=sha(plan_path.read_bytes()), inputs=bindings))
    for label in record['protocol']['labels']:
        source = output/label
        if not (source/'capture.json').exists():
            continue
        capture = json.loads((source/'capture.json').read_bytes())
        for row in capture['calls']:
            if 'before' not in row:
                continue
            for mode in ('captured', 'Reference'):
                folder = output/'replays'/label/f"call-{row['number']}"/mode
                folder.mkdir(parents=True, exist_ok=False)
                command = ['replay-worker', '--plan', str(plan_path), '--source', str(source),
                    '--number', str(row['number']), '--mode', mode, '--output', str(folder)]
                replays.append(dict(label=label, number=row['number'], mode=mode, **execute(command, folder)))
                update_owned(output/'replay-evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=replays))
    for name, digest in bindings.items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError('native capture altered during replay')
    check(record)
    write(output/'evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), native=native, replays=replays,
        production_changed=False))
    print(json.dumps(dict(native=[r['status'] for r in native], replays=[r['status'] for r in replays])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'native-worker', 'replay-worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'native-worker':
            p.add_argument('--index', type=int, choices=range(3), required=True)
        if name == 'replay-worker':
            p.add_argument('--source', type=Path, required=True)
            p.add_argument('--number', type=int, required=True)
            p.add_argument('--mode', choices=('captured', 'Reference'), required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    elif args.command == 'run':
        run(args.plan.resolve(), args.output.resolve())
    elif args.command == 'native-worker':
        native_worker(args.plan.resolve(), args.index, args.output.resolve())
    else:
        replay_worker(args.plan.resolve(), args.source.resolve(), args.number, args.mode, args.output.resolve())
