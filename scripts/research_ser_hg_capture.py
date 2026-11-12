"""Fresh bounded SER-HG captures, never replacements for the lost pytest output.

Six fixed preparations: native unseeded plus Python H seeds 1..5. Defaults and
production safeguards are forwarded unchanged. Private evidence stays in-project.
"""
from __future__ import annotations

import argparse
import ast
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
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed, snapshot_atoms
from scripts.research_thr2_native_capture import capture_context, update_owned

ROOT = Path(__file__).resolve().parents[1]
TEST = ROOT/'tests/test_prepare_broken_geom.py'
CONSTANT = '_SER_INPUT_HG_MISSING_OXT_PRESENT'


def synthetic_input(path=TEST):
    """Extract the tracked literal without importing/executing a pytest module."""
    matches = [node.value for node in ast.parse(path.read_text()).body if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == CONSTANT for target in node.targets)]
    if len(matches) != 1:
        raise ValueError('exact original synthetic input literal is not unique')
    value = ast.literal_eval(matches[0])
    if not isinstance(value, str):
        raise ValueError('synthetic input is not literal text')
    return value.encode()


def measurements(rows):
    import math

    selected = [r for r in rows if r['key'][:3] == ['B', '126', '']]
    points = {r['key'][-1]:r['xyz'] for r in selected}
    if len(points) != len(selected):
        raise ValueError('duplicate target atom identity')
    result = dict(residue=['B', '126', ''], atom_names=sorted(points))
    oxygen = 'OXT' if 'OXT' in points else 'OC1' if 'OC1' in points else None
    if 'OG' in points and 'HG' in points:
        result['og_hg_angstrom'] = math.dist(points['OG'], points['HG'])
        result['satisfies_original_og_hg_assertion'] = result['og_hg_angstrom'] < 1.2
    if oxygen is not None and 'HG' in points:
        result['hg_terminal_oxygen_angstrom'] = math.dist(points['HG'], points[oxygen])
        result['satisfies_original_hg_oxt_assertion'] = result['hg_terminal_oxygen_angstrom'] > 1.
    return result


def findings(topology, positions):
    from dvbfixer.ffutils.geometry import detect_misplaced_hydrogens

    def key(atom):
        return [atom.residue.chain.id, str(atom.residue.id), atom.residue.insertionCode.strip(), atom.residue.name, atom.name]
    return [dict(hydrogen=key(f.h_atom), parent=key(f.parent), reason=f.reason,
        parent_distance_angstrom=10*f.parent_distance_nm, expected_bond_angstrom=10*f.expected_bond_nm,
        coincident_with=key(f.coincident_with) if f.coincident_with is not None else None)
        for f in detect_misplaced_hydrogens(topology, positions)]


def freeze(output):
    import openmm
    import openmm.app.modeller
    import openmm.version

    sources = {str(p):sha(p.read_bytes()) for p in sorted((ROOT/'src/dvbfixer').rglob('*.py'))}
    sources.update({str(p):sha(p.read_bytes()) for p in sorted((ROOT/'scripts').glob('*.py'))})
    inputs = {str(TEST):sha(TEST.read_bytes())}
    for p in (Path(openmm.app.modeller.__file__), ROOT/'docs/research/thr2-hswap-validation.json',
              ROOT/'tmp/pdbfixer-research/thr2-hswap-fast-suite-01.log', ROOT/'tmp/pdbfixer-research/thr2-hswap-fast-suite-01.xml'):
        inputs[str(p)] = sha(p.read_bytes())
    raw = synthetic_input()
    write(output, dict(schema='dvbfixer.ser-hg-capture-plan.v1', python=sys.executable, root=str(ROOT),
        source_hashes=sources, inputs=inputs, synthetic_input_sha256=sha(raw),
        runtime=dict(version=openmm.__version__, full_version=openmm.version.version, revision=openmm.version.git_revision),
        cases=[dict(label='unseeded', seed=None)]+[dict(label=f'seed-{seed}', seed=seed) for seed in range(1, 6)],
        protocol=dict(enrollment='ALL six fixed fresh prepare invocations; no repeat-until-failure or pass filter',
            command='same prepare main([input, -o, output]) as tracked test; no added science flags',
            seed='save/restore Python RNG ONLY around Modeller.addHydrogens for seeds1..5; unseeded natural evolution',
            capture='before/after native H Context System/Integrator/State/checkpoint, topology; H/repair boundaries and full RNG state',
            observation='read-only observations; native minimizer/default settings unchanged; no native reporter or extra energy/force query',
            maximum_h_calls=8, maximum_native_h_calls=8, per_case_seconds=600, total_seconds=3600,
            retention='exclusive project-local directories, intermediate snapshots and all failures/timeouts/nonexecutions; no adaptive expansion',
            limits='fresh synthetic input preparations are not recovered original pytest coordinates/trajectory, independent new structural contexts or production qualification'),
        original_failed_artifacts='unavailable: user copy command reported source path missing; original log/XML retained separately',
        production_changed=False, scientific_thresholds_changed=False))


def check(record):
    import openmm.version

    verify(record)
    if openmm.version.version != record['runtime']['full_version'] or openmm.version.git_revision != record['runtime']['revision']:
        raise ValueError('native runtime provenance differs')


def worker(plan_path, index, output):
    import openmm
    from openmm.app import ForceField, Modeller
    from openmm.unit import kilojoules_per_mole, nanometer

    from dvbfixer.ffutils import geometry
    from dvbfixer.prepare import main

    record = json.loads(plan_path.read_bytes())
    check(record)
    if Path(inspect.getfile(geometry)).resolve() != ROOT/'src/dvbfixer/ffutils/geometry.py':
        raise ValueError('not observing selected current production source')
    case = record['cases'][index]
    original_h, original_create, original_min, original_repair = Modeller.addHydrogens, ForceField.createSystem, openmm.LocalEnergyMinimizer.minimize, geometry.repair_misplaced_hydrogens
    registry, h_calls, native, repairs = {}, [], [], []
    scope = dict(hydrogen=False)

    def progress():
        update_owned(output/'capture.json', dict(case=case, h_calls=h_calls, native=native, repairs=repairs,
            plan_sha256=sha(plan_path.read_bytes())))

    def atoms(prefix, topology, positions):
        rows = snapshot_atoms(topology, positions)
        compressed(output/f'{prefix}-atoms.json.gz', dict(atoms=rows))
        write(output/f'{prefix}-topology.json', dict(bonds=[[a.index, b.index] for a, b in topology.bonds()]))
        return measurements(rows)

    @functools.wraps(original_create)
    def create(ff, topology, *args, **kwargs):
        system = original_create(ff, topology, *args, **kwargs)
        registry[int(system.this)] = (system, topology)
        return system

    @functools.wraps(original_h)
    def hydrogen(model, *args, **kwargs):
        if len(h_calls) >= record['protocol']['maximum_h_calls']:
            raise ValueError('prospective H-call bound exceeded')
        ordinal = len(h_calls)+1
        saved = random.getstate()
        row = dict(number=ordinal, before=atoms(f'h-{ordinal}-before', model.topology, model.positions),
            before_rng=saved, returned=False)
        if case['seed'] is not None:
            random.seed(case['seed'])
        row['initialized_rng'] = random.getstate()
        h_calls.append(row)
        scope['hydrogen'] = True
        progress()
        try:
            result = original_h(model, *args, **kwargs)
            row['returned'] = True
            return result
        finally:
            row['after_rng'] = random.getstate()
            if case['seed'] is not None:
                random.setstate(saved)
            row['restored_rng'] = random.getstate()
            scope['hydrogen'] = False
            row['after'] = atoms(f'h-{ordinal}-after', model.topology, model.positions)
            progress()

    @functools.wraps(original_min)
    def minimize(context, *args, **kwargs):
        if not scope['hydrogen']:
            return original_min(context, *args, **kwargs)
        if len(native) >= record['protocol']['maximum_native_h_calls']:
            raise ValueError('prospective native H-call bound exceeded')
        bound = inspect.signature(original_min).bind(context, *args, **kwargs)
        bound.apply_defaults()
        if bound.arguments['reporter'] is not None:
            raise ValueError('unexpected native reporter')
        system = context.getSystem()
        if int(system.this) in registry:
            topology = registry[int(system.this)][1]
            mapping = 'native ForceField.createSystem topology'
        else:
            # Native no-forcefield fallback constructs this topology locally.
            frame = inspect.currentframe()
            try:
                topology = frame.f_back.f_locals['newTopology']
            finally:
                del frame
            mapping = 'read-only native addHydrogens newTopology frame'
        if topology.getNumAtoms() != system.getNumParticles():
            raise ValueError('native H atom/System inventory differs')
        number = len(native)+1
        tolerance = bound.arguments['tolerance']
        if hasattr(tolerance, 'value_in_unit'):
            tolerance = tolerance.value_in_unit(kilojoules_per_mole/nanometer)
        row = dict(number=number, mapping=mapping, tolerance=float(tolerance),
            max_iterations=bound.arguments['maxIterations'], returned=False)
        native.append(row)
        write(output/f'native-{number}-topology.json', dict(bonds=[[a.index, b.index] for a, b in topology.bonds()]))
        row['before'] = capture_context(context, output, f'native-{number}-before', topology)
        progress()
        try:
            result = original_min(context, *args, **kwargs)
            row['returned'] = True
            return result
        finally:
            row['after'] = capture_context(context, output, f'native-{number}-after', topology)
            progress()

    @functools.wraps(original_repair)
    def repair(topology, positions, *args, **kwargs):
        number = len(repairs)+1
        row = dict(number=number, before=atoms(f'repair-{number}-before', topology, positions),
            findings_before=findings(topology, positions), returned=False)
        repairs.append(row)
        progress()
        try:
            result = original_repair(topology, positions, *args, **kwargs)
            row.update(returned=True, repaired_count=result)
            return result
        finally:
            row['after'] = atoms(f'repair-{number}-after', topology, positions)
            row['findings_after'] = findings(topology, positions)
            progress()

    with patch.object(ForceField, 'createSystem', create), patch.object(Modeller, 'addHydrogens', hydrogen), \
            patch.object(openmm.LocalEnergyMinimizer, 'minimize', minimize), patch.object(geometry, 'repair_misplaced_hydrogens', repair):
        main([str(output/'input.pdb'), '-o', str(output/'prepared.pdb')])
    check(record)
    progress()


def run(plan_path, output):
    record = json.loads(plan_path.read_bytes())
    check(record)
    raw = synthetic_input()
    if sha(raw) != record['synthetic_input_sha256']:
        raise ValueError('frozen synthetic input differs')
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    env = dict(os.environ, PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH', ''),
        PYTHONPATH=os.pathsep.join((str(ROOT/'src'), str(ROOT))))
    started, results = time.monotonic(), []
    for index, case in enumerate(record['cases']):
        folder = output/case['label']
        folder.mkdir()
        (folder/'input.pdb').write_bytes(raw)
        row = dict(case, status='budget_not_executed', returncode=None)
        remaining = record['protocol']['total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path),
                '--index', str(index), '--output', str(folder)]
            with (folder/'prepare.log').open('xb') as log:
                try:
                    process = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
                        timeout=min(remaining, record['protocol']['per_case_seconds']), check=False)
                    row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                except subprocess.TimeoutExpired:
                    row['status'] = 'timeout'
        row['artifacts'] = {p.name:sha(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
        write(folder/'result.json', row)
        results.append(row)
        print(json.dumps({k:v for k, v in row.items() if k != 'artifacts'}), flush=True)
    check(record)
    write(output/'evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=results, production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, choices=range(6), required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    elif args.command == 'run':
        run(args.plan.resolve(), args.output.resolve())
    else:
        worker(args.plan.resolve(), args.index, args.output.resolve())
