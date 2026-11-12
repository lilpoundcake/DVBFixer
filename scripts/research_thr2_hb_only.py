"""Saved-state HB-only initializer intervention, unchanged three native phases."""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import time
from pathlib import Path

from scripts.audit_thr2_native_capture import read, verify_files
from scripts.qualify_hydrogen_postflight import load, reconstruct
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_thr2_hswap import replace_restraint_parameters
from scripts.research_thr2_initializer_frames import capture_data
from scripts.research_thr2_native_capture import check

ROOT = Path(__file__).resolve().parents[1]
HISTORICAL_GEOMETRY = ROOT/'tmp/pdbfixer-research/pre-stereo-baseline-01/src/dvbfixer/ffutils/geometry.py'


def proposed_points(source, initial, arm):
    from scripts.audit_thr2_h_initialization import initializer

    result = [list(p) for p in initial]
    if arm == 'unchanged':
        return result, []
    if arm != 'HB-only-unperturbed':
        raise ValueError('Undeclared initializer arm')
    hb = source['hb']
    base, _ = initializer(source['heavy'][hb['parent']], [source['heavy'][j] for j in hb['neighbors']], hb['draws'])
    result[hb['index']] = base
    return result, [hb['index']]


def transport(requested, actual, mode):
    import math

    import numpy as np

    if mode not in ('captured', 'Reference'):
        raise ValueError('Undeclared numerical mode')
    if len(requested) != len(actual) or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in requested+actual):
        raise ValueError('Transport inventory differs')
    if mode == 'Reference':
        if actual != requested:
            raise ValueError('Reference coordinate transport is not exact')
        return 0.
    maximum = 0.
    for p, q in zip(requested, actual, strict=True):
        for a, b in zip(p, q, strict=True):
            rounded = np.float32(a)
            if not np.isfinite(rounded):
                raise ValueError('Requested coordinate is not finite float32')
            spacing = max(abs(float(np.nextafter(rounded, np.float32('inf')))-float(rounded)),
                          abs(float(np.nextafter(rounded, np.float32('-inf')))-float(rounded)))
            if abs(a-b) > 0.5*spacing+1e-15:
                raise ValueError('Applied coordinate exceeds declared float32 rounding bound')
            maximum = max(maximum, abs(a-b))
    return maximum


def check_plan(record):
    verify(record)
    check(read(Path(record['parent_plan'])))


def freeze(output):
    import openmm
    from openmm.unit import dalton

    from dvbfixer.ffutils.variants import VARIANT_TO_PARENT

    parent_path = ROOT/'docs/research/thr2-native-capture-plan.json'
    parent = read(parent_path)
    check(parent)
    preceding_plans = [ROOT/'docs/research/thr2-initializer-frames-plan.json', ROOT/'docs/research/thr2-hswap-plan.json']
    for preceding in preceding_plans:
        verify(read(preceding))
    root = ROOT/'tmp/pdbfixer-research/thr2-native-capture-run-01'
    execution = read(root/'evidence.json')
    inputs = {str(parent_path):sha(parent_path.read_bytes()), str(root/'evidence.json'):sha((root/'evidence.json').read_bytes()),
              str(HISTORICAL_GEOMETRY):sha(HISTORICAL_GEOMETRY.read_bytes())}
    for preceding in preceding_plans:
        inputs[str(preceding)] = sha(preceding.read_bytes())
        inputs.update(read(preceding)['inputs'])
    for name in ('docs/plans/thr2-hb-only-initializer.md', 'docs/research/thr2-hbridge-admission.json',
                 'docs/research/thr2-initializer-frames-evidence.json'):
        p = ROOT/name
        inputs[str(p)] = sha(p.read_bytes())
    sources = []
    for row in execution['native']:
        folder = root/row['label']
        verify_files(folder, row)
        inputs.update({str(folder/name):digest for name, digest in row['artifacts'].items()})
        calls = read(folder/'capture.json')['calls']
        if [(c['stage'], c['tolerance'], c['max_iterations']) for c in calls] != [
                ('hydrogen', 1., 50), ('whole', 10., 1000), ('whole', 10., 1000)]:
            raise ValueError('Original three native phases required')
        h_atoms = read(folder/f"{calls[0]['before']['prefix']}-atoms.json.gz")['atoms']
        whole_atoms = read(folder/f"{calls[1]['before']['prefix']}-atoms.json.gz")['atoms']
        def canonical(row):
            key = row['key']
            return key[:3]+[VARIANT_TO_PARENT.get(key[3], key[3])]+key[4:], row['element']
        if [canonical(r) for r in h_atoms] != [canonical(r) for r in whole_atoms]:
            raise ValueError('Phase identity differs beyond canonical residue aliases')
        graphs = [sorted(tuple(sorted(pair)) for pair in read(folder/f"call-{i}-topology.json")['bonds']) for i in (1, 2, 3)]
        if graphs[0] != graphs[1] or graphs[1] != graphs[2]:
            raise ValueError('Original phase bond graphs differ')
        for call in calls:
            system = openmm.XmlSerializer.deserialize((folder/f"{call['before']['prefix']}-system.xml").read_text())
            if system.getNumConstraints() or system.usesPeriodicBoundaryConditions() or any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
                raise ValueError('Only original unconstrained nonperiodic systems admitted')
            mobile = [i for i in range(system.getNumParticles()) if system.getParticleMass(i).value_in_unit(dalton) > 0]
            wanted = [i for i, r in enumerate(h_atoms) if r['element'] == 'H'] if call['stage'] == 'hydrogen' else list(range(len(h_atoms)))
            if mobile != wanted:
                raise ValueError('Native mobility mask differs')
        replay = root/'replays'/row['label']/'call-1'/'Reference'/'path.json.gz'
        inputs[str(replay)] = sha(replay.read_bytes())
        capture_data(folder)  # Saved RNG/graph admission only, no intervention.
        sources.append(dict(label=row['label'], source=str(folder), calls=calls, reference_h_control=str(replay)))
    paths = [Path(__file__).resolve(), ROOT/'scripts/audit_thr2_hb_only.py', ROOT/'scripts/research_thr2_initializer_frames.py',
             ROOT/'scripts/audit_thr2_native_capture.py', ROOT/'scripts/research_thr2_hswap.py',
             ROOT/'scripts/qualify_hydrogen_postflight.py', ROOT/'src/dvbfixer/ffutils/variants.py',
             ROOT/'tests/test_thr2_hb_only.py']
    write(output, dict(schema='dvbfixer.thr2-hb-only-plan.v1', parent_plan=str(parent_path),
        inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in paths}, python=parent['python'],
        historical_geometry=str(HISTORICAL_GEOMETRY), sources=sources,
        cases=[dict(label=s['label'], arm=arm, mode=mode) for s in sources
               for arm in ('unchanged', 'HB-only-unperturbed') for mode in ('captured', 'Reference')],
        protocol=dict(maximum_callbacks=4096, per_path_seconds=1200, total_seconds=7200, paths=12, native_calls=36,
            H_native_settings=[1., 50], whole_native_settings=[10., 1000], new_random_sampling=False,
            limits='Six proposals/two numerical modes in one retained poor-heavy context; not independent seeds, native H sampling, original lost run or qualified production initializer'),
        production_changed=False, main_merge_authorized=False))


def make_context(folder, call, mode):
    import openmm

    prefix = call['before']['prefix']
    system = openmm.XmlSerializer.deserialize((folder/f'{prefix}-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((folder/f'{prefix}-integrator.xml').read_text())
    declared = call['before']['platform']
    if mode == 'captured':
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(declared['name']), declared['constructor_properties'])
        context.loadCheckpoint((folder/f'{prefix}.chk').read_bytes())
    else:
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        context.setState(openmm.XmlSerializer.deserialize((folder/f'{prefix}-state.xml').read_text()))
    return system, integrator, context


def worker(plan_path, index, output):
    import openmm
    from openmm.unit import nanometer

    record = read(plan_path)
    check_plan(record)
    case = record['cases'][index]
    source = next(s for s in record['sources'] if s['label'] == case['label'])
    folder = Path(source['source'])
    data = capture_data(folder)
    old_geometry = load('thr2_hb_only_historical_geometry', Path(record['historical_geometry']))
    system, integrator, context = make_context(folder, source['calls'][0], case['mode'])
    original = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    requested, selected = proposed_points(data, original, case['arm'])
    if selected:
        context.setPositions(requested*nanometer)
    phases, bridge = [], None
    try:
        for number, call in enumerate(source['calls'], 1):
            if number == 2:
                returned_positions = context.getState(positions=True).getPositions()
                topology, _ = reconstruct(data['atoms'], read(folder/'call-1-topology.json')['bonds'])
                before = [list(p) for p in returned_positions.value_in_unit(nanometer)]
                repairs = old_geometry.repair_misplaced_hydrogens(topology, returned_positions)
                requested = [list(p) for p in returned_positions.value_in_unit(nanometer)]
                bridge = dict(repairs=repairs, changed_atom_indices=[i for i, (a, b) in enumerate(zip(before, requested, strict=True)) if a != b])
                compressed(output/'historical-bridge.json.gz', dict(before_nm=before, after_nm=requested, **bridge))
                context = None
                system, integrator, context = make_context(folder, call, case['mode'])
                context.setPositions(requested*nanometer)
            elif number == 3:
                target = openmm.XmlSerializer.deserialize((folder/f"{call['before']['prefix']}-system.xml").read_text())
                restraint = replace_restraint_parameters(system, target)
                restraint.updateParametersInContext(context)
                requested = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            phase = output/f'phase-{number}'
            phase.mkdir()
            serialized = openmm.XmlSerializer.serialize(system)
            (phase/'system.xml').write_text(serialized)
            if serialized != (folder/f"{call['before']['prefix']}-system.xml").read_text():
                raise ValueError('Phase System is not the exact original native System')
            initial = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            compressed(phase/'transport.json.gz', dict(requested_nm=requested, actual_nm=initial))
            compressed(phase/'requested.json.gz', dict(positions_nm=requested))
            maximum = transport(requested, initial, case['mode'])
            if number == 1 and any(a != b for i, (a, b) in enumerate(zip(original, initial, strict=True)) if i not in selected):
                raise ValueError('H initial intervention changed unselected atoms')
            reporter = reporter_type(openmm)(len(data['atoms']), record['protocol']['maximum_callbacks'])
            returned = False
            try:
                openmm.LocalEnergyMinimizer.minimize(context, call['tolerance'], call['max_iterations'], reporter)
                returned = True
            finally:
                final = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
                compressed(phase/'path.json.gz', dict(initial_nm=initial, final_nm=final, callbacks=reporter.records, returned=returned))
                (phase/'returned-state.xml').write_text(openmm.XmlSerializer.serialize(context.getState(positions=True, parameters=True)))
                phases.append(dict(number=number, returned=returned, callbacks=len(reporter.records), maximum_transport_component_error_nm=maximum))
                write(output/f'phase-{number}-result.json', dict(case=case, phases=list(phases), bridge=bridge, selected_atom_indices=selected,
                    plan_sha256=sha(plan_path.read_bytes())))
            if number == 1 and any(final[i] != initial[i] for i in data['heavy']):
                raise ValueError('Native H-only stage moved pinned heavy coordinates')
    finally:
        if context is not None:
            del context


def run(plan_path, output):
    record = read(plan_path)
    check_plan(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(ROOT/'src'), str(ROOT))))
    started, results = time.monotonic(), []
    for index, case in enumerate(record['cases']):
        folder = output/f'path-{index:02d}'
        folder.mkdir()
        row = dict(index=index, **case, status='budget_not_executed', returncode=None)
        remaining = record['protocol']['total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path), '--index', str(index), '--output', str(folder)]
            with (folder/'run.log').open('xb') as log:
                process = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=min(remaining, record['protocol']['per_path_seconds']))
                    row.update(status='completed' if code == 0 else 'execution_blocker', returncode=code)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    row['status'] = 'timeout'
        row['artifacts'] = {str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        write(folder/'result.json', row)
        results.append(row)
        print(json.dumps({k:v for k, v in row.items() if k != 'artifacts'}), flush=True)
    check_plan(record)
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
            p.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    elif args.command == 'worker':
        worker(args.plan.resolve(), args.index, args.output.resolve())
    else:
        run(args.plan.resolve(), args.output.resolve())
