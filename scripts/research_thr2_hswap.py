"""Prospectively bounded saved-start H transplant controls; unchanged native FF.

Three baselines + all six ordered donor/recipient pairs x HB/remainder masks,
each once on captured OpenCL and Reference. No new H generation or random seeds.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

from scripts.audit_thr2_native_capture import read, verify_files
from scripts.research_beta_wall import restraint
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_thr2_native_capture import check


def replace_restraint_parameters(system, target):
    """Only apply the exact already-captured production phase-2 parameter set."""
    import openmm

    current, wanted = restraint(system), restraint(target)
    if current.getNumParticles() != wanted.getNumParticles():
        raise ValueError('phase restraint inventory differs')
    for index in range(current.getNumParticles()):
        a, old = current.getParticleParameters(index)
        b, new = wanted.getParticleParameters(index)
        if a != b or list(old[1:]) != list(new[1:]) or new[0] != old[0]/10:
            raise ValueError('not the native same-anchor divide-ten transition')
        current.setParticleParameters(index, b, new)
    if openmm.XmlSerializer.serialize(system) != openmm.XmlSerializer.serialize(target):
        raise ValueError('phase transition changes anything beyond native restraint weights')
    return current


def swap_points(recipient, donor, elements, names, mask):
    if len(recipient) != len(donor) or len(recipient) != len(elements) or len(names) != len(elements):
        raise ValueError('complete swap atom inventory differs')
    hb = [i for i, key in enumerate(names) if key == ['A', '2', '', 'THR', 'HB']]
    if len(hb) != 1 or elements[hb[0]] != 'H':
        raise ValueError('exact target HB not uniquely identified')
    if mask not in ('unchanged', 'HB', 'other-H'):
        raise ValueError('undeclared coordinate mask')
    selected = [i for i, element in enumerate(elements) if element in ('H', 'D')
        and ((mask == 'HB' and i == hb[0]) or (mask == 'other-H' and i != hb[0]))]
    result = [list(p) for p in recipient]
    for i in selected:
        result[i] = list(donor[i])
    return result, selected


def freeze(parent_path, root, angle_path, output):
    import openmm
    from openmm.unit import dalton, nanometer

    parent = read(parent_path)
    check(parent)
    execution = read(root/'evidence.json')
    inputs = {str(parent_path):sha(parent_path.read_bytes()), str(root/'evidence.json'):sha((root/'evidence.json').read_bytes()),
        str(angle_path):sha(angle_path.read_bytes())}
    if execution['plan_sha256'] != sha(parent_path.read_bytes()) or len(execution['native']) != 3:
        raise ValueError('complete three-start source required')
    for row in execution['replays']:
        folder = root/'replays'/row['label']/f"call-{row['number']}"/row['mode']
        verify_files(folder, row)
        inputs.update({str(folder/name):digest for name, digest in row['artifacts'].items()})
    starts = []
    reference = None
    for row in execution['native']:
        if row['status'] != 'completed':
            raise ValueError('incomplete source capture')
        source = root/row['label']
        verify_files(source, row)
        inputs.update({str(source/name):digest for name, digest in row['artifacts'].items()})
        capture = read(source/'capture.json')
        whole = [c for c in capture['calls'] if c['stage'] == 'whole']
        if len(whole) != 2 or any(c['tolerance'] != 10. or c['max_iterations'] != 1000 for c in whole):
            raise ValueError('original two-phase native settings required')
        prefix = whole[0]['before']['prefix']
        system = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
        target = openmm.XmlSerializer.deserialize((source/f"{whole[1]['before']['prefix']}-system.xml").read_text())
        if (system.getNumConstraints() or system.usesPeriodicBoundaryConditions()
                or any(system.isVirtualSite(i) or system.getParticleMass(i).value_in_unit(dalton) <= 0 for i in range(system.getNumParticles()))):
            raise ValueError('unconstrained nonperiodic all-mobile source required')
        replace_restraint_parameters(openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system)), target)
        atoms = read(source/f'{prefix}-atoms.json.gz')['atoms']
        xyz = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        data = dict(system=(source/f'{prefix}-system.xml').read_bytes(),
            phase2=(source/f"{whole[1]['before']['prefix']}-system.xml").read_bytes(),
            keys=[a['key'] for a in atoms], elements=[a['element'] for a in atoms], xyz=xyz)
        if reference is None:
            reference = data
        if (data['system'] != reference['system'] or data['phase2'] != reference['phase2']
                or data['keys'] != reference['keys'] or data['elements'] != reference['elements']
                or any(a != b for a, b, element in zip(xyz, reference['xyz'], data['elements'], strict=True) if element not in ('H', 'D'))):
            raise ValueError('sources are not exact matched-heavy/FF starts')
        starts.append(dict(label=row['label'], source=str(source), calls=whole))
    cases = []
    for recipient in starts:
        arms = [(recipient['label'], 'unchanged')]+[(donor['label'], mask) for donor in starts
            if donor['label'] != recipient['label'] for mask in ('HB', 'other-H')]
        for donor, mask in arms:
            for mode in ('captured', 'Reference'):
                cases.append(dict(recipient=recipient['label'], donor=donor, mask=mask, mode=mode))
    paths = [Path(__file__).resolve(), Path('scripts/audit_thr2_native_capture.py').resolve(),
        Path('scripts/audit_thr2_angle_terms.py').resolve(), Path('scripts/research_beta_wall.py').resolve()]
    write(output, dict(schema='dvbfixer.thr2-hswap-plan.v1', parent_plan=str(parent_path),
        inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in paths}, python=parent['python'],
        starts=starts, cases=cases, policy=dict(paths=30, native_calls=60, maximum_callbacks_per_phase=4096,
            maximum_seconds_per_path=1800, maximum_total_seconds=14400,
            enrollment='three unchanged baselines and ALL six ordered donor/recipient pairs x target HB / complementary H, in both numerical modes, once each',
            intervention='copy named source H coordinates only; exact same initial heavy coordinates, topology, phase1 System and native phase2 transition',
            selection='HB mask selected after saved-only angle decomposition, before any transplant minimization; mechanism cohort, not independent validation',
            physics='native 10/1000 twice, no solvent; no new H, seeds, forces, constraints, repair, reflections or coordinate restoration',
            replay='phase2 chains from its own phase1 result; compare unchanged Reference phase1 against prior saved-start control exactly',
            retention='all outcomes and early failures/budget nonexecutions; no retry, donor/mask/seed/platform/strength search or adaptive extension',
            limits='saved-start counterfactuals on ONE poor-heavy THR2 context, not original lost trajectory, new independent sampling or production qualification'),
        production_changed=False))


def check_plan(record):
    verify(record)
    check(read(Path(record['parent_plan'])))


def worker(plan_path, index, output):
    import openmm
    from openmm.unit import nanometer

    record = read(plan_path)
    check_plan(record)
    case = record['cases'][index]
    recipient = next(c for c in record['starts'] if c['label'] == case['recipient'])
    donor = next(c for c in record['starts'] if c['label'] == case['donor'])
    source, donor_source = Path(recipient['source']), Path(donor['source'])
    prefix = recipient['calls'][0]['before']['prefix']
    donor_prefix = donor['calls'][0]['before']['prefix']
    state = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text())
    donor_state = openmm.XmlSerializer.deserialize((donor_source/f'{donor_prefix}-state.xml').read_text())
    xyz = state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    donor_xyz = donor_state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    atoms = read(source/f'{prefix}-atoms.json.gz')['atoms']
    imposed, indices = swap_points(xyz, donor_xyz, [a['element'] for a in atoms], [a['key'] for a in atoms], case['mask'])
    system = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((source/f'{prefix}-integrator.xml').read_text())
    platform = recipient['calls'][0]['before']['platform']
    if case['mode'] == 'captured':
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(platform['name']), platform['constructor_properties'])
        context.loadCheckpoint((source/f'{prefix}.chk').read_bytes())
    else:
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        context.setState(state)
    if indices:
        context.setPositions(imposed*nanometer)
    phases = []
    try:
        for phase in range(2):
            if phase == 1:
                target = openmm.XmlSerializer.deserialize((source/f"{recipient['calls'][1]['before']['prefix']}-system.xml").read_text())
                positional = replace_restraint_parameters(system, target)
                positional.updateParametersInContext(context)
            folder = output/f'phase-{phase+1}'
            folder.mkdir()
            (folder/'system.xml').write_text(openmm.XmlSerializer.serialize(system))
            initial = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            if phase == 0 and initial != imposed:
                raise ValueError('transport is not exact frozen transplant; retain refusal, do not widen tolerance')
            reporter = reporter_type(openmm)(len(atoms), record['policy']['maximum_callbacks_per_phase'])
            returned = False
            try:
                openmm.LocalEnergyMinimizer.minimize(context, 10., 1000, reporter)
                returned = True
            finally:
                final = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
                compressed(folder/'path.json.gz', dict(initial_nm=initial, final_nm=final, callbacks=reporter.records, returned=returned))
                (folder/'returned-state.xml').write_text(openmm.XmlSerializer.serialize(context.getState(positions=True, parameters=True)))
                phases.append(dict(phase=phase+1, callbacks=len(reporter.records), returned=returned))
                # Each progress record is append-only/exclusive.
                write(output/f'phase-{phase+1}-result.json', dict(case=case, phases=list(phases), selected_atom_indices=indices,
                    plan_sha256=sha(plan_path.read_bytes())))
    finally:
        del context


def run(plan_path, output):
    record = read(plan_path)
    check_plan(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(Path('src').resolve()), str(Path('.').resolve()))))
    results, started = [], time.monotonic()
    for index, case in enumerate(record['cases']):
        folder = output/f'path-{index:02d}'
        folder.mkdir()
        row = dict(index=index, **case, status='budget_not_executed', returncode=None)
        remaining = record['policy']['maximum_total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path),
                '--index', str(index), '--output', str(folder)]
            with (folder/'run.log').open('xb') as log:
                try:
                    process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                        timeout=min(remaining, record['policy']['maximum_seconds_per_path']), check=False)
                    row.update(status='completed' if process.returncode == 0 else 'execution_blocker', returncode=process.returncode)
                except subprocess.TimeoutExpired:
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
    for name in ('parent-plan', 'root', 'angle-audit', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.parent_plan.resolve(), args.root.resolve(), args.angle_audit.resolve(), args.output.resolve())
    elif args.command == 'run':
        run(args.plan.resolve(), args.output.resolve())
    else:
        worker(args.plan.resolve(), args.index, args.output.resolve())
