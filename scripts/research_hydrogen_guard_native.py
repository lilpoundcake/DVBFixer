"""Two prospectively frozen actual-proposal trajectories; no production hook."""
from __future__ import annotations

import argparse
import os
import signal
import subprocess
import time
from pathlib import Path

from scripts.audit_thr2_native_capture import audit_replay, read, verify_files
from scripts.qualify_hydrogen_postflight import load, reconstruct
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_thr2_hb_only import check_plan, make_context, transport
from scripts.research_thr2_hswap import replace_restraint_parameters

ROOT = Path(__file__).resolve().parents[1]


def check(plan):
    verify(plan)
    check_plan(read(Path(plan['control_plan'])))


def freeze(output):
    static_path = ROOT/'docs/research/hydrogen-initialization-guard-static-plan.json'
    static_evidence = ROOT/'docs/research/hydrogen-initialization-guard-static-evidence.json'
    static_audit = ROOT/'docs/research/hydrogen-initialization-guard-static-audit.json'
    control_path = ROOT/'docs/research/thr2-hb-only-plan.json'
    control_evidence = ROOT/'tmp/pdbfixer-research/thr2-hb-only-run-01/evidence.json'
    static, controls = read(static_path), read(control_path)
    verify(static)
    check_plan(controls)
    evidence = read(static_evidence)
    if evidence['plan_sha256'] != sha(static_path.read_bytes()):
        raise ValueError('static binding differs')
    changed = [(c, s) for c in evidence['cases'] for s in c['results'] if s.get('moved')]
    if len(changed) != 1 or changed[0][0]['label'] != 'unseeded' or len(changed[0][1]['moved']) != 1:
        raise ValueError('Expected sole original unseeded-HB activation')
    proposal = ROOT/'tmp/pdbfixer-research/hydrogen-initializer-guard-static-01'/changed[0][1]['proposal_file']
    if sha(proposal.read_bytes()) != evidence['artifacts'][proposal.name]:
        raise ValueError('proposal bytes changed')
    inputs = {**static['inputs'], **static['source_hashes'], **controls['inputs'], **controls['source_hashes']}
    paths = [static_path, static_evidence, static_audit, control_path, control_evidence, proposal,
             ROOT/'docs/plans/hydrogen-initializer-guard-native.md',
             ROOT/'docs/research/thr2-hb-only-evidence.json']
    for p in paths:
        inputs[str(p)] = sha(p.read_bytes())
    reuse = []
    for index, c in enumerate(read(control_evidence)['cases']):
        if c['arm'] != 'unchanged':
            continue
        folder = control_evidence.parent/f'path-{index:02d}'
        verify_files(folder, c)
        inputs.update({str(folder/name):digest for name, digest in c['artifacts'].items()})
        reuse.append(dict(label=c['label'], mode=c['mode'], folder=str(folder),
            role='unchanged-comparator' if c['label'] == 'unseeded' else 'exact-noop-proposal-reuse'))
    source = next(s for s in controls['sources'] if s['label'] == 'unseeded')
    write(output, dict(schema='dvbfixer.private-initializer-native-plan.v1', inputs=inputs,
        source_hashes={str(Path(__file__).resolve()):sha(Path(__file__).read_bytes())},
        control_plan=str(control_path), proposal=str(proposal), source=source, controls=reuse,
        historical_geometry=controls['historical_geometry'], python=controls['python'],
        cases=[dict(label='unseeded', mode=mode) for mode in ('captured', 'Reference')],
        protocol=dict(maximum_callbacks=4096, per_path_seconds=1200, total_seconds=2400,
                      new_paths=2, new_minimizations=6), production_changed=False))


def worker(plan_path, index, output):
    import openmm
    from openmm.unit import nanometer

    plan = read(plan_path)
    check(plan)
    case = plan['cases'][index]
    source = Path(plan['source']['source'])
    calls = plan['source']['calls']
    saved = read(Path(plan['proposal']))
    requested = saved['proposed_nm']
    rows = read(source/'call-1-before-atoms.json.gz')['atoms']
    bonds = read(source/'call-1-topology.json')['bonds']
    historical = load('private_guard_native_old_bridge', Path(plan['historical_geometry']))
    system, integrator, context = make_context(source, calls[0], case['mode'])
    original = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    if original != saved['initial_nm']:
        raise ValueError('Context does not match static proposal source')
    context.setPositions(requested*nanometer)
    try:
        for number, call in enumerate(calls, 1):
            if number == 2:
                positions = context.getState(positions=True).getPositions()
                before = [list(p) for p in positions.value_in_unit(nanometer)]
                topology, _ = reconstruct(rows, bonds)
                repairs = historical.repair_misplaced_hydrogens(topology, positions)
                requested = [list(p) for p in positions.value_in_unit(nanometer)]
                compressed(output/'bridge.json.gz', dict(before_nm=before, after_nm=requested, repairs=repairs))
                context = None
                system, integrator, context = make_context(source, call, case['mode'])
                context.setPositions(requested*nanometer)
            elif number == 3:
                target = openmm.XmlSerializer.deserialize((source/f"{call['before']['prefix']}-system.xml").read_text())
                restraint = replace_restraint_parameters(system, target)
                restraint.updateParametersInContext(context)
                requested = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            phase = output/f'phase-{number}'
            phase.mkdir()
            raw = openmm.XmlSerializer.serialize(system)
            (phase/'system.xml').write_text(raw)
            if raw != (source/f"{call['before']['prefix']}-system.xml").read_text():
                raise ValueError('native System bytes changed')
            initial = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            compressed(phase/'transport.json.gz', dict(requested_nm=requested, actual_nm=initial))
            error = transport(requested, initial, case['mode'])
            if number == 1 and any(a != b for i, (a, b) in enumerate(zip(original, initial, strict=True)) if i not in saved['moved_indices']):
                raise ValueError('unselected initial atoms moved')
            reporter = reporter_type(openmm)(len(rows), plan['protocol']['maximum_callbacks'])
            returned = False
            try:
                openmm.LocalEnergyMinimizer.minimize(context, call['tolerance'], call['max_iterations'], reporter)
                returned = True
            finally:
                final = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
                compressed(phase/'path.json.gz', dict(initial_nm=initial, final_nm=final, callbacks=reporter.records, returned=returned))
                (phase/'returned-state.xml').write_text(openmm.XmlSerializer.serialize(context.getState(positions=True, parameters=True)))
                write(phase/'result.json', dict(returned=returned, callbacks=len(reporter.records), transport_error_nm=error))
    finally:
        if context is not None:
            del context


def run(plan_path, output):
    plan = read(plan_path)
    check(plan)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(ROOT/'src'), str(ROOT))))
    started, cases = time.monotonic(), []
    for index, case in enumerate(plan['cases']):
        folder = output/f'path-{index:02d}'
        folder.mkdir()
        result = dict(index=index, **case, status='not-executed-budget')
        remaining = plan['protocol']['total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            command = [plan['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path), '--output', str(folder), '--index', str(index)]
            with (folder/'run.log').open('xb') as log:
                process = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=min(remaining, plan['protocol']['per_path_seconds']))
                    result.update(status='completed' if code == 0 else 'execution-error', returncode=code)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    result['status'] = 'timeout'
        result['artifacts'] = {str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        cases.append(result)
        write(output/'evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=cases))
    check(plan)


def audit(plan_path, root, output):
    from openmm import Vec3
    from openmm.unit import nanometer

    plan, execution = read(plan_path), read(root/'evidence.json')
    check(plan)
    if (execution['plan_sha256'] != sha(plan_path.read_bytes())
            or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes()
            or [{k:r[k] for k in ('label', 'mode')} for r in execution['cases']] != plan['cases']):
        raise ValueError('native execution enrollment changed')
    saved = read(Path(plan['proposal']))
    source = Path(plan['source']['source'])
    historical = load('private_guard_native_audit_bridge', Path(plan['historical_geometry']))
    rows = read(source/'call-1-before-atoms.json.gz')['atoms']
    heavy = [i for i, r in enumerate(rows) if r['element'] != 'H']
    reports = []
    for result in execution['cases']:
        folder = root/f"path-{result['index']:02d}"
        verify_files(folder, result)
        phases, payloads = [], []
        for number, call in enumerate(plan['source']['calls'], 1):
            phase = folder/f'phase-{number}'
            if not (phase/'path.json.gz').exists():
                continue
            payload = read(phase/'path.json.gz')
            if (phase/'system.xml').read_bytes() != (source/f"{call['before']['prefix']}-system.xml").read_bytes():
                raise ValueError('native phase physics changed')
            t = read(phase/'transport.json.gz')
            if t['actual_nm'] != payload['initial_nm']:
                raise ValueError('transport/path mismatch')
            transport(t['requested_nm'], t['actual_nm'], result['mode'])
            if number == 1:
                if t['requested_nm'] != saved['proposed_nm']:
                    raise ValueError('not the actual static guard proposal')
                if any(p != saved['initial_nm'][i] for i, p in enumerate(t['actual_nm']) if i not in saved['moved_indices']):
                    raise ValueError('unselected initial atom moved')
                points = [payload['initial_nm'], payload['final_nm']]+[
                    [cb['positions_nm'][i:i+3] for i in range(0, len(cb['positions_nm']), 3)] for cb in payload['callbacks']]
                if any(p[i] != saved['initial_nm'][i] for p in points for i in heavy):
                    raise ValueError('pinned heavy moved')
            measured = audit_replay(plan, source, call, phase,
                dict(label=result['label'], mode=result['mode'], status=result['status'], number=number))
            phases.append(measured)
            payloads.append(payload)
        bridge = None
        if len(payloads) >= 2:
            topology, _ = reconstruct(rows, read(source/'call-1-topology.json')['bonds'])
            positions = [Vec3(*p) for p in payloads[0]['final_nm']]*nanometer
            repairs = historical.repair_misplaced_hydrogens(topology, positions)
            expected = [list(p) for p in positions.value_in_unit(nanometer)]
            b = read(folder/'bridge.json.gz')
            if b != dict(before_nm=payloads[0]['final_nm'], after_nm=expected, repairs=repairs):
                raise ValueError('historical bridge mismatch')
            if read(folder/'phase-2/transport.json.gz')['requested_nm'] != expected:
                raise ValueError('phase1 own-return continuity mismatch')
            if any(expected[i] != payloads[0]['final_nm'][i] for i in heavy):
                raise ValueError('bridge moved heavy atoms')
            bridge = dict(repairs=repairs, exact=True)
        if len(payloads) == 3 and payloads[1]['final_nm'] != payloads[2]['initial_nm']:
            raise ValueError('phase2 continuity mismatch')
        reports.append(dict(index=result['index'], label=result['label'], mode=result['mode'],
            status=result['status'], phases=phases, bridge=bridge))
    write(output, dict(schema='dvbfixer.private-initializer-native-audit.v1',
        plan_sha256=sha(plan_path.read_bytes()), execution_sha256=sha((root/'evidence.json').read_bytes()),
        script_sha256=sha(Path(__file__).read_bytes()), cases=reports, reused_controls=plan['controls'],
        production_changed=False, new_h_sampling=False,
        limits='Two changed-proposal paths in one original context; reused controls/noops are not additional sampling or independent activation qualification.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'worker', 'run', 'audit'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--root', type=Path)
    parser.add_argument('--index', type=int)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    elif args.command == 'worker':
        worker(args.plan.resolve(), args.index, args.output.resolve())
    elif args.command == 'run':
        run(args.plan.resolve(), args.output.resolve())
    else:
        audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
