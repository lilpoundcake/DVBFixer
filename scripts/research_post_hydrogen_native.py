"""Post-H-only guard, then two unchanged native whole phases; no H rerun."""
from __future__ import annotations

import argparse
import math
import os
import signal
import subprocess
import time
from pathlib import Path

from scripts.audit_stereo_force_projection import volume_gradient
from scripts.audit_thr2_native_capture import audit_replay, read, verify_files
from scripts.qualify_hydrogen_postflight import load, reconstruct
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type
from scripts.research_legacy_start_controls import compressed
from scripts.research_post_hydrogen_guard import check as check_static
from scripts.research_thr2_hb_only import make_context, transport
from scripts.research_thr2_hswap import replace_restraint_parameters

ROOT = Path(__file__).resolve().parents[1]


def check(plan):
    from scripts.research_ha_minimizer_path import verify

    verify(plan)
    check_static(read(Path(plan['static_plan'])))


def freeze(output):
    static_path = ROOT/'docs/research/post-hydrogen-guard-static-plan.json'
    evidence_path = ROOT/'docs/research/post-hydrogen-guard-static-evidence.json'
    audit_path = ROOT/'docs/research/post-hydrogen-guard-static-audit.json'
    static, evidence = read(static_path), read(evidence_path)
    check_static(static)
    if evidence['plan_sha256'] != sha(static_path.read_bytes()) or len(evidence['cases']) != len(static['cases']):
        raise ValueError('static binding/enrollment differs')
    controls = read(Path(static['control_plan']))
    inputs = {**static['inputs'], **static['source_hashes']}
    for p in (static_path, evidence_path, audit_path):
        inputs[str(p)] = sha(p.read_bytes())
    root = ROOT/'tmp/pdbfixer-research/post-hydrogen-guard-static-01'
    for name, digest in evidence['artifacts'].items():
        if sha((root/name).read_bytes()) != digest:
            raise ValueError('static artifact changed')
        inputs[str(root/name)] = digest
    cases, reuse, unsupported = [], [], []
    for c, r in zip(static['cases'], evidence['cases'], strict=True):
        summary = dict(static_index=r['index'], label=c['label'], mode=c['mode'], static_status=r['status'])
        if c['kind'] != 'original-thr2':
            if r.get('moved') or r['status'] != 'proposed':
                unsupported.append(dict(summary, reason='no captured matching whole-system Context'))
            continue
        if r['status'] == 'proposed' and r['moved']:
            source, = [s for s in controls['sources'] if s['label'] == c['label']]
            cases.append(dict(summary, source=source, control_folder=c['control_folder'],
                              proposal=str(root/r['proposal_file'])))
        else:
            reuse.append(dict(summary, control_folder=c['control_folder'],
                              role='exact-noop-reuse' if r['status'] == 'proposed' else 'refused-not-executed'))
    if not cases or len(cases) > 6:
        raise ValueError('native plan requires one to six changed original paths')
    contract = ROOT/'docs/plans/post-hydrogen-guard.md'
    inputs[str(contract)] = sha(contract.read_bytes())
    files = [Path(__file__).resolve(), ROOT/'tests/test_post_hydrogen_native_accounting.py']
    write(output, dict(schema='dvbfixer.private-post-h-native-plan.v1', static_plan=str(static_path),
        inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in files},
        cases=cases, reused_noops=reuse, independent_cases_without_whole_context=unsupported,
        python=controls['python'], historical_geometry=controls['historical_geometry'],
        protocol=dict(maximum_callbacks=4096, per_path_seconds=1200, total_seconds=7200,
            new_paths=len(cases), new_minimizations=2*len(cases), extra_h_relaxation=False),
        production_changed=False))


def worker(plan_path, index, output):
    import openmm
    from openmm import Vec3
    from openmm.unit import nanometer

    plan = read(plan_path)
    check(plan)
    case = plan['cases'][index]
    source = Path(case['source']['source'])
    saved = read(Path(case['proposal']))
    rows = read(source/'call-1-before-atoms.json.gz')['atoms']
    topology, _ = reconstruct(rows, read(source/'call-1-topology.json')['bonds'])
    geometry = load('post_h_native_historical_bridge', Path(plan['historical_geometry']))
    positions = [Vec3(*p) for p in saved['proposed_nm']]*nanometer
    repairs = geometry.repair_misplaced_hydrogens(topology, positions)
    requested = [list(p) for p in positions.value_in_unit(nanometer)]
    compressed(output/'bridge.json.gz', dict(before_nm=saved['proposed_nm'], after_nm=requested, repairs=repairs))
    calls = case['source']['calls'][1:]
    system, integrator, context = make_context(source, calls[0], case['mode'])
    context.setPositions(requested*nanometer)
    try:
        for number, call in enumerate(calls, 2):
            if number == 3:
                target = openmm.XmlSerializer.deserialize((source/f"{call['before']['prefix']}-system.xml").read_text())
                restraint = replace_restraint_parameters(system, target)
                restraint.updateParametersInContext(context)
                requested = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            phase = output/f'phase-{number}'
            phase.mkdir()
            raw = openmm.XmlSerializer.serialize(system)
            (phase/'system.xml').write_text(raw)
            if raw != (source/f"{call['before']['prefix']}-system.xml").read_text():
                raise ValueError('phase System bytes changed')
            initial = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            compressed(phase/'transport.json.gz', dict(requested_nm=requested, actual_nm=initial))
            error = transport(requested, initial, case['mode'])
            reporter = reporter_type(openmm)(len(rows), plan['protocol']['maximum_callbacks'])
            returned = False
            try:
                openmm.LocalEnergyMinimizer.minimize(context, call['tolerance'], call['max_iterations'], reporter)
                returned = True
            finally:
                final = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
                compressed(phase/'path.json.gz', dict(initial_nm=initial, final_nm=final,
                                                     callbacks=reporter.records, returned=returned))
                (phase/'returned-state.xml').write_text(openmm.XmlSerializer.serialize(context.getState(positions=True, parameters=True)))
                write(phase/'result.json', dict(returned=returned, callbacks=len(reporter.records), transport_error_nm=error))
    finally:
        del context


def invoke(plan, plan_path, index, folder, timeout):
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(ROOT/'src'), str(ROOT))))
    command = [plan['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path),
               '--index', str(index), '--output', str(folder)]
    with (folder/'run.log').open('xb') as log:
        process = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
            return dict(status='completed' if code == 0 else 'execution-error', returncode=code)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            code = process.wait()
            return dict(status='timeout', returncode=code)


def run(plan_path, root):
    plan = read(plan_path)
    check(plan)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    started, records = time.monotonic(), []
    for index, case in enumerate(plan['cases']):
        folder = root/f'path-{index:02d}'
        folder.mkdir()
        result = dict(index=index, label=case['label'], mode=case['mode'],
                      static_index=case['static_index'], status='not-executed-budget', returncode=None)
        remaining = plan['protocol']['total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            try:
                result.update(invoke(plan, plan_path, index, folder, min(remaining, plan['protocol']['per_path_seconds'])))
            except Exception as error:
                result.update(status='orchestration-error', error_type=type(error).__name__, error=str(error))
        result['artifacts'] = {str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        write(root/f'path-{index:02d}-result.json', result)
        records.append(result)
    # Exclusive writer: each per-path record once above, final index ONCE below.
    write(root/'evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=records))
    check(plan)


def audit(plan_path, root, output):
    import openmm
    from openmm import Vec3
    from openmm.unit import nanometer

    plan, execution = read(plan_path), read(root/'evidence.json')
    check(plan)
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('execution binding mismatch')
    if len(execution['cases']) != len(plan['cases']):
        raise ValueError('native case accounting incomplete')
    geometry = load('post_h_native_audit_historical_bridge', Path(plan['historical_geometry']))
    reports = []
    for index, (case, result) in enumerate(zip(plan['cases'], execution['cases'], strict=True)):
        if any(result[k] != case[k] for k in ('label', 'mode', 'static_index')) or result['index'] != index:
            raise ValueError('case order/identity mismatch')
        folder = root/f'path-{index:02d}'
        verify_files(folder, result)
        if read(root/f'path-{index:02d}-result.json') != result:
            raise ValueError('durable per-path record mismatch')
        report = {k:v for k, v in result.items() if k != 'artifacts'}
        saved = read(Path(case['proposal']))
        source = Path(case['source']['source'])
        rows = read(source/'call-1-before-atoms.json.gz')['atoms']
        bonds = read(source/'call-1-topology.json')['bonds']
        topology, _ = reconstruct(rows, bonds)
        h_by_parent = {}
        oh_pairs = []
        for a, b in bonds:
            if rows[a]['element'] == 'H':
                a, b = b, a
            if rows[b]['element'] == 'H' and rows[a]['element'] != 'H':
                h_by_parent.setdefault(a, []).append(b)
                if rows[a]['element'] == 'O':
                    oh_pairs.append((a, b))
        sibling_pairs = [(a, b) for hs in h_by_parent.values() for i, a in enumerate(hs) for b in hs[i+1:]]
        positions = [Vec3(*p) for p in saved['proposed_nm']]*nanometer
        repairs = geometry.repair_misplaced_hydrogens(topology, positions)
        expected = [list(p) for p in positions.value_in_unit(nanometer)]
        bridge_exact = False
        if (folder/'bridge.json.gz').exists():
            bridge_exact = read(folder/'bridge.json.gz') == dict(before_nm=saved['proposed_nm'], after_nm=expected, repairs=repairs)
            if not bridge_exact:
                raise ValueError('historical bridge mismatch')
        unselected = [i for i in range(len(rows)) if i not in saved['moved_indices']]
        if any(expected[i] != saved['initial_nm'][i] for i in unselected):
            raise ValueError('bridge changed atoms outside post-H mask')
        phases, payloads = [], []
        for number, call in enumerate(case['source']['calls'][1:], 2):
            phase = folder/f'phase-{number}'
            if not (phase/'path.json.gz').exists():
                continue
            payload = read(phase/'path.json.gz')
            t = read(phase/'transport.json.gz')
            if (phase/'system.xml').read_bytes() != (source/f"{call['before']['prefix']}-system.xml").read_bytes():
                raise ValueError('phase physics changed')
            if t['actual_nm'] != payload['initial_nm']:
                raise ValueError('applied transport/path mismatch')
            transport(t['requested_nm'], t['actual_nm'], case['mode'])
            state = openmm.XmlSerializer.deserialize((phase/'returned-state.xml').read_text())
            if state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist() != payload['final_nm']:
                raise ValueError('native returned State/path mismatch')
            if number == 2:
                if t['requested_nm'] != expected:
                    raise ValueError('post-H proposal not used after bridge')
                if any(t['actual_nm'][i] != saved['initial_nm'][i] for i in unselected):
                    raise ValueError('copied post-H coordinates changed')
            measured = audit_replay(plan, source, call, phase,
                dict(label=case['label'], mode=case['mode'], status=result['status'], number=number))
            point_sets = [payload['initial_nm']]+[
                [cb['positions_nm'][i:i+3] for i in range(0, len(cb['positions_nm']), 3)] for cb in payload['callbacks']]+[payload['final_nm']]
            names = {r['key'][-1]:i for i, r in enumerate(rows) if r['key'][:4] == ['A', '2', '', 'THR']}
            for node, points in zip(measured['nodes'], point_sets, strict=True):
                v, _, _ = volume_gradient([points[names[n]] for n in ('HB', 'CA', 'OG1', 'CG2')])
                node['hb_four_ligand_volume_angstrom3'] = 1000*v
                node['hb_parent_distance_angstrom'] = 10*math.dist(points[names['HB']], points[names['CB']])
                node['maximum_oh_distance_angstrom'] = max((10*math.dist(points[a], points[b]) for a, b in oh_pairs), default=None)
                node['minimum_same_parent_h_distance_angstrom'] = min((10*math.dist(points[a], points[b]) for a, b in sibling_pairs), default=None)
            phases.append(measured)
            payloads.append(payload)
        if len(payloads) == 2 and payloads[1]['initial_nm'] != payloads[0]['final_nm']:
            raise ValueError('own-return phase continuity mismatch')
        report.update(phases=phases, bridge=dict(repairs=repairs, exact=bridge_exact))
        reports.append(report)
    write(output, dict(schema='dvbfixer.private-post-h-native-audit.v1',
        plan_sha256=sha(plan_path.read_bytes()), execution_sha256=sha((root/'evidence.json').read_bytes()),
        script_sha256=sha(Path(__file__).read_bytes()), cases=reports,
        reused_noops=plan['reused_noops'], independent_cases_without_whole_context=plan['independent_cases_without_whole_context'],
        new_h_only_calls=0, production_changed=False,
        limits='Mode-matched saved H-only returns; only two original whole phases. Local THR2 result, independent real activation lacks a matching whole-system replay.'))


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
