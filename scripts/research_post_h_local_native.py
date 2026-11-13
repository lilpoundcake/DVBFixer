"""Three frozen retained-state paths under one captured 1VCU physics model."""
from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from scripts.audit_post_h_1vcu import audit_whole
from scripts.audit_thr2_native_capture import read, verify_files
from scripts.qualify_glycam_fallback import restore
from scripts.qualify_hydrogen_initialization_guard import atoms_from_rows
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_post_h_1vcu import capture
from scripts.research_post_h_local_noop import check as check_static
from scripts.research_thr2_hb_only import transport
from scripts.research_thr2_hswap import replace_restraint_parameters
from scripts.research_thr2_native_capture import update_owned

from dvbfixer.ffutils.variants import VARIANT_TO_PARENT

ROOT = Path(__file__).resolve().parents[1]


def check(plan):
    verify(plan)
    check_static(read(Path(plan['static_plan'])))


def chemistry(value):
    return [(a.identity, a.element, VARIANT_TO_PARENT.get(a.residue_name, a.residue_name)) for a in value]


def bind(saved, snap):
    """No guessed atom mapping, no canonical-name loss in raw stored artifacts."""
    from scripts.research_post_h_after_repair import input_from_table

    value = input_from_table(saved['input'])
    if chemistry(value.atoms) != chemistry(atoms_from_rows(snap['atoms'])):
        raise ValueError('shared physics exact canonical atom chemistry mismatch')
    if {tuple(sorted(b)) for b in value.bonds} != {tuple(sorted(b)) for b in snap['bonds']}:
        raise ValueError('shared physics bond graph mismatch')
    before = [list(p) for a, p in zip(value.atoms, value.positions_nm, strict=True) if a.element not in ('H', 'D')]
    after = [p for a, p in zip(value.atoms, snap['positions_nm'], strict=True) if a.element not in ('H', 'D')]
    return transport(before, after, 'captured')


def freeze(output):
    import openmm

    static_path = ROOT/'docs/research/post-hydrogen-local-noop-plan.json'
    evidence_path = ROOT/'docs/research/post-hydrogen-local-noop-evidence.json'
    static, evidence = read(static_path), read(evidence_path)
    check_static(static)
    if not evidence['both_states_proposed_and_audited'] or evidence['plan_sha256'] != sha(static_path.read_bytes()):
        raise ValueError('both-state static admission required')
    source = ROOT/'tmp/pdbfixer-research/post-hydrogen-1vcu-public-api-01/case-00'
    original = read(ROOT/'docs/research/post-hydrogen-1vcu-plan.json')
    calls = read(source/'capture.json')['whole_calls']
    if len(calls) != 2 or not all(c['returned'] for c in calls):
        raise ValueError('two retained full-system calls required')
    inputs = {str(p):sha(p.read_bytes()) for p in (static_path, evidence_path,
        ROOT/'docs/plans/post-hydrogen-local-noop-native.md', source/'capture.json')}
    for call in calls:
        prefix = call['before']['prefix']
        for suffix in ('system.xml', 'integrator.xml', 'state.xml', 'atoms.json.gz'):
            p = source/f'{prefix}-{suffix}'
            inputs[str(p)] = sha(p.read_bytes())
    snap = read(source/'whole-1-before-atoms.json.gz')
    local = ROOT/'tmp/pdbfixer-research/post-hydrogen-local-noop-01'
    findings = []
    for i in (0, 1):
        p = local/f'case-{i:02d}-proposal.json.gz'
        digest = sha(p.read_bytes())
        if digest != evidence['private_artifacts'][p.name]:
            raise ValueError('static proposal digest mismatch')
        inputs[str(p)] = digest
        saved = read(p)
        findings.append(bind(saved, snap))
        journal_path = Path(static['cases'][i]['journal'])
        journal = read(journal_path)
        inputs[str(journal_path)] = sha(journal_path.read_bytes())
        raw_meta = {tuple(sorted((a, b))):(t, o) for a, b, t, o in journal['raw_bond_metadata']}
        whole_meta = {tuple(sorted(b)):(m['type'], m['order'])
                      for b, m in zip(snap['bonds'], snap['bond_metadata'], strict=True)}
        if raw_meta != whole_meta:
            raise ValueError('shared physics bond metadata mismatch')
        if i == 0:
            variants = journal['returned']
        elif variants != journal['returned']:
            raise ValueError('raw native selected variants differ')
    system = openmm.XmlSerializer.deserialize((source/'whole-1-before-system.xml').read_text())
    for f in system.getForces():
        if isinstance(f, openmm.CustomExternalForce):
            names = [f.getPerParticleParameterName(i) for i in range(f.getNumPerParticleParameters())]
            if names == ['k', 'x0', 'y0', 'z0']:
                for i in range(f.getNumParticles()):
                    atom, _ = f.getParticleParameters(i)
                    if snap['atoms'][atom]['element'] in ('H', 'D'):
                        raise ValueError('shared anchors include H')
    target = openmm.XmlSerializer.deserialize((source/'whole-2-before-system.xml').read_text())
    replace_restraint_parameters(system, target)
    cases = [dict(label=label, state=i, projected=projected,
        proposal=str(local/f'case-{i:02d}-proposal.json.gz')) for label, i, projected in
        [('baseline-projection', 0, True), ('guarded-unmodified', 1, False), ('guarded-projection', 1, True)]]
    sources = [ROOT/'scripts'/n for n in ('research_post_h_local_native.py', 'research_thr2_hswap.py')]
    sources.append(ROOT/'tests/test_post_h_local_native.py')
    write(output, dict(schema='dvbfixer.private-post-h-local-native-plan.v1',
        static_plan=str(static_path), inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in sources},
        original_policy=original['policy'], source=str(source), calls=calls, cases=cases,
        python=sys.executable, shared_anchor_heavy_rounding_nm=findings,
        policy=dict(per_path_seconds=1800, total_seconds=5400, maximum_callbacks=4096,
            maximum_minimizations=6, new_h_samples=0, new_repairs=0, native_retries=0),
        existing_baseline_control_reused=True, independently_constructed_guarded_forcefield=False,
        production_changed=False))


def worker(plan_path, index, output):
    import openmm
    from openmm.unit import nanometer

    plan = read(plan_path)
    check(plan)
    case, source = plan['cases'][index], Path(plan['source'])
    saved = read(Path(case['proposal']))
    snap = read(source/'whole-1-before-atoms.json.gz')
    bind(saved, snap)
    topology, _ = restore(snap)
    requested = saved['proposal']['result']['positions_nm'] if case['projected'] else saved['input']['positions_nm']
    prefix = plan['calls'][0]['before']['prefix']
    system = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((source/f'{prefix}-integrator.xml').read_text())
    platform = plan['calls'][0]['before']['platform']
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(platform['name']), platform['constructor_properties'])
    calls = []
    try:
        context.setState(openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()))
        context.setPositions(requested*nanometer)
        for call in plan['calls']:
            number = call['number']
            if number == 2:
                target = openmm.XmlSerializer.deserialize((source/f"{call['before']['prefix']}-system.xml").read_text())
                force = replace_restraint_parameters(system, target)
                force.updateParametersInContext(context)
                requested = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            if openmm.XmlSerializer.serialize(system).encode() != (source/f"{call['before']['prefix']}-system.xml").read_bytes():
                raise ValueError('phase System is not exact captured physics')
            row = dict(number=number, returned=False, tolerance=call['tolerance'], max_iterations=call['max_iterations'])
            calls.append(row)
            compressed(output/f'whole-{number}-transport.json.gz', dict(requested_nm=requested))
            row['before'] = capture(context, topology, output, f'whole-{number}-before')
            initial = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            row['transport_error_nm'] = transport(requested, initial, 'captured')
            reporter = reporter_type(openmm)(len(initial), plan['policy']['maximum_callbacks'])
            update_owned(output/'capture.json', dict(calls=calls))
            try:
                openmm.LocalEnergyMinimizer.minimize(context, call['tolerance'], call['max_iterations'], reporter)
                row['returned'] = True
            finally:
                final = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
                compressed(output/f'whole-{number}-path.json.gz', dict(initial_nm=initial, final_nm=final,
                    callbacks=reporter.records, returned=row['returned']))
                row['after'] = capture(context, topology, output, f'whole-{number}-after')
                row['callbacks'] = len(reporter.records)
                update_owned(output/'capture.json', dict(calls=calls))
    finally:
        del context


def run(plan_path, root):
    plan = read(plan_path)
    check(plan)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    rows, start = [], time.monotonic()
    for i, case in enumerate(plan['cases']):
        folder = root/f'path-{i:02d}'
        folder.mkdir()
        row = dict(index=i, label=case['label'], status='not-executed-budget', returncode=None)
        remaining = plan['policy']['total_seconds']-(time.monotonic()-start)
        if remaining > 0:
            env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(ROOT/'src'), str(ROOT))),
                PATH=str(Path(plan['python']).parent)+os.pathsep+os.environ.get('PATH', ''))
            command = [plan['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path),
                '--index', str(i), '--output', str(folder)]
            try:
                with (folder/'run.log').open('xb') as log:
                    process = subprocess.Popen(command, cwd=folder, env=env, stdout=log,
                        stderr=subprocess.STDOUT, start_new_session=True)
                    try:
                        code = process.wait(timeout=min(remaining, plan['policy']['per_path_seconds']))
                        row.update(status='completed' if code == 0 else 'execution-error', returncode=code)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        row.update(status='timeout', returncode=process.wait())
            except Exception as error:
                row.update(status='orchestration-error', error_type=type(error).__name__, error=str(error))
        row['artifacts'] = {str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        write(root/f'path-{i:02d}-result.json', row)
        rows.append(row)
        update_owned(root/'progress.json', dict(cases=rows))
    write(root/'evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=rows))
    check(plan)


def centers(atoms, points):
    import numpy as np

    p = np.asarray(points)
    values = []
    for resid in ('0', '283', '286', '379'):
        names = {a['key'][4]:i for i, a in enumerate(atoms)
                 if a['key'][0:3] == ['A', resid, '']}
        if not {'N', 'CA', 'C', 'CB', 'HA'} <= names.keys():
            raise ValueError('target center missing')
        heavy = p[[names[n] for n in ('N', 'C', 'CB')]]
        hv, fv = heavy-p[names['CA']], heavy-p[names['HA']]
        hvol = float(np.linalg.det(hv))/float(np.prod(np.linalg.norm(hv, axis=1)))
        fvol = float(np.linalg.det(fv))/float(np.prod(np.linalg.norm(fv, axis=1)))
        values.append(dict(resid=resid, heavy_normalized_volume=hvol, four_normalized_volume=fvol,
            sign_conflict=hvol*fvol < 0, degenerate=abs(hvol) <= 1e-4 or abs(fvol) <= 1e-4,
            radius_angstrom=10*float(np.linalg.norm(p[names['HA']]-p[names['CA']]))))
    return values


def audit(plan_path, root, output):
    plan, execution = read(plan_path), read(root/'evidence.json')
    check(plan)
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('native execution binding differs')
    if len(execution['cases']) != 3:
        raise ValueError('native enrollment accounting differs')
    rows = []
    for index, (case, result) in enumerate(zip(plan['cases'], execution['cases'], strict=True)):
        folder = root/f'path-{index:02d}'
        verify_files(folder, result)
        if result['index'] != index or result['label'] != case['label']:
            raise ValueError('native case identity differs')
        row = {k:v for k, v in result.items() if k != 'artifacts'}
        row['phases'] = []
        try:
            saved = read(Path(case['proposal']))
            expected = saved['proposal']['result']['positions_nm'] if case['projected'] else saved['input']['positions_nm']
            calls = read(folder/'capture.json')['calls']
            for call in calls:
                n = call['number']
                if not call.get('after'):
                    continue
                payload = read(folder/f'whole-{n}-path.json.gz')
                requested = read(folder/f'whole-{n}-transport.json.gz')['requested_nm']
                if requested != expected:
                    raise ValueError('proposal/own-return handoff differs')
                transport(requested, payload['initial_nm'], 'captured')
                expected = payload['final_nm']
                if (folder/f'whole-{n}-before-system.xml').read_bytes() != (Path(plan['source'])/f'whole-{n}-before-system.xml').read_bytes():
                    raise ValueError('native System differs')
                measured = audit_whole(dict(policy=plan['original_policy']), folder, call)
                snap = read(folder/f'whole-{n}-before-atoms.json.gz')
                sets = [payload['initial_nm']]+[
                    [c['positions_nm'][i:i+3] for i in range(0, len(c['positions_nm']), 3)]
                    for c in payload['callbacks']]+[payload['final_nm']]
                for node, points in zip(measured['nodes'], sets, strict=True):
                    node['local_centers'] = centers(snap['atoms'], points)
                measured['all_moved_centers_compatible_nondegenerate'] = all(
                    not c['sign_conflict'] and not c['degenerate'] for node in measured['nodes']
                    for c in node['local_centers'] if c['resid'] in ('0', '283', '286'))
                measured['negative_alpha_nodes'] = sum(bool(node['negative_alpha']) for node in measured['nodes'])
                row['phases'].append(measured)
            row['geometry_passed'] = (case['projected'] and len(row['phases']) == 2
                and all(p['returned'] and p['all_moved_centers_compatible_nondegenerate']
                    and not p['negative_alpha_nodes'] for p in row['phases']))
        except Exception as error:
            row.update(audit_error_type=type(error).__name__, audit_error=str(error), geometry_passed=False)
        rows.append(row)
    write(output, dict(schema='dvbfixer.private-post-h-local-native-audit.v1', cases=rows,
        plan_sha256=sha(plan_path.read_bytes()), execution_sha256=sha((root/'evidence.json').read_bytes()),
        production_ready=False, existing_baseline_control_reused=True,
        limits='Shared captured baseline physics, not an independently built guarded FF or live integration. Failed original gradient envelopes and old strict refusals are not erased.'))


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
