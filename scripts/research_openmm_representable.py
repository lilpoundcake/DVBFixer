"""Separate representable-force threshold series; original negative toy is frozen."""
from __future__ import annotations

import argparse
import os
import signal
import subprocess
import time
from pathlib import Path

from scripts import research_openmm_cpu_sign as original

ROOT = Path(__file__).resolve().parents[1]


def freeze(output):
    first_path = ROOT/'docs/research/openmm-cpu-sign-reproducer-plan.json'
    evidence_path = ROOT/'docs/research/openmm-cpu-sign-reproducer-evidence.json'
    first, evidence = original.read(first_path), original.read(evidence_path)
    original.check(first)
    if evidence['synthetic_primary_signal_observed'] or evidence['all_other_controls_passed']:
        raise ValueError('original negative outcomes must be preserved')
    kernel = ROOT/'tmp/pdbfixer-research/minimize-kernel-36a30cb.cc'
    if 'const mm_long limit = 0x4000000000000000' not in kernel.read_text():
        raise ValueError('pinned source detector threshold differs')
    plan = dict(first)
    plan['schema'] = 'dvbfixer.private-representable-force-plan.v1'
    plan['inputs'] = dict(first['inputs'])
    for p in (first_path, evidence_path, kernel, ROOT/'docs/plans/openmm-representable-force-control.md'):
        plan['inputs'][str(p)] = original.sha(p)
    plan['source_hashes'] = dict(first['source_hashes'])
    for p in (Path(__file__).resolve(), ROOT/'tests/test_openmm_representable_control.py'):
        plan['source_hashes'][str(p)] = original.sha(p)
    plan['policy'] = dict(first['policy'], max_iterations=1)
    plan['cases'] = []
    for old in first['cases']:
        trigger = old['stress'] > 1000.
        case = dict(old, stress=5e7 if trigger else 1e7, tier='trigger' if trigger else 'below',
            stress_start=[10., 10., 0.])
        case['label'] = f"{case['backend']}-"+('constrained' if case['constrained'] else 'free')+'-'+case['tier']
        plan['cases'].append(case)
    plan['first_negative_series_rehabilitated'] = False
    original.write(output, plan)


def worker(plan_path, index, output):
    import numpy as np
    import openmm
    from openmm.unit import nanometer

    plan = original.read(plan_path)
    original.check(plan)
    case = plan['cases'][index]
    system = original.make_system(case['stress'], case['constrained'])
    integrator = openmm.VerletIntegrator(0.001)
    integrator.setConstraintTolerance(plan['policy']['integrator_constraint_tolerance'])
    platform = openmm.Platform.getPlatformByName(case['backend'])
    context = openmm.Context(system, integrator, platform, case['properties'])
    positions = [[-0.0500025, 0., 0.], [0.0500025, 0., 0.], case['stress_start']]
    context.setPositions(positions*nanometer)

    def capture(ctx, prefix):
        state = ctx.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
        for label, obj in (('system', system), ('integrator', integrator), ('state', state)):
            with (output/f'{prefix}-{label}.xml').open('x') as stream:
                stream.write(openmm.XmlSerializer.serialize(obj))
        with (output/f'{prefix}.chk').open('xb') as stream:
            stream.write(ctx.createCheckpoint())
        return state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()

    class Recorder(openmm.MinimizationReporter):
        def __init__(self):
            super().__init__()
            self.records = []

        def report(self, iteration, x, grad, args):
            if len(self.records) >= plan['policy']['maximum_callbacks']:
                raise RuntimeError('frozen callback budget exceeded')
            p, g, stats = list(x), list(grad), dict(args)
            if len(p) != 9 or len(g) != 9 or not all(np.isfinite(v) for v in p+g+list(stats.values())):
                raise ValueError('invalid native callback')
            self.records.append(dict(ordinal=len(self.records), iteration=iteration,
                positions_nm=p, gradient=g, statistics=stats))
            return False

    recorder, returned = Recorder(), False
    initial = capture(context, 'before')
    original.write(output/'platform.json', dict(name=platform.getName(), properties={
        name:platform.getPropertyValue(context, name) for name in platform.getPropertyNames()}))
    try:
        openmm.LocalEnergyMinimizer.minimize(context, plan['policy']['tolerance'], plan['policy']['max_iterations'], recorder)
        returned = True
    finally:
        original.write(output/'path.json.gz', dict(requested_nm=positions, initial_nm=initial,
            final_nm=capture(context, 'after'), callbacks=recorder.records, returned=returned))
        del context


def run(plan_path, root):
    plan = original.read(plan_path)
    original.check(plan)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    rows, start = [], time.monotonic()
    for case in plan['cases']:
        folder = root/f"case-{case['index']:02d}"
        folder.mkdir()
        for name in ('cache', 'temp'):
            (folder/name).mkdir()
        row = dict(index=case['index'], label=case['label'], status='not-executed-budget', returncode=None)
        remaining = plan['policy']['total_seconds']-(time.monotonic()-start)
        if remaining > 0:
            env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONPATH=str(ROOT),
                OPENMM_CACHE_DIR=str(folder/'cache'), XDG_CACHE_HOME=str(folder/'cache'),
                TMPDIR=str(folder/'temp'), TMP=str(folder/'temp'), TEMP=str(folder/'temp'))
            command = [plan['runtime']['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path),
                '--index', str(case['index']), '--output', str(folder)]
            try:
                with (folder/'run.log').open('xb') as log:
                    process = subprocess.Popen(command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                    try:
                        code = process.wait(timeout=min(remaining, plan['policy']['per_case_seconds']))
                        row.update(status='completed' if code == 0 else 'execution-error', returncode=code)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        row.update(status='timeout', returncode=process.wait())
            except Exception as error:
                row.update(status='orchestration-error', error_type=type(error).__name__, error=str(error))
        row['artifacts'] = {str(p.relative_to(folder)):original.sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}
        original.write(root/f"case-{case['index']:02d}-result.json", row)
        rows.append(row)
    original.check(plan)
    original.write(root/'execution.json', dict(plan_sha256=original.sha(plan_path), cases=rows))


def audit(plan_path, root, output):
    # Reuse ONLY unchanged numerical/FD comparisons, not the first series' gate.
    numerical_path = root/'original-envelope-readback.json'
    original.audit(plan_path, root, numerical_path)
    numerical = original.read(numerical_path)
    rows = numerical['cases']
    for row in rows:
        row.pop('primary_signal_nodes', None)
        case = row['case']
        row['representable_trigger_nodes'] = []
        raw_path = root/f"case-{case['index']:02d}/path.json.gz"
        if not raw_path.exists():
            continue
        raw = original.read(raw_path)
        for node, callback in zip(row['nodes'], raw['callbacks'], strict=True):
            largest = max(abs(v) for v in callback['gradient'])
            node['maximum_reported_gradient'] = largest
            node['above_detector_below_capacity'] = 2**30 < largest < 2**31
            if (case['backend'] == 'OpenCL' and case['constrained'] and case['tier'] == 'trigger'
                    and node['above_detector_below_capacity'] and node['penalty_gradient_maximum'] > 0
                    and not node['hypotheses']['physical']['allowance_passed']
                    and (node['hypotheses']['source_sign']['strict_passed'] or node['hypotheses']['source_sign_float32']['strict_passed'])
                    and node['finite_difference_passed'] and node['energy_passed'] and node['penalty_energy_passed']):
                row['representable_trigger_nodes'].append(node['ordinal'])
    controls = [r for r in rows if not (r['case']['backend'] == 'OpenCL'
        and r['case']['constrained'] and r['case']['tier'] == 'trigger')]
    original.write(output, dict(schema='dvbfixer.private-representable-force-evidence.v1',
        plan_sha256=numerical['plan_sha256'], execution_sha256=numerical['execution_sha256'], cases=rows,
        representable_primary_signal_observed=any(r.get('representable_trigger_nodes') for r in rows),
        all_eleven_controls_passed=len(controls) == 11 and all(r.get('physical_control_passed', False) for r in controls),
        numerical_readback_sha256=original.sha(numerical_path), first_series_outcome_changed=False,
        installed_dependency_modified=False, internal_fallback_event_directly_observed=False,
        installed_binary_source_identity_attested=False, structural_qualification_rehabilitated=False,
        limits='Separate prospective representable-force threshold control; previous oversized-force failures remain unchanged. Source-consistent public-reporter inconsistency is not correct native physics, direct internal-event evidence or a patched dependency.'))


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
