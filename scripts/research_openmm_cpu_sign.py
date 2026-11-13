"""Isolated three-particle reproducer; no installed OpenMM modification.

Only stdlib, NumPy and public OpenMM APIs. Raw arrays stay in the private run.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    path = Path(path)
    raw = gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()
    return json.loads(raw)


def write(path, value):
    raw = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()
    path = Path(path)
    with path.open('xb') as stream:
        stream.write(gzip.compress(raw, mtime=0) if path.suffix == '.gz' else raw)


def runtime():
    import openmm._openmm
    import openmm.version

    return dict(version=openmm.version.version, revision=openmm.version.git_revision,
        extension_path=openmm._openmm.__file__, python=sys.executable, machine=platform.machine())


def check(plan):
    for name, digest in {**plan['inputs'], **plan['source_hashes']}.items():
        if sha(name) != digest:
            raise ValueError(f'frozen input/source changed: {name}')
    if runtime() != plan['runtime']:
        raise ValueError('runtime metadata changed')


def energy_terms(points, stress, constrained, strength):
    import numpy as np

    p = np.asarray(points, dtype=float)
    r = float(np.linalg.norm(p[1]-p[0]))
    return np.array([0.5*10000.*(r-0.12)**2,
        0.5*stress*float(np.dot(p[2]**2, [1., 3., 1.])),
        0.5*strength*(r-0.1)**2 if constrained else 0.])


def analytic(points, stress, constrained, strength):
    import numpy as np

    p = np.asarray(points, dtype=float)
    delta = p[1]-p[0]
    r = float(np.linalg.norm(delta))
    if r <= 1e-15:
        raise ValueError('degenerate synthetic pair')
    physical, penalty = np.zeros((3, 3)), np.zeros((3, 3))
    bond = 10000.*(r-0.12)*delta/r
    physical[0], physical[1] = -bond, bond
    physical[2] = stress*p[2]*[1., 3., 1.]
    if constrained:
        g = strength*(r-0.1)*delta/r
        penalty[0], penalty[1] = -g, g
    return physical, penalty


def finite_difference(points, stress, constrained, strength):
    import numpy as np

    p = np.asarray(points, dtype=float)
    gradient = np.zeros((3, 3))
    for atom in range(3):
        for axis in range(3):
            left, right = p.copy(), p.copy()
            left[atom, axis] -= 1e-6
            right[atom, axis] += 1e-6
            gradient[atom, axis] = float(np.sum((energy_terms(right, stress, constrained, strength)
                -energy_terms(left, stress, constrained, strength))/(2e-6)))
    return gradient


def make_system(stress, constrained):
    import openmm

    system = openmm.System()
    for _ in range(3):
        system.addParticle(12.)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, 0.12, 10000.)
    system.addForce(bond)
    external = openmm.CustomExternalForce('0.5*ks*(x*x+3*y*y+z*z)')
    external.addGlobalParameter('ks', stress)
    external.addParticle(2, [])
    system.addForce(external)
    if constrained:
        system.addConstraint(0, 1, 0.1)
    return system


def freeze(output):
    import openmm

    files = [Path(__file__).resolve(), ROOT/'tests/test_openmm_cpu_sign_reproducer.py']
    inputs = [ROOT/'docs/plans/openmm-cpu-sign-reproducer.md',
        ROOT/'tmp/pdbfixer-research/common-minimize-36a30cb.cpp',
        ROOT/'docs/research/post-hydrogen-cpu-sign-evidence.json',
        ROOT/'docs/research/post-hydrogen-local-noop-validation.json']
    cases = []
    for backend in ('OpenCL', 'CPU', 'Reference'):
        for constrained in (False, True):
            for stress in (1000., 1e14):
                cases.append(dict(index=len(cases), backend=backend, constrained=constrained, stress=stress,
                    label=f'{backend}-'+('constrained' if constrained else 'free')+('-high' if stress > 1000 else '-low'),
                    properties=dict(DeviceIndex='0', OpenCLPlatformIndex='0', Precision='single',
                        UseCpuPme='false', DisablePmeStream='false') if backend == 'OpenCL' else {}))
    write(output, dict(schema='dvbfixer.private-openmm-sign-reproducer-plan.v1',
        inputs={str(p):sha(p) for p in inputs}, source_hashes={str(p):sha(p) for p in files},
        runtime=runtime(), available_platforms=[openmm.Platform.getPlatform(i).getName()
            for i in range(openmm.Platform.getNumPlatforms())], cases=cases,
        policy=dict(per_case_seconds=30, total_seconds=360, maximum_callbacks=256,
            max_iterations=8, tolerance=10., integrator_constraint_tolerance=1e-5,
            gradient_absolute=1e-3, gradient_relative=2e-5, energy_absolute=1e-3,
            energy_relative=2e-6, penalty_relative=2e-5, constraint_rounding_factor=32,
            finite_difference_absolute=1e-5, finite_difference_relative=1e-6,
            native_calls=12, native_retries=0),
        metadata_only_not_binary_hash_attestation=True,
        preflight_limitations=['Platform.getLoadedPlugins is unavailable',
            'external native-library directory discovery permission denied; not bypassed'],
        installed_dependency_modified=False))


def worker(plan_path, index, output):
    import numpy as np
    import openmm
    from openmm.unit import nanometer

    plan = read(plan_path)
    check(plan)
    case = plan['cases'][index]
    system = make_system(case['stress'], case['constrained'])
    integrator = openmm.VerletIntegrator(0.001)
    integrator.setConstraintTolerance(plan['policy']['integrator_constraint_tolerance'])
    selected = openmm.Platform.getPlatformByName(case['backend'])
    context = openmm.Context(system, integrator, selected, case['properties'])
    points = [[-0.0500025, 0., 0.], [0.0500025, 0., 0.], [1., 1., 0.]]
    context.setPositions(points*nanometer)

    def capture(ctx, prefix):
        state = ctx.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
        for kind, obj in (('system', system), ('integrator', integrator), ('state', state)):
            with (output/f'{prefix}-{kind}.xml').open('x') as stream:
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
                raise RuntimeError('frozen callback bound exceeded, not an early success')
            p, g, stats = list(x), list(grad), dict(args)
            if len(p) != 9 or len(g) != 9 or not all(np.isfinite(v) for v in p+g+list(stats.values())):
                raise ValueError('invalid native callback')
            self.records.append(dict(ordinal=len(self.records), iteration=iteration,
                positions_nm=p, gradient=g, statistics=stats))
            return False

    reporter, returned = Recorder(), False
    initial = capture(context, 'before')
    write(output/'platform.json', dict(name=selected.getName(), properties={
        p:selected.getPropertyValue(context, p) for p in selected.getPropertyNames()}))
    try:
        openmm.LocalEnergyMinimizer.minimize(context, plan['policy']['tolerance'],
            plan['policy']['max_iterations'], reporter)
        returned = True
    finally:
        write(output/'path.json.gz', dict(requested_nm=points, initial_nm=initial,
            final_nm=capture(context, 'after'), callbacks=reporter.records, returned=returned))
        del context


def run(plan_path, root):
    plan = read(plan_path)
    check(plan)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    start, rows = time.monotonic(), []
    for case in plan['cases']:
        folder = root/f"case-{case['index']:02d}"
        folder.mkdir()
        cache, temp = folder/'cache', folder/'temp'
        cache.mkdir()
        temp.mkdir()
        row = dict(index=case['index'], label=case['label'], status='not-executed-budget', returncode=None)
        remaining = plan['policy']['total_seconds']-(time.monotonic()-start)
        if remaining > 0:
            env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OPENMM_CACHE_DIR=str(cache),
                XDG_CACHE_HOME=str(cache), TMPDIR=str(temp), TMP=str(temp), TEMP=str(temp))
            command = [plan['runtime']['python'], str(Path(__file__).resolve()), 'worker',
                '--plan', str(plan_path), '--index', str(case['index']), '--output', str(folder)]
            try:
                with (folder/'run.log').open('xb') as log:
                    process = subprocess.Popen(command, cwd=folder, env=env, stdout=log,
                        stderr=subprocess.STDOUT, start_new_session=True)
                    try:
                        code = process.wait(timeout=min(remaining, plan['policy']['per_case_seconds']))
                        row.update(status='completed' if code == 0 else 'execution-error', returncode=code)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        row.update(status='timeout', returncode=process.wait())
            except Exception as error:
                row.update(status='orchestration-error', error_type=type(error).__name__, error=str(error))
        row['artifacts'] = {str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}
        write(root/f"case-{case['index']:02d}-result.json", row)
        rows.append(row)
    check(plan)
    write(root/'execution.json', dict(plan_sha256=sha(plan_path), cases=rows))


def envelope(expected, actual, policy, allowance):
    import numpy as np

    error = np.abs(expected-actual)
    bound = policy['gradient_absolute']+policy['gradient_relative']*np.maximum(1., np.maximum(np.abs(expected), np.abs(actual)))
    adjusted = bound+allowance
    i, j = np.unravel_index(np.argmax(error-adjusted), error.shape)
    return dict(strict_passed=bool(np.all(error <= bound)), allowance_passed=bool(np.all(error <= adjusted)),
        maximum_residual=float(np.max(error)), failed_components=int(np.sum(error > adjusted)),
        witness=dict(particle=int(i), axis='xyz'[j], expected=float(expected[i, j]),
            reported=float(actual[i, j]), original_bound=float(bound[i, j]),
            constraint_arithmetic_allowance=float(allowance[i, j])))


def audit(plan_path, root, output):
    import numpy as np

    plan, execution = read(plan_path), read(root/'execution.json')
    check(plan)
    if execution['plan_sha256'] != sha(plan_path) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('frozen execution mismatch')
    if len(execution['cases']) != len(plan['cases']):
        raise ValueError('incomplete factorial accounting')
    rows = []
    for case, result in zip(plan['cases'], execution['cases'], strict=True):
        if (result['index'], result['label']) != (case['index'], case['label']):
            raise ValueError('factorial identity differs')
        folder = root/f"case-{case['index']:02d}"
        for name, digest in result['artifacts'].items():
            if sha(folder/name) != digest:
                raise ValueError('private raw artifact changed')
        row = dict(case=case, status=result['status'], returncode=result['returncode'], nodes=[])
        path = folder/'path.json.gz'
        if not path.exists():
            row['missing_native_path'] = True
            rows.append(row)
            continue
        raw = read(path)
        row['native_returned'] = raw['returned']
        for callback in raw['callbacks']:
            p = np.asarray(callback['positions_nm']).reshape((3, 3))
            native = np.asarray(callback['gradient']).reshape((3, 3))
            stats, policy = callback['statistics'], plan['policy']
            strength = stats['restraint strength']
            physical, gp = analytic(p, case['stress'], case['constrained'], strength)
            terms = energy_terms(p, case['stress'], case['constrained'], strength)
            fd = finite_difference(p, case['stress'], case['constrained'], strength)
            correct = physical+gp
            fd_bound = policy['finite_difference_absolute']+policy['finite_difference_relative']*np.maximum(
                1., np.maximum(np.abs(correct), np.abs(fd)))
            allowance = np.zeros((3, 3))
            if case['backend'] == 'OpenCL' and case['constrained']:
                allowance[:2] = (policy['constraint_rounding_factor']*np.finfo(np.float32).eps*abs(strength)
                    *(float(np.linalg.norm(p[1]-p[0]))+0.1))
            hypotheses = dict(physical=correct, source_sign=physical-gp,
                source_sign_float32=(physical-gp).astype(np.float32).astype(np.float64))
            energy, ep = float(np.sum(terms[:2])), float(terms[2])
            row['nodes'].append(dict(ordinal=callback['ordinal'], iteration=callback['iteration'],
                reported_gradient_exceeds_fixed_capacity=bool(np.max(np.abs(native)) > 2**31),
                penalty_gradient_maximum=float(np.max(np.abs(gp))),
                finite_difference_passed=bool(np.all(np.abs(correct-fd) <= fd_bound)),
                finite_difference_maximum_error=float(np.max(np.abs(correct-fd))),
                energy_error=abs(energy-stats['system energy']), penalty_energy_error=abs(ep-stats['restraint energy']),
                energy_passed=abs(energy-stats['system energy']) <= policy['energy_absolute']+policy['energy_relative']*max(1., abs(energy)),
                penalty_energy_passed=abs(ep-stats['restraint energy']) <= policy['energy_absolute']+policy['penalty_relative']*max(1., abs(ep)),
                hypotheses={name:envelope(g, native, policy, allowance) for name, g in hypotheses.items()}))
        row['callbacks'] = len(row['nodes'])
        row['primary_signal_nodes'] = [n['ordinal'] for n in row['nodes']
            if case['backend'] == 'OpenCL' and case['constrained'] and case['stress'] > 1000
            and n['reported_gradient_exceeds_fixed_capacity'] and n['penalty_gradient_maximum'] > 0
            and not n['hypotheses']['physical']['allowance_passed']
            and (n['hypotheses']['source_sign']['allowance_passed'] or n['hypotheses']['source_sign_float32']['allowance_passed'])
            and n['energy_passed'] and n['penalty_energy_passed'] and n['finite_difference_passed']]
        row['physical_control_passed'] = result['status'] == 'completed' and bool(row['nodes']) and raw['returned'] and all(
            n['hypotheses']['physical']['allowance_passed'] and n['finite_difference_passed']
            and n['energy_passed'] and n['penalty_energy_passed'] for n in row['nodes'])
        rows.append(row)
    signal = any(r.get('primary_signal_nodes') for r in rows)
    controls = [r for r in rows if not (r['case']['backend'] == 'OpenCL' and r['case']['constrained'] and r['case']['stress'] > 1000)]
    write(output, dict(schema='dvbfixer.private-openmm-sign-reproducer-evidence.v1',
        plan_sha256=sha(plan_path), execution_sha256=sha(root/'execution.json'), cases=rows,
        synthetic_primary_signal_observed=signal,
        all_other_controls_passed=all(r.get('physical_control_passed', False) for r in controls) and len(controls) == 11,
        internal_cpu_fallback_event_directly_observed=False, installed_binary_source_identity_attested=False,
        installed_dependency_modified=False, old_structural_qualification_rehabilitated=False,
        limits='Synthetic public-reporter gradient observation only. Wrong source-sign agreement is not a correct objective gradient. No protein/native FF patch, main merge, external issue or dependency replacement.'))


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
