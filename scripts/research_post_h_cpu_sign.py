"""Read-only source-sign and float32 transport hypotheses, never a physics fix."""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from scripts.audit_post_h_1vcu import penalty
from scripts.audit_thr2_native_capture import read
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed
from scripts.research_post_h_1vcu import check as check_original

ROOT = Path(__file__).resolve().parents[1]


def predictions(physical, positive_penalty):
    import numpy as np

    physical, positive_penalty = np.asarray(physical), np.asarray(positive_penalty)
    source = physical-positive_penalty
    return dict(physical=physical+positive_penalty, source_sign=source,
                source_sign_float32=source.astype(np.float32).astype(np.float64))


def compare(expected, native, policy, atom_rows, constrained):
    import numpy as np

    native = np.asarray(native).reshape(expected.shape)
    residual = np.abs(expected-native)
    bounds = policy['gradient_absolute_bound']+policy['gradient_relative_bound']*np.maximum(
        1., np.maximum(np.abs(expected), np.abs(native)))
    failed = residual > bounds
    atom, axis = np.unravel_index(np.argmax(residual-bounds), residual.shape)
    return dict(passed=not bool(np.any(failed)), failed_components=int(np.sum(failed)),
        failed_components_on_constrained_atoms=int(np.sum(failed[list(constrained)])) if constrained else 0,
        maximum_absolute_residual=float(np.max(residual)),
        maximum_envelope_excess=float(np.max(residual-bounds)), witness=dict(
            atom_key=atom_rows[atom]['key'], axis='xyz'[axis], expected=float(expected[atom, axis]),
            reported=float(native[atom, axis]), residual=float(residual[atom, axis]),
            original_bound=float(bounds[atom, axis]), constrained=atom in constrained))


def freeze(output):
    old_path = ROOT/'docs/research/post-hydrogen-overflow-readback-plan.json'
    evidence_path = ROOT/'docs/research/post-hydrogen-overflow-readback-evidence.json'
    old = read(old_path)
    verify(old)
    check_original(read(Path(old['parent_plan'])))
    evidence = read(evidence_path)
    if not all(len(b['cases']) == 26 and all(c['status'] == 'readback'
               and not c['original_numerical_envelopes_passed'] for c in b['cases']) for b in evidence['cases']):
        raise ValueError('original complete secondary failures required')
    source_path = ROOT/'tmp/pdbfixer-research/common-minimize-36a30cb.cpp'
    source = source_path.read_text()
    for token in ('hostPositions[indices.y] - hostPositions[indices.x]',
                  'hostGrad[3 * indices.y] -= kdr * delta[0]',
                  'hostGrad[3 * indices.x] += kdr * delta[0]', 'grad.upload(hostGrad, true)'):
        if token not in source:
            raise ValueError('source-sign hypothesis does not match retained exact-revision source')
    inputs = {str(p):sha(p.read_bytes()) for p in (old_path, evidence_path, source_path,
        ROOT/'docs/plans/post-hydrogen-cpu-sign-readback.md')}
    for selected in old['selection']:
        p = Path(old['root'])/f"{selected['call']['before']['prefix']}-atoms.json.gz"
        inputs[str(p)] = sha(p.read_bytes())
    files = [Path(__file__).resolve(), ROOT/'tests/test_post_h_cpu_sign.py']
    write(output, dict(schema='dvbfixer.private-cpu-sign-readback-plan.v1', parent_plan=str(old_path),
        inputs=inputs, source_hashes={str(p):sha(p.read_bytes()) for p in files},
        selection=old['selection'], root=old['root'], policy=dict(maximum_seconds=300,
            contexts=1, queries=26, native_minimizations=0, h_samples=0),
        hypotheses=['physical', 'source_sign', 'source_sign_float32'],
        numerical_policy=old['policy'], physical_qualification_rehabilitated=False))


def check(plan):
    verify(plan)
    old = read(Path(plan['parent_plan']))
    verify(old)
    check_original(read(Path(old['parent_plan'])))


def run(plan_path, root, output):
    import numpy as np
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    plan = read(plan_path)
    check(plan)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    source = Path(plan['root'])
    selected, = plan['selection']
    prefix = selected['call']['before']['prefix']
    system = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((source/f'{prefix}-integrator.xml').read_text())
    if any(system.getParticleMass(i).value_in_unit(dalton) <= 0 for i in range(system.getNumParticles())):
        raise ValueError('all-positive masses required for this fixed source interpretation')
    constraints = []
    for i in range(system.getNumConstraints()):
        x, y, distance = system.getConstraintParameters(i)
        constraints.append((x, y, distance.value_in_unit(nanometer)))
    constrained = {i for x, y, _ in constraints for i in (x, y)}
    atoms = read(source/f'{prefix}-atoms.json.gz')['atoms']
    payload = read(source/f"whole-{selected['number']}-path.json.gz")
    callbacks = [payload['callbacks'][i] for i in selected['ordinals']]
    del payload
    rows, context, start = [], None, time.monotonic()
    try:
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('CPU'))
        context.setState(openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()))
        for ordinal, callback in zip(selected['ordinals'], callbacks, strict=True):
            row = dict(ordinal=ordinal, status='not-executed-budget')
            rows.append(row)
            if time.monotonic()-start > plan['policy']['maximum_seconds']:
                continue
            points = np.asarray(callback['positions_nm']).reshape((-1, 3))
            context.setPositions(points.tolist()*nanometer)
            state = context.getState(energy=True, forces=True, groups=integrator.getIntegrationForceGroups())
            physical = -state.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
            energy = float(state.getPotentialEnergy().value_in_unit(kilojoules_per_mole))
            stats = callback['statistics']
            ep, gp, _, _ = penalty(points, constraints, stats['restraint strength'])
            if not np.all(np.isfinite(physical)) or not np.isfinite(energy):
                raise ValueError('nonfinite CPU readback')
            native = np.asarray(callback['gradient_kj_mol_nm']).reshape((-1, 3))
            hypotheses = predictions(physical, gp)
            compressed(root/f'callback-{ordinal:03d}.json.gz', dict(physical=physical.tolist(),
                positive_penalty_gradient=gp.tolist(), native=native.tolist(), energy=energy,
                penalty_energy=ep))
            row.update(status='readback', energy_error=abs(energy-stats['system energy']),
                penalty_energy_error=abs(ep-stats['restraint energy']), hypotheses={
                    name:compare(value, native, plan['numerical_policy'], atoms, constrained)
                    for name, value in hypotheses.items()})
            bounds = plan['numerical_policy']
            row['original_energy_envelope_passed'] = row['energy_error'] <= (
                bounds['energy_absolute_bound']+bounds['energy_relative_bound']*max(1., abs(energy)))
            row['original_penalty_energy_envelope_passed'] = row['penalty_energy_error'] <= (
                bounds['penalty_energy_absolute_bound']+bounds['penalty_energy_relative_bound']*max(1., abs(stats['restraint energy'])))
    except Exception as error:
        accounted = {r['ordinal'] for r in rows if r['status'] in ('readback', 'not-executed-budget') and 'hypotheses' in r}
        rows = [r for r in rows if r['ordinal'] in accounted]
        rows.extend(dict(ordinal=i, status='readback-error', error_type=type(error).__name__, error=str(error))
                    for i in selected['ordinals'] if i not in accounted)
    finally:
        if context is not None:
            del context
    check(plan)
    write(output, dict(schema='dvbfixer.private-cpu-sign-readback-evidence.v1', cases=rows,
        plan_sha256=sha(plan_path.read_bytes()), original_tolerances_changed=False,
        physical_qualification_rehabilitated=False, installed_binary_attested=False,
        new_minimizations=0, new_h_samples=0,
        private_artifacts={p.name:sha(p.read_bytes()) for p in sorted(root.glob('*')) if p.is_file()},
        limits='Source-sign/quantization attribution is distinct from a physically correct objective gradient. Old strict/native failures remain preserved; no OpenMM patch or public capability.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--root', type=Path)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    else:
        run(args.plan.resolve(), args.root.resolve(), args.output.resolve())
