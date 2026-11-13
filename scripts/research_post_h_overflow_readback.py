"""Read-only CPU/Reference attribution of every failed 1VCU GPU callback."""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from scripts.audit_post_h_1vcu import penalty
from scripts.audit_thr2_native_capture import read
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_post_h_1vcu import check as check_parent

ROOT = Path(__file__).resolve().parents[1]


def compare(physical, energy, points, gradient, statistics, constraints, policy):
    import numpy as np

    ep, gp, _, _ = penalty(points, constraints, statistics['restraint strength'])
    expected = np.asarray(physical)+gp
    native = np.asarray(gradient).reshape(expected.shape)
    error = np.abs(expected-native)
    bound = policy['gradient_absolute_bound']+policy['gradient_relative_bound']*np.maximum(1., np.maximum(np.abs(expected), np.abs(native)))
    ee = abs(energy-statistics['system energy'])
    pe = abs(ep-statistics['restraint energy'])
    passed = (bool(np.all(error <= bound))
        and ee <= policy['energy_absolute_bound']+policy['energy_relative_bound']*max(1., abs(energy))
        and pe <= policy['penalty_energy_absolute_bound']+policy['penalty_energy_relative_bound']*max(1., abs(statistics['restraint energy'])))
    return dict(gradient_maximum_error=float(np.max(error)), system_energy_error=ee,
                penalty_energy_error=pe, original_numerical_envelopes_passed=passed,
                reported_gradient_exceeds_fixed_force_capacity=bool(np.max(np.abs(native)) > 2**31))


def freeze(output):
    parent_path = ROOT/'docs/research/post-hydrogen-1vcu-plan.json'
    parent = read(parent_path)
    check_parent(parent)
    audit_path = ROOT/'docs/research/post-hydrogen-1vcu-audit.json'
    audit = read(audit_path)
    root = ROOT/'tmp/pdbfixer-research/post-hydrogen-1vcu-public-api-01/case-00'
    capture = read(root/'capture.json')
    selection = []
    for phase in audit['cases'][0]['whole_calls']:
        bad = [int(n['stage'].split('-')[1]) for n in phase['nodes'] if n.get('numerical_envelopes_passed') is False]
        if bad:
            selection.append(dict(number=phase['number'], ordinals=bad,
                call=next(c for c in capture['whole_calls'] if c['number'] == phase['number'])))
    if len(selection) != 1 or len(selection[0]['ordinals']) != 26:
        raise ValueError('complete original 26-callback failed-readback selection required')
    inputs = {**parent['inputs'], **parent['source_hashes']}
    files = [parent_path, audit_path, root/'capture.json', ROOT/'docs/plans/post-hydrogen-overflow-readback.md']
    for selected in selection:
        prefix = selected['call']['before']['prefix']
        files += [root/f'{prefix}-{name}.xml' for name in ('system', 'integrator', 'state')]
        files.append(root/f"whole-{selected['number']}-path.json.gz")
    for p in files:
        inputs[str(p)] = sha(p.read_bytes())
    write(output, dict(schema='dvbfixer.private-post-h-overflow-readback-plan.v1', inputs=inputs,
        source_hashes={str(p):sha(p.read_bytes()) for p in
            (Path(__file__).resolve(), ROOT/'tests/test_post_h_overflow_readback.py')},
        parent_plan=str(parent_path), root=str(root), selection=selection, backends=['CPU', 'Reference'],
        policy=dict(parent['policy'], backend_seconds=300, total_readback_seconds=600),
        source_reference='https://github.com/openmm/openmm/blob/36a30cbca54e727b216b606f3c011b67201eb8b4/platforms/common/src/CommonMinimizeKernel.cpp#L575-L637',
        original_gpu_readback_passed=False, new_minimizations=0))


def run(plan_path, output):
    import numpy as np
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    plan = read(plan_path)
    verify(plan)
    check_parent(read(Path(plan['parent_plan'])))
    root = Path(plan['root'])
    selected, = plan['selection']
    prefix = selected['call']['before']['prefix']
    payload = read(root/f"whole-{selected['number']}-path.json.gz")
    callbacks = [payload['callbacks'][i] for i in selected['ordinals']]
    del payload
    reports, started = [], time.monotonic()
    for backend in plan['backends']:
        begin = time.monotonic()
        system = openmm.XmlSerializer.deserialize((root/f'{prefix}-system.xml').read_text())
        integrator = openmm.XmlSerializer.deserialize((root/f'{prefix}-integrator.xml').read_text())
        constraints = []
        for i in range(system.getNumConstraints()):
            a, b, r0 = system.getConstraintParameters(i)
            constraints.append((a, b, r0.value_in_unit(nanometer)))
        if any(system.getParticleMass(i).value_in_unit(dalton) <= 0 for i in range(system.getNumParticles())):
            raise ValueError('this readback requires all-positive whole-System masses')
        context = None
        rows = []
        try:
            platform = openmm.Platform.getPlatformByName(backend)
            context = openmm.Context(system, integrator, platform)
            context.setState(openmm.XmlSerializer.deserialize((root/f'{prefix}-state.xml').read_text()))
            for ordinal, cb in zip(selected['ordinals'], callbacks, strict=True):
                row = dict(ordinal=ordinal, status='not-executed-budget')
                if time.monotonic()-begin <= plan['policy']['backend_seconds'] and time.monotonic()-started <= plan['policy']['total_readback_seconds']:
                    points = np.asarray(cb['positions_nm']).reshape((-1, 3))
                    context.setPositions(points.tolist()*nanometer)
                    state = context.getState(energy=True, forces=True, groups=integrator.getIntegrationForceGroups())
                    forces = state.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
                    energy = float(state.getPotentialEnergy().value_in_unit(kilojoules_per_mole))
                    if not np.all(np.isfinite(forces)) or not np.isfinite(energy):
                        raise ValueError('nonfinite secondary readback')
                    row.update(status='readback', **compare(-forces, energy, points, cb['gradient_kj_mol_nm'], cb['statistics'], constraints, plan['policy']))
                rows.append(row)
        except Exception as error:
            already = {r['ordinal'] for r in rows}
            rows.extend(dict(ordinal=i, status='readback-error', error_type=type(error).__name__, error=str(error))
                        for i in selected['ordinals'] if i not in already)
        finally:
            if context is not None:
                del context
        reports.append(dict(backend=backend, cases=rows))
    verify(plan)
    write(output, dict(schema='dvbfixer.private-post-h-overflow-readback-evidence.v1',
        plan_sha256=sha(plan_path.read_bytes()), cases=reports, original_gpu_readback_passed=False,
        new_minimizations=0, new_h_sampling=0,
        limits='Secondary read-only attribution; original GPU audit failure and guarded refusal remain unchanged. No constrained-stationarity or production-readiness claim.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plan', type=Path)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    else:
        run(args.plan.resolve(), args.output.resolve())
