"""All native whole callbacks: exact replay, physical gradients and terminal state.

Fresh Reference force readback only. Never optimize, repair or silently replace
the returned endpoint with a reported iterate. Report differences explicitly.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_stereo_force_projection import (
    guard_parameters,
    projected_force,
    resolve,
    volume_gradient,
)
from scripts.carbon_H_capture_adapter import read_json
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def finite_points(points, count):
    if len(points) != count or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('native callback positions malformed/nonfinite')


def norm(values):
    return math.sqrt(sum(x*x for x in values))


def centre_reports(case, points, forces, system):
    result = []
    for centre in case['centres']:
        indices = resolve(case['rows'], centre)
        volume, gradient, length = volume_gradient([points[i] for i in indices])
        guard = guard_parameters(system, indices) if centre['names'][0] == 'CA' else None
        if guard is not None:
            q = volume/guard['vref_nm3']
            guard = dict(active=q < guard['qmin'], q=q,
                predicted_projection_kj_mol_nm=guard['k']*max(guard['qmin']-q, 0)/guard['vref_nm3']*length)
        result.append(dict(label=centre['label'], volume_angstrom3=1000*volume,
            native_projection_kj_mol_nm=projected_force(forces, indices, gradient, length), guard=guard))
    return result


def terminal_summary(initial, final, callbacks, nodes):
    previous = [callbacks[-1]['positions_nm'][3*i:3*i+3] for i in range(len(final))] if callbacks else initial
    last = nodes[-2] if callbacks else nodes[0]
    return dict(comparison='last reported iteration' if callbacks else 'initial state; no reported iteration',
        returned_matches_last_reported_nm=final == previous,
        maximum_returned_vs_reported_displacement_angstrom=10*max(math.dist(a, b) for a, b in zip(final, previous, strict=True)),
        returned_minus_reported_energy_kj_mol=nodes[-1]['energy_kj_mol']-last['energy_kj_mol'],
        last_reported_iteration=callbacks[-1]['iteration'] if callbacks else None,
        last_reported_centres=last['centres'], returned_centres=nodes[-1]['centres'])


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    if openmm.__version__ != plan['openmm_version']:
        raise ValueError('native trajectory readback runtime differs')
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['conditioning'], c['arm']) for c in plan['cases']]
    if (run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw
            or [(c['id'], c['seed'], c['conditioning'], c['arm']) for c in run['cases']] != expected):
        raise ValueError('complete matched native path cohort differs')
    results, errors = [], []
    for case, cell in zip(plan['cases'], run['cases'], strict=True):
        folder = root/cell['id']/f"seed-{cell['seed']}"/cell['conditioning']/cell['arm']
        source = Path(case['source'])
        for name, digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('native path artifact changed')
        report = {k:cell[k] for k in ('id', 'seed', 'conditioning', 'arm', 'status', 'returncode')}
        if cell['status'] != 'completed':
            results.append(report)
            continue
        trace = read_json(folder/'trace.json')
        if (trace['case_id'] != list(expected[len(results)]) or trace['plan_sha256'] != sha(raw)
                or trace['reporter_requests_early_exit'] or trace['new_coordinates_or_physics']
                or [p['phase'] for p in trace['phases']] != ['phase1', 'phase2']):
            raise ValueError('native reporter trace contract differs')
        report['phases'] = []
        for phase_record in trace['phases']:
            phase = phase_record['phase']
            path_file = folder/f'{phase}-path.json.gz'
            if sha(path_file.read_bytes()) != phase_record['path_sha256']:
                raise ValueError('native path payload changed')
            path = read_json(path_file)
            before_stage = 'initial' if phase == 'phase1' else 'phase1'
            if (not path['native_call_returned'] or path['initial_nm'] != read_json(source/f'{before_stage}.json.gz')['positions_nm']
                    or path['final_nm'] != read_json(source/f'{phase}.json.gz')['positions_nm']
                    or not phase_record['exact_original_nm_endpoint']):
                raise ValueError('native path start/endpoint is not exact original nm')
            xml = (source/f'{phase}-system.xml').read_text()
            if sha(xml.encode()) != phase_record['system_sha256']:
                raise ValueError('native path force/anchor parameters differ')
            system = openmm.XmlSerializer.deserialize(xml)
            count = system.getNumParticles()
            if system.getNumConstraints() or any(system.getParticleMass(i).value_in_unit(dalton) <= 0 for i in range(count)):
                raise ValueError('native path is not unconstrained/all-mobile')
            callbacks = path['callbacks']
            if len(callbacks) != phase_record['callbacks'] or len(callbacks) > plan['policy']['maximum_callbacks_per_phase']:
                raise ValueError('native callback storage budget differs')
            initial_norm = max(1., math.sqrt(sum(x*x for p in path['initial_nm'] for x in p)/count))
            epsilon = plan['policy']['tolerance_kj_mol_nm']/initial_norm
            integrator = openmm.VerletIntegrator(0.)
            context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
            try:
                initial_state = openmm.XmlSerializer.deserialize(Path(case['state']).read_text())
                if dict(context.getParameters()) != dict(initial_state.getParameters()):
                    raise ValueError('native path global parameters differ from captured State')
                entries = [('initial', path['initial_nm'], None)]
                for i, cb in enumerate(callbacks):
                    if cb['ordinal'] != i or cb['iteration'] != i or len(cb['positions_nm']) != 3*count or len(cb['gradient_kj_mol_nm']) != 3*count:
                        raise ValueError('native unconstrained callback iteration/inventory differs')
                    entries.append((f'callback-{i}', [cb['positions_nm'][3*j:3*j+3] for j in range(count)], cb))
                entries.append(('returned', path['final_nm'], None))
                nodes = []
                for label, points, cb in entries:
                    finite_points(points, count)
                    context.setPositions(points*nanometer)
                    state = context.getState(energy=True, forces=True)
                    energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                    forces = [list(f) for f in state.getForces().value_in_unit(kilojoules_per_mole/nanometer)]
                    force_norm = norm([x for f in forces for x in f])
                    ratio = force_norm/max(1., norm([x for p in points for x in p]))
                    error = None
                    if cb is not None:
                        gradient = cb['gradient_kj_mol_nm']
                        if (any(not math.isfinite(x) for x in gradient)
                                or any(not math.isclose(forces[i][j], -gradient[3*i+j], rel_tol=1e-8, abs_tol=1e-4) for i in range(count) for j in range(3))
                                or not math.isclose(energy, cb['statistics']['system energy'], rel_tol=1e-10, abs_tol=1e-5)
                                or cb['statistics']['restraint energy'] != 0 or cb['statistics']['max constraint error'] != 0):
                            raise ValueError('native whole callback gradient/objective differs from physical readback')
                        error = max(abs(forces[i][j]+gradient[3*i+j]) for i in range(count) for j in range(3))
                        errors.append(error)
                    nodes.append(dict(stage=label, energy_kj_mol=energy,
                        force_l2_norm_kj_mol_nm=force_norm, native_scaled_gradient_ratio=ratio,
                        native_epsilon=epsilon, satisfies_native_gradient_test=ratio <= epsilon,
                        centres=centre_reports(case, points, forces, system),
                        maximum_gradient_readback_error_kj_mol_nm=error))
                report['phases'].append(dict(phase=phase, callbacks=len(callbacks), nodes=nodes,
                    exact_original_nm_start=True, exact_original_nm_endpoint=True,
                    terminal=terminal_summary(path['initial_nm'], path['final_nm'], callbacks, nodes)))
            finally:
                del context
                del integrator
        results.append(report)
    verify(plan)
    phases = [p for c in results if c['status'] == 'completed' for p in c['phases']]
    summary = dict(enrolled=len(results), statuses=dict(Counter(c['status'] for c in results)),
        phases=len(phases), callbacks=sum(p['callbacks'] for p in phases),
        phases_without_reported_iteration=sum(p['callbacks'] == 0 for p in phases),
        returned_differs_from_reported_or_initial_phases=sum(not p['terminal']['returned_matches_last_reported_nm'] for p in phases),
        returned_above_reported_or_initial_energy_phases=sum(p['terminal']['returned_minus_reported_energy_kj_mol'] > 1e-6 for p in phases),
        returned_fails_native_gradient_test=sum(not p['nodes'][-1]['satisfies_native_gradient_test'] for p in phases),
        maximum_gradient_readback_error_kj_mol_nm=max(errors, default=None))
    write(output, dict(schema='dvbfixer.carbon-H-whole-path-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), policy=plan['policy'],
        audit_script_sha256=sha(Path(__file__).read_bytes()), production_changed=False,
        coordinates_restored=False, extra_minimization_executed=False, scientific_qualification=False,
        exact_lbfgs_return_code_observed=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
