"""Retrospective saved-only force attribution on ALL THR24 factorial states.

24 matched paths/72 states, no minimization, integration or force ablation.
Signed-volume projections are measurements, not reconstructed optimizer steps.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

from scripts.audit_stereo_force_projection import (
    focal_centres,
    grouped_system,
    guard_parameters,
    projected_force,
    resolve,
    volume_gradient,
)
from scripts.carbon_H_capture_adapter import read_json
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def finite_gradient_error(points):
    _, analytic, _ = volume_gradient(points)
    errors = []
    for i in range(4):
        for j in range(3):
            positive, negative = [list(p) for p in points], [list(p) for p in points]
            positive[i][j] += 1e-7
            negative[i][j] -= 1e-7
            observed = (volume_gradient(positive)[0]-volume_gradient(negative)[0])/2e-7
            if not math.isclose(observed, analytic[i][j], rel_tol=1e-6, abs_tol=1e-9):
                raise ValueError('Cartesian signed-volume gradient finite difference differs')
            errors.append(abs(observed-analytic[i][j]))
    return max(errors)


def readback(system, points, initial_parameters, saved_energy, centres, rows):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    grouped, groups = grouped_system(system)
    contexts, integrators = [], []
    try:
        for obj in (system, grouped):
            integrator = openmm.VerletIntegrator(0.)
            integrators.append(integrator)
            context = openmm.Context(obj, integrator, openmm.Platform.getPlatformByName('Reference'))
            contexts.append(context)
            if dict(context.getParameters()) != initial_parameters:
                raise ValueError('snapshot global parameters differ from captured native State')
            context.setPositions(points*nanometer)
        native = contexts[0].getState(energy=True, forces=True)
        energy = native.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        original = [list(p) for p in native.getForces().value_in_unit(kilojoules_per_mole/nanometer)]
        if not math.isclose(energy, saved_energy, rel_tol=1e-10, abs_tol=1e-6):
            raise ValueError('original snapshot energy readback differs')
        component_forces, components = [], []
        for group in groups:
            state = contexts[1].getState(energy=True, forces=True, groups=1 << group['index'])
            component_forces.append([list(p) for p in state.getForces().value_in_unit(kilojoules_per_mole/nanometer)])
            components.append(dict(group, energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)))
        errors = []
        for i in range(len(points)):
            for j in range(3):
                summed = sum(f[i][j] for f in component_forces)
                if not math.isclose(summed, original[i][j], rel_tol=1e-8, abs_tol=1e-5):
                    raise ValueError('all-particle grouped force sum differs from original native System')
                errors.append(abs(summed-original[i][j]))
        energy_error = abs(sum(c['energy_kj_mol'] for c in components)-energy)
        if not math.isclose(sum(c['energy_kj_mol'] for c in components), energy, rel_tol=1e-10, abs_tol=1e-6):
            raise ValueError('grouped energy sum differs from original native System')
        measurements = []
        for centre in centres:
            indices = resolve(rows, centre)
            focal = [points[i] for i in indices]
            volume, gradient, norm = volume_gradient(focal)
            finite_error = finite_gradient_error(focal)
            total = projected_force(original, indices, gradient, norm)
            projections = [projected_force(f, indices, gradient, norm) for f in component_forces]
            if not math.isclose(sum(projections), total, rel_tol=1e-8, abs_tol=1e-5):
                raise ValueError('signed-volume group projections do not close')
            guard = guard_parameters(system, indices) if centre['names'][0] == 'CA' else None
            if guard is not None:
                q = volume/guard['vref_nm3']
                deficit = max(guard['qmin']-q, 0.)
                predicted = guard['k']*deficit/guard['vref_nm3']*norm
                if not math.isclose(projections[guard['group']], predicted, rel_tol=1e-8, abs_tol=1e-5):
                    raise ValueError('native CA guard analytic projection differs')
                guard.update(q=q, active=deficit > 0, analytic_projection_kj_mol_nm=predicted)
            measurements.append(dict(centre, indices=indices, volume_angstrom3=1000*volume,
                gradient_nm2=[list(g) for g in gradient], gradient_norm_nm2=norm,
                finite_difference_max_error_nm2=finite_error,
                original_native_projection_kj_mol_nm=total,
                component_projections_kj_mol_nm=projections, guard=guard))
        return dict(original_energy_kj_mol=energy, force_components=components, centres=measurements,
            maximum_native_force_closure_error_kj_mol_nm=max(errors), energy_closure_error_kj_mol=energy_error), dict(
                original_forces_kj_mol_nm=original, component_forces_kj_mol_nm=component_forces)
    finally:
        # Contexts must be released before their integrators.
        context = None
        contexts.clear()
        integrator = None
        integrators.clear()


def plan(whole_plan, whole_root, whole_audit, output):
    raw = whole_plan.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    audit_raw = whole_audit.read_bytes()
    audit = json.loads(audit_raw)
    run_raw = (whole_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if (audit['plan_sha256'] != sha(raw) or audit['run_sha256'] != sha(run_raw)
            or audit['summary']['statuses'] != {'completed':136, 'refused_source':24}
            or audit['summary']['complete_matched_contexts'] != 34):
        raise ValueError('complete whole factorial audit required')
    inputs = {**parent['inputs'], str(whole_plan):sha(raw), str(whole_audit):sha(audit_raw),
        str(whole_root/'evidence.json'):sha(run_raw)}
    selected = [c for c in parent['cases'] if c['id'] == '1SHG-THR24']
    if [(c['seed'], c['conditioning']) for c in selected] != [(s, a) for s in (1, 2) for a in ('raw', 'beta_only', 'full_generated')]:
        raise ValueError('all six original THR24 conditioning contexts required')
    cases = []
    for case in selected:
        for arm in parent['policy']['arms']:
            cell = next(c for c in run['cases'] if (c['id'], c['seed'], c['conditioning'], c['arm'])
                == (case['id'], case['seed'], case['conditioning'], arm))
            if cell['status'] != 'completed':
                raise ValueError('matched THR24 whole path has blocker')
            folder = whole_root/case['id']/f"seed-{case['seed']}"/case['conditioning']/arm
            for name, digest in cell['artifacts'].items():
                if sha((folder/name).read_bytes()) != digest:
                    raise ValueError('original whole factorial artifact changed')
                inputs[str(folder/name)] = digest
            cases.append(dict(id=case['id'], seed=case['seed'], conditioning=case['conditioning'], arm=arm,
                source=str(folder), state=case['whole']['state'], rows=case['whole']['rows'],
                centres=focal_centres(case['whole'])))
    sources = {name:sha(Path(name).read_bytes()) for name in (
        'scripts/attribute_carbon_H_whole.py', 'scripts/audit_stereo_force_projection.py',
        'scripts/research_dunbrack_prior.py', 'scripts/research_ha_minimizer_path.py')}
    write(output, dict(schema='dvbfixer.carbon-H-whole-attribution-plan.v1', cases=cases,
        inputs=inputs, source_hashes=sources, openmm_version=parent['openmm_version'],
        policy=dict(timing='retrospective to known near-planar THR24 signs; frozen before new force-contributor inspection',
            enrollment='all six THR24 conditioning contexts/two existing seeds, all four arms, all three states;24 paths/72 states and216 focal measurements',
            centres='target CA, target THR-beta and context GLN16 CA, exact case-sensitive identities',
            method='fresh original and measurement-group-cloned native Reference Contexts; parameter-exact restoration of force grouping, full force/energy sums and Cartesian finite-difference signed-volume gradients',
            interpretation='unit-volume-gradient projections are not actual quasi-Newton steps, integrated work, a captured crossing or unique pair causality',
            limits='saved-only; no minimization/integration/new H, force ablation, seed/phase/strength/iteration search, native-accuracy score, repair or production qualification'),
        production_changed=False, minimization_executed=False))


def run(plan_path, output):
    import openmm

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('native readback runtime differs')
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    cases = []
    for ordinal, case in enumerate(record['cases']):
        state = openmm.XmlSerializer.deserialize(Path(case['state']).read_text())
        source = Path(case['source'])
        frames = []
        for stage in ('initial', 'phase1', 'phase2'):
            data = read_json(source/f'{stage}.json.gz')
            system = openmm.XmlSerializer.deserialize((source/f'{stage}-system.xml').read_text())
            report, forces = readback(system, data['positions_nm'], dict(state.getParameters()),
                data['energy_kj_mol'], case['centres'], case['rows'])
            artifact = output/f'case-{ordinal+1:02d}-{stage}-forces.json.gz'
            with artifact.open('xb') as stream:
                stream.write(gzip.compress(json.dumps(forces, sort_keys=True, allow_nan=False).encode(), mtime=0))
            frames.append(dict(stage=stage, raw_artifact=artifact.name, raw_sha256=sha(artifact.read_bytes()), **report))
        cases.append(dict(case_id=[case['id'], case['seed'], case['conditioning'], case['arm']], frames=frames))
    verify(record)
    nodes = [f for c in cases for f in c['frames']]
    summary = dict(paths=len(cases), states=len(nodes), focal_measurements=sum(len(n['centres']) for n in nodes),
        maximum_native_force_closure_error_kj_mol_nm=max(n['maximum_native_force_closure_error_kj_mol_nm'] for n in nodes),
        maximum_energy_closure_error_kj_mol=max(n['energy_closure_error_kj_mol'] for n in nodes),
        maximum_finite_difference_error_nm2=max(c['finite_difference_max_error_nm2'] for n in nodes for c in n['centres']))
    write(output/'evidence.json', dict(schema='dvbfixer.carbon-H-whole-attribution.v1', cases=cases,
        summary=summary, plan_sha256=sha(raw), policy=record['policy'], production_changed=False,
        minimization_executed=False, scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    for name in ('whole-plan', 'whole-root', 'whole-audit', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.whole_plan, args.whole_root, args.whole_audit, args.output)
    else:
        run(args.plan, args.output)
