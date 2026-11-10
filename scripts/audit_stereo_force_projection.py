"""Retrospective force readback at ALL saved cohort phase boundaries.

No minimization, integration, new H placement, native-reference scoring or repair.
Old scientific source bindings resolve against an explicit preserved baseline;
new analysis binds the active source after the diagnostic-only cherry-picks.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path

from scripts.audit_legacy_start_controls import load_atoms
from scripts.research_dunbrack_prior import sha, write

from dvbfixer.diagnose.geometry import stereocentre_geometry


def historical_verify(record, root):
    for name, digest in {**record['inputs'], **record['source_hashes']}.items():
        path = Path(name)
        if sha((path if path.is_absolute() else root / path).read_bytes()) != digest:
            raise ValueError(f'historical binding changed: {name}')


def verify(record):
    for name, digest in {**record['inputs'], **record['source_hashes']}.items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError(f'analysis binding changed: {name}')


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def volume_gradient(points):
    """Triple product and Cartesian gradient; points ordered centre, A, B, C.

    Coordinates in nm produce volume in nm³ and gradients in nm². A unit-gradient
    force projection is in kJ/mol/nm; it is not a constrained minimizer step,
    mass-weighted acceleration, integrated work or evidence of a whole trajectory.
    """
    if len(points) != 4 or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('four finite three-dimensional points required')
    a, b, c = [tuple(x-y for x, y in zip(p, points[0])) for p in points[1:]]
    ga, gb, gc = cross(b, c), cross(c, a), cross(a, b)
    gradient = [tuple(-sum(x) for x in zip(ga, gb, gc)), ga, gb, gc]
    volume = sum(x*y for x, y in zip(a, ga))
    norm = math.sqrt(sum(x*x for g in gradient for x in g))
    if not math.isfinite(volume) or not math.isfinite(norm) or norm <= 1e-18:
        raise ValueError('undefined volume gradient')
    return volume, gradient, norm


def projected_force(forces, indices, gradient, norm):
    value = sum(sum(x*y for x, y in zip(forces[i], g)) for i, g in zip(indices, gradient)) / norm
    if not math.isfinite(value):
        raise ValueError('nonfinite force projection')
    return value


def focal_centres(case):
    target = case['target']
    identity = [target['chain_id'], target['sequence_number'], target['insertion_code'], case['residue']]
    result = [dict(label='target_ca', residue=identity, names=['CA', 'N', 'C', 'CB'])]
    if case['residue'] == 'THR':
        result.append(dict(label='target_beta', residue=identity, names=['CB', 'CA', 'OG1', 'CG2']))
    result.append(dict(label='context_ca', residue=['A', '16', '', 'GLN'], names=['CA', 'N', 'C', 'CB']))
    return result


def resolve(rows, centre):
    mapped = {tuple(row['key']): i for i, row in enumerate(rows)}
    if len(mapped) != len(rows):
        raise ValueError('duplicate saved atom identity')
    keys = [tuple(centre['residue'] + [name]) for name in centre['names']]
    if any(key not in mapped for key in keys):
        raise ValueError('focal centre missing an exact atom')
    return [mapped[key] for key in keys]


def grouped_system(system):
    """Only measurement groups change; require exact restoration of serialized FF."""
    import openmm

    if system.getNumForces() > 32:
        raise ValueError('force inventory exceeds available measurement groups')
    original = openmm.XmlSerializer.serialize(system)
    grouped = openmm.XmlSerializer.deserialize(original)
    groups = []
    for i, force in enumerate(grouped.getForces()):
        row = dict(index=i, type=type(force).__name__, name=force.getName(), original_group=force.getForceGroup())
        if isinstance(force, openmm.NonbondedForce):
            row['original_reciprocal_group'] = force.getReciprocalSpaceForceGroup()
            force.setReciprocalSpaceForceGroup(i)
        force.setForceGroup(i)
        groups.append(row)
    restored = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(grouped))
    for force, row in zip(restored.getForces(), groups, strict=True):
        force.setForceGroup(row['original_group'])
        if 'original_reciprocal_group' in row:
            force.setReciprocalSpaceForceGroup(row['original_reciprocal_group'])
    if openmm.XmlSerializer.serialize(restored) != original:
        raise ValueError('measurement grouping changed more than group metadata')
    return grouped, groups


def guard_parameters(system, indices):
    import openmm

    matches = []
    for i, force in enumerate(system.getForces()):
        if not isinstance(force, openmm.CustomCompoundBondForce):
            continue
        if force.getNumParticlesPerBond() != 4 or force.getNumPerBondParameters() != 1:
            continue
        if force.getPerBondParameterName(0) != 'vref' or 'q=v/vref;' not in force.getEnergyFunction():
            continue
        globals_ = {force.getGlobalParameterName(j): force.getGlobalParameterDefaultValue(j)
                    for j in range(force.getNumGlobalParameters())}
        if globals_.keys() != {'k', 'qmin'}:
            raise ValueError('unrecognized CA guard global parameters')
        for j in range(force.getNumBonds()):
            particles, parameters = force.getBondParameters(j)
            if list(particles) == [indices[1], indices[0], indices[2], indices[3]]:
                matches.append(dict(group=i, vref_nm3=parameters[0], **globals_))
    if len(matches) > 1:
        raise ValueError('ambiguous C-alpha guard reference')
    return matches[0] if matches else None


def phase_readback(folder, phase, case):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    before_xml = (folder/f'phase-{phase}-before-system.xml').read_text()
    if before_xml != (folder/f'phase-{phase}-after-system.xml').read_text():
        raise ValueError('system parameters changed within captured native phase')
    original = openmm.XmlSerializer.deserialize(before_xml)
    system, groups = grouped_system(original)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        # The original pipeline has no setParameter calls; bind the recorded
        # initial globals and require equality to each serialized phase's defaults.
        initial = openmm.XmlSerializer.deserialize((folder/'initial-state.xml').read_text())
        if dict(initial.getParameters()) != dict(context.getParameters()):
            raise ValueError('recorded globals differ from serialized defaults')
        results = []
        for side in ('before', 'after'):
            rows = load_atoms(folder/f'phase-{phase}-{side}.json.gz')
            if len(rows) != system.getNumParticles():
                raise ValueError('system and saved particle inventory differ')
            points = [tuple(x/10 for x in row['xyz']) for row in rows]
            context.setPositions(points * nanometer)
            total = context.getState(energy=True, forces=True)
            total_energy = total.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            total_forces = total.getForces().value_in_unit(kilojoules_per_mole/nanometer)
            component_forces = []
            components = []
            for group in groups:
                state = context.getState(energy=True, forces=True, groups={group['index']})
                component_forces.append(state.getForces().value_in_unit(kilojoules_per_mole/nanometer))
                components.append(dict(group, energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)))
            if not math.isfinite(total_energy) or not math.isclose(sum(c['energy_kj_mol'] for c in components), total_energy, rel_tol=1e-10, abs_tol=1e-6):
                raise ValueError('group energy sum differs from total')
            centres = []
            for centre in focal_centres(case):
                indices = resolve(rows, centre)
                volume, gradient, norm = volume_gradient([points[i] for i in indices])
                projections = [projected_force(f, indices, gradient, norm) for f in component_forces]
                combined = projected_force(total_forces, indices, gradient, norm)
                if not math.isclose(sum(projections), combined, rel_tol=1e-10, abs_tol=1e-6):
                    raise ValueError('projected component sum differs from total')
                measurement = stereocentre_geometry(*(rows[i]['xyz'] for i in indices))
                guard = guard_parameters(system, indices) if centre['names'][0] == 'CA' else None
                if guard is not None:
                    q = volume / guard['vref_nm3']
                    deficit = max(guard['qmin']-q, 0)
                    predicted = guard['k'] * deficit / guard['vref_nm3'] * norm
                    if not math.isclose(projections[guard['group']], predicted, rel_tol=1e-8, abs_tol=1e-6):
                        raise ValueError('C-alpha guard force differs from analytic derivative')
                    guard.update(q=q, active=deficit > 0, energy_kj_mol=0.5*guard['k']*deficit**2,
                                 analytic_projection_kj_mol_nm=predicted)
                centres.append(dict(centre, indices=indices, volume_angstrom3=volume*1000,
                    angle_degrees=measurement.out_of_plane_degrees if measurement else None,
                    total_projection_kj_mol_nm=combined, component_projections_kj_mol_nm=projections, guard=guard))
            results.append(dict(phase=phase, side=side, total_energy_kj_mol=total_energy,
                                force_components=components, centres=centres))
        return results
    finally:
        del context
        del integrator


def plan(parent_path, run_root, historical_root, output):
    import openmm

    parent_raw = parent_path.read_bytes()
    parent = json.loads(parent_raw)
    historical_verify(parent, historical_root)
    run_raw = (run_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(parent_raw) or (run_root/'executed-plan.json').read_bytes() != parent_raw:
        raise ValueError('cohort execution binding differs')
    expected = [(c['id'], s, a) for c in parent['cases'] for s in parent['protocol']['seeds'] for a in parent['protocol']['arms']]
    if [(r['id'], r['seed'], r['arm']) for r in run['cases']] != expected:
        raise ValueError('incomplete original cohort enrollment')
    inputs = {str(parent_path): sha(parent_raw), str(run_root/'evidence.json'): sha(run_raw)}
    for row in run['cases']:
        folder = run_root/row['id']/f"seed-{row['seed']}"/row['arm']
        for name, digest in row['artifacts'].items():
            path = folder/name
            if sha(path.read_bytes()) != digest:
                raise ValueError('saved artifact changed before enrollment')
            inputs[str(path)] = digest
    sources = ['scripts/audit_stereo_force_projection.py', 'scripts/audit_legacy_start_controls.py',
               'scripts/research_dunbrack_prior.py', 'src/dvbfixer/diagnose/geometry.py',
               'src/dvbfixer/diagnose/chemistry.py', 'src/dvbfixer/minimize/pipeline.py']
    write(output, dict(schema='dvbfixer.stereo-force-projection-plan.v1', cases=parent['cases'],
        executions=run['cases'], run_root=str(run_root), historical_root=str(historical_root),
        parent_plan=str(parent_path), inputs=inputs, source_hashes={p: sha(Path(p).read_bytes()) for p in sources},
        openmm_version=openmm.__version__,
        analysis_git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        historical_git_head=subprocess.check_output(['git', '-C', str(historical_root), 'rev-parse', 'HEAD'], text=True).strip(),
        policy=dict(timing='retrospective; known endpoint failures, frozen before force evaluation',
            enrollment='ALL24 saved cohort cells, both sides of both phases; target CA, THR beta only, and known GLN16 context CA',
            platform='Reference static evaluations only, no integration/minimization/constraints projection',
            projection='Cartesian force dot unit volume gradient; positive increases ordered signed volume, not a trajectory or causal proof',
            groups='clone only; unique force groups with exact serialization restoration and additive force/energy checks',
            parameters='recorded initial globals must equal each phase system defaults; original pipeline contains no setParameter calls',
            limits='no new seeds/weights/H/native scoring/profile choice/coordinate repair; preserve all failures'),
        production_changed=False, native_coordinates_read=False))


def run(plan_path, output):
    import openmm

    record = json.loads(plan_path.read_bytes())
    verify(record)
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('OpenMM version changed since analysis enrollment')
    historical_verify(json.loads(Path(record['parent_plan']).read_bytes()), Path(record['historical_root']))
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    cases = {c['id']: c for c in record['cases']}
    results = []
    for execution in record['executions']:
        row = {k: execution[k] for k in ('id', 'seed', 'arm', 'status')}
        folder = Path(record['run_root'])/row['id']/f"seed-{row['seed']}"/row['arm']
        if row['status'] == 'completed':
            row['states'] = []
            try:
                for phase in (1, 2):
                    row['states'].extend(phase_readback(folder, phase, cases[row['id']]))
                row['status'] = 'evaluated'
            except Exception as error:
                row.update(status='readback_blocker', error_type=type(error).__name__, error=str(error))
        results.append(row)
        write(output/f"cell-{len(results):02d}.json", row)
    verify(record)
    report = dict(schema='dvbfixer.stereo-force-projection.v1', cases=results,
        plan_sha256=sha(plan_path.read_bytes()), static_evaluation_only=True, production_changed=False,
        native_coordinates_read=False, scientific_qualification=False,
        summary=dict(enrolled=len(results), evaluated=sum(r['status']=='evaluated' for r in results),
            readback_blockers=sum(r['status']=='readback_blocker' for r in results),
            states=sum(len(r.get('states', [])) for r in results),
            centres=sum(len(s['centres']) for r in results for s in r.get('states', []))))
    write(output/'evidence.json', report)
    print(json.dumps(report['summary'], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--parent', type=Path, required=True)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--historical-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.parent, args.run, args.historical_root, args.output)
    else:
        run(args.plan, args.output)
