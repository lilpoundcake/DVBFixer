"""Independent grid/geometry/native-gradient audit of paired methyl-H starts."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_ha_minimizer_path import outward_frame, projection, sphere_tangent
from scripts.audit_ha_start_control import ha_geometry
from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.audit_stereo_force_projection import grouped_system
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def rotation_invariants(rows, before, after, degrees):
    mapped = {tuple(r['key']): i for i, r in enumerate(rows)}
    if len(mapped) != len(rows) or len(before) != len(rows) or len(after) != len(rows):
        raise ValueError('methyl rotation inventory differs')
    required = [('A', '24', '', 'THR', n) for n in ('CG2', 'CB', 'HG21', 'HG22', 'HG23')]
    if any(key not in mapped for key in required):
        raise ValueError('methyl identity missing')
    indices = [mapped[key] for key in required]
    moved = set(indices[2:])
    if any(before[i] != after[i] for i in range(len(rows)) if i not in moved):
        raise ValueError('methyl rotation moved another starting coordinate')
    origin, end = (before[i] for i in indices[:2])
    axis = [x-y for x, y in zip(end, origin)]
    norm = math.hypot(*axis)
    if norm <= 1e-12:
        raise ValueError('undefined methyl axis')
    axis = [x/norm for x in axis]
    for i in moved:
        first, second = ([x-y for x, y in zip(points[i], origin)] for points in (before, after))
        if not math.isclose(math.hypot(*first), math.hypot(*second), rel_tol=1e-11, abs_tol=1e-13):
            raise ValueError('methyl rotation changed C-H distance')
        dots = [sum(x*y for x, y in zip(v, axis)) for v in (first, second)]
        if not math.isclose(*dots, rel_tol=1e-11, abs_tol=1e-13):
            raise ValueError('methyl rotation changed axial coordinate/bond angle')
        a, b = [[x-d*y for x, y in zip(v, axis)] for v, d in zip((first, second), dots)]
        cross = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
        angle = math.degrees(math.atan2(sum(x*y for x, y in zip(axis, cross)), sum(x*y for x, y in zip(a, b))))
        if abs((angle-degrees+180)%360-180) > 1e-8:
            raise ValueError('methyl orientation differs from declared angle')
        for j in moved:
            if not math.isclose(math.dist(before[i], before[j]), math.dist(after[i], after[j]), rel_tol=1e-11, abs_tol=1e-13):
                raise ValueError('methyl rotation changed internal H geometry')
    return indices


def choose_independently(candidates):
    if [c['angle_degrees'] for c in candidates] != list(range(0, 360, 15)):
        raise ValueError('incomplete or reordered initial energy grid')
    if any(not math.isfinite(c['energy_kj_mol']) for c in candidates):
        raise ValueError('nonfinite grid energy')
    best = min(c['energy_kj_mol'] for c in candidates)
    eligible = [c['angle_degrees'] for c in candidates if c['energy_kj_mol']-best <= 1e-8]
    def signed(angle):
        return angle if angle <= 180 else angle-360

    return min(eligible, key=lambda a: (abs(signed(a)), signed(a)))


def h_bonds(system, rows, points_nm):
    import openmm
    from openmm.unit import nanometer

    values = []
    for force in system.getForces():
        if not isinstance(force, openmm.HarmonicBondForce):
            continue
        for j in range(force.getNumBonds()):
            a, b, reference, _ = force.getBondParameters(j)
            if rows[a]['element'] not in {'H', 'D'} and rows[b]['element'] not in {'H', 'D'}:
                continue
            observed = math.dist(points_nm[a], points_nm[b])
            deviation = 10*(observed-reference.value_in_unit(nanometer))
            values.append(dict(atoms=[rows[a]['key'], rows[b]['key']], length_angstrom=10*observed,
                deviation_angstrom=deviation))
    return dict(count=len(values), maximum_absolute_deviation_angstrom=max((abs(v['deviation_angstrom']) for v in values), default=None),
        rms_deviation_angstrom=math.sqrt(sum(v['deviation_angstrom']**2 for v in values)/len(values)) if values else None,
        largest_deviations=sorted(values, key=lambda v: -abs(v['deviation_angstrom']))[:3])


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw:
        raise ValueError('methyl-start executed plan differs')
    expected = [(c['id'], c['seed'], c['source_arm'], a) for c in plan['cases'] for a in plan['policy']['arms']]
    if [(r['id'], r['seed'], r['source_arm'], r['arm']) for r in run['cases']] != expected:
        raise ValueError('methyl-start cohort differs')
    declarations = {(c['id'], c['seed'], c['source_arm']): c for c in plan['cases']}
    results, gradient_errors = [], []
    for record in run['cases']:
        case = declarations[(record['id'], record['seed'], record['source_arm'])]
        folder = root/case['id']/f"seed-{case['seed']}"/case['source_arm']/record['arm']
        for name, digest in record['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('retained methyl-start artifact changed')
        report = {k: record[k] for k in ('id', 'seed', 'source_arm', 'arm', 'status', 'returncode')}
        if record['status'] != 'completed':
            results.append(report)
            continue
        source = Path(case['source'])
        template = load_atoms(source/'h-internal-before.json.gz')
        original = json.loads(gzip.decompress(Path(case['path']).read_bytes()))['initial_nm']
        data = json.loads(gzip.decompress((folder/'path.json.gz').read_bytes()))
        trace = json.loads((folder/'trace.json').read_bytes())
        if trace['plan_sha256'] != sha(raw) or trace['case'] != case or trace['arm'] != record['arm'] or trace['whole_system_minimization']:
            raise ValueError('methyl-start trace contract differs')
        system = openmm.XmlSerializer.deserialize((source/'h-internal-system.xml').read_text())
        if sha(openmm.XmlSerializer.serialize(system).encode()) != trace['serialized_system_sha256']:
            raise ValueError('methyl-start force/parameter/mass identity differs')
        choice = choose_independently(trace['candidates'])
        if choice != trace['selected_degrees'] or trace['applied_degrees'] != (0 if record['arm'] == 'unchanged' else choice):
            raise ValueError('methyl selection differs from frozen initial-energy rule')
        indices = rotation_invariants(template, original, data['initial_nm'], trace['applied_degrees'])
        applied_candidate = next(c for c in trace['candidates'] if c['angle_degrees'] == trace['applied_degrees'])
        if [data['initial_nm'][i] for i in indices[2:]] != applied_candidate['H_points_nm']:
            raise ValueError('applied methyl coordinates differ from the recorded selected grid candidate')
        for label in ('initial', 'final'):
            if len(data[label+'_atoms']) != len(template) or any(row['key'] != expected['key'] or row['element'] != expected['element'] or row['xyz'] != [10*x for x in point]
                    for row, expected, point in zip(data[label+'_atoms'], template, data[label+'_nm'], strict=True)):
                raise ValueError('methyl-control nm/atom endpoint inventories differ')
        if not data['callbacks'] or len(data['callbacks']) != trace['callbacks'] or len(data['callbacks']) > plan['policy']['maximum_callbacks'] or data['callbacks'][-1]['positions_nm'] != [x for p in data['final_nm'] for x in p]:
            raise ValueError('methyl native callback/final endpoint differs')
        mobile = [i for i in range(len(template)) if system.getParticleMass(i).value_in_unit(dalton) > 0]
        heavy = [i for i, r in enumerate(template) if r['element'] not in {'H', 'D'}]
        if any(i in mobile for i in heavy):
            raise ValueError('H-only system has mobile heavy particles')
        grouped, groups = grouped_system(system)
        nb_groups = [g['index'] for g in groups if g['type'] == 'NonbondedForce']
        if len(nb_groups) != 1:
            raise ValueError('unique nonbonded group missing')
        integrator, grouped_integrator = openmm.VerletIntegrator(0.), openmm.VerletIntegrator(0.)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        grouped_context = openmm.Context(grouped, grouped_integrator, openmm.Platform.getPlatformByName('Reference'))
        try:
            for candidate in trace['candidates']:
                points = [list(p) for p in original]
                for i, p in zip(indices[2:], candidate['H_points_nm'], strict=True):
                    points[i] = p
                rotation_invariants(template, original, points, candidate['angle_degrees'])
                context.setPositions(points*nanometer)
                energy = context.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                if not math.isclose(energy, candidate['energy_kj_mol'], rel_tol=1e-12, abs_tol=1e-8):
                    raise ValueError('recorded candidate energy differs from independent native readback')
            centre, outward = outward_frame(template)
            ha_index = next(i for i, r in enumerate(template) if r['key'] == ['A', '16', '', 'GLN', 'HA'])
            nodes = [('initial', None, data['initial_nm'], None)]
            for ordinal, cb in enumerate(data['callbacks']):
                if cb['ordinal'] != ordinal or len(cb['positions_nm']) != 3*len(template):
                    raise ValueError('native callback inventory differs')
                points = [cb['positions_nm'][3*i:3*i+3] for i in range(len(template))]
                nodes.append(('callback', ordinal, points, cb))
            nodes.append(('final', None, data['final_nm'], None))
            published = []
            for label, ordinal, points, cb in nodes:
                if any(points[i] != original[i] for i in heavy):
                    raise ValueError('H-stage rotated/moved an existing heavy observation')
                rows = [dict(r, xyz=[10*x for x in p]) for r, p in zip(template, points, strict=True)]
                context.setPositions(points*nanometer)
                state = context.getState(energy=True, forces=True)
                energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                forces = state.getForces().value_in_unit(kilojoules_per_mole/nanometer)
                error = None
                if cb is not None:
                    gradient = cb['gradient_kj_mol_nm']
                    if len(gradient) != 3*len(template):
                        raise ValueError('native gradient inventory differs')
                    differences = [abs(forces[i][j]+gradient[3*i+j]) for i in mobile for j in range(3)]
                    if any(not math.isclose(forces[i][j], -gradient[3*i+j], rel_tol=1e-7, abs_tol=1e-4) for i in mobile for j in range(3)):
                        raise ValueError('native mobile-atom gradient differs from physical force')
                    if cb['statistics']['restraint energy'] != 0 or cb['statistics']['max constraint error'] != 0 or not math.isclose(cb['statistics']['system energy'], energy, rel_tol=1e-8, abs_tol=1e-4):
                        raise ValueError('native callback objective differs from physical energy')
                    error = max(differences, default=0.)
                    gradient_errors.append(error)
                cosine, _ = sphere_tangent(centre, rows[ha_index]['xyz'], outward)
                item = dict(label=label, ordinal=ordinal, geometry=ha_geometry(rows), outward_cosine=cosine,
                    energy_kj_mol=energy, HA_outward_force_kj_mol_nm=projection(forces[ha_index], outward),
                    maximum_mobile_gradient_error_kj_mol_nm=error)
                if label != 'callback':
                    grouped_context.setPositions(points*nanometer)
                    nb_force = grouped_context.getState(forces=True, groups=1 << nb_groups[0]).getForces().value_in_unit(kilojoules_per_mole/nanometer)[ha_index]
                    item.update(HA_nonbonded_outward_force_kj_mol_nm=projection(nb_force, outward),
                        methyl_HA_distances_angstrom={template[i]['key'][-1]: math.dist(rows[i]['xyz'], rows[ha_index]['xyz']) for i in indices[2:]},
                        H_bonds=h_bonds(system, template, points))
                published.append(item)
            baseline = coordinate_delta(load_atoms(source/'h-internal-after.json.gz'), data['final_atoms'])
            if record['arm'] == 'unchanged' and any(not baseline[k]['exact'] for k in ('heavy', 'hydrogen')):
                raise ValueError('unchanged methyl control not exactly reproduced')
            report.update(applied_degrees=trace['applied_degrees'], selected_degrees=choice,
                candidates=[{k: c[k] for k in ('angle_degrees', 'energy_kj_mol')} for c in trace['candidates']],
                callbacks=len(data['callbacks']), nodes=published, exact_nonmethyl_start=True,
                exact_heavy_path=True, baseline_endpoint_difference=baseline,
                any_inward_callback=any(s['outward_cosine'] < 0 for s in published if s['label'] == 'callback'),
                any_all_acute_callback=any(max(s['geometry']['ha_angles_degrees'].values()) < 90 for s in published if s['label'] == 'callback'))
        finally:
            del context
            del grouped_context
            del integrator
            del grouped_integrator
        results.append(report)
    completed = [r for r in results if r['status'] == 'completed']
    summary = dict(enrolled=len(results), statuses=dict(Counter(r['status'] for r in results)), arms={})
    for arm in plan['policy']['arms']:
        cells = [r for r in completed if r['arm'] == arm]
        summary['arms'][arm] = dict(completed=len(cells), any_inward_callback=sum(r['any_inward_callback'] for r in cells),
            any_all_acute_callback=sum(r['any_all_acute_callback'] for r in cells),
            final_all_acute=sum(max(r['nodes'][-1]['geometry']['ha_angles_degrees'].values()) < 90 for r in cells))
        if arm == 'unchanged':
            summary['arms'][arm]['exact_original_endpoints'] = sum(all(r['baseline_endpoint_difference'][k]['exact'] for k in ('heavy', 'hydrogen')) for r in cells)
    summary['maximum_mobile_gradient_error_kj_mol_nm'] = max(gradient_errors, default=None)
    write(output, dict(schema='dvbfixer.methyl-start-control-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        driver_sha256=plan['source_hashes'][str(Path('scripts/research_methyl_start_control.py').resolve())],
        protocol=plan['policy'], openmm_version=plan['openmm_version'], production_changed=False,
        native_coordinates_read=False, scientific_qualification=False, whole_system_minimization=False,
        original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
