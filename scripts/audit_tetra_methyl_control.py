"""Independent canonical-angle/grid/gradient audit of all tetra-methyl controls."""
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
from scripts.audit_methyl_start_control import choose_independently, h_bonds
from scripts.audit_stereo_force_projection import grouped_system
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def normalize(vector):
    length = math.hypot(*vector)
    if length <= 1e-12 or not all(math.isfinite(x) for x in vector):
        raise ValueError('undefined tetrahedral frame/vector')
    return [x/length for x in vector]


def tetrahedral_invariants(rows, before, after, phase):
    mapped = {tuple(r['key']): i for i, r in enumerate(rows)}
    required = [('A', '24', '', 'THR', n) for n in ('CG2', 'CB', 'HG21', 'HG22', 'HG23', 'CA')]
    if len(mapped) != len(rows) or any(k not in mapped for k in required) or len(before) != len(rows) or len(after) != len(rows):
        raise ValueError('tetrahedral atom inventory differs')
    indices = [mapped[k] for k in required]
    owned = set(indices[2:5])
    if any(before[i] != after[i] for i in range(len(rows)) if i not in owned):
        raise ValueError('tetrahedral start changed an unowned coordinate')
    if phase is None:
        if before != after:
            raise ValueError('unchanged tetrahedral control moved a coordinate')
        return indices
    if phase not in range(0, 360, 15):
        raise ValueError('undeclared tetrahedral phase')
    origin = before[indices[0]]
    axis = normalize([x-y for x, y in zip(before[indices[1]], origin)])
    ca = [x-y for x, y in zip(before[indices[-1]], origin)]
    normal = normalize(cross(axis, ca))
    first = normalize(cross(normal, axis))
    directions = []
    for ordinal, i in enumerate(indices[2:5]):
        if not math.isclose(math.dist(before[i], origin), math.dist(after[i], origin), rel_tol=1e-11, abs_tol=1e-13):
            raise ValueError('tetrahedral start changed original C-H radius')
        direction = normalize([x-y for x, y in zip(after[i], origin)])
        if not math.isclose(sum(x*y for x, y in zip(direction, axis)), -1/3, abs_tol=1e-10):
            raise ValueError('methyl H-C-CB angle is not canonical tetrahedral')
        azimuth = math.degrees(math.atan2(sum(x*y for x, y in zip(direction, normal)), sum(x*y for x, y in zip(direction, first))))
        if abs((azimuth-phase-120*ordinal+180)%360-180) > 1e-8:
            raise ValueError('canonical methyl azimuth differs from declared phase')
        directions.append(direction)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        if not math.isclose(sum(x*y for x, y in zip(directions[i], directions[j])), -1/3, abs_tol=1e-10):
            raise ValueError('methyl H-C-H angle is not canonical tetrahedral')
    return indices


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['source_arm'], a) for c in plan['cases'] for a in plan['policy']['arms']]
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw or [(c['id'], c['seed'], c['source_arm'], c['arm']) for c in run['cases']] != expected:
        raise ValueError('tetra-methyl executed enrollment differs')
    by_identity = {(c['id'], c['seed'], c['source_arm']): c for c in plan['cases']}
    results, errors = [], []
    for record in run['cases']:
        case = by_identity[(record['id'], record['seed'], record['source_arm'])]
        folder = root/case['id']/f"seed-{case['seed']}"/case['source_arm']/record['arm']
        for name, digest in record['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('tetra-methyl retained artifact changed')
        report = {k: record[k] for k in ('id', 'seed', 'source_arm', 'arm', 'status', 'returncode')}
        if record['status'] != 'completed':
            results.append(report)
            continue
        source = Path(case['source'])
        rows = load_atoms(source/'h-internal-before.json.gz')
        original = json.loads(gzip.decompress(Path(case['path']).read_bytes()))['initial_nm']
        path = json.loads(gzip.decompress((folder/'path.json.gz').read_bytes()))
        trace = json.loads((folder/'trace.json').read_bytes())
        choice = choose_independently(trace['candidates'])
        phase = None if record['arm'] == 'unchanged' else 0 if record['arm'] == 'tetrahedral_zero' else choice
        if trace['plan_sha256'] != sha(raw) or trace['case'] != case or trace['arm'] != record['arm'] or trace['selected_degrees'] != choice or trace['applied_degrees'] != phase or trace['whole_system_minimization']:
            raise ValueError('tetra-methyl selection/trace contract differs')
        indices = tetrahedral_invariants(rows, original, path['initial_nm'], phase)
        original_directions = [normalize([x-y for x, y in zip(original[i], original[indices[0]])]) for i in indices[2:5]]
        original_angles = [math.degrees(math.acos(max(-1., min(1., sum(x*y for x, y in zip(original_directions[i], original_directions[j])))))) for i, j in ((0, 1), (0, 2), (1, 2))]
        if any(not math.isclose(a, b, abs_tol=1e-8) for a, b in zip(original_angles, case['original_methyl_HCH_degrees'], strict=True)) or any(not math.isclose(a, math.degrees(math.acos(-1/3)), abs_tol=1e-8) for a in case['canonical_methyl_HCH_degrees']):
            raise ValueError('predeclared methyl angular inventories differ from source/canonical geometry')
        if phase is not None:
            candidate = next(c for c in trace['candidates'] if c['angle_degrees'] == phase)
            if [path['initial_nm'][i] for i in indices[2:5]] != candidate['H_points_nm']:
                raise ValueError('applied canonical coordinates differ from logged candidate')
        for label in ('initial', 'final'):
            if any(a['key'] != r['key'] or a['element'] != r['element'] or a['xyz'] != [10*x for x in p] for a, r, p in zip(path[label+'_atoms'], rows, path[label+'_nm'], strict=True)):
                raise ValueError('canonical nm/atom endpoints differ')
        if not path['callbacks'] or len(path['callbacks']) != trace['callbacks'] or len(path['callbacks']) > 512 or path['callbacks'][-1]['positions_nm'] != [x for p in path['final_nm'] for x in p]:
            raise ValueError('canonical native callback inventory differs')
        system = openmm.XmlSerializer.deserialize((source/'h-internal-system.xml').read_text())
        if sha(openmm.XmlSerializer.serialize(system).encode()) != trace['serialized_system_sha256']:
            raise ValueError('tetrahedral start changed original System/parameters/masses')
        mobile = [i for i in range(len(rows)) if system.getParticleMass(i).value_in_unit(dalton) > 0]
        heavy = [i for i, r in enumerate(rows) if r['element'] not in {'H', 'D'}]
        if set(mobile).intersection(heavy):
            raise ValueError('canonical control has mobile heavy particles')
        grouped, groups = grouped_system(system)
        nb = [g['index'] for g in groups if g['type'] == 'NonbondedForce']
        if len(nb) != 1:
            raise ValueError('unique native nonbonded force missing')
        integrator, grouped_integrator = openmm.VerletIntegrator(0.), openmm.VerletIntegrator(0.)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        grouped_context = openmm.Context(grouped, grouped_integrator, openmm.Platform.getPlatformByName('Reference'))
        try:
            context.setPositions(original*nanometer)
            original_energy = context.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            if not math.isclose(original_energy, trace['original_energy_kj_mol'], rel_tol=1e-12, abs_tol=1e-8):
                raise ValueError('canonical original-energy readback differs')
            for candidate in trace['candidates']:
                points = [list(p) for p in original]
                for i, p in zip(indices[2:5], candidate['H_points_nm'], strict=True):
                    points[i] = p
                tetrahedral_invariants(rows, original, points, candidate['angle_degrees'])
                context.setPositions(points*nanometer)
                energy = context.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                if not math.isclose(energy, candidate['energy_kj_mol'], rel_tol=1e-12, abs_tol=1e-8):
                    raise ValueError('canonical candidate energy differs from independent native readback')
            centre, outward = outward_frame(rows)
            ha = next(i for i, r in enumerate(rows) if r['key'] == ['A', '16', '', 'GLN', 'HA'])
            nodes = [('initial', None, path['initial_nm'], None)]
            for ordinal, cb in enumerate(path['callbacks']):
                if cb['ordinal'] != ordinal or len(cb['positions_nm']) != 3*len(rows) or len(cb['gradient_kj_mol_nm']) != 3*len(rows):
                    raise ValueError('canonical callback particle inventory differs')
                points = [cb['positions_nm'][3*i:3*i+3] for i in range(len(rows))]
                nodes.append(('callback', ordinal, points, cb))
            nodes.append(('final', None, path['final_nm'], None))
            published = []
            for label, ordinal, points, cb in nodes:
                if any(points[i] != original[i] for i in heavy):
                    raise ValueError('canonical H-stage moved an original heavy coordinate')
                context.setPositions(points*nanometer)
                state = context.getState(energy=True, forces=True)
                energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                forces = state.getForces().value_in_unit(kilojoules_per_mole/nanometer)
                error = None
                if cb is not None:
                    gradient = cb['gradient_kj_mol_nm']
                    difference = [abs(forces[i][j]+gradient[3*i+j]) for i in mobile for j in range(3)]
                    if any(not math.isclose(forces[i][j], -gradient[3*i+j], rel_tol=1e-7, abs_tol=1e-4) for i in mobile for j in range(3)) or cb['statistics']['restraint energy'] != 0 or cb['statistics']['max constraint error'] != 0 or not math.isclose(cb['statistics']['system energy'], energy, rel_tol=1e-8, abs_tol=1e-4):
                        raise ValueError('canonical callback gradient/objective differs from native physical readback')
                    error = max(difference, default=0.)
                    errors.append(error)
                material = [dict(r, xyz=[10*x for x in p]) for r, p in zip(rows, points, strict=True)]
                cosine, _ = sphere_tangent(centre, material[ha]['xyz'], outward)
                item = dict(label=label, ordinal=ordinal, geometry=ha_geometry(material), outward_cosine=cosine,
                    energy_kj_mol=energy, HA_outward_force_kj_mol_nm=projection(forces[ha], outward),
                    maximum_mobile_gradient_error_kj_mol_nm=error)
                if cb is None:
                    grouped_context.setPositions(points*nanometer)
                    nb_force = grouped_context.getState(forces=True, groups=1 << nb[0]).getForces().value_in_unit(kilojoules_per_mole/nanometer)[ha]
                    item.update(HA_nonbonded_outward_force_kj_mol_nm=projection(nb_force, outward),
                        methyl_HA_distances_angstrom={rows[i]['key'][-1]: math.dist(material[i]['xyz'], material[ha]['xyz']) for i in indices[2:5]},
                        H_bonds=h_bonds(system, rows, points))
                published.append(item)
            baseline = coordinate_delta(load_atoms(source/'h-internal-after.json.gz'), path['final_atoms'])
            if record['arm'] == 'unchanged' and any(not baseline[k]['exact'] for k in ('heavy', 'hydrogen')):
                raise ValueError('canonical unchanged control is not exact')
            report.update(selected_degrees=choice, applied_degrees=phase,
                original_methyl_HCH_degrees=case['original_methyl_HCH_degrees'],
                canonical_methyl_HCH_degrees=case['canonical_methyl_HCH_degrees'],
                candidates=[{k: c[k] for k in ('angle_degrees', 'energy_kj_mol')} for c in trace['candidates']],
                callbacks=len(path['callbacks']), nodes=published, exact_nonmethyl_start=True,
                exact_heavy_path=True, baseline_endpoint_difference=baseline,
                any_inward_callback=any(s['outward_cosine'] < 0 for s in published if s['label'] == 'callback'),
                any_all_acute_callback=any(max(s['geometry']['ha_angles_degrees'].values()) < 90 for s in published if s['label'] == 'callback'))
        finally:
            del context
            del grouped_context
            del integrator
            del grouped_integrator
        results.append(report)
    summary = dict(enrolled=len(results), statuses=dict(Counter(r['status'] for r in results)), arms={})
    for arm in plan['policy']['arms']:
        cells = [r for r in results if r['status'] == 'completed' and r['arm'] == arm]
        summary['arms'][arm] = dict(completed=len(cells), any_inward_callback=sum(r['any_inward_callback'] for r in cells),
            any_all_acute_callback=sum(r['any_all_acute_callback'] for r in cells),
            final_all_acute=sum(max(r['nodes'][-1]['geometry']['ha_angles_degrees'].values()) < 90 for r in cells))
        if arm == 'unchanged':
            summary['arms'][arm]['exact_original_endpoints'] = sum(all(r['baseline_endpoint_difference'][k]['exact'] for k in ('heavy', 'hydrogen')) for r in cells)
    summary['maximum_mobile_gradient_error_kj_mol_nm'] = max(errors, default=None)
    write(output, dict(schema='dvbfixer.tetra-methyl-start-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        driver_sha256=plan['source_hashes'][str(Path('scripts/research_tetra_methyl_control.py').resolve())],
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
