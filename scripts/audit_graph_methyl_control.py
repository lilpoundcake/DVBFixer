"""Independent native-energy/geometry/gradient audit of every graph methyl arm."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_methyl_start_control import h_bonds
from scripts.audit_methyl_whole_control import hydrogen_angles
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def unit(v):
    radius = math.hypot(*v)
    if radius <= 1e-12 or not math.isfinite(radius):
        raise ValueError('undefined graph-audit frame')
    return [x/radius for x in v]


def independently_place(original, group, phase):
    origin = original[group['parent']]
    axis = unit([x-y for x, y in zip(original[group['axis']], origin)])
    reference = [x-y for x, y in zip(original[group['frame']], origin)]
    normal = unit(cross(axis, reference))
    first = unit(cross(normal, axis))
    owned = []
    for ordinal, i in enumerate(group['hydrogens']):
        radius = math.dist(original[i], origin)
        azimuth = math.radians(phase+120*ordinal)
        direction = [-a/3+math.sqrt(8/9)*(math.cos(azimuth)*u+math.sin(azimuth)*v) for a, u, v in zip(axis, first, normal)]
        owned.append([o+radius*d for o, d in zip(origin, direction)])
    return owned


def select_grid(grid):
    if [r['phase'] for r in grid] != list(range(0, 360, 15)) or any(not math.isfinite(r['energy_kj_mol']) for r in grid):
        raise ValueError('incomplete/nonfinite/reordered graph energy grid')
    minimum = min(r['energy_kj_mol'] for r in grid)
    eligible = [r['phase'] for r in grid if r['energy_kj_mol'] <= minimum+1e-8]
    signed = [a if a <= 180 else a-360 for a in eligible]
    return min(signed, key=lambda a: (abs(a), a))%360


def ha_geometry(rows, bonds, points):
    adjacent = [set() for _ in rows]
    for a, b in bonds:
        adjacent[a].add(b)
        adjacent[b].add(a)
    observations = []
    for i, row in enumerate(rows):
        if row['key'][-1] != 'CA' or row['element'] != 'C':
            continue
        heavy = [j for j in adjacent[i] if rows[j]['element'] not in {'H', 'D'}]
        hydrogens = [j for j in adjacent[i] if rows[j]['element'] in {'H', 'D'}]
        for h in hydrogens:
            direction = unit([x-y for x, y in zip(points[h], points[i])])
            angles = []
            for j in heavy:
                vector = unit([x-y for x, y in zip(points[j], points[i])])
                cosine = sum(x*y for x, y in zip(direction, vector))
                angles.append(math.degrees(math.acos(max(-1., min(1., cosine)))))
            observations.append(dict(atom=rows[h]['key'], radius_angstrom=10*math.dist(points[h], points[i]),
                heavy_neighbor_angles_degrees=angles,
                all_acute_three_heavy=len(heavy) == 3 and max(angles) < 90))
    return dict(HA_count=len(observations), all_acute=[v for v in observations if v['all_acute_three_heavy']],
        minimum_CA_H_angstrom=min((v['radius_angstrom'] for v in observations), default=None),
        maximum_CA_H_angstrom=max((v['radius_angstrom'] for v in observations), default=None))


def h_close_pairs(rows, bonds, points):
    """Descriptive <1 Å H-containing contacts, excluding bonded/angle pairs."""
    from scipy.spatial import cKDTree

    adjacent = [set() for _ in rows]
    for a, b in bonds:
        adjacent[a].add(b)
        adjacent[b].add(a)
    pairs = []
    for a, b in sorted(cKDTree(points).query_pairs(.1)):
        if (rows[a]['element'] not in {'H', 'D'} and rows[b]['element'] not in {'H', 'D'}) or b in adjacent[a] or adjacent[a] & adjacent[b]:
            continue
        pairs.append(dict(atoms=[rows[a]['key'], rows[b]['key']], distance_angstrom=10*math.dist(points[a], points[b])))
    return dict(below_1_angstrom=len(pairs), below_0_1_angstrom=sum(p['distance_angstrom'] < .1 for p in pairs),
        closest_pairs=sorted(pairs, key=lambda p: p['distance_angstrom'])[:5],
        scope='geometric descriptor only, not a new acceptance gate; directly bonded and common-neighbor angle pairs excluded')


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], arm) for c in plan['cases'] for arm in plan['policy']['arms']]
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw or [(c['id'], c['seed'], c['arm']) for c in run['cases']] != expected:
        raise ValueError('complete graph study enrollment differs')
    cases = {(c['id'], c['seed']): c for c in plan['cases']}
    results, grid_errors, gradient_errors = [], [], []
    for cell in run['cases']:
        source_case = cases[(cell['id'], cell['seed'])]
        source = Path(source_case['source'])
        folder = root/cell['id']/f"seed-{cell['seed']}"/cell['arm']
        for name, digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('graph study retained artifact changed')
        report = {k: cell[k] for k in ('id', 'seed', 'arm', 'status', 'returncode')}
        if cell['status'] != 'completed':
            results.append(report)
            continue
        trace = json.loads((folder/'trace.json').read_bytes())
        path = json.loads(gzip.decompress((folder/'path.json.gz').read_bytes()))
        original = json.loads(gzip.decompress((source/'H-path.json.gz').read_bytes()))
        ownership = json.loads((source/'H-ownership.json').read_bytes())
        rows, bonds, groups = ownership['rows'], ownership['bonds'], ownership['graph']['groups']
        if trace['plan_sha256'] != sha(raw) or trace['case'] != source_case or trace['arm'] != cell['arm'] or trace['groups'] != len(groups):
            raise ValueError('graph trace contract differs')
        original_points = original['initial_nm']
        system = openmm.XmlSerializer.deserialize((source/'H-system.xml').read_text())
        if sha(openmm.XmlSerializer.serialize(system).encode()) != trace['system_sha256']:
            raise ValueError('graph control System changed')
        mobile = [i for i in range(len(rows)) if system.getParticleMass(i).value_in_unit(dalton) > 0]
        heavy = [i for i, r in enumerate(rows) if r['element'] not in {'H', 'D'}]
        if set(mobile) & set(heavy):
            raise ValueError('mobile heavy atom in H System')
        owned = {i for g in groups for i in g['hydrogens']}
        if any(path['initial_nm'][i] != original_points[i] for i in range(len(rows)) if i not in owned):
            raise ValueError('graph initialization changed unowned coordinates')
        integrator = openmm.VerletIntegrator(0.)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        try:
            independently_expected = [p[:] for p in original_points]
            if cell['arm'] != 'unchanged':
                for group in groups:
                    for i, p in zip(group['hydrogens'], independently_place(original_points, group, 0), strict=True):
                        independently_expected[i] = p
            selected = cell['arm'] == 'graph_tetrahedral_energy_selected'
            if len(trace['grids']) != (len(groups) if selected else 0):
                raise ValueError('graph energy pass count differs')
            for group, grid in zip(groups if selected else [], trace['grids'], strict=True):
                chosen = select_grid(grid['energies'])
                if grid['parent_key'] != group['parent_key'] or chosen != grid['selected_phase']:
                    raise ValueError('graph selection/tie/parent order differs')
                for candidate in grid['energies']:
                    trial = [p[:] for p in independently_expected]
                    for i, p in zip(group['hydrogens'], independently_place(original_points, group, candidate['phase']), strict=True):
                        trial[i] = p
                    context.setPositions(trial*nanometer)
                    actual = context.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                    error = abs(actual-candidate['energy_kj_mol'])
                    grid_errors.append(error)
                    if not math.isclose(actual, candidate['energy_kj_mol'], rel_tol=1e-10, abs_tol=1e-5):
                        raise ValueError('independent original initial-energy grid differs')
                for i, p in zip(group['hydrogens'], independently_place(original_points, group, chosen), strict=True):
                    independently_expected[i] = p
            if len(path['initial_nm']) != len(rows) or any(math.dist(a, b) > 1e-12 for a, b in zip(path['initial_nm'], independently_expected, strict=True)):
                raise ValueError('applied graph geometry differs from independent canonical/selection replay')
            nodes = [('initial', path['initial_nm'], None)]
            for ordinal, cb in enumerate(path['callbacks']):
                if cb['ordinal'] != ordinal or len(cb['positions_nm']) != 3*len(rows) or len(cb['gradient_kj_mol_nm']) != 3*len(rows):
                    raise ValueError('graph native callback inventory differs')
                nodes.append((f'callback-{ordinal}', [cb['positions_nm'][3*i:3*i+3] for i in range(len(rows))], cb))
            nodes.append(('final', path['final_nm'], None))
            if len(path['callbacks']) != trace['callbacks'] or len(path['callbacks']) > plan['policy']['maximum_callbacks']:
                raise ValueError('graph callback count differs')
            descriptions = []
            for label, points, cb in nodes:
                if len(points) != len(rows) or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
                    raise ValueError('nonfinite/invalid graph control coordinates')
                if any(points[i] != original_points[i] for i in heavy):
                    raise ValueError('graph H path moved heavy observations')
                context.setPositions(points*nanometer)
                state = context.getState(energy=True, forces=True)
                energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                forces = state.getForces().value_in_unit(kilojoules_per_mole/nanometer)
                max_error = None
                if cb is not None:
                    gradient = cb['gradient_kj_mol_nm']
                    if any(not math.isclose(forces[i][j], -gradient[3*i+j], rel_tol=1e-7, abs_tol=1e-4) for i in mobile for j in range(3)):
                        raise ValueError('graph native gradient differs from physical force')
                    if cb['statistics']['restraint energy'] != 0 or cb['statistics']['max constraint error'] != 0 or not math.isclose(energy, cb['statistics']['system energy'], rel_tol=1e-8, abs_tol=1e-4):
                        raise ValueError('graph native objective differs')
                    max_error = max((abs(forces[i][j]+gradient[3*i+j]) for i in mobile for j in range(3)), default=0.)
                    gradient_errors.append(max_error)
                item = dict(stage=label, HA=ha_geometry(rows, bonds, points), energy_kj_mol=energy,
                    maximum_mobile_gradient_error_kj_mol_nm=max_error)
                if cb is None:
                    item.update(H_bonds=h_bonds(system, rows, points), H_angles=hydrogen_angles(system, rows, points),
                        close_H_pairs=h_close_pairs(rows, bonds, points))
                descriptions.append(item)
            exact = path['final_nm'] == original['final_nm']
            if exact != trace['exact_original_endpoint'] or (cell['arm'] == 'unchanged' and not exact):
                raise ValueError('graph unchanged nm endpoint not exact')
            report.update(groups=len(groups), callbacks=len(path['callbacks']), nodes=descriptions,
                exact_heavy_path=True, exact_unowned_start=True, exact_original_endpoint=exact,
                grids=trace['grids'], any_all_acute_callback=any(n['HA']['all_acute'] for n in descriptions[1:-1]))
        finally:
            del context
            del integrator
        results.append(report)
    summary = dict(enrolled=len(results), statuses=dict(Counter(c['status'] for c in results)), arms={})
    for arm in plan['policy']['arms']:
        cells = [c for c in results if c['status'] == 'completed' and c['arm'] == arm]
        summary['arms'][arm] = dict(completed=len(cells), groups=sum(c['groups'] for c in cells),
            any_all_acute_callback_cells=sum(c['any_all_acute_callback'] for c in cells),
            final_all_acute_HA=sum(len(c['nodes'][-1]['HA']['all_acute']) for c in cells),
            exact_original_endpoints=sum(c['exact_original_endpoint'] for c in cells))
    summary.update(maximum_grid_energy_readback_error_kj_mol=max(grid_errors, default=None),
        maximum_mobile_gradient_error_kj_mol_nm=max(gradient_errors, default=None))
    write(output, dict(schema='dvbfixer.graph-methyl-control-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        protocol=plan['policy'], production_changed=False, whole_minimization_executed=False,
        scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
