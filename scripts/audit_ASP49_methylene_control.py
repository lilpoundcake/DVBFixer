"""Independent ownership/slot/native-gradient audit of all prospective CH2 cells."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.attribute_graph_methyl_ha import frame, projection
from scripts.audit_graph_methyl_control import h_close_pairs, ha_geometry
from scripts.audit_methyl_start_control import h_bonds
from scripts.audit_methyl_whole_control import hydrogen_angles
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def check_initial(rows, before, after, groups):
    owned = {i for group in groups for i in group['hydrogens']}
    if len(before) != len(rows) or len(after) != len(rows) or any(before[i] != after[i] for i in range(len(rows)) if i not in owned):
        raise ValueError('CH2 initialization moved an unowned coordinate')
    for group in groups:
        centre = before[group['parent']]
        vectors = [[x-y for x, y in zip(before[i], centre)] for i in group['heavy']]
        normal = [vectors[0][1]*vectors[1][2]-vectors[0][2]*vectors[1][1], vectors[0][2]*vectors[1][0]-vectors[0][0]*vectors[1][2], vectors[0][0]*vectors[1][1]-vectors[0][1]*vectors[1][0]]
        if math.hypot(*normal) <= 1e-12:
            raise ValueError('undefined independent CH2 plane')
        if group['hydrogens'] != sorted(group['hydrogens'], key=lambda i: tuple(rows[i]['key'])) or group['heavy'] != sorted(group['heavy'], key=lambda i: tuple(rows[i]['key'])):
            raise ValueError('CH2 slot identity ordering differs')
        for sign, i in zip((1., -1.), group['hydrogens'], strict=True):
            vector = [x-y for x, y in zip(after[i], centre)]
            radius = math.hypot(*vector)
            if not math.isclose(radius, math.dist(before[i], centre), rel_tol=1e-11, abs_tol=1e-13):
                raise ValueError('CH2 initialization changed original H radius')
            for heavy in vectors:
                cosine = sum(x*y for x, y in zip(vector, heavy))/(radius*math.hypot(*heavy))
                if not math.isclose(cosine, -1/3, abs_tol=1e-10):
                    raise ValueError('CH2 H-heavy slot is not heavy-compatible tetrahedral')
            if sign*sum(x*y for x, y in zip(vector, normal)) <= 0:
                raise ValueError('CH2 deterministic H slot sign differs')


def focal(rows, points):
    target = ['A', '49', '', 'ASP', 'HA']
    ca, ha, outward = frame(rows, points, target)
    vector = [x-y for x, y in zip(points[ha], points[ca])]
    norm = math.hypot(*vector)
    angles = {}
    for name in ('N', 'C', 'CB'):
        i = next(i for i, r in enumerate(rows) if r['key'] == [*target[:4], name])
        direction = [x-y for x, y in zip(points[i], points[ca])]
        cosine = sum(x*y for x, y in zip(vector, direction))/(norm*math.hypot(*direction))
        angles[name] = math.degrees(math.acos(max(-1., min(1., cosine))))
    distances = {}
    for name in ('HB2', 'HB3'):
        i = next(i for i, r in enumerate(rows) if r['key'] == ['A', '48', '', 'LEU', name])
        distances[name] = 10*math.dist(points[ha], points[i])
    return dict(outward_cosine=projection(vector, outward)/norm, CA_H_angstrom=10*norm,
        HA_heavy_angles_degrees=angles, LEU48_H_distances_angstrom=distances)


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['source_arm'], arm) for c in plan['cases'] for arm in plan['policy']['arms']]
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw or [(c['id'], c['seed'], c['source_arm'], c['arm']) for c in run['cases']] != expected:
        raise ValueError('complete CH2 cohort differs')
    declarations = {(c['id'], c['seed'], c['source_arm']): c for c in plan['cases']}
    results, errors = [], []
    for cell in run['cases']:
        case = declarations[(cell['id'], cell['seed'], cell['source_arm'])]
        source = Path(case['source'])
        folder = root/cell['id']/f"seed-{cell['seed']}"/cell['source_arm']/cell['arm']
        for name, digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('CH2 artifact changed')
        result = {k: cell[k] for k in ('id', 'seed', 'source_arm', 'arm', 'status', 'returncode')}
        if cell['status'] != 'completed':
            results.append(result)
            continue
        trace = json.loads((folder/'trace.json').read_bytes())
        path = json.loads(gzip.decompress((folder/'path.json.gz').read_bytes()))
        original = json.loads(gzip.decompress(Path(case['path']).read_bytes()))
        ownership = json.loads((source/'H-ownership.json').read_bytes())
        rows, bonds = ownership['rows'], ownership['bonds']
        groups = [] if cell['arm'] == 'unchanged' else case['local'] if cell['arm'] == 'LEU48_CH2' else case['groups']
        if trace['case'] != case or trace['arm'] != cell['arm'] or trace['plan_sha256'] != sha(raw) or trace['groups'] != groups:
            raise ValueError('CH2 applied intervention differs')
        adjacent = [set() for _ in rows]
        for a, b in bonds:
            adjacent[a].add(b)
            adjacent[b].add(a)
        for group in groups:
            parent = group['parent']
            if rows[parent]['element'] != 'C' or adjacent[parent] != set(group['heavy']+group['hydrogens']) or not set(group['hydrogens']) <= set(ownership['new_h']) & set(ownership['mobile']):
                raise ValueError('CH2 native graph/ownership differs')
            if any(rows[i]['element'] not in {'H', 'D'} for i in group['hydrogens']) or any(rows[i]['element'] in {'H', 'D'} for i in group['heavy']):
                raise ValueError('CH2 native element inventory differs')
        check_initial(rows, original['initial_nm'], path['initial_nm'], groups)
        system = openmm.XmlSerializer.deserialize((source/'H-system.xml').read_text())
        if sha(openmm.XmlSerializer.serialize(system).encode()) != trace['system_sha256']:
            raise ValueError('CH2 System changed')
        mobile = [i for i in range(len(rows)) if system.getParticleMass(i).value_in_unit(dalton) > 0]
        fixed = set(range(len(rows)))-set(mobile)
        integrator = openmm.VerletIntegrator(0.)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        try:
            nodes = [('initial', path['initial_nm'], None)]
            for i, cb in enumerate(path['callbacks']):
                if cb['ordinal'] != i or len(cb['positions_nm']) != 3*len(rows) or len(cb['gradient_kj_mol_nm']) != 3*len(rows):
                    raise ValueError('CH2 callback inventory differs')
                nodes.append((f'callback-{i}', [cb['positions_nm'][3*j:3*j+3] for j in range(len(rows))], cb))
            nodes.append(('final', path['final_nm'], None))
            if len(path['callbacks']) != trace['callbacks'] or len(path['callbacks']) > plan['policy']['maximum_callbacks']:
                raise ValueError('CH2 callback count differs')
            measured = []
            for label, points, cb in nodes:
                if len(points) != len(rows) or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
                    raise ValueError('CH2 nonfinite/malformed coordinates')
                if any(points[i] != original['initial_nm'][i] for i in fixed):
                    raise ValueError('CH2 H-only path moved a fixed observation')
                context.setPositions(points*nanometer)
                state = context.getState(energy=True, forces=True)
                energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                max_error = None
                if cb is not None:
                    forces = state.getForces().value_in_unit(kilojoules_per_mole/nanometer)
                    gradient = cb['gradient_kj_mol_nm']
                    if any(not math.isclose(forces[i][j], -gradient[3*i+j], rel_tol=1e-7, abs_tol=1e-4) for i in mobile for j in range(3)) or cb['statistics']['restraint energy'] != 0 or cb['statistics']['max constraint error'] != 0 or not math.isclose(energy, cb['statistics']['system energy'], rel_tol=1e-8, abs_tol=1e-4):
                        raise ValueError('CH2 physical force/objective differs from native callback')
                    max_error = max((abs(forces[i][j]+gradient[3*i+j]) for i in mobile for j in range(3)), default=0.)
                    errors.append(max_error)
                node = dict(stage=label, focal=focal(rows, points), global_HA=ha_geometry(rows, bonds, points),
                    energy_kj_mol=energy, maximum_mobile_gradient_error_kj_mol_nm=max_error)
                if cb is None:
                    node.update(H_bonds=h_bonds(system, rows, points), H_angles=hydrogen_angles(system, rows, points),
                        close_H_pairs=h_close_pairs(rows, bonds, points))
                measured.append(node)
            exact = path['final_nm'] == original['final_nm']
            if exact != trace['exact_original_nm_endpoint'] or (cell['arm'] == 'unchanged' and not exact):
                raise ValueError('CH2 unchanged endpoint not exact')
            result.update(changed_groups=len(groups), callbacks=len(path['callbacks']), exact_original_nm_endpoint=exact,
                exact_unowned_start=True, exact_fixed_path=True, nodes=measured)
        finally:
            del context
            del integrator
        results.append(result)
    summary = dict(enrolled=len(results), statuses=dict(Counter(c['status'] for c in results)), arms={})
    for arm in plan['policy']['arms']:
        cells = [c for c in results if c['status'] == 'completed' and c['arm'] == arm]
        summary['arms'][arm] = dict(completed=len(cells), exact_original_nm_endpoints=sum(c['exact_original_nm_endpoint'] for c in cells),
            any_focal_inward_callback=sum(any(n['focal']['outward_cosine'] < 0 for n in c['nodes'][1:-1]) for c in cells),
            final_focal_inward=sum(c['nodes'][-1]['focal']['outward_cosine'] < 0 for c in cells),
            final_global_all_acute_HA=sum(len(c['nodes'][-1]['global_HA']['all_acute']) for c in cells))
    summary['maximum_mobile_gradient_error_kj_mol_nm'] = max(errors, default=None)
    write(output, dict(schema='dvbfixer.ASP49-methylene-control-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), policy=plan['policy'],
        audit_script_sha256=sha(Path(__file__).read_bytes()), production_changed=False,
        whole_minimization_executed=False, scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
