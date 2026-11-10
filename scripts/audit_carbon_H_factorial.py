"""Independent complete factorial ownership, geometry and native-force readback.

No minimization/optimization. Preserve stress/natural strata and known replays.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_ASP49_methylene_control import check_initial
from scripts.audit_graph_methyl_control import h_close_pairs, ha_geometry
from scripts.audit_methyl_start_control import h_bonds
from scripts.audit_methyl_whole_control import hydrogen_angles
from scripts.carbon_H_capture_adapter import adapt
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def check_methyl(rows, before, after, groups):
    owned = {i for g in groups for i in g['hydrogens']}
    if len(after) != len(before) or any(before[i] != after[i] for i in range(len(before)) if i not in owned):
        raise ValueError('methyl initializer moved an unowned coordinate')
    def unit(vector):
        norm = math.hypot(*vector)
        if norm <= 1e-12:
            raise ValueError('undefined independent methyl frame')
        return [x/norm for x in vector]
    for g in groups:
        origin = before[g['parent']]
        axis = unit([x-y for x, y in zip(before[g['axis']], origin)])
        ref = [x-y for x, y in zip(before[g['frame']], origin)]
        dot = sum(x*y for x, y in zip(ref, axis))
        first = unit([x-dot*y for x, y in zip(ref, axis)])
        second = [axis[1]*first[2]-axis[2]*first[1], axis[2]*first[0]-axis[0]*first[2], axis[0]*first[1]-axis[1]*first[0]]
        if g['hydrogens'] != sorted(g['hydrogens'], key=lambda i: tuple(rows[i]['key'])):
            raise ValueError('methyl name-ordered H slots differ')
        directions = []
        for ordinal, i in enumerate(g['hydrogens']):
            vector = [x-y for x, y in zip(after[i], origin)]
            radius = math.hypot(*vector)
            if not math.isclose(radius, math.dist(before[i], origin), rel_tol=1e-11, abs_tol=1e-13):
                raise ValueError('methyl original radius changed')
            direction = unit(vector)
            expected = (-1/3, math.sqrt(8/9)*math.cos(2*math.pi*ordinal/3), math.sqrt(8/9)*math.sin(2*math.pi*ordinal/3))
            actual = [sum(x*y for x, y in zip(direction, basis)) for basis in (axis, first, second)]
            if any(not math.isclose(a, b, abs_tol=1e-10) for a, b in zip(actual, expected)):
                raise ValueError('fixed-zero methyl tetrahedral slot differs')
            directions.append(direction)
        for i, j in ((0, 1), (0, 2), (1, 2)):
            if not math.isclose(sum(x*y for x, y in zip(directions[i], directions[j])), -1/3, abs_tol=1e-10):
                raise ValueError('methyl H-H geometry differs')


def check_factorial_initial(rows, before, after, case, arm):
    if len(before) != len(rows) or len(after) != len(rows) or any(
            len(p) != 3 or not all(math.isfinite(x) for x in p) for p in [*before, *after]):
        raise ValueError('factorial initial coordinates malformed')
    methyl = case['methyl']['groups'] if arm in {'all_methyl_zero', 'combined_zero'} else []
    methylene = case['methylene']['groups'] if arm in {'all_methylene', 'combined_zero'} else []
    first = {i for g in methyl for i in g['hydrogens']}
    second = {i for g in methylene for i in g['hydrogens']}
    if first & second or any(before[i] != after[i] for i in range(len(rows)) if i not in first | second):
        raise ValueError('factorial overlap or nonowned initial change')
    # Isolate disjoint coordinates for each independent geometry checker.
    check_methyl(rows, before, [after[i] if i in first else p for i, p in enumerate(before)], methyl)
    check_initial(rows, before, [after[i] if i in second else p for i, p in enumerate(before)], methylene)


def known_replay(case, arm):
    root = Path('tmp/pdbfixer-research')
    if case['context_kind'] != 'prepared-protein-view':
        return None
    if arm == 'all_methyl_zero':
        run_root = root/'graph-methyl-control-run-01'
        audit_path = Path('docs/research/graph-methyl-control-evidence.json')
        folder = run_root/case['id']/f"seed-{case['seed']}"/'graph_tetrahedral_zero'
    elif case['id'] == '1TIM' and arm in {'all_methylene', 'combined_zero'}:
        source_arm = 'unchanged' if arm == 'all_methylene' else 'graph_tetrahedral_zero'
        run_root = root/'ASP49-methylene-control-run-01'
        audit_path = Path('docs/research/ASP49-methylene-control-evidence.json')
        folder = run_root/'1TIM'/f"seed-{case['seed']}"/source_arm/'all_carbon_CH2'
    else:
        return None
    manifest_raw = (folder/'result.json').read_bytes()
    manifest = json.loads(manifest_raw)
    run_raw = (run_root/'evidence.json').read_bytes()
    audit_raw = audit_path.read_bytes()
    if (json.loads(audit_raw)['run_sha256'] != sha(run_raw)
            or manifest not in json.loads(run_raw)['cases']):
        raise ValueError('known replay does not match previously published run binding')
    raw = (folder/'path.json.gz').read_bytes()
    if manifest['status'] != 'completed' or manifest['artifacts']['path.json.gz'] != sha(raw):
        raise ValueError('known historical replay artifact differs')
    return dict(path=str(folder/'path.json.gz'), manifest_sha256=sha(manifest_raw),
        path_sha256=sha(raw), published_audit_sha256=sha(audit_raw),
        run_sha256=sha(run_raw), data=json.loads(gzip.decompress(raw)))


def summarize(cases):
    arms = ('unchanged', 'all_methyl_zero', 'all_methylene', 'combined_zero')
    result = dict(enrolled=len(cases), statuses=dict(Counter(c['status'] for c in cases)), strata={})
    complete = {}
    for c in cases:
        complete.setdefault((c['id'], c['seed'], c['conditioning']), set())
        if c['status'] == 'completed':
            complete[(c['id'], c['seed'], c['conditioning'])].add(c['arm'])
    matched = {k for k, v in complete.items() if v == set(arms)}
    result['complete_matched_contexts'] = len(matched)
    for kind in ('geometry-reconstruction-stress', 'prepared-protein-view'):
        rows = [c for c in cases if c['context_kind'] == kind]
        stratum = dict(enrolled=len(rows), statuses=dict(Counter(c['status'] for c in rows)), arms={})
        for arm in arms:
            cells = [c for c in rows if c['arm'] == arm and c['status'] == 'completed']
            stratum['arms'][arm] = dict(completed=len(cells),
                matched_contexts=sum((c['id'], c['seed'], c['conditioning']) in matched for c in cells),
                exact_native_endpoints=sum(c['exact_original_nm_endpoint'] for c in cells),
                exact_known_replays=sum(c.get('known_replay') is not None for c in cells),
                any_all_acute_callback_cells=sum(any(n['HA']['all_acute'] for n in c['nodes'][1:-1]) for c in cells),
                final_all_acute_HA=sum(len(c['nodes'][-1]['HA']['all_acute']) for c in cells))
        result['strata'][kind] = stratum
    return result


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    audit_sources = {name: sha(Path(name).read_bytes()) for name in (
        'scripts/audit_carbon_H_factorial.py', 'scripts/audit_ASP49_methylene_control.py',
        'scripts/audit_graph_methyl_control.py', 'scripts/audit_methyl_start_control.py',
        'scripts/audit_methyl_whole_control.py')}
    if openmm.__version__ != plan['openmm_version']:
        raise ValueError('factorial native readback runtime differs')
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['conditioning'], arm) for c in plan['cases'] for arm in plan['policy']['arms']]
    if (run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw
            or [(c['id'], c['seed'], c['conditioning'], c['arm']) for c in run['cases']] != expected):
        raise ValueError('complete factorial execution enrollment differs')
    declarations = {(c['id'], c['seed'], c['conditioning']): c for c in plan['cases']}
    results, errors = [], []
    for cell in run['cases']:
        case = declarations[(cell['id'], cell['seed'], cell['conditioning'])]
        folder = root/cell['id']/f"seed-{cell['seed']}"/cell['conditioning']/cell['arm']
        for name, digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('factorial artifact changed')
        report = {k: cell[k] for k in ('id', 'seed', 'conditioning', 'arm', 'status', 'returncode')}
        report.update(context_kind=case['context_kind'], replay_labels=case.get('replay_labels', []))
        if cell['status'] != 'completed':
            if case['status'] == 'refused_source' and cell['status'] != 'refused_source':
                raise ValueError('inherited source refusal was changed')
            results.append(report)
            continue
        if case['status'] != 'admitted':
            raise ValueError('refused capture executed')
        data = adapt(case)
        rows, bonds = data['rows'], data['bonds']
        path = json.loads(gzip.decompress((folder/'path.json.gz').read_bytes()))
        trace = json.loads((folder/'trace.json').read_bytes())
        if (trace['plan_sha256'] != sha(raw) or trace['case_id'] != [case['id'], case['seed'], case['conditioning']]
                or trace['arm'] != cell['arm'] or trace['whole_minimization_executed'] or trace['production_changed']):
            raise ValueError('factorial worker trace differs')
        check_factorial_initial(rows, data['initial_nm'], path['initial_nm'], case, cell['arm'])
        system = data['system']
        if trace['system_sha256'] != sha(openmm.XmlSerializer.serialize(system).encode()):
            raise ValueError('factorial System changed')
        exact = path['final_nm'] == data['final_nm']
        if trace['exact_original_nm_endpoint'] != exact or (cell['arm'] == 'unchanged' and not exact):
            raise ValueError('unchanged factorial endpoint not exact')
        replay = known_replay(case, cell['arm'])
        if replay is not None:
            prior = replay.pop('data')
            if path['initial_nm'] != prior['initial_nm'] or path['final_nm'] != prior['final_nm']:
                raise ValueError('known computational replay start/endpoint differs')
        integrator = openmm.VerletIntegrator(0.)
        context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
        try:
            nodes = [('initial', path['initial_nm'], None)]
            for i, cb in enumerate(path['callbacks']):
                if cb['ordinal'] != i or len(cb['positions_nm']) != 3*len(rows) or len(cb['gradient_kj_mol_nm']) != 3*len(rows):
                    raise ValueError('factorial callback inventory differs')
                nodes.append((f'callback-{i}', [cb['positions_nm'][3*j:3*j+3] for j in range(len(rows))], cb))
            nodes.append(('final', path['final_nm'], None))
            if len(path['callbacks']) != trace['callbacks'] or len(path['callbacks']) > plan['policy']['maximum_callbacks']:
                raise ValueError('factorial callback budget differs')
            measured = []
            fixed = set(range(len(rows)))-set(data['mobile'])
            for label, points, cb in nodes:
                if len(points) != len(rows) or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
                    raise ValueError('factorial malformed/nonfinite path')
                if any(points[i] != data['initial_nm'][i] for i in fixed):
                    raise ValueError('factorial H path moved fixed observations')
                context.setPositions(points*nanometer)
                state = context.getState(energy=True, forces=True)
                energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                error = None
                if cb is not None:
                    forces = state.getForces().value_in_unit(kilojoules_per_mole/nanometer)
                    gradient = cb['gradient_kj_mol_nm']
                    if (any(not math.isclose(forces[i][j], -gradient[3*i+j], rel_tol=1e-7, abs_tol=1e-4)
                            for i in data['mobile'] for j in range(3))
                            or cb['statistics']['restraint energy'] != 0 or cb['statistics']['max constraint error'] != 0
                            or not math.isclose(energy, cb['statistics']['system energy'], rel_tol=1e-8, abs_tol=1e-4)):
                        raise ValueError('factorial native mobile-gradient/objective closure differs')
                    error = max((abs(forces[i][j]+gradient[3*i+j]) for i in data['mobile'] for j in range(3)), default=0.)
                    errors.append(error)
                node = dict(stage=label, energy_kj_mol=energy, HA=ha_geometry(rows, bonds, points),
                    H_bonds=h_bonds(system, rows, points), H_angles=hydrogen_angles(system, rows, points),
                    maximum_mobile_gradient_error_kj_mol_nm=error)
                if cb is None:
                    node['close_H_pairs'] = h_close_pairs(rows, bonds, points)
                measured.append(node)
            report.update(nodes=measured, callbacks=len(path['callbacks']), exact_original_nm_endpoint=exact,
                exact_nonowned_start=True, exact_fixed_path=True, known_replay=replay)
        finally:
            del context
            del integrator
        results.append(report)
    verify(plan)
    verify(dict(inputs={}, source_hashes=audit_sources))
    summary = summarize(results)
    summary.update(callbacks=sum(c.get('callbacks', 0) for c in results),
        nodes=sum(len(c.get('nodes', [])) for c in results),
        maximum_mobile_gradient_error_kj_mol_nm=max(errors, default=None))
    write(output, dict(schema='dvbfixer.carbon-H-factorial-audit.v1', cases=results,
        summary=summary, plan_sha256=sha(raw), run_sha256=sha(run_raw), policy=plan['policy'],
        audit_script_sha256=sha(Path(__file__).read_bytes()), production_changed=False,
        audit_source_hashes=audit_sources,
        whole_minimization_executed=False, scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
