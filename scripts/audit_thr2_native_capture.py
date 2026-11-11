"""Independent saved-only THR2 geometry/gradient audit; never runs an optimizer."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_stereo_force_projection import volume_gradient
from scripts.branched_stereo_oracle import rdkit_labels
from scripts.research_dunbrack_prior import sha, write
from scripts.research_legacy_start_controls import snapshot_atoms
from scripts.research_thr2_native_capture import check


def read(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def verify_files(folder, row):
    for name, digest in row.get('artifacts', {}).items():
        if sha((folder/name).read_bytes()) != digest:
            raise ValueError(f'captured artifact changed: {folder/name}')


def delta(a, b, factor=10.):
    if len(a) != len(b) or not a:
        raise ValueError('coordinate inventory mismatch')
    distances = [factor*math.dist(x, y) for x, y in zip(a, b, strict=True)]
    return dict(exact=a == b, maximum_angstrom=max(distances),
        rms_angstrom=math.sqrt(sum(x*x for x in distances)/len(distances)))


def geometry(rows, points):
    if len(rows) != len(points):
        raise ValueError('atom/coordinate inventory mismatch')
    mapping = {tuple(row['key']):i for i, row in enumerate(rows)}
    if len(mapping) != len(rows):
        raise ValueError('duplicate atom identity')
    residues = {key[:4] for key in mapping}
    alpha = []
    for residue in sorted(residues):
        keys = [(*residue, name) for name in ('CA', 'N', 'C', 'CB')]
        if all(key in mapping for key in keys):
            volume, _, _ = volume_gradient([points[mapping[key]] for key in keys])
            alpha.append(dict(residue=list(residue), volume_angstrom3=1000*volume))
    target = ('A', '2', '', 'THR')
    names = ('CB', 'CA', 'OG1', 'CG2')
    volume, _, _ = volume_gradient([points[mapping[(*target, name)]] for name in names])
    xyz = {name:tuple(10*x for x in points[i]) for key, i in mapping.items()
        if key[:4] == target for name in (key[-1],)}
    cip = rdkit_labels(xyz, 'THR')
    return dict(target_ca_volume_angstrom3=next(r['volume_angstrom3'] for r in alpha if r['residue'] == list(target)),
        target_cb_volume_angstrom3=1000*volume, target_cip=cip,
        negative_alpha=[r for r in alpha if r['volume_angstrom3'] < 0],
        minimum_alpha_volume_angstrom3=min(r['volume_angstrom3'] for r in alpha), alpha_count=len(alpha))


def pdb_geometry(path):
    from openmm.app import PDBFile

    pdb = PDBFile(str(path))
    rows = snapshot_atoms(pdb.topology, pdb.positions)
    return geometry(rows, [[x/10 for x in row['xyz']] for row in rows])


def states(source, call):
    import openmm
    from openmm.unit import nanometer

    result = {}
    for side in ('before', 'after'):
        if side not in call:
            continue
        prefix = call[side]['prefix']
        state = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text())
        points = state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        rows = read(source/f'{prefix}-atoms.json.gz')['atoms']
        result[side] = dict(geometry=geometry(rows, points), rows=rows, points=points)
    return result


def audit_replay(plan, source, call, folder, result):
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    verify_files(folder, result)
    report = {key:result[key] for key in ('label', 'number', 'mode', 'status')}
    if not (folder/'path.json.gz').exists():
        return report
    payload = read(folder/'path.json.gz')
    prefix = call['before']['prefix']
    rows = read(source/f'{prefix}-atoms.json.gz')['atoms']
    system = openmm.XmlSerializer.deserialize((source/f'{prefix}-system.xml').read_text())
    count = system.getNumParticles()
    if count != len(rows) or system.getNumConstraints() or any(system.isVirtualSite(i) for i in range(count)):
        raise ValueError('gradient audit requires exact inventory, no constraints/virtual sites')
    mobile = [i for i in range(count) if system.getParticleMass(i).value_in_unit(dalton) > 0]
    integrator = openmm.XmlSerializer.deserialize((source/f'{prefix}-integrator.xml').read_text())
    declared = call['before']['platform']
    if result['mode'] == 'Reference':
        platform = openmm.Platform.getPlatformByName('Reference')
        properties = {}
    else:
        platform = openmm.Platform.getPlatformByName(declared['name'])
        properties = declared['constructor_properties']
    context = openmm.Context(system, integrator, platform, properties)
    initial_state = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text())
    context.setState(initial_state)
    before = initial_state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    after_prefix = call.get('after', {}).get('prefix')
    native_after = None
    if after_prefix:
        native_after = openmm.XmlSerializer.deserialize((source/f'{after_prefix}-state.xml').read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    callbacks = payload['callbacks']
    if len(callbacks) > plan['protocol']['maximum_callbacks']:
        raise ValueError('callback budget exceeded')
    entries = [('initial', payload['initial_nm'], None)]
    for ordinal, cb in enumerate(callbacks):
        if cb['ordinal'] != ordinal or len(cb['positions_nm']) != 3*count or len(cb['gradient_kj_mol_nm']) != 3*count:
            raise ValueError('callback identity/inventory mismatch')
        entries.append((f'callback-{ordinal}', [cb['positions_nm'][3*i:3*i+3] for i in range(count)], cb))
    entries.append(('returned', payload['final_nm'], None))
    nodes, gradient_errors, objective_errors = [], [], []
    initial_norm = max(1., math.sqrt(sum(x*x for p in payload['initial_nm'] for x in p)/count))
    epsilon = call['tolerance']/initial_norm
    try:
        for label, points, cb in entries:
            if len(points) != count or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
                raise ValueError('nonfinite/malformed saved coordinates')
            context.setPositions(points*nanometer)
            state = context.getState(energy=True, forces=True)
            energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            forces = state.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
            if not math.isfinite(energy) or not all(math.isfinite(forces[i][j]) for i in mobile for j in range(3)):
                raise ValueError('nonfinite saved-state physical readback')
            node = dict(stage=label, geometry=geometry(rows, points), energy_kj_mol=energy,
                mobile_force_l2_kj_mol_nm=math.sqrt(sum(forces[i][j]**2 for i in mobile for j in range(3))))
            ratio = node['mobile_force_l2_kj_mol_nm']/max(1., math.sqrt(sum(x*x for p in points for x in p)))
            node.update(native_scaled_gradient_ratio=ratio, native_epsilon=epsilon,
                satisfies_native_gradient_test=ratio <= epsilon)
            if cb is not None:
                gradient = cb['gradient_kj_mol_nm']
                error = max(abs(forces[i][j]+gradient[3*i+j]) for i in mobile for j in range(3))
                objective = abs(energy-cb['statistics']['system energy'])
                node.update(maximum_mobile_gradient_readback_error_kj_mol_nm=float(error),
                    objective_readback_error_kj_mol=objective, native_statistics=cb['statistics'])
                gradient_errors.append(float(error))
                objective_errors.append(objective)
                if cb['statistics']['restraint energy'] != 0 or cb['statistics']['max constraint error'] != 0:
                    raise ValueError('unexpected constraint penalty in unconstrained replay')
            nodes.append(node)
    finally:
        del context
    last_points = entries[-2][1]
    crossing = []
    for previous, current in zip(nodes, nodes[1:]):
        a = previous['geometry']['target_cb_volume_angstrom3']
        b = current['geometry']['target_cb_volume_angstrom3']
        if a*b < 0:
            crossing.append(dict(before_stage=previous['stage'], after_stage=current['stage'],
                before_volume_angstrom3=a, after_volume_angstrom3=b,
                energy_change_kj_mol=current['energy_kj_mol']-previous['energy_kj_mol']))
    report.update(stage=call['stage'], callbacks=len(callbacks), returned=payload['returned'],
        mobile_particles=len(mobile), particles=count, start_delta=delta(before, payload['initial_nm']),
        native_endpoint_delta=delta(native_after, payload['final_nm']) if native_after is not None else None,
        terminal=dict(returned_vs_last_reported=delta(last_points, payload['final_nm']),
            returned_minus_last_reported_energy_kj_mol=nodes[-1]['energy_kj_mol']-nodes[-2]['energy_kj_mol'],
            last_reported_geometry=nodes[-2]['geometry'], returned_geometry=nodes[-1]['geometry'],
            returned_satisfies_native_gradient_test=nodes[-1]['satisfies_native_gradient_test'],
            last_reported_iteration=callbacks[-1]['iteration'] if callbacks else None,
            comparison='last reported callback' if callbacks else 'initial (no callbacks)'),
        maximum_mobile_gradient_readback_error_kj_mol_nm=max(gradient_errors, default=None),
        maximum_objective_readback_error_kj_mol=max(objective_errors, default=None),
        target_beta_sign_crossings=crossing,
        negative_target_beta_callbacks=sum(n['geometry']['target_cb_volume_angstrom3'] < 0 for n in nodes[1:-1]),
        negative_alpha_nodes=sum(bool(n['geometry']['negative_alpha']) for n in nodes), nodes=nodes)
    return report


def audit(plan_path, root, private_output, public_output):
    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    check(plan)
    run = read(root/'evidence.json')
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw:
        raise ValueError('execution plan differs')
    if [r['label'] for r in run['native']] != plan['protocol']['labels']:
        raise ValueError('native enrollment differs')
    for name, digest in read(root/'replay-inputs.json')['inputs'].items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError('replay input binding changed')
    native, replay_rows, expected_replays = [], [], []
    for index, result in enumerate(run['native']):
        source = root/result['label']
        verify_files(source, result)
        report = dict(label=result['label'], status=result['status'], calls=[], h_calls=[])
        if (source/'capture.json').exists():
            capture = read(source/'capture.json')
            if capture['plan_sha256'] != sha(raw) or capture['label'] != result['label']:
                raise ValueError('native capture provenance differs')
            for ordinal, h in enumerate(capture['h_calls'], 1):
                seed = plan['protocol']['h_seeds'][index]
                expected_rng = h['after_rng'] if seed is None else h['before_rng']
                if h['seed'] != seed or h['restored_rng'] != expected_rng:
                    raise ValueError('H RNG intervention/restoration differs')
                first = read(source/f'h-{ordinal}-before.json.gz')['atoms']
                second = read(source/f'h-{ordinal}-after.json.gz')['atoms']
                a = {tuple(r['key']):r['xyz'] for r in first if r['element'] not in ('H', 'D')}
                b = {tuple(r['key']):r['xyz'] for r in second if r['element'] not in ('H', 'D')}
                if a.keys() != b.keys():
                    raise ValueError('H stage lost/added heavy identity')
                report['h_calls'].append(dict(returned=h['returned'], seed=seed,
                    heavy_delta=delta(list(a.values()), [b[k] for k in a], factor=1.)))
            for call in capture['calls']:
                measured = states(source, call)
                report['calls'].append(dict(number=call['number'], stage=call['stage'], returned=call['returned'],
                    tolerance=call['tolerance'], max_iterations=call['max_iterations'],
                    before=measured.get('before', {}).get('geometry'), after=measured.get('after', {}).get('geometry'),
                    platform=call.get('before', {}).get('platform')))
                if 'before' in call:
                    for mode in ('captured', 'Reference'):
                        expected_replays.append((result['label'], call['number'], mode))
            if (source/'minimized.pdb').exists():
                report['final_pdb'] = pdb_geometry(source/'minimized.pdb')
        native.append(report)
    if [(r['label'], r['number'], r['mode']) for r in run['replays']] != expected_replays:
        raise ValueError('complete saved-call replay enrollment differs')
    for result in run['replays']:
        source = root/result['label']
        call = next(c for c in read(source/'capture.json')['calls'] if c['number'] == result['number'])
        folder = root/'replays'/result['label']/f"call-{result['number']}"/result['mode']
        replay_rows.append(audit_replay(plan, source, call, folder, result))
    summary = dict(native_statuses=dict(Counter(r['status'] for r in native)),
        replay_statuses=dict(Counter(r['status'] for r in replay_rows)),
        native_calls=sum(len(r['calls']) for r in native),
        callbacks=sum(r.get('callbacks', 0) for r in replay_rows),
        final_target_cb_cip={r['label']:r.get('final_pdb', {}).get('target_cip', {}).get('CB') for r in native},
        by_mode={})
    for mode in ('captured', 'Reference'):
        rows = [r for r in replay_rows if r['mode'] == mode and 'terminal' in r]
        summary['by_mode'][mode] = dict(calls=len(rows), callbacks=sum(r['callbacks'] for r in rows),
            exact_starts=sum(r['start_delta']['exact'] for r in rows),
            exact_native_endpoints=sum(r['native_endpoint_delta'] is not None and r['native_endpoint_delta']['exact'] for r in rows),
            returned_differs_last_reported=sum(not r['terminal']['returned_vs_last_reported']['exact'] for r in rows),
            returned_higher_energy=sum(r['terminal']['returned_minus_last_reported_energy_kj_mol'] > 1e-6 for r in rows),
            maximum_mobile_gradient_readback_error_kj_mol_nm=max((r['maximum_mobile_gradient_readback_error_kj_mol_nm'] for r in rows if r['maximum_mobile_gradient_readback_error_kj_mol_nm'] is not None), default=None),
            maximum_objective_readback_error_kj_mol=max((r['maximum_objective_readback_error_kj_mol'] for r in rows if r['maximum_objective_readback_error_kj_mol'] is not None), default=None))
    historical_pdb = next(Path(name) for name in plan['inputs'] if name.endswith('/minimized.pdb'))
    common = dict(schema='dvbfixer.thr2-native-audit.v1', plan_sha256=sha(raw),
        execution_sha256=sha((root/'evidence.json').read_bytes()), audit_script_sha256=sha(Path(__file__).read_bytes()),
        summary=summary, native=native, historical_final_pdb=pdb_geometry(historical_pdb),
        limits='OpenCL readback deltas are measured, not asserted bit equality; Reference is a numerical control, not historical causality',
        production_changed=False, additional_minimization=False)
    write(private_output, dict(common, replays=replay_rows))
    public_rows = [{k:v for k, v in row.items() if k != 'nodes'} for row in replay_rows]
    write(public_output, dict(common, replays=public_rows, private_audit_sha256=sha(private_output.read_bytes())))
    check(plan)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'private-output', 'public-output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.private_output.resolve(), args.public_output.resolve())
