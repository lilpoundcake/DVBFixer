"""Independent mask/transport/physics/stereo audit of the complete H-swap cohort."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.audit_stereo_force_projection import projected_force, volume_gradient
from scripts.audit_thr2_native_capture import audit_replay, delta, read, verify_files
from scripts.research_dunbrack_prior import sha, write
from scripts.research_thr2_hswap import check_plan


def independent_mask_check(atoms, recipient, donor, actual, mask):
    """Deliberately does not reuse the experiment's swap_points implementation."""
    if len(atoms) != len(recipient) or len(donor) != len(actual) or len(atoms) != len(actual):
        raise ValueError('complete start inventory differs')
    hb = [i for i, a in enumerate(atoms) if a['key'] == ['A', '2', '', 'THR', 'HB']]
    if len(hb) != 1 or atoms[hb[0]]['element'] != 'H':
        raise ValueError('target HB identity differs')
    selected = []
    for i, atom in enumerate(atoms):
        hydrogen = atom['element'] in ('H', 'D')
        choose = (mask == 'HB' and i == hb[0]) or (mask == 'other-H' and hydrogen and i != hb[0])
        if choose:
            selected.append(i)
        if actual[i] != (donor[i] if choose else recipient[i]):
            raise ValueError('exact selected/unselected start ownership violated')
    if mask not in ('unchanged', 'HB', 'other-H'):
        raise ValueError('undeclared atom mask')
    return selected


def initial_projection(source, native, folder, mode, points, atoms):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    system = openmm.XmlSerializer.deserialize((folder/'phase-1/system.xml').read_text())
    integrator = openmm.VerletIntegrator(0.)
    declared = native['before']['platform']
    name = 'Reference' if mode == 'Reference' else declared['name']
    props = {} if mode == 'Reference' else declared['constructor_properties']
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(name), props)
    try:
        context.setState(openmm.XmlSerializer.deserialize((source/f"{native['before']['prefix']}-state.xml").read_text()))
        context.setPositions(points*nanometer)
        forces = context.getState(forces=True).getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
        keys = {tuple(row['key']):i for i, row in enumerate(atoms)}
        ids = [keys[('A', '2', '', 'THR', name)] for name in ('CB', 'CA', 'OG1', 'CG2')]
        _, gradient, length = volume_gradient([points[i] for i in ids])
        return projected_force(forces, ids, gradient, length)
    finally:
        del context


def audit(plan_path, root, private_output, public_output):
    import openmm
    from openmm.unit import nanometer

    record = read(plan_path)
    check_plan(record)
    execution = read(root/'evidence.json')
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('executed H-swap plan differs')
    if [{k:r[k] for k in ('recipient', 'donor', 'mask', 'mode')} for r in execution['cases']] != record['cases']:
        raise ValueError('complete prospective 30-path enrollment differs')
    starts = {c['label']:c for c in record['starts']}
    reports, errors, groups = [], [], {}
    replay_plan = dict(protocol=dict(maximum_callbacks=record['policy']['maximum_callbacks_per_phase']))
    for result in execution['cases']:
        folder = root/f"path-{result['index']:02d}"
        verify_files(folder, result)
        report = {k:v for k, v in result.items() if k not in ('artifacts',)}
        report['phases'] = []
        receiver, donor = starts[result['recipient']], starts[result['donor']]
        source, donor_source = Path(receiver['source']), Path(donor['source'])
        prefix = receiver['calls'][0]['before']['prefix']
        atoms = read(source/f'{prefix}-atoms.json.gz')['atoms']
        recipient_xyz = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        donor_prefix = donor['calls'][0]['before']['prefix']
        donor_xyz = openmm.XmlSerializer.deserialize((donor_source/f'{donor_prefix}-state.xml').read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        previous, group_key, endpoints = None, None, []
        for phase, native in enumerate(receiver['calls'], 1):
            path_folder = folder/f'phase-{phase}'
            if not (path_folder/'path.json.gz').exists():
                continue
            payload = read(path_folder/'path.json.gz')
            original_system = source/f"{native['before']['prefix']}-system.xml"
            if (path_folder/'system.xml').read_bytes() != original_system.read_bytes():
                raise ValueError('H-swap native System/phase parameters differ')
            if phase == 1:
                selected = independent_mask_check(atoms, recipient_xyz, donor_xyz, payload['initial_nm'], result['mask'])
                if selected != read(folder/'phase-1-result.json')['selected_atom_indices']:
                    raise ValueError('transported mask differs from declared observer record')
                report['selected_atoms'] = len(selected)
                report['exact_start_ownership'] = True
                report['initial_cb_projection_kj_mol_nm'] = initial_projection(source, native, folder, result['mode'], payload['initial_nm'], atoms)
                group_key = sha(json.dumps(payload['initial_nm'], separators=(',', ':'), allow_nan=False).encode())
                report['initial_coordinates_sha256'] = group_key
                if group_key not in groups:
                    groups[group_key] = dict(initial=payload['initial_nm'], members=[])
                elif groups[group_key]['initial'] != payload['initial_nm']:
                    raise ValueError('initial-coordinate hash collision')
                if result['mask'] == 'unchanged' and result['mode'] == 'Reference':
                    prior = source.parent/'replays'/result['recipient']/f"call-{native['number']}"/'Reference/path.json.gz'
                    control = read(prior)
                    report['exact_prior_reference_phase1_start'] = payload['initial_nm'] == control['initial_nm']
                    report['exact_prior_reference_phase1_endpoint'] = payload['final_nm'] == control['final_nm']
                    if not report['exact_prior_reference_phase1_start'] or not report['exact_prior_reference_phase1_endpoint']:
                        errors.append(dict(index=result['index'], problem='unchanged Reference control differs; no retry'))
            elif payload['initial_nm'] != previous:
                raise ValueError('phase2 is not chained from this path phase1 return')
            previous = payload['final_nm']
            endpoints.append(payload['final_nm'])
            measured = audit_replay(replay_plan, source, native, path_folder,
                dict(label=result['recipient'], number=native['number'], mode=result['mode'], status=result['status']))
            measured['phase'] = phase
            # These recipient-native deltas are expected under an intervention,
            # not mismatched-start admission failures or controls.
            measured['delta_from_recipient_native_start'] = measured.pop('start_delta')
            measured['delta_from_recipient_native_endpoint'] = measured.pop('native_endpoint_delta')
            report['phases'].append(measured)
        if report['phases']:
            report['final_geometry'] = report['phases'][-1]['terminal']['returned_geometry']
        if group_key is not None:
            groups[group_key]['members'].append(dict(index=result['index'], mode=result['mode'], endpoints=endpoints))
        reports.append(report)
    baselines = {(r['recipient'], r['mode']):r['final_geometry']['target_cip'].get('CB')
        for r in reports if r['mask'] == 'unchanged' and 'final_geometry' in r}
    for row in reports:
        hb_source = row['donor'] if row['mask'] == 'HB' else row['recipient']
        rest_source = row['donor'] if row['mask'] == 'other-H' else row['recipient']
        expected = baselines.get((hb_source, row['mode']))
        observed = row.get('final_geometry', {}).get('target_cip', {}).get('CB')
        row.update(hb_source=hb_source, other_h_source=rest_source,
            expected_cb_under_hb_localization=expected,
            matches_hb_localization=expected is not None and observed == expected)
    duplicates = []
    for digest, group in groups.items():
        reference = [m for m in group['members'] if m['mode'] == 'Reference']
        captured = [m for m in group['members'] if m['mode'] == 'captured']
        reference_equal = all(m['endpoints'] == reference[0]['endpoints'] for m in reference[1:]) if len(reference) > 1 else None
        native_deltas = [delta(a, b)['maximum_angstrom'] for m in captured[1:]
            for a, b in zip(captured[0]['endpoints'], m['endpoints'], strict=True)] if len(captured) > 1 else []
        duplicates.append(dict(initial_coordinates_sha256=digest, indices=[m['index'] for m in group['members']],
            reference_repeated_endpoints_exact=reference_equal,
            maximum_captured_repeated_endpoint_delta_angstrom=max(native_deltas, default=None)))
    phases = [p for r in reports for p in r['phases']]
    summary = dict(paths=len(reports), statuses=dict(Counter(r['status'] for r in reports)),
        phases=len(phases), callbacks=sum(p['callbacks'] for p in phases),
        exact_start_ownership=sum(r.get('exact_start_ownership', False) for r in reports),
        distinct_initial_coordinate_sets=len(groups),
        matches_hb_localization=sum(r['matches_hb_localization'] for r in reports),
        reference_reciprocal_duplicate_groups=sum(g['reference_repeated_endpoints_exact'] is not None for g in duplicates),
        reference_reciprocal_duplicate_failures=sum(g['reference_repeated_endpoints_exact'] is False for g in duplicates),
        prior_reference_controls=[dict(index=r['index'], start=r['exact_prior_reference_phase1_start'], endpoint=r['exact_prior_reference_phase1_endpoint'])
            for r in reports if 'exact_prior_reference_phase1_start' in r],
        negative_alpha_nodes=sum(p['negative_alpha_nodes'] for p in phases),
        returned_differs_last_reported=sum(not p['terminal']['returned_vs_last_reported']['exact'] for p in phases),
        maximum_mobile_gradient_readback_error_kj_mol_nm=max((p['maximum_mobile_gradient_readback_error_kj_mol_nm'] for p in phases if p['maximum_mobile_gradient_readback_error_kj_mol_nm'] is not None), default=None),
        maximum_objective_readback_error_kj_mol=max((p['maximum_objective_readback_error_kj_mol'] for p in phases if p['maximum_objective_readback_error_kj_mol'] is not None), default=None),
        final_cb_by_mode_mask={mode:{mask:dict(Counter(r.get('final_geometry', {}).get('target_cip', {}).get('CB', 'not_measured')
            for r in reports if r['mode'] == mode and r['mask'] == mask)) for mask in ('unchanged', 'HB', 'other-H')}
            for mode in ('captured', 'Reference')})
    common = dict(schema='dvbfixer.thr2-hswap-audit.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()), audit_script_sha256=sha(Path(__file__).read_bytes()),
        summary=summary, exact_input_groups=duplicates, errors=errors, production_changed=False,
        extra_minimization=False, limits=record['policy']['limits'])
    write(private_output, dict(common, cases=reports))
    public = []
    for row in reports:
        public.append(dict(row, phases=[{k:v for k, v in p.items() if k != 'nodes'} for p in row['phases']]))
    write(public_output, dict(common, cases=public, private_audit_sha256=sha(private_output.read_bytes())))
    check_plan(record)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'private-output', 'public-output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.private_output.resolve(), args.public_output.resolve())
