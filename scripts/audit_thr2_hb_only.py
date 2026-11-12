"""Independent saved native-path gradient/geometry and HB/bridge identity audit."""
from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path

from scripts.audit_thr2_native_capture import audit_replay, read, verify_files
from scripts.qualify_hydrogen_postflight import load, reconstruct
from scripts.research_dunbrack_prior import sha, write
from scripts.research_thr2_hb_only import check_plan, proposed_points, transport
from scripts.research_thr2_initializer_frames import capture_data


def audit(plan_path, root, output):
    import openmm
    from openmm import Vec3
    from openmm.unit import nanometer

    plan, execution = read(plan_path), read(root/'evidence.json')
    check_plan(plan)
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('Frozen execution plan differs')
    if [{k:c[k] for k in ('label', 'arm', 'mode')} for c in execution['cases']] != plan['cases']:
        raise ValueError('Case enrollment differs')
    old_geometry = load('thr2_hb_only_audit_historical', Path(plan['historical_geometry']))
    reports = []
    for case in execution['cases']:
        folder = root/f"path-{case['index']:02d}"
        verify_files(folder, case)
        source = next(s for s in plan['sources'] if s['label'] == case['label'])
        source_folder = Path(source['source'])
        data = capture_data(source_folder)
        sibling_pairs = [(a['index'], b['index']) for i, a in enumerate(data['specs'])
                         for b in data['specs'][i+1:] if a['parent'] == b['parent']]
        old_initial = openmm.XmlSerializer.deserialize((source_folder/f"{source['calls'][0]['before']['prefix']}-state.xml").read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        wanted, selected = proposed_points(data, old_initial, case['arm'])
        report = {k:v for k, v in case.items() if k != 'artifacts'}
        phase_reports, paths = [], []
        for number, call in enumerate(source['calls'], 1):
            phase = folder/f'phase-{number}'
            if not (phase/'path.json.gz').exists():
                continue
            if (phase/'system.xml').read_bytes() != (source_folder/f"{call['before']['prefix']}-system.xml").read_bytes():
                raise ValueError('Native phase System bytes differ')
            path = read(phase/'path.json.gz')
            paths.append(path)
            requested = read(phase/'requested.json.gz')['positions_nm']
            transport(requested, path['initial_nm'], case['mode'])
            physical = audit_replay(plan, source_folder, call, phase,
                dict(label=case['label'], number=number, mode=case['mode'], status=case['status']))
            physical['exact_phase_system'] = True
            # The four-ligand volume is separate from the physical oracle's
            # heavy-center/CIP classification and is reported at every node.
            rows = read(source_folder/f"{call['before']['prefix']}-atoms.json.gz")['atoms']
            mapping = {tuple(r['key']):i for i, r in enumerate(rows)}
            ids = [mapping[('A', '2', '', 'THR', name)] for name in ('HB', 'CA', 'OG1', 'CG2')]
            from scripts.audit_stereo_force_projection import volume_gradient
            points = [path['initial_nm']]+[
                [cb['positions_nm'][i:i+3] for i in range(0, len(cb['positions_nm']), 3)] for cb in path['callbacks']]+[path['final_nm']]
            for node, coordinates in zip(physical['nodes'], points, strict=True):
                v, _, _ = volume_gradient([coordinates[i] for i in ids])
                node['hb_four_ligand_volume_angstrom3'] = 1000*v
                node['hb_heavy_sign_conflict'] = v*node['geometry']['target_cb_volume_angstrom3'] < 0
                oxygen_h = [10*math.dist(coordinates[s['index']], coordinates[s['parent']])
                            for s in data['specs'] if data['atoms'][s['parent']]['element'] == 'O']
                node['maximum_oh_distance_angstrom'] = max(oxygen_h, default=None)
                node['minimum_same_parent_h_distance_angstrom'] = min(
                    (10*math.dist(coordinates[a], coordinates[b]) for a, b in sibling_pairs), default=None)
            phase_reports.append(physical)
        report['phases'] = phase_reports
        report['exact_initial_requested_mask'] = bool(paths) and read(folder/'phase-1/requested.json.gz')['positions_nm'] == wanted
        report['initial_unselected_exact'] = bool(paths) and all(p == q for i, (p, q) in enumerate(zip(old_initial, paths[0]['initial_nm'], strict=True)) if i not in selected)
        if paths:
            h_points = [paths[0]['initial_nm'], paths[0]['final_nm']]+[
                [cb['positions_nm'][i:i+3] for i in range(0, len(cb['positions_nm']), 3)] for cb in paths[0]['callbacks']]
            report['h_stage_heavy_exact'] = all(p[i] == paths[0]['initial_nm'][i] for p in h_points for i in data['heavy'])
        report['reference_original_h_control_exact'] = None
        if case['arm'] == 'unchanged' and case['mode'] == 'Reference' and paths:
            report['reference_original_h_control_exact'] = paths[0]['final_nm'] == read(Path(source['reference_h_control']))['final_nm']
        report['bridge'] = None
        if len(paths) >= 2:
            saved = read(folder/'historical-bridge.json.gz')
            topology, _ = reconstruct(data['atoms'], read(source_folder/'call-1-topology.json')['bonds'])
            positions = [Vec3(*p) for p in paths[0]['final_nm']]*nanometer
            count = old_geometry.repair_misplaced_hydrogens(topology, positions)
            expected = [list(p) for p in positions.value_in_unit(nanometer)]
            report['bridge'] = dict(repairs=count,
                exact=saved['before_nm'] == paths[0]['final_nm'] and saved['after_nm'] == expected
                    and count == saved['repairs'] and read(folder/'phase-2/requested.json.gz')['positions_nm'] == expected,
                heavy_exact=all(expected[i] == paths[0]['final_nm'][i] for i in data['heavy']))
        report['own_phase2_continuity_exact'] = len(paths) == 3 and paths[2]['initial_nm'] == paths[1]['final_nm']
        report['audit_passed'] = (case['status'] == 'completed' and len(paths) == 3
            and report['exact_initial_requested_mask'] and report['initial_unselected_exact']
            and report['h_stage_heavy_exact'] and report['reference_original_h_control_exact'] is not False
            and report['bridge']['exact'] and report['bridge']['heavy_exact']
            and report['own_phase2_continuity_exact'] and all(p['returned'] for p in phase_reports))
        reports.append(report)
    phases = [p for r in reports for p in r['phases']]
    summary = dict(paths=len(reports), statuses=dict(Counter(r['status'] for r in reports)),
        native_phases=len(phases), callbacks=sum(p['callbacks'] for p in phases),
        accounting_and_identity_passed=all(r['audit_passed'] for r in reports),
        exact_reference_h_controls=sum(r['reference_original_h_control_exact'] is True for r in reports),
        negative_alpha_nodes=sum(p['negative_alpha_nodes'] for p in phases),
        terminal_position_mismatches=sum(not p['terminal']['returned_vs_last_reported']['exact'] for p in phases),
        maximum_mobile_gradient_readback_error_kj_mol_nm=max((p['maximum_mobile_gradient_readback_error_kj_mol_nm'] or 0.) for p in phases),
        maximum_objective_readback_error_kj_mol=max((p['maximum_objective_readback_error_kj_mol'] or 0.) for p in phases))
    write(output, dict(schema='dvbfixer.thr2-hb-only-audit.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()), script_sha256=sha(Path(__file__).read_bytes()),
        summary=summary, cases=reports, production_changed=False, new_random_sampling=False,
        limits=plan['protocol']['limits']))
    print(summary)
    for row in reports:
        print(row['label'], row['arm'], row['mode'], [p['terminal']['returned_geometry']['target_cip'] for p in row['phases']])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
