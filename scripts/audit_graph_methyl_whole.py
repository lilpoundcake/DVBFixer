"""Saved-only native whole-phase readback; no candidate/complete-geometry gate."""
from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

from scripts.audit_graph_methyl_control import h_close_pairs, ha_geometry
from scripts.audit_legacy_start_controls import coordinate_delta
from scripts.audit_legacy_stereo_trace import determinant
from scripts.audit_methyl_start_control import h_bonds
from scripts.audit_methyl_whole_control import hydrogen_angles
from scripts.preflight_graph_methyl_transport import mapped_start
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify

from dvbfixer.ffutils.dat import DatRecord


def global_shape(rows, rounded=False):
    groups = {}
    for row in rows:
        if row['element'] in {'H', 'D', None}:
            continue
        key = tuple(row['key'][:4])
        point = [float(f'{x:.3f}') for x in row['xyz']] if rounded else row['xyz']
        if row['key'][-1] in groups.setdefault(key, {}):
            raise ValueError('ambiguous whole-phase heavy identity')
        groups[key][row['key'][-1]] = point
    ca, beta = [], []
    for key, points in groups.items():
        if {'CA', 'N', 'C', 'CB'} <= points.keys():
            ca.append(dict(residue=key, volume_angstrom3=determinant(points, 'CA', 'N', 'C', 'CB')))
        if key[-1] == 'THR' and {'CB', 'CA', 'OG1', 'CG2'} <= points.keys():
            beta.append(dict(residue=key, volume_angstrom3=determinant(points, 'CB', 'CA', 'OG1', 'CG2')))
    return dict(CA_centres=len(ca), THR_beta_centres=len(beta),
        CA_signed_volumes=ca, THR_beta_signed_volumes=beta,
        CA_at_or_below_0_5=[c for c in ca if c['volume_angstrom3'] <= .5],
        THR_beta_at_or_below_0_5=[c for c in beta if c['volume_angstrom3'] <= .5],
        negative_CA=[c for c in ca if c['volume_angstrom3'] < 0],
        negative_THR_beta=[c for c in beta if c['volume_angstrom3'] < 0],
        scope='unchanged 0.5A3 global shape descriptor, NOT target/comprehensive chemistry/clash acceptance')


def audit(plan_path, root, baseline_root, output):
    import openmm

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    baseline_raw = (baseline_root/'evidence.json').read_bytes()
    run, baseline = json.loads(run_raw), json.loads(baseline_raw)
    if (run['plan_sha256'] != sha(raw) or baseline['plan_sha256'] != sha(raw)
            or run['baseline_run_sha256'] != sha(baseline_raw) or run['mode'] != 'continuations'
            or baseline['mode'] != 'native_baseline' or (root/'executed-plan.json').read_bytes() != raw
            or (baseline_root/'executed-plan.json').read_bytes() != raw):
        raise ValueError('whole baseline/continuation binding differs')
    expected = [(c['id'], c['seed'], arm) for c in plan['cases'] for arm in plan['policy']['arms']]
    if [(c['id'], c['seed'], c['arm']) for c in run['cases']] != expected:
        raise ValueError('whole continuation enrollment differs')
    if [(c['id'], c['seed']) for c in baseline['cases']] != [(c['id'], c['seed']) for c in plan['cases']]:
        raise ValueError('whole baseline enrollment differs')
    for cell in baseline['cases']:
        folder = baseline_root/cell['id']/f"seed-{cell['seed']}"
        for name, digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('baseline artifact changed')
    declarations = {(c['id'], c['seed']): c for c in plan['cases']}
    results = []
    for cell in run['cases']:
        case = declarations[(cell['id'], cell['seed'])]
        source = Path(case['source'])
        folder = root/cell['id']/f"seed-{cell['seed']}"/cell['arm']
        baseline_folder = baseline_root/cell['id']/f"seed-{cell['seed']}"
        for name, digest in cell['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('whole-phase artifact changed')
        report = {k: cell[k] for k in ('id', 'seed', 'arm', 'status', 'returncode')}
        if cell['status'] != 'completed':
            results.append(report)
            continue
        trace = json.loads((folder/'trace.json').read_bytes())
        if trace['case'] != case or trace['plan_sha256'] != sha(raw) or trace['arm'] != cell['arm'] or trace['emergency_reflection'] or trace['observations_restored']:
            raise ValueError('whole-phase trace contract differs')
        if [f['stage'] for f in trace['frames']] != ['initial', 'phase1', 'phase2']:
            raise ValueError('native phase inventory differs')
        topology = json.loads((source/'whole-topology.json').read_bytes())
        h_topology = json.loads((source/'H-topology.json').read_bytes())
        h_folder = Path(plan['h_root'])/case['id']/f"seed-{case['seed']}"/cell['arm']
        h_path = json.loads(gzip.decompress((h_folder/'path.json.gz').read_bytes()))
        expected_initial, repairs = mapped_start(h_topology['atoms'], h_topology['bonds'], h_path['final_nm'], topology['H_to_whole_indices'])
        actual_initial = json.loads(gzip.decompress((folder/'initial.json.gz').read_bytes()))['positions_nm']
        if expected_initial != actual_initial or repairs != trace['post_H_repairs']:
            raise ValueError('exact H-to-whole repair/lineage bridge differs')
        rows = topology['atoms']
        initial_rows = [dict(r, xyz=[10*x for x in p]) for r, p in zip(rows, actual_initial, strict=True)]
        added = DatRecord.load(source/'prepared.dat').added_keys()
        observed_indices = [i for i, r in enumerate(rows) if r['element'] not in {'H', 'D', None}
            and (r['key'][0], r['key'][1], r['key'][2], r['key'][-1]) not in added]
        report.update(exact_H_bridge=True, frames=[])
        for frame in trace['frames']:
            stage = frame['stage']
            frame_raw = (folder/f'{stage}.json.gz').read_bytes()
            system_raw = (folder/f'{stage}-system.xml').read_bytes()
            if sha(frame_raw) != frame['coordinate_sha256'] or sha(system_raw) != frame['system_sha256']:
                raise ValueError('whole frame snapshot hash differs')
            if system_raw != (baseline_folder/f'{stage}-system.xml').read_bytes():
                raise ValueError('whole forces/anchors differ from native baseline')
            data = json.loads(gzip.decompress(frame_raw))
            points = data['positions_nm']
            baseline_points = json.loads(gzip.decompress((baseline_folder/f'{stage}.json.gz').read_bytes()))['positions_nm']
            exact = points == baseline_points
            if exact != frame['exact_baseline_nm_state'] or (cell['arm'] == 'unchanged' and not exact):
                raise ValueError('unchanged native whole baseline not exact')
            material = [dict(r, xyz=[10*x for x in p]) for r, p in zip(rows, points, strict=True)]
            system = openmm.XmlSerializer.deserialize(system_raw.decode())
            report['frames'].append(dict(stage=stage, energy_kj_mol=data['energy_kj_mol'],
                raw_global_shape=global_shape(material), rounded_global_shape=global_shape(material, rounded=True),
                HA=ha_geometry(rows, topology['bonds'], points), H_bonds=h_bonds(system, rows, points),
                H_angles=hydrogen_angles(system, rows, points), close_H_pairs=h_close_pairs(rows, topology['bonds'], points),
                motion_from_full_start=coordinate_delta(initial_rows, material),
                original_heavy_motion=coordinate_delta([initial_rows[i] for i in observed_indices], [material[i] for i in observed_indices]),
                exact_native_baseline_nm_state=exact))
        results.append(report)
    summary = dict(enrolled=len(results), statuses=dict(Counter(c['status'] for c in results)), arms={})
    for arm in plan['policy']['arms']:
        cells = [c for c in results if c['status'] == 'completed' and c['arm'] == arm]
        final = [c['frames'][-1] for c in cells]
        summary['arms'][arm] = dict(completed=len(cells), final_all_acute_HA=sum(len(f['HA']['all_acute']) for f in final),
            final_negative_CA=sum(len(f['raw_global_shape']['negative_CA']) for f in final),
            final_CA_at_or_below_0_5=sum(len(f['raw_global_shape']['CA_at_or_below_0_5']) for f in final),
            final_THR_beta_at_or_below_0_5=sum(len(f['raw_global_shape']['THR_beta_at_or_below_0_5']) for f in final),
            exact_native_baseline_cells=sum(all(f['exact_native_baseline_nm_state'] for f in c['frames']) for c in cells))
    write(output, dict(schema='dvbfixer.graph-methyl-whole-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), baseline_run_sha256=sha(baseline_raw),
        audit_script_sha256=sha(Path(__file__).read_bytes()), policy=plan['policy'], production_changed=False,
        safety_net_executed=False, scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--baseline-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.baseline_root, args.output)
