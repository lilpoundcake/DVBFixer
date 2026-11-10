"""Independent saved-only audit of every declared one-HA paired replay.

No minimization, H regeneration, native-reference scoring or failure deletion.
Use unchanged actual-frame raw/rounded target and whole-context heavy gates.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.audit_legacy_stereo_trace import determinant, raw_rows
from scripts.audit_start_geometry import geometry, heavy_map
from scripts.audit_stereo_force_projection import verify
from scripts.research_dunbrack_prior import sha, write

HA = ('A', '16', '', 'GLN', 'HA')


def inventory(rows):
    mapped = {tuple(r['key']): r for r in rows}
    if len(mapped) != len(rows):
        raise ValueError('duplicate replay atom identity')
    if any(len(r['xyz']) != 3 or not all(math.isfinite(x) for x in r['xyz']) for r in rows):
        raise ValueError('nonfinite replay coordinate')
    return mapped


def ha_geometry(rows):
    mapped = inventory(rows)
    points = {name: mapped[('A', '16', '', 'GLN', name)]['xyz'] for name in ('CA', 'N', 'C', 'CB', 'HA')}
    vectors = {name: [x-y for x, y in zip(points[name], points['CA'])] for name in ('N', 'C', 'CB', 'HA')}
    lengths = {name: math.hypot(*v) for name, v in vectors.items()}
    if min(lengths.values()) <= 1e-12:
        raise ValueError('undefined replay GLN16 angle')
    units = {name: [x/lengths[name] for x in vector] for name, vector in vectors.items()}
    angles = {name: math.degrees(math.acos(max(-1., min(1., sum(x*y for x, y in zip(units[name], units['HA']))))))
              for name in ('N', 'C', 'CB')}
    return dict(ha_angles_degrees=angles, ca_ha_distance_angstrom=lengths['HA'],
        ca_volume_angstrom3=determinant(points, 'CA', 'N', 'C', 'CB'))


def intervention(before, after, arm):
    a, b = inventory(before), inventory(after)
    if a.keys() != b.keys() or any(a[k]['element'] != b[k]['element'] for k in a):
        raise ValueError('replay initial inventory differs')
    changed = [k for k in a if a[k]['xyz'] != b[k]['xyz']]
    if any(k != HA for k in changed) or (arm == 'unchanged' and changed):
        raise ValueError('undeclared replay initial coordinate change')
    start, end = ha_geometry(before), ha_geometry(after)
    if not math.isclose(start['ca_ha_distance_angstrom'], end['ca_ha_distance_angstrom'], rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError('one-HA replay changed CA-HA bond distance')
    if arm == 'outward_HA':
        ca = b[('A', '16', '', 'GLN', 'CA')]['xyz']
        unit_sum = [0., 0., 0.]
        for name in ('N', 'C', 'CB'):
            point = b[('A', '16', '', 'GLN', name)]['xyz']
            norm = math.dist(point, ca)
            unit_sum = [total+(x-y)/norm for total, x, y in zip(unit_sum, point, ca)]
        norm = math.hypot(*unit_sum)
        actual = [(x-y)/end['ca_ha_distance_angstrom'] for x, y in zip(b[HA]['xyz'], ca)]
        if norm <= 1e-12 or any(abs(x+y/norm) > 1e-10 for x, y in zip(actual, unit_sum)):
            raise ValueError('HA position differs from declared negative-unit-sum construction')
        if min(end['ha_angles_degrees'].values()) <= 90:
            raise ValueError('constructed HA is not outward to all heavy neighbours')
    return dict(changed_atoms=changed, initial_HA_movement_angstrom=math.dist(a[HA]['xyz'], b[HA]['xyz']),
        non_HA_coordinates_exact=True, before=start, after=end)


def observed_motion(source, rows):
    a, b = heavy_map(raw_rows(source)), heavy_map(rows)
    if a.keys()-b.keys():
        raise ValueError('original heavy observation missing from replay')
    distances = [math.dist(r['xyz'], b[key]['xyz']) for key, r in a.items()]
    return dict(observed_heavy_atoms=len(distances), maximum_angstrom=max(distances, default=None),
        rms_angstrom=math.sqrt(sum(d*d for d in distances)/len(distances)) if distances else None,
        exact=all(d == 0 for d in distances), observations_restored=False)


def audit(plan_path, root, output):
    import openmm

    plan_raw = plan_path.read_bytes()
    plan = json.loads(plan_raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(plan_raw) or (root/'executed-plan.json').read_bytes() != plan_raw:
        raise ValueError('replay executed plan differs')
    expected = [(c['id'], c['seed'], c['source_arm'], arm) for c in plan['cases'] for arm in plan['policy']['arms']]
    if [(r['id'], r['seed'], r['source_arm'], r['arm']) for r in run['cases']] != expected:
        raise ValueError('replay cohort omitted, duplicated or reordered cells')
    enrolled = {(c['id'], c['seed'], c['source_arm']): c for c in plan['cases']}
    results = []
    for row in run['cases']:
        case = enrolled[(row['id'], row['seed'], row['source_arm'])]
        source = Path(case['source'])
        folder = root/row['id']/f"source-seed-{row['seed']}"/row['source_arm']/row['arm']
        for name, digest in row['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('replay retained artifact changed')
        report = {k: row[k] for k in ('id', 'seed', 'source_arm', 'arm', 'status')}
        if row['status'] != 'completed':
            report['returncode'] = row['returncode']
            results.append(report)
            continue
        trace = json.loads((folder/'trace.json').read_bytes())
        if trace['plan_sha256'] != sha(plan_raw) or trace['arm'] != row['arm'] or plan['cases'][trace['case_index']] != case:
            raise ValueError('replay trace identity differs')
        if trace['platform'] != 'Reference' or not trace['no_integration'] or not trace['no_H_regeneration'] or trace['emergency_reflection']:
            raise ValueError('replay physics/operation trace differs')
        if [f['stage'] for f in trace['frames']] != ['initial', 'phase1', 'phase2']:
            raise ValueError('replay omitted/reordered saved frames')
        before = load_atoms(source/'phase-1-before.json.gz')
        initial = load_atoms(folder/'initial.json.gz')
        report['intervention'] = intervention(before, initial, row['arm'])
        change = json.loads((folder/'input-change.json').read_bytes())
        expected_allowed = [] if row['arm'] == 'unchanged' else [list(HA)]
        if change['arm'] != row['arm'] or change['allowed_atoms'] != expected_allowed or change['changed_atoms'] != [list(k) for k in report['intervention']['changed_atoms']]:
            raise ValueError('replay intervention declaration differs from independent coordinate audit')
        if change['original_state_sha256'] != sha((source/'initial-state.xml').read_bytes()) or tuple(before[change['HA_index']]['key']) != HA:
            raise ValueError('replay source-State digest or HA index differs')
        bonds = json.loads((source/'topology.json').read_bytes())['bonds']
        frames = []
        for frame in trace['frames']:
            stage = frame['stage']
            path = folder/f'{stage}.json.gz'
            system_path = folder/f'{stage}-system.xml'
            if sha(path.read_bytes()) != frame['sha256'] or sha(system_path.read_bytes()) != frame['system_sha256']:
                raise ValueError('replay frame/force digest differs')
            expected_system = source/('phase-2-before-system.xml' if stage == 'phase2' else 'phase-1-before-system.xml')
            serialize = openmm.XmlSerializer.serialize
            deserialize = openmm.XmlSerializer.deserialize
            if serialize(deserialize(system_path.read_text())) != serialize(deserialize(expected_system.read_text())):
                raise ValueError('replay force/anchor/reference parameters changed')
            rows = load_atoms(path)
            delta = coordinate_delta(initial, rows)
            material = dict(stage=stage, raw=geometry(rows, bonds, case['target'], case['residue']),
                rounded=geometry(rows, bonds, case['target'], case['residue'], rounded=True),
                gln16=ha_geometry(rows), movement_from_replay_start=delta,
                original_observation_motion=observed_motion(Path(case['original_masked_source']).read_bytes(), rows))
            if stage != 'initial':
                baseline = load_atoms(source/f"phase-{1 if stage == 'phase1' else 2}-after.json.gz")
                material['saved_endpoint_difference'] = coordinate_delta(baseline, rows)
            frames.append(material)
        report['frames'] = frames
        results.append(report)
    completed = [r for r in results if r['status'] == 'completed']
    summary = dict(enrolled=len(results), statuses=dict(Counter(r['status'] for r in results)), arms={})
    for arm in plan['policy']['arms']:
        cells = [r for r in completed if r['arm'] == arm]
        summary['arms'][arm] = dict(completed=len(cells),
            final_raw_passes=sum(r['frames'][-1]['raw']['accepted'] for r in cells),
            final_rounded_passes=sum(r['frames'][-1]['rounded']['accepted'] for r in cells),
            final_gln16_nearplanar_or_non_L=sum(r['frames'][-1]['gln16']['ca_volume_angstrom3'] <= .5 for r in cells))
        if arm == 'unchanged':
            summary['arms'][arm]['saved_endpoint_exact_cells'] = sum(all(
                f['saved_endpoint_difference'][kind]['exact'] for f in r['frames'][1:] for kind in ('heavy', 'hydrogen')) for r in cells)
            summary['arms'][arm]['maximum_saved_endpoint_difference_angstrom'] = max((
                f['saved_endpoint_difference'][kind]['maximum_angstrom'] for r in cells for f in r['frames'][1:]
                for kind in ('heavy', 'hydrogen')), default=None)
    pairs = []
    by_identity = {(r['id'], r['seed'], r['source_arm'], r['arm']): r for r in completed}
    for case in plan['cases']:
        arms = [by_identity.get((case['id'], case['seed'], case['source_arm'], arm)) for arm in ('unchanged', 'outward_HA')]
        if any(r is None for r in arms):
            continue
        first, second = [r['frames'][-1] for r in arms]
        pairs.append(dict(id=case['id'], seed=case['seed'], source_arm=case['source_arm'],
            unchanged_raw_accepted=first['raw']['accepted'], outward_raw_accepted=second['raw']['accepted'],
            unchanged_rounded_accepted=first['rounded']['accepted'], outward_rounded_accepted=second['rounded']['accepted'],
            unchanged_gln16_volume_angstrom3=first['gln16']['ca_volume_angstrom3'],
            outward_gln16_volume_angstrom3=second['gln16']['ca_volume_angstrom3']))
    summary['complete_pairs'] = len(pairs)
    write(output, dict(schema='dvbfixer.generated-ha-start-control-audit.v1', cases=results, summary=summary,
        pairs=pairs,
        plan_sha256=sha(plan_raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        driver_sha256=plan['source_hashes']['scripts/research_ha_start_control.py'],
        openmm_version=plan['openmm_version'], protocol=plan['policy'],
        production_changed=False, native_coordinates_read=False, scientific_qualification=False,
        original_thr2_causality_resolved=False, source_observations_restored=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
