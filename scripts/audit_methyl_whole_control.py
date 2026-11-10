"""Saved-only exact bridge/force/geometry audit of all methyl whole-system arms."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_ha_start_control import ha_geometry, observed_motion
from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.audit_methyl_start_control import h_bonds
from scripts.audit_start_geometry import geometry
from scripts.research_dunbrack_prior import sha, write
from scripts.research_h_stage_transport import post_h_repair
from scripts.research_ha_minimizer_path import verify


def hydrogen_angles(system, rows, points_nm):
    import openmm
    from openmm.unit import radian

    values = []
    for force in system.getForces():
        if not isinstance(force, openmm.HarmonicAngleForce):
            continue
        for index in range(force.getNumAngles()):
            a, b, c, reference, _ = force.getAngleParameters(index)
            if all(rows[i]['element'] not in {'H', 'D'} for i in (a, b, c)):
                continue
            first, second = ([x-y for x, y in zip(points_nm[i], points_nm[b])] for i in (a, c))
            norm = math.hypot(*first)*math.hypot(*second)
            if norm <= 1e-24 or not all(math.isfinite(x) for x in first+second):
                raise ValueError('undefined final hydrogen-containing angle')
            observed = math.degrees(math.acos(max(-1., min(1., sum(x*y for x, y in zip(first, second))/norm))))
            expected = math.degrees(reference.value_in_unit(radian))
            values.append(dict(atoms=[rows[i]['key'] for i in (a, b, c)], observed_degrees=observed,
                reference_degrees=expected, deviation_degrees=observed-expected))
    return dict(count=len(values), maximum_absolute_deviation_degrees=max((abs(v['deviation_degrees']) for v in values), default=None),
        rms_deviation_degrees=math.sqrt(sum(v['deviation_degrees']**2 for v in values)/len(values)) if values else None,
        largest_deviations=sorted(values, key=lambda v: -abs(v['deviation_degrees']))[:3])


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import nanometer

    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['source_arm'], arm) for c in plan['cases'] for arm in plan['arms']]
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw or [(r['id'], r['seed'], r['source_arm'], r['arm']) for r in run['cases']] != expected:
        raise ValueError('whole methyl cohort omitted/reordered cells or changed plan')
    declarations = {(c['id'], c['seed'], c['source_arm']): c for c in plan['cases']}
    results = []
    for record in run['cases']:
        case = declarations[(record['id'], record['seed'], record['source_arm'])]
        folder = root/case['id']/f"seed-{case['seed']}"/case['source_arm']/record['arm']
        for name, digest in record['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('retained whole methyl artifact changed')
        report = {k: record[k] for k in ('id', 'seed', 'source_arm', 'arm', 'status', 'returncode')}
        if record['status'] != 'completed':
            results.append(report)
            continue
        source, original = Path(case['source']), Path(case['original_whole'])
        trace = json.loads((folder/'trace.json').read_bytes())
        if trace['plan_sha256'] != sha(raw) or trace['case'] != case or trace['arm'] != record['arm'] or trace['emergency_reflection'] or trace['restored_observations']:
            raise ValueError('whole methyl trace contract differs')
        if [f['stage'] for f in trace['frames']] != ['initial', 'phase1', 'phase2']:
            raise ValueError('whole methyl native phase inventory differs')
        template = load_atoms(source/'h-internal-before.json.gz')
        bonds = json.loads((source/'topology.json').read_bytes())['bonds']
        h_folder = Path(plan['h_run'])/case['id']/f"seed-{case['seed']}"/case['source_arm']/record['arm']
        h_data = json.loads(gzip.decompress((h_folder/'path.json.gz').read_bytes()))
        points, repairs = post_h_repair(template, bonds, h_data['final_nm'])
        if repairs != trace['post_H_repairs']:
            raise ValueError('whole methyl repair bridge inventory differs')
        initial = load_atoms(folder/'initial.json.gz')
        expected_start = [dict(r, xyz=[10*x for x in p]) for r, p in zip(template, points, strict=True)]
        if any(not coordinate_delta(expected_start, initial)[k]['exact'] for k in ('heavy', 'hydrogen')):
            raise ValueError('whole-system initial coordinates differ from independently replayed post-H repair')
        if any(not coordinate_delta([r for r in template if r['element'] not in {'H', 'D'}],
                                    [r for r in initial if r['element'] not in {'H', 'D'}])[k]['exact'] for k in ('heavy',)):
            raise ValueError('H bridge moved an original heavy coordinate')
        if record['arm'] == 'unchanged':
            state = openmm.XmlSerializer.deserialize((original/'initial-state.xml').read_text())
            if points != [list(p) for p in state.getPositions().value_in_unit(nanometer)]:
                raise ValueError('whole unchanged nm start is not exact')
        report.update(post_H_repairs=repairs, exact_H_bridge=True, exact_heavy_bridge=True, frames=[])
        for frame in trace['frames']:
            stage = frame['stage']
            path, system_path = folder/f'{stage}.json.gz', folder/f'{stage}-system.xml'
            if sha(path.read_bytes()) != frame['atom_sha256'] or sha(system_path.read_bytes()) != frame['system_sha256']:
                raise ValueError('whole frame atom/System digest differs')
            system = openmm.XmlSerializer.deserialize(system_path.read_text())
            expected_path = original/('phase-2-before-system.xml' if stage == 'phase2' else 'phase-1-before-system.xml')
            if openmm.XmlSerializer.serialize(system) != openmm.XmlSerializer.serialize(openmm.XmlSerializer.deserialize(expected_path.read_text())):
                raise ValueError('whole native force/anchor/reference parameters differ')
            rows = load_atoms(path)
            reference = load_atoms(original/('phase-1-before.json.gz' if stage == 'initial' else f"phase-{1 if stage == 'phase1' else 2}-after.json.gz"))
            baseline = coordinate_delta(reference, rows)
            if baseline != frame['baseline_difference'] or (record['arm'] == 'unchanged' and any(not baseline[k]['exact'] for k in ('heavy', 'hydrogen'))):
                raise ValueError('whole baseline difference/unchanged endpoint not independently reproduced')
            points_nm = [[x/10 for x in r['xyz']] for r in rows]
            report['frames'].append(dict(stage=stage, energy_kj_mol=frame['energy_kj_mol'],
                raw=geometry(rows, bonds, case['target'], case['residue']),
                rounded=geometry(rows, bonds, case['target'], case['residue'], rounded=True),
                gln16=ha_geometry(rows), H_bonds=h_bonds(system, rows, points_nm),
                H_angles=hydrogen_angles(system, rows, points_nm),
                movement_from_whole_start=coordinate_delta(initial, rows),
                original_observation_motion=observed_motion(Path(case['observed_source']).read_bytes(), rows),
                baseline_endpoint_difference=baseline))
        results.append(report)
    completed = [r for r in results if r['status'] == 'completed']
    summary = dict(enrolled=len(results), statuses=dict(Counter(r['status'] for r in results)), arms={})
    for arm in plan['arms']:
        cells = [r for r in completed if r['arm'] == arm]
        summary['arms'][arm] = dict(completed=len(cells),
            final_raw_heavy_passes=sum(c['frames'][-1]['raw']['accepted'] for c in cells),
            final_rounded_heavy_passes=sum(c['frames'][-1]['rounded']['accepted'] for c in cells),
            final_gln16_nearplanar_or_non_L=sum(c['frames'][-1]['gln16']['ca_volume_angstrom3'] <= .5 for c in cells))
        if arm == 'unchanged':
            summary['arms'][arm]['exact_saved_endpoint_cells'] = sum(all(f['baseline_endpoint_difference'][k]['exact'] for f in c['frames'] for k in ('heavy', 'hydrogen')) for c in cells)
    write(output, dict(schema='dvbfixer.methyl-whole-control-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        driver_sha256=plan['source_hashes'][str(Path('scripts/research_methyl_whole_control.py').resolve())],
        transport_script_sha256=plan['source_hashes'][str(Path('scripts/research_h_stage_transport.py').resolve())],
        protocol=plan['policy'], openmm_version=plan['openmm_version'], production_changed=False,
        native_coordinates_read=False, scientific_qualification=False, source_observations_restored=False,
        original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
