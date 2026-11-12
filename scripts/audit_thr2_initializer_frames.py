"""Independent NumPy replay of every frozen formula frame, no Context/optimizer."""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from scripts.audit_thr2_native_capture import read
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_thr2_initializer_frames import capture_data


def audit(plan_path, root, output):
    import numpy as np

    plan, execution = read(plan_path), read(root/'evidence.json')
    verify(plan)
    if execution['plan_sha256'] != sha(plan_path.read_bytes()):
        raise ValueError('Execution plan binding differs')
    for name, digest in execution['artifacts'].items():
        if sha((root/name).read_bytes()) != digest:
            raise ValueError('Private formula-frame artifact changed')
    if [dict(label=c['label'], rotation=c['rotation'], translation=c['translation']) for c in execution['cases']] != plan['cases']:
        raise ValueError('Case enrollment differs')
    data = {s['label']:capture_data(Path(s['source'])) for s in plan['sources']}
    reports, inverse_frames, maxima = [], {}, Counter()
    tolerance = plan['policy']['algebra_error_bound_nm']
    vbound = plan['policy']['volume_degenerate_bound_angstrom3']
    for case in execution['cases']:
        if case['status'] != 'completed':
            reports.append(dict(case, passed=False))
            continue
        source = data[case['label']]
        matrix = np.array(plan['rotations'][case['rotation']], dtype=float)
        shift = np.array(plan['translations_nm'][case['translation']], dtype=float)
        if not np.array_equal(matrix@matrix.T, np.eye(3)) or round(float(np.linalg.det(matrix))) != 1:
            raise ValueError('Not an exact proper signed-permutation rotation')
        heavy = {i:matrix@np.array(p)+shift for i, p in source['heavy'].items()}
        arm_reports = []
        if [r['arm'] for r in case['arms']] != plan['arms']:
            raise ValueError('Formula-arm enrollment differs')
        for arm in case['arms']:
            result = dict(arm)
            if arm['status'] != 'completed':
                result['passed'] = False
                arm_reports.append(result)
                continue
            payload = read(root/f"frame-{case['index']:03d}-{arm['arm']}.json.gz")
            poses = np.asarray(payload['positions_nm'])
            expected = np.zeros_like(poses)
            for index, point in heavy.items():
                expected[index] = point
            for spec in source['specs']:
                parent = heavy[spec['parent']]
                raw = np.sum([parent-heavy[j] for j in spec['neighbors']], axis=0)
                if np.linalg.norm(raw) <= 1e-12:
                    raise ValueError('Unadmitted singular independent formula')
                offset = 0.05*np.asarray(spec['draws'])
                if arm['arm'] == 'covariant-offset':
                    offset = matrix@offset
                elif arm['arm'] == 'zero-offset':
                    offset = np.zeros(3)
                delta = 0.1*raw/np.linalg.norm(raw)+offset
                expected[spec['index']] = parent+0.1*delta/np.linalg.norm(delta)
            error = float(np.max(np.linalg.norm(poses-expected, axis=1)))
            maxima['independent_formula_error_nm'] = max(maxima['independent_formula_error_nm'], error)
            restored = (poses-shift)@matrix
            inverse_frames[case['label'], case['rotation'], case['translation'], arm['arm']] = restored
            conflicts, degenerates, target, volume_error, heavy_invariant_error = [], [], None, 0., 0.
            for center in source['centers']:
                p, h, ns = center['parent'], center['index'], center['neighbors']
                hv = float(1000*np.linalg.det(poses[ns]-poses[p]))
                fv = float(1000*np.linalg.det(poses[ns]-poses[h]))
                original_heavy = float(1000*np.linalg.det(np.array([source['heavy'][j] for j in ns])-np.array(source['heavy'][p])))
                heavy_invariant_error = max(heavy_invariant_error, abs(hv-original_heavy))
                degenerate = abs(hv) <= vbound or abs(fv) <= vbound
                conflict = None if degenerate else hv*fv < 0
                if conflict:
                    conflicts.append(source['atoms'][h]['key'])
                if degenerate:
                    degenerates.append(source['atoms'][h]['key'])
                if h == source['hb']['index']:
                    target = dict(heavy_volume_angstrom3=hv, four_ligand_volume_angstrom3=fv,
                                  sign_conflict=conflict, degenerate=degenerate)
            if target is None:
                raise ValueError('Target was not independently measured')
            volume_error = max(abs(target[k]-arm['target'][k]) for k in ('heavy_volume_angstrom3', 'four_ligand_volume_angstrom3'))
            result.update(independent_formula_maximum_error_nm=error,
                independent_target_volume_error_angstrom3=volume_error,
                heavy_volume_invariance_error_angstrom3=heavy_invariant_error,
                passed=error <= tolerance and volume_error <= vbound and heavy_invariant_error <= vbound
                    and target['sign_conflict'] == arm['target']['sign_conflict']
                    and conflicts == arm['conflicting_carbon_hydrogens']
                    and degenerates == arm['degenerate_carbon_hydrogens']
                    and arm['maximum_heavy_backtransform_error_nm'] <= tolerance
                    and (arm['arm'] == 'fixed-lab' or arm['maximum_full_h_backtransform_error_nm'] <= tolerance)
                    and (arm['arm'] == 'zero-offset' or arm['exact_identity_formula'] is not False))
            arm_reports.append(result)
        reports.append(dict(index=case['index'], label=case['label'], rotation=case['rotation'],
                            translation=case['translation'], arms=arm_reports,
                            passed=all(a['passed'] for a in arm_reports)))
    translation_errors = []
    identity_rotation = plan['rotations'].index([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    covariance_checks = []
    for row in reports:
        for arm in row.get('arms', []):
            if arm['status'] != 'completed':
                continue
            current = inverse_frames[row['label'], row['rotation'], row['translation'], arm['arm']]
            reference_arm = 'zero-offset' if arm['arm'] == 'zero-offset' else 'fixed-lab'
            reference = inverse_frames.get((row['label'], identity_rotation, 0, reference_arm))
            if reference is None:
                arm['passed'] = False
                covariance_checks.append(dict(label=row['label'], rotation=row['rotation'],
                    translation=row['translation'], arm=arm['arm'], passed=False))
                continue
            error = float(np.max(np.linalg.norm(current-reference, axis=1)))
            passed = arm['arm'] == 'fixed-lab' or error <= tolerance
            arm['passed'] = arm['passed'] and passed
            covariance_checks.append(dict(label=row['label'], rotation=row['rotation'],
                translation=row['translation'], arm=arm['arm'], independent_maximum_error_nm=error,
                required_covariant=arm['arm'] != 'fixed-lab', passed=passed))
        if row.get('arms'):
            row['passed'] = all(arm['passed'] for arm in row['arms'])
    for label in data:
        for rotation in range(24):
            for arm in plan['arms']:
                key = label, rotation, 0, arm
                other = label, rotation, 1, arm
                if key not in inverse_frames or other not in inverse_frames:
                    translation_errors.append(dict(label=label, rotation=rotation, arm=arm, passed=False))
                    continue
                delta = float(np.max(np.linalg.norm(inverse_frames[key]-inverse_frames[other], axis=1)))
                translation_errors.append(dict(label=label, rotation=rotation, arm=arm,
                                               maximum_error_nm=delta, passed=delta <= tolerance))
    counts = []
    for label in data:
        for arm in plan['arms']:
            selected = [a for row in reports if row['label'] == label and row['translation'] == 0
                        for a in row.get('arms', []) if a['arm'] == arm and a['status'] == 'completed']
            counts.append(dict(label=label, arm=arm, rotations=len(selected),
                target_sign_conflicts=sum(a['target']['sign_conflict'] is True for a in selected),
                target_degenerate=sum(a['target']['degenerate'] for a in selected),
                carbon_hydrogen_conflicts=dict(Counter('/'.join(key) for a in selected for key in a['conflicting_carbon_hydrogens']))))
    summary = dict(frames=len(reports), formula_arm_records=sum(len(r.get('arms', [])) for r in reports),
        admitted_sources=len(data), algebra_controls_passed=all(r['passed'] for r in reports)
            and all(r['passed'] for r in translation_errors), native_calls=0, new_random_sampling=False,
        maximum_independent_formula_error_nm=maxima['independent_formula_error_nm'],
        maximum_translation_pair_error_nm=max((r.get('maximum_error_nm', 0.) for r in translation_errors), default=0.))
    write(output, dict(schema='dvbfixer.thr2-initializer-frames-audit.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()), audit_script_sha256=sha(Path(__file__).read_bytes()),
        summary=summary, orientation_counts=counts, cases=reports, translation_checks=translation_errors,
        independent_covariance_checks=covariance_checks,
        initialization_checks=execution['initialization_checks'], production_changed=False,
        limits=plan['policy']['limits']))
    print(summary)
    print(counts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
