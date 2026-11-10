"""Independent readback of saved stage coordinates and template/hybrid geometry.

No PDBFixer instance, coordinate generator, minimizer, reflection or dynamics is
called. Matrix determinants/triple products are computed directly from saved data.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

import pdbfixer.pdbfixer as pf
from scripts.branched_stereo_oracle import rdkit_labels
from scripts.branched_template_baseline import preflight
from scripts.research_rosetta_joint import sha
from scripts.research_rosetta_repack import write_json


def volume(points):
    n, ca, c, cb = (points[name] for name in ('N', 'CA', 'C', 'CB'))
    u, v, w = ([a-b for a, b in zip(p, ca, strict=True)] for p in (n, c, cb))
    return ((u[1]*v[2]-u[2]*v[1])*w[0]+(u[2]*v[0]-u[0]*v[2])*w[1]
            +(u[0]*v[1]-u[1]*v[0])*w[2])


def determinant(matrix):
    a, b, c = matrix
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def read_pdb(source):
    atoms = []
    for line in source.decode('ascii').splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            atoms.append(dict(key=[line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()],
                              xyz=[float(line[i:i+8]) for i in (30, 38, 46)]))
    return atoms


def points_equal(left, right, tolerance=1e-10):
    return left.keys() == right.keys() and all(math.dist(left[k], right[k]) <= tolerance for k in left)


def audit_points(points, recorded_volume, residue_name=None, recorded_cip=None):
    measured = volume(points)
    if not math.isclose(measured, recorded_volume, rel_tol=1e-12, abs_tol=1e-11):
        raise ValueError('independent signed volume differs')
    if residue_name is not None and rdkit_labels(points, residue_name) != recorded_cip:
        raise ValueError('independent target CIP differs')
    return measured


def audit(plan_path, run_root, output):
    plan_raw = plan_path.read_bytes()
    plan = json.loads(plan_raw)
    for filename, digest in plan['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('enrolled experiment implementation changed')
    if preflight() != plan['environment'] or sha(Path(pf.__file__).read_bytes()) != plan['installed_pdbfixer_sha256']:
        raise ValueError('installed scientific source/resources changed')
    if (run_root/'executed-plan.json').read_bytes() != plan_raw:
        raise ValueError('executed plan changed')
    raw = (run_root/'evidence.json').read_bytes()
    run = json.loads(raw)
    if run['plan_sha256'] != sha(plan_raw) or not run['installed_source_unchanged']:
        raise ValueError('run plan/source binding differs')
    expected_names = [case['id']+'-'+arm['overlay']+'-'+arm['relaxation'] for case in plan['cases'] for arm in plan['arms']]
    if set(expected_names) != set(run['results']):
        raise ValueError('factorial cell omitted or added')
    summaries, case_roles, seed_count, event_count, iteration_count = {}, {}, 0, 0, 0
    template_evidence = {}
    for declaration in plan['cases']:
        identifier, residue_name = declaration['id'], declaration['residue_name']
        source = Path(declaration['input_path']).read_bytes()
        if sha(source) != declaration['input_sha256'] or sha(Path(declaration['prior_proposal_path']).read_bytes()) != declaration['prior_proposal_sha256']:
            raise ValueError('original source/prior proposal changed')
        target = declaration['target']
        target_key = target['chain_id'], target['sequence_number'], target['insertion_code']
        observed = {tuple(a['key']): a['xyz'] for a in read_pdb(source)}
        target_observed = {key[3]: p for key, p in observed.items() if key[:3] == target_key}
        filenames = [name for name in plan['environment']['resource_sha256'] if Path(name).stem.upper() == residue_name]
        if len(filenames) != 1:
            raise ValueError('exact canonical template resource absent/ambiguous')
        template_path = Path(pf.__file__).parent/filenames[0]
        template_raw = template_path.read_bytes()
        if sha(template_raw) != plan['environment']['resource_sha256'][filenames[0]]:
            raise ValueError('original template resource changed')
        template_atoms = read_pdb(template_raw)
        template = {a['key'][3]: a['xyz'] for a in template_atoms}
        template_evidence[residue_name] = dict(resource=filenames[0], sha256=sha(template_raw),
                                             points_angstrom=template, ca_volume_angstrom3=volume(template))
        expected_q = [[x/10. for x in a['xyz']] for a in template_atoms if a['key'][3] in target_observed]
        prior = json.loads(Path(declaration['prior_proposal_path']).read_bytes())
        case_roles[identifier] = declaration['scope']
        for arm in plan['arms']:
            name = identifier+'-'+arm['overlay']+'-'+arm['relaxation']
            root = run_root/name
            cell_raw = (root/'evidence.json').read_bytes()
            cell = json.loads(cell_raw)
            if sha(cell_raw) != run['results'][name]['evidence_sha256'] or cell['arm'] != arm or not cell['hooks_restored']:
                raise ValueError('stage trace/arm/hook binding differs')
            if (root/'input.pdb').read_bytes() != source:
                raise ValueError('arm source bytes changed')
            seed_events = {}
            returned = []
            for event in cell['events']:
                coordinate_raw = (root/event['coordinates']).read_bytes()
                if sha(coordinate_raw) != event['coordinate_sha256']:
                    raise ValueError('stage coordinate digest differs')
                payload = json.loads(gzip.decompress(coordinate_raw))['atoms']
                coords = {tuple(a['key']): a['xyz'] for a in payload}
                if len(coords) != len(payload) or not observed.keys() <= coords.keys():
                    raise ValueError('stage atom identities lost/duplicated')
                movement = max(math.dist(coords[key], p) for key, p in observed.items())
                if movement > 1e-6 or not math.isclose(movement, event['observed_max_movement_angstrom'], abs_tol=1e-12):
                    raise ValueError('stage observation movement differs')
                points = {key[3]: p for key, p in coords.items() if key[:3] == target_key}
                if points != event['target_points_angstrom']:
                    raise ValueError('target projection differs from full stage')
                audit_points(points, event['ca_volume_angstrom3'], residue_name, event['target_cip'])
                seed_events.setdefault(event['seed'], []).append(event)
                if event['stage'] == 'returned-seed':
                    returned.append((event, payload))
                event_count += 1
            if list(seed_events) != list(range(1, len(seed_events)+1)) or len(seed_events) > 5:
                raise ValueError('helper seed inventory differs')
            seed_count += len(seed_events)
            if len(returned) != len(seed_events) or len(returned) != run['results'][name]['seed_count']:
                raise ValueError('missing returned seed state')
            if arm['relaxation'] == 'none' and cell['minimizer_iterations']:
                raise ValueError('no-relaxation control executed minimizer iterations')
            for step in cell['minimizer_iterations']:
                audit_points(step['target_points_angstrom'], step['ca_volume_angstrom3'])
                if any(math.dist(step['target_points_angstrom'][atom], xyz) > 1e-6 for atom, xyz in target_observed.items()):
                    raise ValueError('minimizer iterate moved observed target backbone')
                if step['observed_max_movement_angstrom'] > 1e-6:
                    raise ValueError('minimizer reporter moved observations')
                iteration_count += 1
            overlays = cell['overlays']
            if len(overlays) != 2*len(seed_events):
                raise ValueError('unexpected extra/missing construction/overlay call')
            hybrid_checks = []
            for overlay, construction in zip(overlays, [e for e in cell['events'] if e['stage'].startswith('construction-')], strict=True):
                matrix = overlay['rotation']
                det = determinant(matrix)
                if not math.isclose(det, overlay['determinant'], abs_tol=1e-12) or not math.isclose(abs(det), 1., abs_tol=1e-12):
                    raise ValueError('overlay determinant differs')
                if arm['overlay'] == 'proper' and det < 0:
                    raise ValueError('proper intervention still reflects')
                q = overlay['template_points_nm']
                if len(q) != len(expected_q) or any(math.dist(a, b) > 1e-12 for a, b in zip(q, expected_q, strict=True)):
                    raise ValueError('recorded anchors differ from original template resource/order')
                if overlay['seed'] != construction['seed']:
                    raise ValueError('overlay/seed construction order differs')
                transformed = {}
                for atom, p in template.items():
                    centered = [p[k]/10.+overlay['translate_template_nm'][k] for k in range(3)]
                    transformed[atom] = [10.*(sum(matrix[k][j]*centered[j] for j in range(3))+overlay['translate_reference_nm'][k]) for k in range(3)]
                hybrid = dict(transformed)
                hybrid.update(target_observed)
                actual = construction['target_points_angstrom']
                if not points_equal({k: hybrid[k] for k in actual}, actual):
                    raise ValueError('template+unchanged-observation hybrid does not reproduce construction')
                full_volume = volume(transformed)
                if not math.isclose(full_volume, det*volume(template), rel_tol=1e-11, abs_tol=1e-11):
                    raise ValueError('full-template volume fails independent determinant law')
                hybrid_checks.append(dict(seed=overlay['seed'], stage=construction['stage'], determinant=det,
                                           fully_transformed_template_volume_angstrom3=full_volume,
                                           actual_hybrid_volume_angstrom3=volume(actual),
                                           fully_transformed_template_cip=rdkit_labels(transformed, residue_name),
                                           actual_hybrid_cip=rdkit_labels(actual, residue_name)))
            if arm == dict(overlay='stock', relaxation='normal'):
                if not cell['stock_normal_exact_seed_reproduction'] or len(returned) != len(prior['snapshots']):
                    raise ValueError('untraced seed baseline not reproduced')
                for (event, payload), snapshot in zip(returned, prior['snapshots'], strict=True):
                    old = {(a['component']['chain_id'], a['component']['sequence_number'], a['component']['insertion_code'], a['name']): (element, p)
                           for a, element, p in snapshot['atoms']}
                    if {tuple(a['key']): (a['element'], a['xyz']) for a in payload} != old or event['seed'] != snapshot['seed']:
                        raise ValueError('original untraced seed coordinates differ')
            by_seed = []
            for seed, events in seed_events.items():
                heavy = next(e for e in events if e['stage'] == 'construction-heavy')
                final = next(e for e in events if e['stage'] == 'returned-seed')
                before = [e for e in events if e['stage'] == 'before-minimize']
                after = [e for e in events if e['stage'] == 'after-minimize']
                if len(before) != len(after):
                    raise ValueError('minimization boundary missing')
                for start, end in zip(before, after, strict=True):
                    if arm['relaxation'] == 'none' and start['target_points_angstrom'] != end['target_points_angstrom']:
                        raise ValueError('bypassed minimizer moved generated target')
                if arm['relaxation'] == 'none' and heavy['target_points_angstrom'] != final['target_points_angstrom']:
                    raise ValueError('no-relaxation final differs from actual construction')
                by_seed.append(dict(seed=seed, initial_volume=heavy['ca_volume_angstrom3'], initial_cip=heavy['target_cip'],
                                    returned_volume=final['ca_volume_angstrom3'], returned_cip=final['target_cip'],
                                    minimize_calls=len(before), dynamics_calls=sum(e['stage'] == 'before-dynamics' for e in events),
                                    minimization_boundaries=[dict(before_volume=a['ca_volume_angstrom3'], after_volume=b['ca_volume_angstrom3'])
                                                             for a, b in zip(before, after, strict=True)]))
            summaries[name] = dict(case_id=identifier, arm=arm, seed_count=len(seed_events), seeds=by_seed,
                                    hybrid_checks=hybrid_checks,
                                    identical_returned_seed_payloads=len({e['coordinate_sha256'] for e, _ in returned}) == 1,
                                    stock_normal_exact_seed_reproduction=cell['stock_normal_exact_seed_reproduction'],
                                    arm_evidence_sha256=sha(cell_raw))
    aggregate = dict(cases=len(plan['cases']), arms=len(summaries), returned_seeds=seed_count,
                     independently_checked_stage_snapshots=event_count, independently_checked_minimizer_iterates=iteration_count,
                     dynamics_step_calls=sum(s['dynamics_calls'] for c in summaries.values() for s in c['seeds']),
                     exact_untraced_baseline_arms=sum(c['stock_normal_exact_seed_reproduction'] is True for c in summaries.values()),
                     exact_untraced_baseline_seed_states=sum(c['seed_count'] for c in summaries.values() if c['stock_normal_exact_seed_reproduction'] is True))
    write_json(output, dict(schema='dvbfixer.template-causality-audit.v1', aggregate=aggregate, cases=summaries,
                            case_roles=case_roles, original_template_evidence=template_evidence,
                            experiment_plan=plan, experiment_plan_sha256=sha(plan_raw),
                            run_path=str(run_root), run_evidence_sha256=sha(raw),
                            audit_script_sha256=sha(Path(__file__).read_bytes()),
                            production_changed=False, scientific_qualification=False,
                            scope='retrospective case-specific initialization/relaxation ablations; not a production repair or guarantee'))
    print(json.dumps(aggregate, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
