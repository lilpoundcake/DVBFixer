"""Independent used-anchor/stage/common-geometry replay; never rebuild or relax.

The legacy trace records supplied anchors, not used anchors. Verify their separate
transport, proper-fit optimum, source observations and full-template/hybrid split.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

import numpy as np
import pdbfixer.pdbfixer as pf
from scripts.audit_template_causality import (
    audit_points,
    determinant,
    points_equal,
    read_pdb,
    volume,
)
from scripts.branched_stereo_oracle import rdkit_labels
from scripts.branched_template_baseline import preflight
from scripts.compare_branched_completion import measure
from scripts.research_rosetta_joint import sha
from scripts.research_rosetta_repack import write_json

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def fit_check(call, overlay, expected_names, reference, template, anchor_set):
    if anchor_set not in {'all', 'ncac'}:
        raise ValueError('undeclared anchor policy')
    used_indices = list(range(4)) if anchor_set == 'all' else [i for i, n in enumerate(expected_names) if n != 'O']
    if (call['supplied_names'] != expected_names or call['used_indices'] != used_indices
            or call['used_names'] != [expected_names[i] for i in used_indices]):
        raise ValueError('actual used-anchor identity/order differs from declaration')
    if (call['supplied_reference_points_nm'] != overlay['reference_points_nm']
            or call['supplied_template_points_nm'] != overlay['template_points_nm']):
        raise ValueError('supplied-anchor transport differs from tracer')
    for field in ('rotation', 'translate_template_nm', 'translate_reference_nm'):
        if call[field] != overlay[field]:
            raise ValueError('used-anchor transform differs from actual upstream transform')
    p = np.asarray(call['supplied_reference_points_nm'])
    q = np.asarray(call['supplied_template_points_nm'])
    if p.shape != (4, 3) or q.shape != (4, 3) or not np.allclose(p, reference, rtol=0., atol=1e-12) or not np.allclose(q, template, rtol=0., atol=1e-12):
        raise ValueError('supplied anchors do not match original observed/resource points')
    used_p, used_q = p[used_indices], q[used_indices]
    if call['used_reference_points_nm'] != used_p.tolist() or call['used_template_points_nm'] != used_q.tolist():
        raise ValueError('used coordinates differ from declared subset')
    matrix = np.asarray(call['rotation'])
    det = determinant(matrix.tolist())
    if not np.allclose(matrix @ matrix.T, np.eye(3), rtol=0., atol=1e-10) or not math.isclose(det, 1., abs_tol=1e-10):
        raise ValueError('anchor fit is not a proper orthogonal rotation')
    if not np.allclose(call['translate_template_nm'], -used_q.mean(axis=0), rtol=0., atol=1e-12) or not np.allclose(call['translate_reference_nm'], used_p.mean(axis=0), rtol=0., atol=1e-12):
        raise ValueError('centering does not use declared anchors')
    transformed_q = (matrix @ (used_q+np.asarray(call['translate_template_nm'])).T).T+call['translate_reference_nm']
    residual = float(np.sum((transformed_q-used_p)**2))
    pc, qc = used_p-used_p.mean(axis=0), used_q-used_q.mean(axis=0)
    covariance = qc.T @ pc
    singular = np.linalg.svd(covariance, compute_uv=False)
    optimum = float(np.sum(pc**2)+np.sum(qc**2)-2.*(singular[0]+singular[1]+(1. if np.linalg.det(covariance) >= 0 else -1.)*singular[2]))
    if not math.isclose(residual, optimum, rel_tol=1e-9, abs_tol=1e-11):
        raise ValueError('recorded rotation is not optimal for declared used anchors')
    return dict(used_names=call['used_names'], determinant=det, used_anchor_rmsd_angstrom=10.*math.sqrt(max(0., residual)/len(used_indices)),
                optimum_residual_nm2=optimum)


def audit(plan_path, run_root, output):
    plan_raw = plan_path.read_bytes()
    plan = json.loads(plan_raw)
    for filename, digest in plan['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('frozen anchor experiment implementation changed')
    if preflight() != plan['environment'] or sha(Path(pf.__file__).read_bytes()) != plan['installed_pdbfixer_sha256']:
        raise ValueError('scientific source/resources changed')
    if (run_root/'executed-plan.json').read_bytes() != plan_raw:
        raise ValueError('executed anchor plan differs')
    run_raw = (run_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    parent_raw = (Path(plan['parent_run_path'])/'evidence.json').read_bytes()
    if sha(Path(plan['parent_plan_path']).read_bytes()) != plan['parent_plan_sha256'] or sha(parent_raw) != plan['parent_run_sha256']:
        raise ValueError('parent experiment binding changed')
    parent = json.loads(parent_raw)
    if run['plan_sha256'] != sha(plan_raw) or not run['installed_source_unchanged']:
        raise ValueError('run plan/installed source binding differs')
    names = [c['id']+'-'+a['anchor_set']+'-'+a['relaxation'] for c in plan['cases'] for a in plan['arms']]
    if set(run['results']) != set(names):
        raise ValueError('factorial arm omitted or added')
    cases, full_snapshots, iterates, fits, seeds, repeated_arms, repeated_seeds = {}, 0, 0, 0, 0, 0, 0
    for declaration in plan['cases']:
        identifier = declaration['id']
        target = declaration['target']
        target_key = target['chain_id'], target['sequence_number'], target['insertion_code']
        source = Path(declaration['input_path']).read_bytes()
        if sha(source) != declaration['input_sha256']:
            raise ValueError('original masked source changed')
        observed = {tuple(a['key']): a['xyz'] for a in read_pdb(source)}
        target_observed = {key[3]: xyz for key, xyz in observed.items() if key[:3] == target_key}
        request = plan_branched_source(source, ComponentInstanceRef(**target)).request
        data = declaration['anchor_data']
        resource_raw = (Path(pf.__file__).parent/data['resource']).read_bytes()
        if sha(resource_raw) != data['resource_sha256']:
            raise ValueError('canonical template resource changed')
        resource_atoms = read_pdb(resource_raw)
        template = {a['key'][3]: a['xyz'] for a in resource_atoms}
        supplied_names = [a['key'][3] for a in resource_atoms if a['key'][3] in target_observed]
        expected_p = [[x/10. for x in target_observed[n]] for n in supplied_names]
        expected_q = [[x/10. for x in template[n]] for n in supplied_names]
        if data['names'] != supplied_names or data['reference_points_nm'] != expected_p or data['template_points_nm'] != expected_q:
            raise ValueError('enrolled anchor mapping differs from original input/template')
        for arm in plan['arms']:
            name = identifier+'-'+arm['anchor_set']+'-'+arm['relaxation']
            root = run_root/name
            cell_raw = (root/'evidence.json').read_bytes()
            transport_raw = (root/'anchor-policy.json').read_bytes()
            cell, transport = json.loads(cell_raw), json.loads(transport_raw)
            index = run['results'][name]
            if (sha(cell_raw) != index['trace_sha256'] or sha(transport_raw) != index['transport_sha256']
                    or transport['arm'] != arm or not transport['research_hook_restored'] or not cell['hooks_restored']
                    or cell['arm'] != dict(overlay='proper', relaxation=arm['relaxation'])):
                raise ValueError('trace/used-anchor transport/hook metadata differs')
            if (root/'input.pdb').read_bytes() != source:
                raise ValueError('arm input bytes changed')
            by_seed, payloads = {}, {}
            for event in cell['events']:
                raw = (root/event['coordinates']).read_bytes()
                if sha(raw) != event['coordinate_sha256']:
                    raise ValueError('stage coordinate digest differs')
                payload = json.loads(gzip.decompress(raw))['atoms']
                coords = {tuple(a['key']): a['xyz'] for a in payload}
                if len(coords) != len(payload) or not observed.keys() <= coords.keys():
                    raise ValueError('stage atom identities lost or duplicated')
                movement = max(math.dist(coords[k], p) for k, p in observed.items())
                if movement > 1e-6 or not math.isclose(movement, event['observed_max_movement_angstrom'], abs_tol=1e-12):
                    raise ValueError('observed source atom movement differs')
                points = {k[3]: p for k, p in coords.items() if k[:3] == target_key}
                if points != event['target_points_angstrom']:
                    raise ValueError('target points differ from complete stage')
                audit_points(points, event['ca_volume_angstrom3'], declaration['residue_name'], event['target_cip'])
                by_seed.setdefault(event['seed'], []).append(event)
                payloads[event['coordinates']] = raw
                full_snapshots += 1
            if list(by_seed) != list(range(1, len(by_seed)+1)) or len(by_seed) > 5 or len(by_seed) != index['seed_count']:
                raise ValueError('helper seed inventory differs')
            seeds += len(by_seed)
            if arm['relaxation'] == 'none' and cell['minimizer_iterations']:
                raise ValueError('no-relaxation control executed minimizer iterations')
            for step in cell['minimizer_iterations']:
                audit_points(step['target_points_angstrom'], step['ca_volume_angstrom3'])
                if step['observed_max_movement_angstrom'] > 1e-6 or any(math.dist(step['target_points_angstrom'][atom], p) > 1e-6 for atom, p in target_observed.items()):
                    raise ValueError('minimizer iterate moved observed atoms')
                iterates += 1
            constructions = [e for e in cell['events'] if e['stage'].startswith('construction-')]
            if len(constructions) != 2*len(by_seed) or len(transport['calls']) != len(constructions) or len(cell['overlays']) != len(constructions):
                raise ValueError('used-anchor construction transport omitted or added a fit')
            fit_reports = []
            for supplied, used, stage in zip(cell['overlays'], transport['calls'], constructions, strict=True):
                checked = fit_check(used, supplied, supplied_names, expected_p, expected_q, arm['anchor_set'])
                matrix = used['rotation']
                transformed = {}
                for atom, point in template.items():
                    centered = [point[k]/10.+used['translate_template_nm'][k] for k in range(3)]
                    transformed[atom] = [10.*(sum(matrix[k][j]*centered[j] for j in range(3))+used['translate_reference_nm'][k]) for k in range(3)]
                hybrid = dict(transformed)
                hybrid.update(target_observed)  # O stays an observation even when not used in fitting.
                actual = stage['target_points_angstrom']
                if not points_equal({atom: hybrid[atom] for atom in actual}, actual):
                    raise ValueError('used-anchor transform does not reproduce actual hybrid')
                if not math.isclose(volume(transformed), determinant(matrix)*volume(template), rel_tol=1e-11, abs_tol=1e-11):
                    raise ValueError('fully transformed template volume law differs')
                fit_reports.append(dict(checked, seed=stage['seed'], stage=stage['stage'],
                                        fully_transformed_template_volume=volume(transformed),
                                        actual_hybrid_volume=volume(actual), fully_transformed_template_cip=rdkit_labels(transformed, declaration['residue_name']),
                                        actual_hybrid_cip=stage['target_cip']))
                fits += 1
            if arm['anchor_set'] == 'all':
                parent_name = identifier+'-proper-'+arm['relaxation']
                parent_root = Path(plan['parent_run_path'])/parent_name
                if sha((parent_root/'evidence.json').read_bytes()) != parent['results'][parent_name]['evidence_sha256']:
                    raise ValueError('original parent trace digest differs')
                if not transport['all_anchor_parent_trace_exact'] or cell_raw != (parent_root/'evidence.json').read_bytes():
                    raise ValueError('complete four-anchor parent trace not reproduced exactly')
                for filename, payload in payloads.items():
                    if payload != (parent_root/filename).read_bytes():
                        raise ValueError('four-anchor parent full-stage payload differs')
                repeated_arms += 1
                repeated_seeds += len(by_seed)
            returns = []
            for seed, events in by_seed.items():
                initial = next(e for e in events if e['stage'] == 'construction-heavy')
                final = next(e for e in events if e['stage'] == 'returned-seed')
                before = [e for e in events if e['stage'] == 'before-minimize']
                after = [e for e in events if e['stage'] == 'after-minimize']
                if len(before) != len(after):
                    raise ValueError('missing minimization boundary')
                if arm['relaxation'] == 'none' and initial['target_points_angstrom'] != final['target_points_angstrom']:
                    raise ValueError('no-relaxation arm changed constructed target')
                points = final['target_points_angstrom']
                added = tuple(Addition(ExactAtomRef(request.residue, atom), atom[0], tuple(points[atom])) for atom in request.missing_atoms)
                measured, _ = measure(source, request, added, anchored=False)
                returns.append(dict(seed=seed, initial_volume=initial['ca_volume_angstrom3'], initial_cip=initial['target_cip'],
                                    returned_volume=final['ca_volume_angstrom3'], returned_cip=final['target_cip'],
                                    common_audit=dict(accepted=measured.accepted, raw=measured.raw, rounded=measured.rounded,
                                                      global_stereo_findings=measured.global_stereo_findings, rdkit_target_cip=measured.rdkit_target_cip,
                                                      observed_preserved=measured.observed_preserved),
                                    minimize_calls=len(before), dynamics_calls=sum(e['stage'] == 'before-dynamics' for e in events),
                                    minimization_boundaries=[dict(before_volume=a['ca_volume_angstrom3'], after_volume=b['ca_volume_angstrom3'],
                                                                  before_cip=a['target_cip'], after_cip=b['target_cip'])
                                                             for a, b in zip(before, after, strict=True)]))
            cases[name] = dict(case_id=identifier, arm=arm, seed_count=len(by_seed), seeds=returns, fit_checks=fit_reports,
                              identical_returned_seed_payloads=len({e['coordinate_sha256'] for e in cell['events'] if e['stage'] == 'returned-seed'}) == 1,
                              all_anchor_parent_trace_exact=transport['all_anchor_parent_trace_exact'],
                              trace_sha256=sha(cell_raw), transport_sha256=sha(transport_raw))
    normal_subset = [c for c in cases.values() if c['arm'] == dict(anchor_set='ncac', relaxation='normal')]
    aggregate = dict(cases=len(plan['cases']), arms=len(cases), returned_seeds=seeds, full_stage_snapshots=full_snapshots,
                     target_minimizer_iterates=iterates, used_anchor_fits=fits, exact_parent_arms=repeated_arms, exact_parent_seed_states=repeated_seeds,
                     dynamics_step_calls=sum(s['dynamics_calls'] for c in cases.values() for s in c['seeds']),
                     ncac_normal_initial_ca_S=sum(c['seeds'][0]['initial_cip'].get('CA') == 'S' for c in normal_subset),
                     ncac_normal_final_ca_S=sum(c['seeds'][-1]['returned_cip'].get('CA') == 'S' for c in normal_subset),
                     ncac_normal_final_common_passes=sum(c['seeds'][-1]['common_audit']['accepted'] for c in normal_subset))
    write_json(output, dict(schema='dvbfixer.template-anchor-audit.v1', aggregate=aggregate, cases=cases,
                            experiment_plan=plan, experiment_plan_sha256=sha(plan_raw), run_path=str(run_root), run_evidence_sha256=sha(run_raw),
                            audit_script_sha256=sha(Path(__file__).read_bytes()), production_changed=False, scientific_qualification=False,
                            scope='retrospective four-case anchor/relaxation controls, not a production remedy or independent accuracy cohort'))
    print(json.dumps(aggregate, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
