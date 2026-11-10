"""Independent saved-candidate readback for the separately declared nine-chi grid.

No coordinate generator, minimizer, packer or truth-based selector is called.
Frozen dense geometry/materialization is replayed; original consumers unchanged.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import math
import statistics
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from scripts.audit_rosetta_joint import norm, report_coordinates, restore_candidate
from scripts.research_rosetta_dense import CHIS, audit_target, joint_audits, materialize
from scripts.research_rosetta_joint import context_signature, merge_views, sha, verify
from scripts.research_rosetta_repack import identity, observed_atoms, write_json
from scripts.rosetta_dense_diagnostics import diagnose
from scripts.rosetta_fixed_contract import choose_scored_own, import_audit
from scripts.summarize_rosetta_fixed_scoring import source_records_preserved

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.protein_completion.cohort import _enumerate_joint
from dvbfixer.protein_completion.contracts import Attempt, CompletionResult
from dvbfixer.protein_completion.mixed_cohort import plan_mixed_cohort
from dvbfixer.protein_completion.source import _insert_candidates


def replay(source, case, saved):
    plan = plan_mixed_cohort(source, tuple(ComponentInstanceRef(**t['residue']) for t in case['targets']))
    results = []
    for target, raw in zip(plan.targets, saved['targets'], strict=True):
        expected = 1 if target.mode == 'anchored' else len(CHIS)
        if len(raw['attempts']) != expected or raw['request_digest'] != target.request.digest:
            raise ValueError('target attempt count or source binding differs')
        attempts = []
        for i, row in enumerate(raw['attempts']):
            candidate = restore_candidate(row['candidate'])
            if target.mode != 'anchored' and candidate.torsion_degrees != CHIS[i]:
                raise ValueError('dense grid order/membership differs')
            audit = audit_target(target, candidate)
            if norm(asdict(audit)) != row['audit']:
                raise ValueError('saved candidate audit differs')
            attempts.append(Attempt(candidate, audit))
        passing = [i for i, a in enumerate(attempts) if a.audit.accepted]
        selected = max(passing, key=lambda i: min(0., attempts[i].audit.minimum_nonbonded_clearance_angstrom or 0.)) if passing else None
        if selected != raw['selected_index']:
            raise ValueError('independent clearance winner differs')
        results.append(CompletionResult(**{k: v for k, v in raw.items() if k != 'attempts'}, attempts=tuple(attempts)))
    result = _enumerate_joint(tuple(results), lambda candidates: joint_audits(plan, candidates), plan.digest,
                              max_combinations=81, completed_reason='private_nine_chi_joint')
    if norm(asdict(result)) != saved:
        raise ValueError('complete dense joint enumeration differs')
    pools = [[i for i, a in enumerate(r.attempts) if a.audit.accepted] for r in results]
    if [tuple(a.indices) for a in result.attempts] != list(itertools.product(*pools)):
        raise ValueError('incomplete Cartesian pool')
    return plan, result


def check_scored_cell(plan_root, score_root, case, cell, environment):
    folder = score_root/case['id']/cell['id']
    if json.loads((folder/'score.json').read_bytes()) != cell:
        raise ValueError('cell record differs from evidence index')
    pdb = (plan_root/cell['pdb']).read_bytes()
    records, snapshots = [], []
    for i, view in enumerate(case['views']):
        root = folder/f'view-{i}'
        record = json.loads((root/'score.json').read_bytes())
        payload = (root/'pose-snapshot.json.gz').read_bytes()
        if sha(payload) != record['snapshot_sha256']:
            raise ValueError('snapshot digest differs')
        snapshot = json.loads(gzip.decompress(payload))
        audit = import_audit(observed_atoms(pdb), snapshot['heavy_atoms'], view['expected_external'], snapshot['external'],
                             view['target_internal'], snapshot['target_internal'], snapshot['allowed_extra'])
        if norm(audit) != record['import_audit'] or not audit['accepted']:
            raise ValueError('exact physical-heavy/graph audit failed')
        termini = {(r['chain'], r['resid'], r['icode']) for r in snapshot['residues'] if 'CtermProteinFull' in r['name']}
        outgoing = {tuple(a if a[-1] == 'C' else b)[:4] for a, b in view['expected_external'] if {a[-1], b[-1]} == {'C', 'N'}}
        if any(k[-1] != 'OXT' or tuple(k[:3]) not in termini or tuple(k[:4]) in outgoing for k in snapshot['allowed_extra']):
            raise ValueError('physical-heavy exception is not terminal OXT')
        if any(not a['virtual_by_rosetta_type'] for a in snapshot['virtual_atoms']):
            raise ValueError('virtual atom lacks type evidence')
        if {identity(a) for a in snapshot['heavy_atoms']} & {identity(a) for a in snapshot['virtual_atoms']}:
            raise ValueError('physical and virtual inventories overlap')
        terms = record['score_terms']
        if {k: t['weight'] for k, t in terms.items()} != environment['weights']:
            raise ValueError('score weights changed')
        if not math.isclose(sum(t['raw']*t['weight'] for t in terms.values()), record['total_score'], abs_tol=1e-8, rel_tol=1e-10):
            raise ValueError('weighted score sum differs')
        records.append(record)
        snapshots.append(snapshot)
    merged = merge_views(records, snapshots)
    if records != cell['views']:
        raise ValueError('target-view score records differ from cell index')
    if merged['score_status'] != 'scored' or merged['total_score'] != cell['total_score'] or cell['score_status'] != 'scored':
        raise ValueError('full-pose two-target score merge failed')
    signature = context_signature(snapshots[0], case['views'])
    if signature != cell['context_signature']:
        raise ValueError('fixed context signature differs')
    return signature


def audit(plan_root, score_root, coarse_score_root, output):
    raw_plan = (plan_root/'plan.json').read_bytes()
    record = json.loads(raw_plan)
    verify(record, plan_root)
    enrollment = record['enrollment']
    executed_enrollment = (plan_root/'executed-enrollment.json').read_bytes()
    if sha(executed_enrollment) != record['enrollment_sha256'] or json.loads(executed_enrollment) != enrollment:
        raise ValueError('executed enrollment binding differs')
    parent_root = Path(enrollment['parent_plan_path'])
    parent_raw = (parent_root/'plan.json').read_bytes()
    if sha(parent_raw) != enrollment['parent_plan_sha256']:
        raise ValueError('parent plan changed')
    parent = json.loads(parent_raw)
    verify(parent, parent_root)
    coarse_raw = (coarse_score_root/'evidence.json').read_bytes()
    coarse_scores = json.loads(coarse_raw)
    if coarse_scores['plan_sha256'] != sha(parent_raw) or (coarse_score_root/'executed-plan.json').read_bytes() != parent_raw:
        raise ValueError('coarse score execution differs from unchanged parent plan')
    if [c['id'] for c in record['cases']] != enrollment['case_ids'] or enrollment['chis'] != list(CHIS):
        raise ValueError('case list or grid differs from enrollment')
    raw_scores = (score_root/'evidence.json').read_bytes()
    scores = json.loads(raw_scores)
    if set(scores['results']) != set(enrollment['case_ids']):
        raise ValueError('score evidence omitted or added a frozen case')
    environment_raw = (score_root/'environment.json').read_bytes()
    environment = json.loads(environment_raw)
    if (score_root/'executed-plan.json').read_bytes() != raw_plan or scores['plan_sha256'] != sha(raw_plan):
        raise ValueError('executed plan binding differs')
    if enrollment['score_policy']['expected_rosetta_revision'] not in environment['version']:
        raise ValueError('Rosetta version changed')
    cases, counts, coupling_counts, changes, factorized_changes = {}, Counter(), Counter(), [], []
    coarse_errors, dense_errors = [], []
    for case, old_case in zip(record['cases'], parent['cases'], strict=True):
        name = case['id']
        for field in ('id', 'targets', 'views', 'generated', 'selection_geometry'):
            if case[field] != old_case[field]:
                raise ValueError('frozen source/target view differs')
        root = plan_root/name
        for filename in ('source.pdb', 'reference.pdb', 'native-control.pdb'):
            if (root/filename).read_bytes() != (parent_root/name/filename).read_bytes():
                raise ValueError('source/reference bytes changed')
        source, reference = (root/'source.pdb').read_bytes(), (root/'reference.pdb').read_bytes()
        if case['status'] != 'planned':
            cases[name] = dict(status=case['status'], reason=case['reason'])
            continue
        saved = json.loads((root/'reconstruction.json').read_bytes())
        plan, result = replay(source, case, saved)
        if norm(asdict(plan)) != json.loads((root/'source-plan.json').read_bytes()):
            raise ValueError('source plan readback differs')
        old_result = json.loads((parent_root/name/'reconstruction.json').read_bytes())
        for new, old in zip(result.targets, old_result['targets'], strict=True):
            if norm([asdict(a) for a in new.attempts[:len(old['attempts'])]]) != old['attempts']:
                raise ValueError('exact coarse subset differs')
        row = scores['results'][name]
        if [r['id'] for r in row['cells']] != [c['id'] for c in case['cells']]:
            raise ValueError('score table omitted a candidate or native control')
        choice = choose_scored_own([r for r in row['cells'] if r['arm'] == 'own'])
        if choice != row['selection'] or choice != json.loads((score_root/name/'selection-before-native.json').read_bytes()):
            raise ValueError('selection differs from complete feasible pool')
        reports, signatures = {}, []
        for cell, declared in zip(row['cells'], case['cells'], strict=True):
            if any(cell[k] != v for k, v in declared.items()):
                raise ValueError('scored candidate declaration changed')
            pdb = (plan_root/cell['pdb']).read_bytes()
            if not source_records_preserved(source, pdb, case['generated']):
                raise ValueError('original source records changed')
            if cell['arm'] == 'own':
                attempt = result.attempts[cell['index']]
                candidates = tuple(t.attempts[i].candidate for t, i in zip(result.targets, attempt.indices, strict=True))
                if _insert_candidates(source, tuple(t.request for t in plan.targets), candidates) != pdb:
                    raise ValueError('scored coordinates differ from saved proposals')
                rounded_ok = False
                if attempt.accepted:
                    try:
                        checked, dat = materialize(source, plan, candidates)
                    except ValueError:
                        counts['rounded_rejections'] += 1
                    else:
                        rounded_ok = True
                        loaded = DatRecord.load((plan_root/cell['pdb']).with_suffix('.dat'))
                        if checked != pdb or loaded.to_dict() != dat.to_dict():
                            raise ValueError('canonical private generated-only materialization differs')
                        expected = {(a[0], str(a[1]), a[2].strip(), a[4]) for a in case['generated']}
                        if loaded.added_keys() != expected or len(loaded.added_atoms) != len(expected):
                            raise ValueError('dat owns observations or omits generated atoms')
                        counts['dat_records'] += 1
                else:
                    counts['raw_joint_rejections'] += 1
                if cell['old_source_accepted'] != attempt.accepted or cell['common_accepted'] != rounded_ok:
                    raise ValueError('raw/rounded eligibility differs')
                if not rounded_ok and (plan_root/cell['pdb']).with_suffix('.dat').exists():
                    raise ValueError('rejected candidate owns a sidecar')
                counts['own_combinations'] += 1
            else:
                counts['native_controls'] += 1
            signatures.append(check_scored_cell(plan_root, score_root, case, cell, environment))
            if cell['arm'] == 'native_control' or all(i < 3 for i in cell['indices']):
                previous = next(c for c in coarse_scores['results'][name]['cells']
                                if c['arm'] == cell['arm'] and c['indices'] == cell['indices'])
                if (parent_root/previous['pdb']).read_bytes() != pdb or previous['total_score'] != cell['total_score']:
                    raise ValueError('unchanged coarse/native pose or full score does not reproduce')
                for i in range(2):
                    old_snapshot = json.loads(gzip.decompress((coarse_score_root/name/previous['id']/f'view-{i}'/'pose-snapshot.json.gz').read_bytes()))
                    new_snapshot = json.loads(gzip.decompress((score_root/name/cell['id']/f'view-{i}'/'pose-snapshot.json.gz').read_bytes()))
                    if old_snapshot != new_snapshot:
                        raise ValueError('unchanged coarse/native import context does not reproduce')
                counts['coarse_native_exact_repeats'] += 1
            reports[cell['id']] = dict(index=cell['index'], indices=cell['indices'], arm=cell['arm'],
                                      eligible=cell['common_accepted'], score=cell['total_score'],
                                      coordinates=report_coordinates(plan, pdb, reference),
                                      score_record_sha256=sha((score_root/name/cell['id']/'score.json').read_bytes()))
            counts['scored_poses'] += 1
        if len(set(signatures)) != 1:
            raise ValueError('fixed context differs across candidates')
        diagnostic = diagnose(row['cells'])
        coupling_counts[diagnostic['coupling']['status']] += 1
        by_index = {r['index']: r for r in reports.values() if r['arm'] == 'own'}
        reference_indices = diagnostic['coupling'].get('reference_indices')
        diagnostic['reference_joint_eligible'] = next((r['eligible'] for r in by_index.values() if r['indices'] == reference_indices), None)
        selected = choice['selected_index']
        if selected is not None:
            counts['selected_cases'] += 1
            for suffix in ('.pdb', '.dat'):
                if (score_root/name/('selected-candidate'+suffix)).read_bytes() != (root/f'joint-{selected}{suffix}').read_bytes():
                    raise ValueError('selected private artifact differs')
        elif any((score_root/name/('selected-candidate'+s)).exists() for s in ('.pdb', '.dat')):
            raise ValueError('unselected case has selected output')
        coarse = diagnostic['coarse_selection']['selected_index']
        if coarse is not None and selected is not None:
            a = by_index[coarse]['coordinates']['pooled_exact_name_rmsd']
            b = by_index[selected]['coordinates']['pooled_exact_name_rmsd']
            coarse_errors.append(a)
            dense_errors.append(b)
            if coarse != selected:
                changes.append(dict(case=name, coarse_index=coarse, dense_index=selected, coarse_rmsd=a, dense_rmsd=b))
        factorized = diagnostic.get('feasible_factorized_index')
        if factorized is not None and factorized != selected:
            factorized_changes.append(dict(case=name, factorized_index=factorized, dense_index=selected,
                                           full_score_gain=diagnostic['dense_vs_factorized_score_gain'],
                                           factorized_rmsd=by_index[factorized]['coordinates']['pooled_exact_name_rmsd'],
                                           dense_rmsd=by_index[selected]['coordinates']['pooled_exact_name_rmsd']))
        independent = [t.selected_index for t in result.targets]
        independent_index = next((i for i, a in enumerate(result.attempts) if list(a.indices) == independent), None)
        independent_eligible = by_index[independent_index]['eligible'] if independent_index is not None else None
        counts['independent_clearance_joint_failures'] += int(independent_eligible is False)
        counts['conditional_independent_joint_failures'] += int(diagnostic.get('independent_conditional_joint_eligible') is False)
        cases[name] = dict(status='evaluated', source_feasible_counts=[sum(a.audit.accepted for a in t.attempts) for t in result.targets],
                           source_attempt_counts=[len(t.attempts) for t in result.targets], diagnostic=diagnostic,
                           independent_clearance_indices=independent, independent_clearance_eligible=independent_eligible,
                           cells=reports)
    aggregate = dict(counts, declared_cases=len(record['cases']), target_view_imports=2*counts['scored_poses'],
                     rounded_rejections=counts['rounded_rejections'],
                     coupling_statuses=dict(coupling_counts), coarse_to_dense_changes=changes, factorized_to_joint_changes=factorized_changes,
                     matched_coarse_dense_cases=len(coarse_errors),
                     median_coarse_rmsd=statistics.median(coarse_errors) if coarse_errors else None,
                     median_dense_rmsd=statistics.median(dense_errors) if dense_errors else None)
    excerpt = {k: v for k, v in enrollment.items() if k != 'score_policy'}
    write_json(output, dict(schema='dvbfixer.private-nine-chi-audit.v1', aggregate=aggregate, cases=cases,
                            plan_path=str(plan_root/'plan.json'), plan_sha256=sha(raw_plan),
                            raw_score_path=str(score_root/'evidence.json'), raw_score_sha256=sha(raw_scores),
                            coarse_score_path=str(coarse_score_root/'evidence.json'), coarse_score_sha256=sha(coarse_raw),
                            enrollment_sha256=record['enrollment_sha256'], enrollment_protocol_excerpt=excerpt,
                            environment=environment, environment_sha256=sha(environment_raw),
                            audit_script_sha256=sha(Path(__file__).read_bytes()), public_capability=False,
                            scientific_qualification=False, scope='conditional small-grid controls in already studied correlated cases, not public completion or independent-family qualification'))
    print(json.dumps(aggregate, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--coarse-scores', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.scores, args.coarse_scores, args.output)
