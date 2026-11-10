"""Independent saved-candidate exact-joint replay; never regenerates coordinates."""
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

from scripts.branched_stereo_oracle import rdkit_labels
from scripts.research_rosetta_joint import (
    context_signature,
    coupling_diagnostic,
    merge_views,
    sha,
    verify,
)
from scripts.research_rosetta_repack import identity, observed_atoms, write_json
from scripts.rosetta_fixed_contract import choose_scored_own, import_audit
from scripts.summarize_rosetta_fixed_scoring import source_records_preserved
from scripts.unfiltered_branched_comparison import reporting_scores

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.branched_source import _thr_beta_checks
from dvbfixer.protein_completion.cohort import _enumerate_joint
from dvbfixer.protein_completion.contracts import Addition, Attempt, Candidate, CompletionResult
from dvbfixer.protein_completion.mixed_cohort import (
    _audit,
    _joint,
    materialize_mixed_cohort,
    plan_mixed_cohort,
)
from dvbfixer.protein_completion.source import _global_ca, _insert_candidates, _layout
from dvbfixer.protein_completion.validation import sidechain_bonds


def norm(value):
    return json.loads(json.dumps(value, allow_nan=False))


def restore_candidate(raw: dict) -> Candidate:
    if raw['rotamer_evidence'] is not None:
        raise ValueError('joint control does not accept an empirical authority')
    additions = tuple(Addition(ExactAtomRef(ComponentInstanceRef(**a['atom']['component']), a['atom']['name']),
                               a['element'], tuple(a['xyz'])) for a in raw['additions'])
    return Candidate(raw['request_digest'], additions, raw['torsion_degrees'])


def replay(source: bytes, case: dict, saved: dict):
    targets = tuple(ComponentInstanceRef(**t['residue']) for t in case['targets'])
    plan = plan_mixed_cohort(source, targets)
    rebuilt = []
    for target, raw in zip(plan.targets, saved['targets'], strict=True):
        if target.request.digest != raw['request_digest']:
            raise ValueError('source-bound request changed')
        attempts = []
        for row in raw['attempts']:
            candidate = restore_candidate(row['candidate'])
            audit = _audit(target, candidate)
            if norm(asdict(audit)) != row['audit']:
                raise ValueError('individual saved geometry differs')
            attempts.append(Attempt(candidate, audit))
        passing = [i for i, a in enumerate(attempts) if a.audit.accepted]
        selected = max(passing, key=lambda i: min(0., attempts[i].audit.minimum_nonbonded_clearance_angstrom
                                                if attempts[i].audit.minimum_nonbonded_clearance_angstrom is not None else 0.)) if passing else None
        if selected != raw['selected_index']:
            raise ValueError('independent per-target clearance winner changed')
        rebuilt.append(CompletionResult(**{k: v for k, v in raw.items() if k != 'attempts'}, attempts=tuple(attempts)))
    result = _enumerate_joint(tuple(rebuilt), lambda candidates: _joint(plan, candidates), plan.digest,
                              max_combinations=9, completed_reason='mixed_joint_geometry_control')
    if norm(asdict(result)) != saved:
        raise ValueError('exact canonical joint enumeration differs from saved result')
    pools = [[i for i, a in enumerate(t.attempts) if a.audit.accepted] for t in result.targets]
    expected = list(itertools.product(*pools))
    if [tuple(a.indices) for a in result.attempts] != expected:
        raise ValueError('joint search is not the complete Cartesian product')
    return plan, result


def report_coordinates(plan, pdb: bytes, reference: bytes) -> dict:
    atoms = _parse_atom_lines(pdb.decode().splitlines(keepends=True))
    truth = _parse_atom_lines(reference.decode().splitlines(keepends=True))
    _, layout, _ = _layout(pdb)
    _global_ca(layout, allow_missing_cb=None)
    _thr_beta_checks(layout, pending_target=None)
    values = []
    for target in plan.targets:
        request = target.request
        points = {a.name: a.position for a in atoms if _component_ref(a) == request.residue}
        native = {a.name: a.position for a in truth if _component_ref(a) == request.residue}
        labels = rdkit_labels(points, request.residue_name)
        expected = {'CA': 'S', 'CB': 'R'} if request.residue_name == 'THR' else {'CA': 'S'}
        if labels != expected:
            raise ValueError('independent target stereo mismatch')
        additions = tuple(Addition(ExactAtomRef(request.residue, n), n[0], points[n]) for n in request.missing_atoms)
        values.append(dict(residue=asdict(request.residue), residue_name=request.residue_name,
                           cip=labels, scores=reporting_scores(request, additions, native)))
    count = sum(v['scores']['atom_count'] for v in values)
    return dict(targets=values, global_stereo_pass=True, generated_atom_count=count,
                pooled_exact_name_rmsd=math.sqrt(sum(v['scores']['exact_name_rmsd_angstrom']**2*v['scores']['atom_count'] for v in values)/count),
                pooled_generated_only_val_symmetry_rmsd=math.sqrt(sum(v['scores']['generated_only_val_symmetry_rmsd_angstrom']**2*v['scores']['atom_count'] for v in values)/count),
                reporting_only=True)


def audit(plan_root: Path, score_root: Path, output: Path) -> None:
    raw_plan = (plan_root/'plan.json').read_bytes()
    record = json.loads(raw_plan)
    verify(record, plan_root)
    raw_scores = (score_root/'evidence.json').read_bytes()
    scores = json.loads(raw_scores)
    environment_raw = (score_root/'environment.json').read_bytes()
    environment = json.loads(environment_raw)
    if record['enrollment']['score_policy']['expected_rosetta_revision'] not in environment['version']:
        raise ValueError('executed Rosetta build differs from enrollment')
    if scores['plan_sha256'] != sha(raw_plan) or (score_root/'executed-plan.json').read_bytes() != raw_plan:
        raise ValueError('executed joint plan changed')
    if sha((plan_root/'executed-enrollment.json').read_bytes()) != record['enrollment_sha256']:
        raise ValueError('enrollment binding changed')
    cases, changes = {}, []
    cell_statuses, coupling_statuses = Counter(), Counter()
    raw_rejected, rounded_rejected, eligible_count, dat_count, view_count = 0, 0, 0, 0, 0
    baseline_errors, scored_errors = [], []
    for case in record['cases']:
        name = case['id']
        row = scores['results'][name]
        if case['status'] != 'planned':
            cases[name] = dict(status=case['status'], reason=case['reason'])
            continue
        source = (plan_root/name/'source.pdb').read_bytes()
        reference = (plan_root/name/'reference.pdb').read_bytes()
        saved = json.loads((plan_root/name/'reconstruction.json').read_bytes())
        plan, replayed = replay(source, case, saved)
        if norm(asdict(plan)) != json.loads((plan_root/name/'source-plan.json').read_bytes()):
            raise ValueError('mixed plan differs on source replay')
        native_rows = observed_atoms((plan_root/name/'native-control.pdb').read_bytes())
        lookup = {(a['chain'], str(a['resid']), a['icode'].strip(), a['atom']): identity(a) for a in native_rows}

        def key(atom, lookup=lookup):
            c = atom.component
            return lookup[c.chain_id, c.sequence_number, c.insertion_code, atom.name]

        if len(case['views']) != len(plan.targets):
            raise ValueError('a target graph view is absent')
        for target, view in zip(plan.targets, case['views'], strict=True):
            request = target.request
            expected_internal = [(key(ExactAtomRef(request.residue, a)), key(ExactAtomRef(request.residue, b)))
                                 for a, b in sidechain_bonds(request.residue_name)
                                 if (request.residue.chain_id, request.residue.sequence_number, request.residue.insertion_code, a) in lookup
                                 and (request.residue.chain_id, request.residue.sequence_number, request.residue.insertion_code, b) in lookup]
            expected_external = [(key(a), key(b)) for a, b in request.links]
            if norm(view) != norm(dict(target_key=key(ExactAtomRef(request.residue, 'CA'))[:4],
                                       target_internal=expected_internal, expected_external=expected_external)):
                raise ValueError('target graph transport differs from canonical source-derived graph')
        independent = [t.selected_index for t in replayed.targets]
        independent_index = next((i for i, a in enumerate(replayed.attempts) if list(a.indices) == independent), None)
        if independent_index != case['independent_joint_index'] or replayed.selected_index != case['raw_clearance_joint_index']:
            raise ValueError('clearance baseline transport differs')
        own_rows = [c for c in row['cells'] if c['arm'] == 'own']
        if len(own_rows) != len(replayed.attempts):
            raise ValueError('scored search omitted a joint combination')
        selection = choose_scored_own(own_rows)
        if selection != row['selection'] or selection != json.loads((score_root/name/'selection-before-native.json').read_bytes()):
            raise ValueError('score-selected winner differs from complete feasible pool')
        cells, signatures = [], []
        for cell in row['cells']:
            path = score_root/name/cell['id']
            if json.loads((path/'score.json').read_bytes()) != cell:
                raise ValueError('cell record differs from root evidence')
            pdb = (plan_root/cell['pdb']).read_bytes()
            if not source_records_preserved(source, pdb, case['generated']):
                raise ValueError('original source records changed')
            if cell['arm'] == 'native_control':
                actual = {identity(a): a['xyz'] for a in observed_atoms(pdb)}
                truth = {identity(a): a['xyz'] for a in observed_atoms(reference)}
                if any(actual[tuple(k)] != truth[tuple(k)] for k in case['generated']):
                    raise ValueError('native control differs from withheld coordinates')
            if cell['arm'] == 'own':
                attempt = replayed.attempts[cell['index']]
                candidates = tuple(t.attempts[j].candidate for t, j in zip(replayed.targets, attempt.indices, strict=True))
                if _insert_candidates(source, tuple(t.request for t in plan.targets), candidates) != pdb:
                    raise ValueError('scored coordinates differ from saved proposals')
                raw_ok = attempt.accepted
                rounded_ok = False
                if raw_ok:
                    try:
                        materialized = materialize_mixed_cohort(source, plan, candidates)
                    except ValueError:
                        rounded_rejected += 1
                    else:
                        rounded_ok = True
                        sidecar = DatRecord.load((plan_root/cell['pdb']).with_suffix('.dat'))
                        if sidecar.to_dict() != materialized.dat_record().to_dict():
                            raise ValueError('canonical generated-only dat differs')
                        generated = {(a[0], str(a[1]), a[2].strip(), a[4]) for a in case['generated']}
                        if sidecar.added_keys() != generated or len(sidecar.added_atoms) != len(generated):
                            raise ValueError('dat ownership is not exactly the new atoms')
                        dat_count += 1
                else:
                    raw_rejected += 1
                if cell['old_source_accepted'] != raw_ok or cell['common_accepted'] != rounded_ok:
                    raise ValueError('joint raw/rounded eligibility differs')
                eligible_count += int(rounded_ok)
            records, snapshots = [], []
            for i, view in enumerate(case['views']):
                folder = path/f'view-{i}'
                score = json.loads((folder/'score.json').read_bytes())
                payload = (folder/'pose-snapshot.json.gz').read_bytes()
                if sha(payload) != score['snapshot_sha256']:
                    raise ValueError('snapshot digest changed')
                snapshot = json.loads(gzip.decompress(payload))
                replayed_import = import_audit(observed_atoms(pdb), snapshot['heavy_atoms'], view['expected_external'], snapshot['external'],
                                              view['target_internal'], snapshot['target_internal'], snapshot['allowed_extra'])
                if norm(replayed_import) != score['import_audit'] or not replayed_import['accepted']:
                    raise ValueError('independent exact physical/topology import audit failed')
                # Bounded physical extras must genuinely be terminal OXT, not arbitrary names.
                termini = {(r['chain'], r['resid'], r['icode']) for r in snapshot['residues'] if 'CtermProteinFull' in r['name']}
                if any(k[-1] != 'OXT' or tuple(k[:3]) not in termini for k in snapshot['allowed_extra']):
                    raise ValueError('unexpected physical-heavy exception')
                outgoing = {tuple(a if a[-1] == 'C' else b)[:4] for a, b in view['expected_external'] if {a[-1], b[-1]} == {'C', 'N'}}
                if any(tuple(k[:4]) in outgoing for k in snapshot['allowed_extra']):
                    raise ValueError('extra OXT is not a terminal atom')
                if any(not v['virtual_by_rosetta_type'] for v in snapshot['virtual_atoms']):
                    raise ValueError('virtual atom lacks type classification')
                if {identity(a) for a in snapshot['heavy_atoms']} & {identity(a) for a in snapshot['virtual_atoms']}:
                    raise ValueError('physical and virtual inventories overlap')
                total = sum(t['raw']*t['weight'] for t in score['score_terms'].values())
                if {name: term['weight'] for name, term in score['score_terms'].items()} != environment['weights']:
                    raise ValueError('score weights differ from execution environment')
                if not math.isclose(total, score['total_score'], abs_tol=1e-8, rel_tol=1e-10):
                    raise ValueError('weighted score-term sum differs')
                records.append(score)
                snapshots.append(snapshot)
                view_count += 1
            merged = merge_views(records, snapshots)
            if merged['score_status'] != cell['score_status'] or merged['total_score'] != cell['total_score']:
                raise ValueError('two-view pose/score merge differs')
            signature = context_signature(snapshots[0], case['views'])
            if signature != cell['context_signature']:
                raise ValueError('fixed context signature differs')
            signatures.append(signature)
            report = report_coordinates(plan, pdb, reference)
            cell_statuses[cell['score_status']] += 1
            cells.append(dict(id=cell['id'], arm=cell['arm'], index=cell['index'], indices=cell['indices'],
                              score=cell['total_score'], eligible=cell['common_accepted'],
                              coordinate_report=report, score_record_sha256=sha((path/'score.json').read_bytes())))
        if len(set(signatures)) > 1:
            raise ValueError('fixed context differs across the joint candidate matrix')
        coupling = coupling_diagnostic(row['cells'])
        coupling_statuses[coupling['status']] += 1
        own_by_index = {c['index']: c for c in cells if c['arm'] == 'own'}
        baseline = own_by_index.get(case['raw_clearance_joint_index'])
        independent = own_by_index.get(case['independent_joint_index'])
        selected = own_by_index.get(selection['selected_index'])
        if selected is None and any((score_root/name/('selected-candidate'+suffix)).exists() for suffix in ('.pdb', '.dat')):
            raise ValueError('an empty feasible pool has a selected artifact')
        if selected is not None:
            if not selected['eligible']:
                raise ValueError('score selection bypassed joint eligibility')
            original = plan_root/case['id']/f"joint-{selected['index']}"
            for suffix in ('.pdb', '.dat'):
                if original.with_suffix(suffix).read_bytes() != (score_root/name/('selected-candidate'+suffix)).read_bytes():
                    raise ValueError('private selected structure/dat differs from selected feasible cell')
            if baseline is not None and baseline['eligible']:
                a, b = baseline['coordinate_report']['pooled_exact_name_rmsd'], selected['coordinate_report']['pooled_exact_name_rmsd']
                baseline_errors.append(a)
                scored_errors.append(b)
                if baseline['index'] != selected['index']:
                    changes.append(dict(id=name, previous_index=baseline['index'], selected_index=selected['index'],
                                        previous_rmsd=a, selected_rmsd=b))
        cases[name] = dict(status='evaluated', individual_feasible_counts=case['individual_feasible_counts'],
                           observed_ca_distance_angstrom=math.dist(
                               next(o.xyz for o in plan.targets[0].request.observations if o.atom.component == plan.targets[0].request.residue and o.atom.name == 'CA'),
                               next(o.xyz for o in plan.targets[1].request.observations if o.atom.component == plan.targets[1].request.residue and o.atom.name == 'CA')),
                           independent_joint_index=case['independent_joint_index'],
                           independent_joint_eligible=independent['eligible'] if independent else None,
                           raw_clearance_joint_index=case['raw_clearance_joint_index'], selection=selection,
                           coupling=coupling, cells=cells)
    aggregate = dict(declared_cases=len(record['cases']), cell_statuses=dict(cell_statuses),
                     target_view_imports=view_count, own_joint_combinations=sum(c['arm'] == 'own' for r in cases.values() for c in r.get('cells', [])),
                     native_controls=sum(c['arm'] == 'native_control' for r in cases.values() for c in r.get('cells', [])),
                     eligible_joint_combinations=eligible_count, raw_joint_rejections=raw_rejected, rounded_rejections=rounded_rejected,
                     verified_dat_sidecars=dat_count,
                     selection_statuses=dict(Counter(r['selection']['status'] for r in cases.values() if r['status'] == 'evaluated')),
                     coupling_statuses=dict(coupling_statuses), selection_changes=changes,
                     independent_choice_joint_rejections=sum(r.get('independent_joint_eligible') is False for r in cases.values()),
                     same_case_median_baseline_rmsd=statistics.median(baseline_errors) if baseline_errors else None,
                     same_case_median_scored_rmsd=statistics.median(scored_errors) if scored_errors else None,
                     matched_cases=len(baseline_errors))
    enrollment_excerpt = {k: v for k, v in record['enrollment'].items() if k != 'score_policy'}
    enrollment_excerpt['score_policy'] = {k: record['enrollment']['score_policy'][k] for k in (
        'score_function', 'seed', 'init_flags', 'hydrogen_protocol', 'virtual_atom_policy', 'extra_heavy_policy',
        'expected_rosetta_revision', 'selection_policy', 'native_policy', 'topology_scope')}
    write_json(output, dict(schema='dvbfixer.exact-joint-audit.v1', aggregate=aggregate, cases=cases,
                            plan_path=str(plan_root/'plan.json'), plan_sha256=sha(raw_plan),
                            raw_score_path=str(score_root/'evidence.json'), raw_score_sha256=sha(raw_scores),
                            environment=environment, environment_sha256=sha(environment_raw),
                            audit_script_sha256=sha(Path(__file__).read_bytes()),
                            enrollment_protocol_excerpt=enrollment_excerpt,
                            original_enrollment_sha256=record['enrollment_sha256'], scientific_qualification=False,
                            scope='exhaustive small matrices on previously studied systems; not interaction-graph approximation, global conformational optimum or independent-family validation'))
    print(json.dumps(aggregate, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.scores, args.output)


if __name__ == '__main__':
    main()
