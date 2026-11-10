"""Retrospective exhaustive reference sensitivity on saved, already audited scores.

Stdlib only: no coordinate generation, geometry changes, native-based selection,
Rosetta import or energy calculation. Keep every reference, including rejected
ones, and separate comparisons using entirely feasible score terms.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path

REPORTING_TOLERANCE = 1e-9  # Reporting only; never changes a selector or geometry.


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_exclusive(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        stream.write((json.dumps(record, sort_keys=True, indent=2, allow_nan=False)+'\n').encode())


def analyze_matrix(rows: list[dict]) -> dict:
    """Consume only index/indices/score/eligibility, never coordinates or RMSD."""
    matrix, indices = {}, set()
    for row in rows:
        key = tuple(row['indices'])
        if (len(key) != 2 or any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < 9 for i in key)
                or isinstance(row['index'], bool) or not isinstance(row['index'], int) or not 0 <= row['index'] < 81
                or not isinstance(row['eligible'], bool) or isinstance(row['score'], bool)
                or not isinstance(row['score'], (int, float)) or not math.isfinite(row['score'])):
            raise ValueError('invalid bounded matrix cell')
        if key in matrix or row['index'] in indices:
            raise ValueError('duplicate matrix indices')
        # Deliberately strip all extraneous fields before any selection.
        matrix[key] = dict(index=row['index'], indices=list(key), score=row['score'], eligible=row['eligible'])
        indices.add(row['index'])
    if not matrix:
        return dict(status='empty_matrix', dimensions=[0, 0], references=[], feasible_rectangles=[],
                    rectangle_count=0, feasible_rectangle_count=0, max_feasible_absolute_delta=None)
    left = sorted({i for i, _ in matrix})
    right = sorted({j for _, j in matrix})
    if len(matrix) != len(left)*len(right):
        raise ValueError('incomplete source-feasible Cartesian score matrix')
    feasible = [r for r in matrix.values() if r['eligible']]
    best = min(feasible, key=lambda r: (r['score'], r['index'])) if feasible else None
    rectangles, rejected_rectangles = [], 0
    for i, k in itertools.combinations(left, 2):
        for j, l in itertools.combinations(right, 2):
            corners = [matrix[i, j], matrix[i, l], matrix[k, j], matrix[k, l]]
            if not all(c['eligible'] for c in corners):
                rejected_rectangles += 1
                continue
            delta = math.fsum([corners[0]['score'], corners[3]['score'], -corners[1]['score'], -corners[2]['score']])
            rectangles.append(dict(rows=[i, k], columns=[j, l], cell_indices=[c['index'] for c in corners],
                                   delta=delta, above_reporting_tolerance=abs(delta) > REPORTING_TOLERANCE))
    references = []
    for a, b in itertools.product(left, right):
        reference = matrix[a, b]

        def surrogate(row, a=a, b=b, reference=reference):
            i, j = row['indices']
            return math.fsum([matrix[i, b]['score'], matrix[a, j]['score'], -reference['score']])

        # All reference-axis terms, not merely the anchor itself, must be valid
        # before calling this a feasible-term global comparison.
        used = {(a, b)} | {(r['indices'][0], b) for r in feasible} | {(a, r['indices'][1]) for r in feasible}
        rejected_terms = sorted(matrix[key]['index'] for key in used if not matrix[key]['eligible'])
        supported = [r for r in feasible if reference['eligible']
                     and matrix[r['indices'][0], b]['eligible'] and matrix[a, r['indices'][1]]['eligible']]
        factorized = min(feasible, key=lambda r: (surrogate(r), r['index'])) if feasible else None
        local_full = min(supported, key=lambda r: (r['score'], r['index'])) if supported else None
        local_factorized = min(supported, key=lambda r: (surrogate(r), r['index'])) if supported else None
        gap = factorized['score']-best['score'] if best else None
        local_gap = local_factorized['score']-local_full['score'] if local_full else None
        classification = ('rejected_reference' if not reference['eligible'] else
                          'feasible_reference_with_rejected_terms' if rejected_terms else 'all_global_terms_feasible')
        references.append(dict(reference_indices=[a, b], reference_index=reference['index'],
                               reference_eligible=reference['eligible'], classification=classification,
                               rejected_term_indices=rejected_terms,
                               global_factorized_index=factorized['index'] if factorized else None,
                               global_full_index=best['index'] if best else None,
                               global_selection_differs=factorized['index'] != best['index'] if best else None,
                               global_full_score_gap=gap,
                               global_score_gap_above_tolerance=gap > REPORTING_TOLERANCE if best else None,
                               supported_pool_indices=sorted(r['index'] for r in supported),
                               supported_full_index=local_full['index'] if local_full else None,
                               supported_factorized_index=local_factorized['index'] if local_factorized else None,
                               supported_contains_global_winner=any(r['index'] == best['index'] for r in supported) if best else None,
                               supported_selection_differs=local_factorized['index'] != local_full['index'] if local_full else None,
                               supported_full_score_gap=local_gap,
                               supported_score_gap_above_tolerance=local_gap > REPORTING_TOLERANCE if local_full else None))
    by_class = {}
    for classification in ('rejected_reference', 'feasible_reference_with_rejected_terms', 'all_global_terms_feasible'):
        group = [r for r in references if r['classification'] == classification]
        by_class[classification] = dict(references=len(group),
                                        differing_global_choices=sum(r['global_selection_differs'] is True for r in group),
                                        positive_global_score_gaps=sum(r['global_score_gap_above_tolerance'] is True for r in group))
    max_delta = max((abs(r['delta']) for r in rectangles), default=None)
    return dict(status='two_axis_matrix' if min(len(left), len(right)) > 1 else 'uninformative_single_option_axis',
                dimensions=[len(left), len(right)], source_feasible_cells=len(matrix), joint_feasible_cells=len(feasible),
                global_full_index=best['index'] if best else None, reference_classes=by_class,
                references=references, feasible_rectangles=rectangles,
                rectangle_count=len(rectangles)+rejected_rectangles, feasible_rectangle_count=len(rectangles),
                rejected_rectangle_count=rejected_rectangles, max_feasible_absolute_delta=max_delta,
                feasible_rectangles_above_tolerance=sum(r['above_reporting_tolerance'] for r in rectangles),
                supported_comparisons=sum(r['supported_full_index'] is not None for r in references),
                supported_global_winner_retained=sum(r['supported_contains_global_winner'] is True for r in references),
                supported_choice_disagreements=sum(r['supported_selection_differs'] is True for r in references),
                supported_positive_score_gaps=sum(r['supported_score_gap_above_tolerance'] is True for r in references))


def plan(audit_path: Path, output: Path) -> None:
    raw = audit_path.read_bytes()
    source = json.loads(raw)
    if source['schema'] != 'dvbfixer.private-nine-chi-audit.v1':
        raise ValueError('expected original dense audit')
    dependencies = {str(Path(__file__).resolve().relative_to(Path.cwd())): sha(Path(__file__).read_bytes()),
                    'scripts/audit_rosetta_dense.py': source['audit_script_sha256'],
                    **source['enrollment_protocol_excerpt']['source_hashes']}
    for filename, digest in dependencies.items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('audited implementation changed')
    inputs = {str(audit_path): sha(raw), source['plan_path']: source['plan_sha256'],
              source['raw_score_path']: source['raw_score_sha256'], source['coarse_score_path']: source['coarse_score_sha256']}
    write_exclusive(output, dict(schema='dvbfixer.reference-sensitivity-plan.v1', input_audit=str(audit_path),
                                 inputs=inputs, source_hashes=dependencies,
                                 case_ids=list(source['cases']), retrospective=True, new_sampling=False,
                                 reporting_tolerance=REPORTING_TOLERANCE,
                                 protocol='all Cartesian anchors, no favorable reference selection; same global feasible set per anchor; separately compare on supported feasible-term subset; enumerate every distinct fully feasible 2x2 rectangle',
                                 interpretation='descriptive analysis of already seen matrices, not new independent prospective qualification',
                                 truth_used_for_selection=False, public_capability=False))


def analyze(plan_path: Path, output: Path) -> None:
    plan_raw = plan_path.read_bytes()
    declaration = json.loads(plan_raw)
    for filename, digest in {**declaration['source_hashes'], **declaration['inputs']}.items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('frozen analysis input or implementation changed')
    source = json.loads(Path(declaration['input_audit']).read_bytes())
    if list(source['cases']) != declaration['case_ids']:
        raise ValueError('saved case inventory changed')
    cases, counts, statuses, classes = {}, Counter(), Counter(), {}
    for name, case in source['cases'].items():
        if case['status'] != 'evaluated':
            raise ValueError('this descriptive plan requires each original evaluated matrix, not omission')
        rows = [cell for cell in case['cells'].values() if cell['arm'] == 'own']
        result = analyze_matrix(rows)
        if result['dimensions'] != case['source_feasible_counts']:
            raise ValueError('source-feasible axis inventory differs')
        if result['global_full_index'] != case['diagnostic']['dense_selection']['selected_index']:
            raise ValueError('analysis changed frozen full-score selection')
        cases[name] = result
        statuses[result['status']] += 1
        for field in ('source_feasible_cells', 'joint_feasible_cells', 'rectangle_count', 'feasible_rectangle_count',
                      'rejected_rectangle_count', 'feasible_rectangles_above_tolerance', 'supported_comparisons',
                      'supported_global_winner_retained', 'supported_choice_disagreements', 'supported_positive_score_gaps'):
            counts[field] += result[field]
        counts['cases_with_feasible_rectangles'] += int(result['feasible_rectangle_count'] > 0)
        counts['cases_with_nonadditive_feasible_rectangles'] += int(result['feasible_rectangles_above_tolerance'] > 0)
        counts['cases_with_supported_positive_gaps'] += int(result['supported_positive_score_gaps'] > 0)
        for label, summary in result['reference_classes'].items():
            classes.setdefault(label, Counter()).update(summary)
    aggregate = dict(counts, cases=len(cases), matrix_statuses=dict(statuses),
                     reference_classes={k: dict(v) for k, v in classes.items()})
    write_exclusive(output, dict(schema='dvbfixer.reference-sensitivity-evidence.v1', aggregate=aggregate, cases=cases,
                                 analysis_plan=declaration, analysis_plan_sha256=sha(plan_raw),
                                 script_sha256=sha(Path(__file__).read_bytes()), retrospective=True,
                                 new_sampling=False, new_scoring=False, selected_artifacts_changed=False,
                                 public_capability=False, scientific_qualification=False))
    print(json.dumps(aggregate, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--audit', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('analyze')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.audit, args.output)
    else:
        analyze(args.plan, args.output)
