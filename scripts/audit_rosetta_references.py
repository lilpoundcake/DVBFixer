"""Independent stdlib-only replay of every saved-reference comparison.

Use original raw scored cells rather than the analysis's audit-derived matrices.
Enumerate feasible rectangles as sets of four cells, not Cartesian axis pairs.
No imports from the analyzer, geometry stack or scoring runtime.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


def verify_matrix(rows, reported):
    cells = {tuple(r['indices']): r for r in rows}
    if len(cells) != len(rows) or len({r['index'] for r in rows}) != len(rows):
        raise ValueError('duplicate raw cell')
    left_values = {i for i, _ in cells}
    right_values = {j for _, j in cells}
    if not cells or set(cells) != set(itertools.product(left_values, right_values)):
        raise ValueError('complete nonempty raw Cartesian matrix required')
    feasible = [r for r in rows if r['eligible']]
    rectangles = {}
    for group in itertools.combinations(feasible, 4):
        keys = {tuple(r['indices']) for r in group}
        left = {key[0] for key in keys}
        right = {key[1] for key in keys}
        if len(left) != 2 or len(right) != 2 or keys != {(i, j) for i in left for j in right}:
            continue
        i, k = sorted(left)
        j, l = sorted(right)
        delta = math.fsum([cells[i, j]['score'], cells[k, l]['score'], -cells[i, l]['score'], -cells[k, j]['score']])
        rectangles[i, k, j, l] = delta
    listed = {(r['rows'][0], r['rows'][1], r['columns'][0], r['columns'][1]): r for r in reported['feasible_rectangles']}
    if len(listed) != len(reported['feasible_rectangles']) or set(listed) != set(rectangles):
        raise ValueError('feasible rectangle inventory differs')
    for key, delta in rectangles.items():
        item = listed[key]
        i, k, j, l = key
        indices = [cells[p]['index'] for p in ((i, j), (i, l), (k, j), (k, l))]
        if item['cell_indices'] != indices or item['delta'] != delta or item['above_reporting_tolerance'] != (abs(delta) > 1e-9):
            raise ValueError('feasible rectangle values differ')
    anchors = {tuple(r['reference_indices']): r for r in reported['references']}
    if len(anchors) != len(reported['references']) or set(anchors) != set(cells):
        raise ValueError('all-reference inventory differs')
    best = sorted(feasible, key=lambda r: (r['score'], r['index']))[0] if feasible else None
    classes = {key: dict(references=0, differing_global_choices=0, positive_global_score_gaps=0) for key in
               ('rejected_reference', 'feasible_reference_with_rejected_terms', 'all_global_terms_feasible')}
    supported_count = retained = support_disagrees = support_gaps = 0
    for key, anchor in cells.items():
        a, b = key
        reported_anchor = anchors[key]
        candidates = []
        supported = []
        terms = {key}
        for cell in feasible:
            i, j = cell['indices']
            terms.update(((i, b), (a, j)))
            value = math.fsum([cells[i, b]['score'], cells[a, j]['score'], -anchor['score']])
            candidates.append((value, cell['index'], cell))
            if all(cells[p]['eligible'] for p in (key, (i, j), (i, b), (a, j))):
                supported.append((value, cell['index'], cell))
        factorized = min(candidates, key=lambda t: t[:2])[2] if candidates else None
        subset_full = min((t[2] for t in supported), key=lambda r: (r['score'], r['index'])) if supported else None
        subset_factorized = min(supported, key=lambda t: t[:2])[2] if supported else None
        rejected = sorted(cells[p]['index'] for p in terms if not cells[p]['eligible'])
        label = ('rejected_reference' if not anchor['eligible'] else
                 'feasible_reference_with_rejected_terms' if rejected else 'all_global_terms_feasible')
        expected = dict(reference_indices=list(key), reference_index=anchor['index'], reference_eligible=anchor['eligible'],
                        classification=label, rejected_term_indices=rejected,
                        global_factorized_index=factorized['index'] if factorized else None,
                        global_full_index=best['index'] if best else None,
                        global_selection_differs=factorized['index'] != best['index'] if best else None,
                        global_full_score_gap=factorized['score']-best['score'] if best else None,
                        global_score_gap_above_tolerance=factorized['score']-best['score'] > 1e-9 if best else None,
                        supported_pool_indices=sorted(t[1] for t in supported),
                        supported_full_index=subset_full['index'] if subset_full else None,
                        supported_factorized_index=subset_factorized['index'] if subset_factorized else None,
                        supported_contains_global_winner=any(t[1] == best['index'] for t in supported) if best else None,
                        supported_selection_differs=subset_factorized['index'] != subset_full['index'] if subset_full else None,
                        supported_full_score_gap=subset_factorized['score']-subset_full['score'] if subset_full else None,
                        supported_score_gap_above_tolerance=subset_factorized['score']-subset_full['score'] > 1e-9 if subset_full else None)
        if expected != reported_anchor:
            raise ValueError('reference/supported-pool comparison differs')
        classes[label]['references'] += 1
        classes[label]['differing_global_choices'] += int(expected['global_selection_differs'] is True)
        classes[label]['positive_global_score_gaps'] += int(expected['global_score_gap_above_tolerance'] is True)
        supported_count += int(subset_full is not None)
        retained += int(expected['supported_contains_global_winner'] is True)
        support_disagrees += int(expected['supported_selection_differs'] is True)
        support_gaps += int(expected['supported_score_gap_above_tolerance'] is True)
    nleft = len({i for i, _ in cells})
    nright = len({j for _, j in cells})
    total_rectangles = math.comb(nleft, 2)*math.comb(nright, 2)
    expected_summary = dict(
        status='two_axis_matrix' if min(nleft, nright) > 1 else 'uninformative_single_option_axis',
        dimensions=[nleft, nright], source_feasible_cells=len(cells), joint_feasible_cells=len(feasible),
        global_full_index=best['index'] if best else None, reference_classes=classes,
        rectangle_count=total_rectangles, feasible_rectangle_count=len(rectangles),
        rejected_rectangle_count=total_rectangles-len(rectangles),
        max_feasible_absolute_delta=max((abs(v) for v in rectangles.values()), default=None),
        feasible_rectangles_above_tolerance=sum(abs(v) > 1e-9 for v in rectangles.values()),
        supported_comparisons=supported_count, supported_global_winner_retained=retained,
        supported_choice_disagreements=support_disagrees, supported_positive_score_gaps=support_gaps)
    if any(reported[k] != v for k, v in expected_summary.items()):
        raise ValueError('matrix summary differs from raw-cell replay')
    return expected_summary


def audit(plan_path, analysis_path, output):
    plan_raw, analysis_raw = plan_path.read_bytes(), analysis_path.read_bytes()
    declaration, analysis = json.loads(plan_raw), json.loads(analysis_raw)
    if (analysis.get('retrospective') is not True or analysis.get('new_sampling') is not False
            or analysis.get('new_scoring') is not False or analysis.get('selected_artifacts_changed') is not False):
        raise ValueError('analysis provenance flags differ')
    if analysis['analysis_plan_sha256'] != sha(plan_raw) or analysis['analysis_plan'] != declaration:
        raise ValueError('analysis plan binding differs')
    for filename, digest in {**declaration['source_hashes'], **declaration['inputs']}.items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('analysis dependency changed')
    if analysis['script_sha256'] != sha(Path('scripts/analyze_rosetta_references.py').read_bytes()):
        raise ValueError('analyzer implementation changed')
    original = json.loads(Path(declaration['input_audit']).read_bytes())
    raw_scores = json.loads(Path(original['raw_score_path']).read_bytes())
    if set(analysis['cases']) != set(declaration['case_ids']) or set(raw_scores['results']) != set(declaration['case_ids']):
        raise ValueError('analysis omitted or added a saved case')
    summaries, summed, statuses, classes = {}, {}, {}, {}
    for name in declaration['case_ids']:
        rows = [dict(index=r['index'], indices=r['indices'], score=r['total_score'],
                     eligible=r['old_source_accepted'] and r['common_accepted'])
                for r in raw_scores['results'][name]['cells'] if r['arm'] == 'own']
        summary = verify_matrix(rows, analysis['cases'][name])
        summaries[name] = summary
        statuses[summary['status']] = statuses.get(summary['status'], 0)+1
        for field in ('source_feasible_cells', 'joint_feasible_cells', 'rectangle_count', 'feasible_rectangle_count',
                      'rejected_rectangle_count', 'feasible_rectangles_above_tolerance', 'supported_comparisons',
                      'supported_global_winner_retained', 'supported_choice_disagreements', 'supported_positive_score_gaps'):
            summed[field] = summed.get(field, 0)+summary[field]
        for name_field, value in (
            ('cases_with_feasible_rectangles', summary['feasible_rectangle_count'] > 0),
            ('cases_with_nonadditive_feasible_rectangles', summary['feasible_rectangles_above_tolerance'] > 0),
            ('cases_with_supported_positive_gaps', summary['supported_positive_score_gaps'] > 0)):
            summed[name_field] = summed.get(name_field, 0)+int(value)
        for label, group in summary['reference_classes'].items():
            for key, value in group.items():
                classes.setdefault(label, {})[key] = classes.setdefault(label, {}).get(key, 0)+value
    expected = dict(summed, cases=len(summaries), matrix_statuses=statuses, reference_classes=classes)
    if expected != analysis['aggregate']:
        raise ValueError('global aggregate differs')
    result = dict(schema='dvbfixer.reference-sensitivity-audit.v1', aggregate=expected, cases=summaries,
                  analysis_sha256=sha(analysis_raw), plan_sha256=sha(plan_raw),
                  audit_script_sha256=sha(Path(__file__).read_bytes()),
                  independent_rectangle_method='four eligible cells, not analyzer Cartesian-axis enumeration',
                  independent_data_source='original raw saved score cells, not analyzer matrix or coordinates',
                  no_resampling=True, no_new_scoring=True, all_references_retained=True,
                  retrospective=True, scientific_qualification=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('xb') as stream:
        stream.write((json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n').encode())
    print(json.dumps(expected, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.analysis, args.output)
