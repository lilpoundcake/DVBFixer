"""Pure exhaustive saved-score diagnostics; no scientific imports or native truth."""
from __future__ import annotations

import copy
import itertools
import math
import sys

import pytest
from scripts import analyze_rosetta_references as module


def matrix(scores, rejected=()):
    return [dict(index=i*len(scores[0])+j, indices=[i, j], score=value, eligible=(i, j) not in rejected)
            for i, row in enumerate(scores) for j, value in enumerate(row)]


def test_additive_matrix_has_zero_feasible_rectangle_deltas_for_every_anchor():
    result = module.analyze_matrix(matrix([[3., 5., 8.], [5., 7., 10.]]))
    assert len(result['references']) == 6 and result['feasible_rectangle_count'] == 3
    assert result['max_feasible_absolute_delta'] == 0.
    assert all(not r['global_selection_differs'] and not r['supported_selection_differs'] for r in result['references'])
    assert result['reference_classes']['all_global_terms_feasible']['references'] == 6


def test_genuine_coupling_remains_detectable_using_only_feasible_corners():
    result = module.analyze_matrix(matrix([[0., -2.], [-1., 10.]]))
    assert result['feasible_rectangle_count'] == 1 and result['max_feasible_absolute_delta'] == 13.
    assert result['supported_positive_score_gaps'] > 0
    assert all(r['classification'] == 'all_global_terms_feasible' for r in result['references'])
    # Every reference is retained, not just the one with the largest difference.
    assert {tuple(r['reference_indices']) for r in result['references']} == set(itertools.product(range(2), repeat=2))


def test_feasible_reference_does_not_imply_feasible_surrogate_terms_or_rectangle():
    result = module.analyze_matrix(matrix([[0., -2.], [-1., 100.]], rejected={(1, 1)}))
    assert result['feasible_rectangle_count'] == 0 and result['max_feasible_absolute_delta'] is None
    assert result['rejected_rectangle_count'] == 1
    assert result['reference_classes']['rejected_reference']['references'] == 1
    assert result['reference_classes']['feasible_reference_with_rejected_terms']['references'] == 2
    assert result['reference_classes']['all_global_terms_feasible']['references'] == 1
    assert result['supported_positive_score_gaps'] == 0
    for reference in result['references']:
        if reference['reference_indices'] == [1, 1]:
            assert reference['supported_pool_indices'] == []
            assert reference['supported_full_index'] is None


def test_supported_pool_comparison_never_pretends_excluded_global_winner_is_retained():
    result = module.analyze_matrix(matrix([[0., -10.], [-1., 5.]], rejected={(1, 1)}))
    reference = next(r for r in result['references'] if r['reference_indices'] == [1, 0])
    assert not reference['supported_contains_global_winner']
    assert reference['supported_full_index'] != result['global_full_index']


def test_single_axis_empty_and_no_joint_feasible_cases_are_distinct():
    single = module.analyze_matrix(matrix([[4., 2., 3.]]))
    assert single['status'] == 'uninformative_single_option_axis'
    assert single['feasible_rectangle_count'] == 0 and single['max_feasible_absolute_delta'] is None
    assert module.analyze_matrix([])['status'] == 'empty_matrix'
    none = module.analyze_matrix(matrix([[0., 1.], [2., 3.]], rejected=set(itertools.product(range(2), repeat=2))))
    assert none['global_full_index'] is None and none['supported_comparisons'] == 0


@pytest.mark.parametrize('mutation', ['duplicate', 'missing', 'nan', 'bool_score', 'bool_index', 'wrong_eligibility', 'off_budget'])
def test_incomplete_duplicate_and_malformed_matrices_fail_closed(mutation):
    rows = matrix([[0., 1.], [2., 3.]])
    if mutation == 'duplicate':
        rows.append(copy.deepcopy(rows[0]))
    elif mutation == 'missing':
        rows.pop()
    elif mutation == 'nan':
        rows[0]['score'] = float('nan')
    elif mutation == 'bool_score':
        rows[0]['score'] = True
    elif mutation == 'bool_index':
        rows[0]['indices'][0] = False
    elif mutation == 'wrong_eligibility':
        rows[0]['eligible'] = 1
    else:
        rows[0]['indices'][0] = 9
    with pytest.raises(ValueError):
        module.analyze_matrix(rows)


def test_reporting_tolerance_does_not_rewrite_selector_or_hide_numerical_ties():
    result = module.analyze_matrix(matrix([[0., 0.], [0., 5e-10]]))
    assert result['feasible_rectangles'][0]['delta'] == 5e-10
    assert not result['feasible_rectangles'][0]['above_reporting_tolerance']
    assert result['global_full_index'] == 0


def test_truth_and_coordinates_are_ignored_and_inputs_are_not_mutated():
    rows = matrix([[0., -2.], [-1., 10.]])
    expected = module.analyze_matrix(rows)
    for row in rows:
        row['native_rmsd'] = -999 if row['index'] == 3 else 999
        row['coordinates'] = {'bogus': True}
    before = copy.deepcopy(rows)
    assert module.analyze_matrix(rows) == expected and rows == before


def test_feasible_rectangle_deltas_invariant_to_additive_row_column_offsets():
    rows = matrix([[0., -2., 8.], [-1., 10., 4.], [5., 2., -5.]])
    first = module.analyze_matrix(rows)
    shifted = copy.deepcopy(rows)
    for row in shifted:
        i, j = row['indices']
        row['score'] += 100*i-20*j+7
    second = module.analyze_matrix(shifted)
    assert [r['delta'] for r in first['feasible_rectangles']] == [r['delta'] for r in second['feasible_rectangles']]


def test_rectangle_inventory_and_all_anchor_choices_match_independent_brute_force():
    rows = matrix([[0., -2., 8.], [-1., 10., 4.], [5., 2., -5.]], rejected={(0, 2), (1, 1)})
    result = module.analyze_matrix(rows)
    lookup = {tuple(r['indices']): r for r in rows}
    expected = []
    for corners in itertools.combinations(rows, 4):
        pairs = {tuple(r['indices']) for r in corners}
        left, right = {i for i, _ in pairs}, {j for _, j in pairs}
        if len(left) == len(right) == 2 and pairs == set(itertools.product(left, right)) and all(r['eligible'] for r in corners):
            i, k = sorted(left)
            j, l = sorted(right)
            expected.append(math.fsum([lookup[i, j]['score'], lookup[k, l]['score'], -lookup[i, l]['score'], -lookup[k, j]['score']]))
    assert sorted(expected) == sorted(r['delta'] for r in result['feasible_rectangles'])
    for reference in result['references']:
        a, b = reference['reference_indices']
        feasible = [r for r in rows if r['eligible']]
        brute = sorted(feasible, key=lambda r: (math.fsum([lookup[r['indices'][0], b]['score'], lookup[a, r['indices'][1]]['score'], -lookup[a, b]['score']]), r['index']))[0]
        assert reference['global_factorized_index'] == brute['index']


def test_module_has_no_science_stack_import_requirement():
    assert module.__doc__ and 'pyrosetta' not in module.__dict__
    assert not any(name.startswith('dvbfixer') for name in module.__dict__)
    assert sys.modules[module.__name__] is module


def test_independent_raw_cell_auditor_matches_full_report():
    from scripts.audit_rosetta_references import verify_matrix

    rows = matrix([[0., -2., 8.], [-1., 10., 4.], [5., 2., -5.]], rejected={(0, 2), (1, 1)})
    report = module.analyze_matrix(rows)
    replayed = verify_matrix(rows, report)
    assert replayed['reference_classes'] == report['reference_classes']
    assert replayed['feasible_rectangle_count'] == report['feasible_rectangle_count']


@pytest.mark.parametrize('tamper', ['rectangle', 'reference', 'supported_pool', 'aggregate', 'rejected_terms'])
def test_independent_auditor_rejects_omitted_or_altered_analysis(tamper):
    from scripts.audit_rosetta_references import verify_matrix

    rows = matrix([[0., -2., 8.], [-1., 10., 4.], [5., 2., -5.]], rejected={(0, 2), (1, 1)})
    report = module.analyze_matrix(rows)
    if tamper == 'rectangle':
        report['feasible_rectangles'].pop()
    elif tamper == 'reference':
        report['references'].pop()
    elif tamper == 'supported_pool':
        report['references'][0]['supported_pool_indices'].append(999)
    elif tamper == 'aggregate':
        report['supported_positive_score_gaps'] += 1
    else:
        report['references'][0]['rejected_term_indices'].append(999)
    with pytest.raises(ValueError):
        verify_matrix(rows, report)
