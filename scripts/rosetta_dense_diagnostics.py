"""Frozen, truth-free dense/coarse and conditional factorized comparisons."""
from __future__ import annotations

import math

from scripts.research_rosetta_joint import coupling_diagnostic
from scripts.rosetta_fixed_contract import choose_scored_own


def diagnose(rows: list[dict]) -> dict:
    own = [r for r in rows if r['arm'] == 'own']
    dense = choose_scored_own(own)
    coarse = choose_scored_own([r for r in own if all(i < 3 for i in r['indices'])])
    coupling = coupling_diagnostic(own)
    result = dict(dense_selection=dense, coarse_selection=coarse, coupling=coupling)
    if not own:
        return dict(result, factorized_status='empty_matrix')
    if any(r['score_status'] != 'scored' or isinstance(r.get('total_score'), bool)
           or not isinstance(r.get('total_score'), (float, int)) or not math.isfinite(r['total_score']) for r in own):
        return dict(result, factorized_status='incomplete_scores')
    matrix = {tuple(r['indices']): r for r in own}
    if len(matrix) != len(own):
        raise ValueError('duplicate joint indices')
    left = sorted({i for i, _ in matrix})
    right = sorted({j for _, j in matrix})
    if len(matrix) != len(left)*len(right):
        return dict(result, factorized_status='incomplete_matrix')
    a, b = left[0], right[0]
    i = min(left, key=lambda i: (matrix[i, b]['total_score'], i))
    j = min(right, key=lambda j: (matrix[a, j]['total_score'], j))
    independent = matrix[i, j]
    feasible = [r for r in own if r['old_source_accepted'] and r['common_accepted']]

    def factorized(row):
        i, j = row['indices']
        return matrix[i, b]['total_score']+matrix[a, j]['total_score']-matrix[a, b]['total_score']

    winner = min(feasible, key=lambda r: (factorized(r), r['index'])) if feasible else None
    dense_row = next((r for r in own if r['index'] == dense['selected_index']), None)
    return dict(result, factorized_status='complete', reference_indices=[a, b],
                independent_conditional_indices=[i, j], independent_conditional_index=independent['index'],
                independent_conditional_joint_eligible=independent['old_source_accepted'] and independent['common_accepted'],
                feasible_factorized_index=winner['index'] if winner else None,
                factorized_selected_full_score=winner['total_score'] if winner else None,
                dense_vs_factorized_score_gain=winner['total_score']-dense_row['total_score'] if winner and dense_row else None,
                interpretation='conditional reference-axis full-score surrogate, not a physical independent-residue energy decomposition')
