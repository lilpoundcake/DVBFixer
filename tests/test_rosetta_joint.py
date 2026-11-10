"""Exact two-target mechanics and fail-closed multi-view scoring, without PyRosetta."""
from __future__ import annotations

import copy
import json
import sys
from dataclasses import asdict

import pytest
from scripts import research_rosetta_joint as module


def snapshot():
    return dict(heavy_atoms=[], hydrogens=[], virtual_atoms=[], residues=[], external=[], allowed_extra=[])


def scored(value=1.):
    return dict(score_status='scored', total_score=value, import_audit={'accepted': True})


def test_both_target_views_required_and_both_import_audits_must_pass():
    with pytest.raises(ValueError, match='exactly two'):
        module.merge_views([scored()], [snapshot()])
    bad = scored()
    bad['import_audit']['accepted'] = False
    result = module.merge_views([scored(), bad], [snapshot(), snapshot()])
    assert result['score_status'] == 'target_view_import_blocker' and result['total_score'] is None


@pytest.mark.parametrize('field', ['heavy_atoms', 'hydrogens', 'virtual_atoms', 'residues', 'external', 'allowed_extra'])
def test_target_views_cannot_be_different_poses(field):
    first, second = snapshot(), snapshot()
    second[field] = ['different']
    assert module.merge_views([scored(), scored()], [first, second])['score_status'] == 'inconsistent_target_view_pose'


@pytest.mark.parametrize('value', [True, float('nan'), float('inf'), 1.000001])
def test_scores_cannot_be_averaged_or_nonfinite(value):
    assert module.merge_views([scored(), scored(value)], [snapshot(), snapshot()])['score_status'] == 'inconsistent_target_view_score'


def test_different_target_graph_views_can_describe_one_identical_molecule():
    first, second = snapshot(), snapshot()
    first['target_internal'], second['target_internal'] = ['target-one-graph'], ['target-two-graph']
    assert module.merge_views([scored(), scored()], [first, second]) == dict(score_status='scored', total_score=1.)


def test_target_hydrogens_may_vary_but_fixed_context_and_virtuals_may_not():
    view = [{'target_key': ['A', 1, ' ', 'VAL']}, {'target_key': ['A', 2, ' ', 'THR']}]
    first = snapshot()
    first['hydrogens'] = [dict(chain='A', resid=i, icode=' ', atom='H', xyz=[float(i), 0., 0.]) for i in (1, 2, 3)]
    second = copy.deepcopy(first)
    second['hydrogens'][0]['xyz'][0] += 1
    assert module.context_signature(first, view) == module.context_signature(second, view)
    second['hydrogens'][2]['xyz'][0] += 1
    assert module.context_signature(first, view) != module.context_signature(second, view)
    third = copy.deepcopy(first)
    third['virtual_atoms'] = ['changed']
    assert module.context_signature(first, view) != module.context_signature(third, view)


def matrix_row(i, j, value):
    return dict(arm='own', indices=[i, j], score_status='scored', total_score=value)


def test_full_matrix_nonadditivity_uses_first_indices_not_best_score_or_truth():
    rows = [matrix_row(i, j, i+2*j+(5 if i == j == 1 else 0)) for i in (0, 1) for j in (0, 1)]
    result = module.coupling_diagnostic(rows)
    assert result['reference_indices'] == [0, 0] and result['max_absolute_delta'] == 5.
    assert result['status'] == 'informative_two_axis_matrix'
    for r in rows:
        r['native_rmsd'] = 0 if r['indices'] == [1, 1] else 100
    assert module.coupling_diagnostic(rows) == result


def test_incomplete_or_degenerate_matrices_are_not_claimed_as_absent_interactions():
    rows = [matrix_row(0, j, float(j)) for j in (0, 1)]
    assert module.coupling_diagnostic(rows)['status'] == 'uninformative_single_option_axis'
    assert module.coupling_diagnostic([])['status'] == 'empty_source_feasible_matrix'
    assert module.coupling_diagnostic(rows+[matrix_row(1, 1, 5.)])['status'] == 'incomplete_score_matrix'
    with pytest.raises(ValueError, match='duplicate'):
        module.coupling_diagnostic(rows+rows)


def test_enrollment_keeps_all_twelve_cases_and_never_runs_admission_or_reconstruction(tmp_path, monkeypatch):
    from dvbfixer.protein_completion import mixed_cohort

    def never(*args, **kwargs):
        raise AssertionError('enrollment must not preview admission/reconstruction')

    monkeypatch.setattr(mixed_cohort, 'plan_mixed_cohort', never)
    monkeypatch.setattr(mixed_cohort, 'complete_mixed_cohort', never)
    imported = 'pyrosetta' in sys.modules
    module.enroll(tmp_path/'enrollment')
    result = json.loads((tmp_path/'enrollment/enrollment.json').read_bytes())
    assert len(result['cases']) == 12 and result['max_combinations'] == 9
    assert {c['mask_class'] for c in result['cases']} == {'distal-distal', 'whole-whole', 'anchored-distal'}
    assert sum(c['system'] == '1CRN' for c in result['cases']) == 6
    assert any(c['targets'][0]['residue']['sequence_number'] == '8' and c['mask_class'] == 'anchored-distal' for c in result['cases'])
    assert ('pyrosetta' in sys.modules) == imported


def test_enrolled_implementation_and_input_hashes_checked_before_use(tmp_path):
    path = tmp_path/'source.pdb'
    path.write_bytes(b'original')
    module.verify(dict(source_hashes={}, artifacts={'source.pdb': module.sha(path.read_bytes())}), tmp_path)
    with pytest.raises(ValueError, match='artifact changed'):
        module.verify(dict(source_hashes={}, artifacts={'source.pdb': '0'*64}), tmp_path)


@pytest.fixture
def saved_joint(tmp_path):
    from dvbfixer.domain.structure_identity import ComponentInstanceRef
    from dvbfixer.protein_completion.mixed_cohort import complete_mixed_cohort, plan_mixed_cohort

    module.enroll(tmp_path/'enrolled')
    enrollment = json.loads((tmp_path/'enrolled/enrollment.json').read_bytes())
    case = enrollment['cases'][0]
    source = (tmp_path/'enrolled'/case['source']).read_bytes()
    plan = plan_mixed_cohort(source, tuple(ComponentInstanceRef(**t['residue']) for t in case['targets']))
    result = complete_mixed_cohort(source, plan)
    return source, case, json.loads(json.dumps(asdict(result)))


def test_independent_joint_replay_does_not_generate_coordinates(saved_joint, monkeypatch):
    from scripts.audit_rosetta_joint import replay

    from dvbfixer.protein_completion import mixed_cohort

    def never(*args, **kwargs):
        raise AssertionError('readback must not call a coordinate generator')

    monkeypatch.setattr(mixed_cohort, 'complete_sidechain', never)
    monkeypatch.setattr(mixed_cohort, 'complete_anchored_branch', never)
    _, result = replay(*saved_joint)
    assert len(result.attempts) == len(saved_joint[2]['attempts'])


def test_missing_combination_and_changed_candidate_are_not_accepted_as_replay(saved_joint):
    from scripts.audit_rosetta_joint import replay

    source, case, saved = saved_joint
    omitted = copy.deepcopy(saved)
    omitted['attempts'].pop()
    with pytest.raises(ValueError, match='canonical joint enumeration'):
        replay(source, case, omitted)
    changed = copy.deepcopy(saved)
    changed['targets'][0]['attempts'][0]['candidate']['additions'][0]['xyz'][0] += 3.
    with pytest.raises(ValueError, match='saved geometry differs'):
        replay(source, case, changed)


def test_no_unverified_rotamer_authority_can_enter_joint_replay():
    from scripts.audit_rosetta_joint import restore_candidate

    with pytest.raises(ValueError, match='empirical authority'):
        restore_candidate({'rotamer_evidence': {}})
