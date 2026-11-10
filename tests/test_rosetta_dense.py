"""Synthetic pre-execution controls; no dense study case is sampled here."""
import copy
import json
from dataclasses import asdict, replace
from types import SimpleNamespace

import pytest
from scripts import research_rosetta_dense as module
from scripts.rosetta_dense_diagnostics import diagnose

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.mixed_cohort import MixedTarget
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.validation import audit_candidate


def target(name):
    ref = ComponentInstanceRef(1, 'a', '82', 'A')
    points = {'N': (-1.46, 0., 0.), 'CA': (0., 0., 0.), 'C': (.55, 1.43, 0.), 'O': (1.76, 1.66, 0.)}
    branch = module.BRANCHES[name]
    request = CompletionRequest('a'*64, ref, name, tuple(Observation(ExactAtomRef(ref, n), n[0], xyz) for n, xyz in points.items()),
                                ('CB', branch.first, branch.second))
    return MixedTarget(request, 'paired')


@pytest.mark.parametrize('name', ['VAL', 'THR'])
def test_dense_grid_contains_exact_coarse_proposals_but_does_not_change_old_admission(name):
    t = target(name)
    before = t.request.digest
    dense = module.propose(t)
    coarse = complete_sidechain(t.request)
    assert dense.attempts[:3] == coarse.attempts and len(dense.attempts) == 9
    assert any(a.audit.accepted for a in dense.attempts[3:])
    assert all(f.startswith('generated_clash:') for a in dense.attempts for f in a.audit.findings)
    assert t.request.digest == before
    for a in dense.attempts[3:]:
        assert not audit_candidate(t.request, a.candidate).accepted


def test_dense_gate_rejects_undeclared_angle_wrong_metadata_and_moved_generated_atom():
    t = target('THR')
    candidate = module.propose(t).attempts[3].candidate
    assert module.audit_target(t, replace(candidate, torsion_degrees=-79.)).findings == ('dense_grid_or_authority_mismatch',)
    assert not module.audit_target(t, replace(candidate, torsion_degrees=60.)).accepted
    moved = replace(candidate, additions=tuple(replace(a, xyz=(999., 999., 999.)) for a in candidate.additions))
    assert not module.audit_target(t, moved).accepted


def test_joint_audit_never_omits_a_target():
    t = target('VAL')
    with pytest.raises(ValueError, match='exact immutable'):
        module.joint_audits(SimpleNamespace(targets=(t, t)), (module.propose(t).attempts[0].candidate,))


def row(i, j, score, eligible=True):
    return dict(arm='own', index=10*i+j, indices=[i, j], old_source_accepted=eligible,
                common_accepted=eligible, score_status='scored', total_score=score)


def test_conditional_independent_choice_can_fail_joint_gate_without_truth_selection():
    rows = [row(0, 0, 0.), row(0, 1, -2.), row(1, 0, -1.), row(1, 1, 100., False)]
    result = diagnose(rows)
    assert result['independent_conditional_indices'] == [1, 1]
    assert not result['independent_conditional_joint_eligible']
    assert result['dense_selection']['selected_index'] == 1
    assert result['feasible_factorized_index'] == 1
    for r in rows:
        r['native_rmsd'] = 0 if r['index'] == 11 else 99
    assert diagnose(rows) == result


def test_factorized_vs_joint_difference_is_measured_with_identical_geometry_constraints():
    rows = [row(0, 0, 0.), row(0, 1, -2.), row(1, 0, -1.), row(1, 1, 10.)]
    result = diagnose(rows)
    assert result['feasible_factorized_index'] == 11
    assert result['dense_selection']['selected_index'] == 1
    assert result['dense_vs_factorized_score_gain'] == 12.
    assert result['coupling']['max_absolute_delta'] == 13.


def test_factorized_diagnostic_requires_all_source_feasible_scores_even_if_geometry_rejected():
    rows = [row(0, 0, 0.), row(0, 1, -2.), row(1, 0, -1.), row(1, 1, 10., False)]
    rows[-1]['score_status'] = 'execution_blocker'
    assert diagnose(rows)['factorized_status'] == 'incomplete_scores'
    assert diagnose([])['factorized_status'] == 'empty_matrix'


def test_anchored_mode_stays_one_exact_solution_without_rotating_the_observed_partner():
    from dvbfixer.protein_completion.anchored_branch import complete_anchored_branch

    t = target('THR')
    complete = complete_sidechain(t.request).attempts[1].candidate
    fixed = tuple(Observation(a.atom, a.element, a.xyz) for a in complete.additions if a.atom.name != 'OG1')
    request = replace(t.request, observations=t.request.observations+fixed, missing_atoms=('OG1',))
    partner = next(o for o in fixed if o.atom.name == 'CG2')
    anchored = MixedTarget(request, 'anchored', partner)
    assert module.propose(anchored) == complete_anchored_branch(request)
    assert len(module.propose(anchored).attempts) == 1


def test_dense_threonine_reflection_is_rejected_with_all_observations_unchanged():
    from dvbfixer.protein_completion.geometry import add, cross, dot, scaled, sub, unit

    t = target('THR')
    candidate = module.propose(t).attempts[3].candidate
    points = {o.atom.name: o.xyz for o in t.request.observations} | {a.atom.name: a.xyz for a in candidate.additions}
    normal = unit(cross(sub(points['CA'], points['CB']), sub(points['OG1'], points['CB'])))
    reflected = add(points['CG2'], scaled(normal, -2*dot(sub(points['CG2'], points['CB']), normal)))
    changed = replace(candidate, additions=tuple(replace(a, xyz=reflected) if a.atom.name == 'CG2' else a for a in candidate.additions))
    before = t.request.digest
    audit = module.audit_target(t, changed)
    assert not audit.accepted and 'thr_cb_not_R_or_nearplanar' in audit.findings
    assert t.request.digest == before


@pytest.fixture
def saved_synthetic_dense():
    lines, refs = [], []
    for j, name in enumerate(('VAL', 'THR')):
        req = target(name).request
        ref = replace(req.residue, chain_id='A' if j == 0 else 'a')
        refs.append(ref)
        for i, observation in enumerate(req.observations):
            line = _format_generated_atom(serial=20*j+i+1, name=observation.atom.name, element=observation.element,
                                          position=tuple(x+(25*j if k == 0 else 0) for k, x in enumerate(observation.xyz)),
                                          residue_name=name, chain_id=ref.chain_id, sequence_number=ref.sequence_number, insertion_code=ref.insertion_code)
            lines.append('ATOM  '+line[6:])
        lines.append('TER\n')
    source = (''.join(lines)+'END\n').encode()
    plan = module.plan_mixed_cohort(source, tuple(refs))
    targets = tuple(module.propose(t) for t in plan.targets)
    result = module._enumerate_joint(targets, lambda c: module.joint_audits(plan, c), plan.digest,
                                     max_combinations=81, completed_reason='private_nine_chi_joint')
    case = dict(targets=[dict(residue=asdict(r)) for r in refs])
    return source, case, json.loads(json.dumps(asdict(result)))


def test_saved_dense_replay_never_generates_coordinates_and_materialization_owns_only_additions(saved_synthetic_dense, monkeypatch):
    from scripts.audit_rosetta_dense import replay

    def never(*args, **kwargs):
        raise AssertionError('audit must not generate proposals')

    monkeypatch.setattr(module, 'propose', never)
    source, case, saved = saved_synthetic_dense
    plan, result = replay(source, case, saved)
    chosen = result.attempts[result.selected_index]
    candidates = tuple(r.attempts[i].candidate for r, i in zip(result.targets, chosen.indices, strict=True))
    pdb, dat = module.materialize(source, plan, candidates)
    assert pdb and len(dat.added_atoms) == 6
    assert all(key[-1] in {'CB', 'CG1', 'CG2', 'OG1'} for key in dat.added_keys())


def test_saved_dense_replay_rejects_omission_grid_tampering_and_changed_coordinate(saved_synthetic_dense):
    from scripts.audit_rosetta_dense import replay

    source, case, saved = saved_synthetic_dense
    missing = copy.deepcopy(saved)
    missing['attempts'].pop()
    with pytest.raises(ValueError, match='enumeration differs'):
        replay(source, case, missing)
    grid = copy.deepcopy(saved)
    grid['targets'][0]['attempts'][3]['candidate']['torsion_degrees'] = -79.
    with pytest.raises(ValueError, match='grid order/membership'):
        replay(source, case, grid)
    moved = copy.deepcopy(saved)
    moved['targets'][0]['attempts'][3]['candidate']['additions'][0]['xyz'][0] += 10.
    with pytest.raises(ValueError, match='candidate audit differs'):
        replay(source, case, moved)
