"""Nonempirical graph fragment search; references never enter proposal/ranking."""
from __future__ import annotations

import math
from dataclasses import replace
from itertools import combinations

import pytest
from scripts.chemical_reconstruction_all20_controls import (
    HEAVY_POLICY,
    PAIRS,
    sidechain_heavy_mask,
    synthetic_input,
)

from dvbfixer.chemical_reconstruction import conformer_search
from dvbfixer.chemical_reconstruction.conformer_search import (
    audit_search_result,
    search_conformers,
    selected_attempt_index,
)
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Observation,
    Origin,
    Ownership,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan
from dvbfixer.chemical_reconstruction.proposals import CandidateAudit
from dvbfixer.chemical_reconstruction.steric_ranking import compile_steric_context, score_sterics
from dvbfixer.chemical_reconstruction.torsion_geometry import (
    ConformerError,
    FragmentTorsion,
    audit_torsion_offsets,
    compile_fragment_torsions,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def source_case(source='SER'):
    state = dict(PAIRS)[source]
    full, authority, _ = synthetic_input(source, state, {})
    generated = frozenset(o.identity for o in full.observations if o.identity.component == authority.identity
                          and o.identity.name in sidechain_heavy_mask(authority))
    request = replace(full, observations=tuple(o for o in full.observations if o.identity not in generated), ownership=Ownership(generated))
    plan = compile_peptide_request_plan(request, HEAVY_POLICY)
    roots = tuple(sorted(o.identity for o in request.observations if o.identity.component == authority.identity))
    return request, plan, authority, roots


def test_single_axis_grid_trace_and_exact_zero_baseline():
    request, plan, authority, roots = source_case()
    result = search_conformers(request, plan, authority, roots)
    assert result.plan is not None and len(result.plan.topology.rotors) == 1
    assert len(result.attempts) == 3
    assert [a.offsets for a in result.attempts] == [(0.,), (-120.,), (120.,)]
    assert result.attempts[0].raw.additions is result.baseline.raw.additions
    assert result.attempts[0].selectable
    assert all(a.raw is not None and a.raw_trace.passed for a in result.attempts)
    assert result.selected_index is not None and not audit_search_result(request, authority, result)


def test_proper_whole_branch_motion_preserves_all_internal_distances_and_stereo():
    request, plan, authority, roots = source_case('THR')
    result = search_conformers(request, plan, authority, roots)
    rotor, = result.plan.topology.rotors
    assert {a.name for a in rotor.moving} == {'OG1', 'CG2'}
    original = {o.identity: o.position for o in result.baseline.raw.additions}
    for attempt in result.attempts:
        assert attempt.raw is not None and attempt.raw_trace.passed
        points = {o.identity: o.position for o in attempt.raw.additions}
        for a, b in combinations(rotor.moving, 2):
            assert math.dist(points[a], points[b]) == pytest.approx(math.dist(original[a], original[b]), abs=1e-10)
    assert not audit_search_result(request, authority, result)


def test_aromatic_ring_is_rigid_with_two_acyclic_axes_not_ring_bond_sampling():
    request, plan, authority, roots = source_case('PHE')
    result = search_conformers(request, plan, authority, roots)
    assert len(result.plan.topology.rotors) == 2 and len(result.attempts) == 9
    ring = {a.name for a in authority.state.graph.atoms if a.aromatic}
    original = {o.identity.name: o.position for o in result.baseline.raw.additions}
    for attempt in result.attempts:
        assert attempt.raw is not None and attempt.raw_trace.passed
        points = {o.identity.name: o.position for o in attempt.raw.additions}
        for a, b in combinations(sorted(ring), 2):
            assert math.dist(points[a], points[b]) == pytest.approx(math.dist(original[a], original[b]), abs=1e-10)


def test_independent_audit_cannot_call_proposal_or_rotator(monkeypatch):
    from rdkit.Chem import rdMolTransforms
    request, plan, authority, roots = source_case()
    result = search_conformers(request, plan, authority, roots)
    def forbidden(*args, **kwargs):
        pytest.fail('measurement audit called a constructor')
    monkeypatch.setattr(conformer_search, 'propose', forbidden)
    monkeypatch.setattr(conformer_search, 'rotate_fragments', forbidden)
    monkeypatch.setattr(rdMolTransforms, 'SetDihedralDeg', forbidden)
    assert not audit_search_result(request, authority, result)


@pytest.mark.parametrize('source,rotors', [('ALA', 0), ('ARG', 4), ('ASN', 2), ('ASP', 2), ('CYS', 1),
                                        ('GLN', 3), ('GLU', 3), ('HIS', 2), ('ILE', 2), ('LEU', 2),
                                        ('LYS', 4), ('MET', 3), ('PHE', 2), ('PRO', 0), ('SER', 1),
                                        ('THR', 1), ('TRP', 2), ('TYR', 2), ('VAL', 1)])
def test_source_graph_axes_are_not_a_per_residue_chi_constructor(source, rotors):
    request, _, authority, roots = source_case(source)
    topology = compile_fragment_torsions(request, authority, roots)
    assert len(topology.rotors) == rotors
    assert all(set(r.moving) <= request.ownership.generated for r in topology.rotors)
    if source == 'ARG':
        assert not any({a.name for a in r.anchors[1:3]} == {'NE', 'CZ'} for r in topology.rotors)
        assert any(e.reason == 'conjugated-or-unqualified-hybridization' for e in topology.excluded)
    if source == 'PRO':
        assert any(e.reason == 'cyclic-axis' for e in topology.excluded)


def test_multi_chi_cartesian_product_is_complete_and_never_truncated():
    request, plan, authority, roots = source_case('LYS')
    result = search_conformers(request, plan, authority, roots)
    assert len(result.plan.topology.rotors) == 4 and len(result.attempts) == 81
    assert [a.index for a in result.attempts] == list(range(81))
    assert len({a.offsets for a in result.attempts}) == 81
    assert all(a.raw is not None and a.raw_trace.passed for a in result.attempts)
    assert not audit_search_result(request, authority, result)


def test_excess_grid_is_refused_before_enumeration_not_cropped():
    request, plan, authority, roots = source_case()
    result = search_conformers(request, plan, authority, roots)
    target = authority.identity
    # Contract-only fictitious axes, never passed to the rotator or source compiler.
    rotors = tuple(FragmentTorsion(tuple(ExactAtomRef(target, f'{prefix}{i}') for prefix in ('a', 'u', 'v', 'd')),
                                   (ExactAtomRef(target, f'd{i}'),), i) for i in range(6))
    topology = replace(result.plan.topology, rotors=rotors)
    with pytest.raises(ConformerError, match='grid-budget'):
        replace(result.plan, topology=topology)


def test_observed_distal_anchor_blocks_whole_fragment_motion():
    request, plan, authority, roots = source_case()
    base = propose(request, plan)
    distal = next(o for o in base.raw.additions if o.identity.name == 'OG')
    request = replace(request, observations=(*request.observations, distal),
                      ownership=Ownership(request.ownership.generated-{distal.identity}))
    plan = compile_peptide_request_plan(request, HEAVY_POLICY)
    result = search_conformers(request, plan, authority, roots)
    assert result.plan is not None and not result.plan.topology.rotors
    assert any(e.reason == 'observed-downstream-anchor' for e in result.plan.topology.excluded)
    assert len(result.attempts) == 1 and result.attempts[0].raw.additions == result.baseline.raw.additions


def test_frozen_baseline_mismatch_stops_before_axes_or_rotations(monkeypatch):
    request, plan, authority, roots = source_case()
    monkeypatch.setattr(conformer_search, 'compile_fragment_torsions', lambda *a: pytest.fail('compiled axes after mismatch'))
    monkeypatch.setattr(conformer_search, 'rotate_fragments', lambda *a: pytest.fail('rotated after mismatch'))
    result = search_conformers(request, plan, authority, roots, expected_baseline_signature='0'*64)
    assert result.status == 'baseline-refused' and result.plan is None and not result.attempts


def test_ordinary_constructor_refuses_search_strategy_instead_of_ignoring_grid():
    request, plan, authority, roots = source_case()
    result = search_conformers(request, plan, authority, roots)
    legacy = propose(request, result.plan)
    assert legacy.reason == 'strategy-needs-dedicated-constructor' and legacy.raw is None
    assert not result.md_ready and not result.native_conformation_qualified


def test_grid_is_proper_rigid_equivariant_without_moving_source_observations():
    request, plan, authority, roots = source_case('THR')
    result = search_conformers(request, plan, authority, roots)
    def transform(p):
        return (-p[1]+5., p[0]-2., p[2]+7.)
    other_request = replace(request, observations=tuple(replace(o, position=transform(o.position)) for o in request.observations))
    before = other_request.observations
    other = search_conformers(other_request, replace(plan, request_digest=other_request.digest), authority, roots)
    assert other_request.observations is before
    for a, b in zip(result.attempts, other.attempts, strict=True):
        assert a.offsets == b.offsets and b.raw_trace.passed
        for left, right in zip(a.raw.additions, b.raw.additions, strict=True):
            assert left.identity == right.identity and right.position == pytest.approx(transform(left.position), abs=1e-9)


def test_context_obstacle_changes_steric_selection_without_a_reference():
    request, plan, authority, roots = source_case()
    baseline = propose(request, plan)
    positions = {o.identity.name: o.position for o in baseline.raw.additions}
    p, cb = positions['OG'], positions['CB']
    length = math.dist(p, cb)
    obstacle = tuple(p[i]+.9*(p[i]-cb[i])/length for i in range(3))
    identity = ExactAtomRef(ComponentInstanceRef(1, 'L', '900'), 'Q')
    observation = Observation(identity, 'C', obstacle, Unit.ANGSTROM)
    pin = digest_record(('SYNTHETIC contact obstacle, no native holdout', request.digest, observation))
    request = replace(request, source_sha256=pin, observations=(*request.observations, observation),
                      components=(*request.components, Component(identity.component, ChemistryKind.UNKNOWN)),
                      origin=Origin('synthetic-grid-obstacle-v1', pin))
    result = search_conformers(request, replace(plan, request_digest=request.digest), authority, roots)
    assert result.selected_index is not None and result.selected_index != 0
    assert result.attempts[result.selected_index].rounded_sterics.sum_squared_overlap_angstrom2 < result.attempts[0].rounded_sterics.sum_squared_overlap_angstrom2


def test_steric_pair_accounting_uses_canonical_radii_and_bond_distance_exclusions():
    from dvbfixer.diagnose.geometry import VDW_RADII_ANGSTROM
    request, plan, _, _ = source_case()
    base = propose(request, plan)
    context = compile_steric_context(request)
    score = score_sterics(request, plan, base.raw, context)
    positions = {o.identity: o.angstrom for o in request.observations} | {o.identity: o.position for o in base.raw.additions}
    elements = {o.identity: o.element for o in (*request.observations, *base.raw.additions)}
    pairs = {(a, b) for a, b, _ in context.pairs}
    assert len(pairs) == len(context.pairs) == score.pair_count
    assert all({a, b} & request.ownership.generated for a, b in pairs)
    assert not any({a.name, b.name} == {'CB', 'OG'} for a, b in pairs)
    expected = []
    for a, b, radius in context.pairs:
        assert radius == VDW_RADII_ANGSTROM[elements[a]]+VDW_RADII_ANGSTROM[elements[b]]
        expected.append(max(0., radius-math.dist(positions[a], positions[b])))
    assert score.sum_squared_overlap_angstrom2 == math.fsum(x*x for x in expected)


def test_unknown_radii_and_pair_budget_refuse_not_skip_context(monkeypatch):
    from dvbfixer.chemical_reconstruction import steric_ranking
    request, _, _, _ = source_case()
    outside = next(o for o in request.observations if o.identity.component not in request.targets)
    bad = replace(request, observations=tuple(replace(o, element='XE') if o.identity == outside.identity else o for o in request.observations))
    with pytest.raises(ConformerError, match='unknown-context-radii'):
        compile_steric_context(bad)
    monkeypatch.setattr(steric_ranking, 'MAX_PAIR_EVALUATIONS', 1)
    with pytest.raises(ConformerError, match='pair-budget'):
        compile_steric_context(request)


def test_fractional_context_is_not_downweighted_by_ranker():
    request, plan, _, _ = source_case()
    base = propose(request, plan)
    old = score_sterics(request, plan, base.raw, compile_steric_context(request))
    outside = next(o for o in request.observations if o.identity.component not in request.targets)
    changed = replace(request, observations=tuple(replace(o, occupancy=.3) if o.identity == outside.identity else o for o in request.observations))
    other_plan = replace(plan, request_digest=changed.digest)
    other = replace(base.raw, request_digest=changed.digest, plan_digest=other_plan.digest)
    assert score_sterics(changed, other_plan, other, compile_steric_context(changed)) == old


def test_trace_checks_unsampled_bonds_and_rejects_wrong_declared_grid():
    request, plan, authority, roots = source_case('ARG')
    result = search_conformers(request, plan, authority, roots)
    selected = result.attempts[1]
    assert selected.raw_trace.passed
    wrong = audit_torsion_offsets(request, result.plan.base_candidate, selected.raw, result.plan.topology,
                                  (0.,)*len(result.plan.topology.rotors), rounded=False)
    assert not wrong.passed and wrong.checked_paths > len(result.plan.topology.rotors)


@pytest.mark.parametrize('mutation', ['selection', 'offset', 'missing-attempt', 'raw-position', 'rounded-score'])
def test_saved_search_auditor_rejects_tampered_selection_trace_and_scores(mutation):
    request, plan, authority, roots = source_case()
    result = search_conformers(request, plan, authority, roots)
    attempt = result.attempts[0]
    if mutation == 'selection':
        result = replace(result, selected_index=999)
    elif mutation == 'missing-attempt':
        result = replace(result, attempts=result.attempts[:-1])
    elif mutation == 'offset':
        result = replace(result, attempts=(replace(attempt, offsets=(120.,)), *result.attempts[1:]))
    elif mutation == 'raw-position':
        atom = attempt.raw.additions[-1]
        raw = replace(attempt.raw, additions=(*attempt.raw.additions[:-1], replace(atom, position=(atom.position[0]+.2, atom.position[1], atom.position[2]))))
        result = replace(result, attempts=(replace(attempt, raw=raw), *result.attempts[1:]))
    else:
        score = replace(attempt.rounded_sterics, sum_squared_overlap_angstrom2=999.)
        result = replace(result, attempts=(replace(attempt, rounded_sterics=score), *result.attempts[1:]))
    assert audit_search_result(request, authority, result)


def test_tie_break_is_stable_index_and_failed_geometry_is_not_selectable():
    request, plan, authority, roots = source_case()
    result = search_conformers(request, plan, authority, roots)
    first = result.attempts[0]
    second = replace(first, index=1)
    assert selected_attempt_index((second, first)) == 0
    failed = replace(first, raw_audit=CandidateAudit(False, ('forced-unit-negative',), (), ()))
    assert not failed.selectable and selected_attempt_index((failed, second)) == 1


def test_grid_policy_rejects_bool_or_unversioned_offsets():
    request, plan, authority, roots = source_case()
    result = search_conformers(request, plan, authority, roots)
    with pytest.raises(ConformerError):
        replace(result.plan, grid_offsets=(False, -120., 120.))
    with pytest.raises(ConformerError):
        replace(result.plan, grid_offsets=(0., -60., 60.))


def test_twenty_frozen_entries_keep_baseline_signatures_and_negative_native_outcome(tmp_path, monkeypatch):
    import inspect

    from scripts import chemical_reconstruction_conformer_controls as controls
    assert not {'reference', 'heldout', 'rmsd'} & inspect.signature(search_conformers).parameters.keys()
    original_search = controls.search_conformers
    original_reference = controls.heldout_after_proposal
    original_score = controls.score_heavy_holdout
    state = {'finished': False, 'references': 0}
    def wrapped_search(*args, **kwargs):
        state['finished'] = False
        result = original_search(*args, **kwargs)
        state['finished'] = True
        return result
    def wrapped_reference(*args, **kwargs):
        assert state['finished']
        state['references'] += 1
        return original_reference(*args, **kwargs)
    def wrapped_score(*args, **kwargs):
        assert state['finished']
        return original_score(*args, **kwargs)
    monkeypatch.setattr(controls, 'search_conformers', wrapped_search)
    monkeypatch.setattr(controls, 'heldout_after_proposal', wrapped_reference)
    monkeypatch.setattr(controls, 'score_heavy_holdout', wrapped_score)
    result = controls.run_controls(tmp_path/'controls')
    assert result['declared'] == result['executed_entries'] == 20
    assert sum(c['status'] == 'selected' for c in result['controls']) == 18
    assert sum(c['status'] == 'not-applicable-no-sidechain-heavy-atoms' for c in result['controls']) == 1
    assert sum(c['status'] == 'refused' for c in result['controls']) == 1
    assert sum(c.get('attempts', 0) for c in result['controls']) == 302
    assert state['references'] == 18
    import json
    asp = json.loads((tmp_path/'controls/observed-ASP-ASH/evidence.json').read_text())
    scores = asp['postselection_native_metrics']
    assert scores['selected']['rounded']['rmsd_angstrom'] > scores['baseline']['rounded']['rmsd_angstrom']
    with pytest.raises(FileExistsError):
        controls.run_controls(tmp_path/'controls')
