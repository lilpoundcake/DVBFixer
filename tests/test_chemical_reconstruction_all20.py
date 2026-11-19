"""All20 bounded source coverage, fixed-frame scoring and reference isolation."""
from __future__ import annotations

from collections import Counter
from dataclasses import replace

import pytest
from scripts import chemical_reconstruction_all20_controls as controls

from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import Ownership, Unit
from dvbfixer.chemical_reconstruction.heldout_scoring import (
    HoldoutScoringError,
    score_heavy_holdout,
)
from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan
from dvbfixer.chemical_reconstruction.proposals import Strategy


def heavy_case(source='THR', state='THR'):
    full, authority, _ = controls.synthetic_input(source, state, {})
    generated = frozenset(o.identity for o in full.observations if o.identity.component == authority.identity
                          and o.identity.name in controls.sidechain_heavy_mask(authority))
    request = replace(full, observations=tuple(o for o in full.observations if o.identity not in generated), ownership=Ownership(generated))
    plan = compile_peptide_request_plan(request, controls.HEAVY_POLICY)
    result = propose(request, plan)
    reference = tuple(o for o in full.observations if o.identity in generated)
    return request, plan, result, reference


@pytest.mark.parametrize('source,state', controls.PAIRS)
def test_twenty_source_state_chemical_proofs_use_unchanged_source_graph(source, state):
    evidence = {}
    request, authority, _ = controls.synthetic_input(source, state, evidence)
    assert authority.inventory.state_name == state and authority.parent.component.component_id == source
    graph = authority.state.graph
    old = {a.name: a for a in authority.parent.state.graph.atoms}
    assert all(a == old[a.name] for a in graph.atoms)
    assert all(b in authority.parent.state.graph.bonds for b in graph.bonds)
    assert request.ownership.generated == frozenset()
    assert all(o.element not in {'H', 'D'} for o in request.observations)
    assert authority.removed_atoms == ({'OXT', 'HXT', 'H'} if source == 'PRO' else {'OXT', 'HXT', 'H2'})
    assert controls.sidechain_heavy_mask(authority).isdisjoint({'N', 'CA', 'C', 'O'})
    if state in {'ASH', 'GLH', 'HIP'}:
        assert source != state  # Explicit caller choice, never source PDB renaming.


def test_exact_unchanged_frame_score_has_no_alignment_or_acceptance_gate():
    request, plan, result, reference = heavy_case()
    zero = score_heavy_holdout(request, plan, result.raw, reference)
    assert zero.atom_count == 3 and zero.rmsd_angstrom == pytest.approx(0., abs=1e-12)
    displaced = replace(result.raw, additions=tuple(replace(o, position=(o.position[0]+3., o.position[1], o.position[2]))
                                                   for o in result.raw.additions))
    score = score_heavy_holdout(request, plan, displaced, reference)
    assert score.rmsd_angstrom == pytest.approx(3.) and score.maximum_angstrom == pytest.approx(3.)
    assert score.scope == 'exact-name-fixed-frame-heavy-only-no-symmetry-fit-or-acceptance-gate'
    assert not hasattr(score, 'passed')
    reverse = score_heavy_holdout(request, plan, displaced, tuple(reversed(reference)))
    assert reverse == score


@pytest.mark.parametrize('mutation', ['missing-reference', 'missing-candidate', 'duplicate', 'wrong-chain-case', 'unit',
                                     'occupancy', 'element', 'candidate-unit', 'candidate-element', 'plan-binding', 'strategy'])
def test_incomplete_or_incompatible_scores_refuse_instead_of_partial_metrics(mutation):
    request, plan, result, reference = heavy_case()
    candidate = result.raw
    first = reference[0]
    if mutation == 'missing-reference':
        reference = reference[1:]
    elif mutation == 'missing-candidate':
        candidate = replace(candidate, additions=candidate.additions[1:])
    elif mutation == 'duplicate':
        reference = (reference[1], *reference[1:])
    elif mutation == 'wrong-chain-case':
        identity = replace(first.identity, component=replace(first.identity.component, chain_id=first.identity.component.chain_id.upper()))
        reference = (replace(first, identity=identity), *reference[1:])
    elif mutation == 'unit':
        reference = (replace(first, unit=Unit.NANOMETER), *reference[1:])
    elif mutation == 'occupancy':
        reference = (replace(first, occupancy=.5), *reference[1:])
    elif mutation == 'element':
        reference = (replace(first, element='H'), *reference[1:])
    elif mutation == 'candidate-unit':
        candidate = replace(candidate, additions=(replace(candidate.additions[0], unit=Unit.NANOMETER), *candidate.additions[1:]))
    elif mutation == 'candidate-element':
        candidate = replace(candidate, additions=(replace(candidate.additions[0], element='H'), *candidate.additions[1:]))
    elif mutation == 'plan-binding':
        candidate = replace(candidate, plan_digest='0'*64)
    else:
        candidate = replace(candidate, strategy=Strategy.TEMPLATE)
    with pytest.raises(HoldoutScoringError):
        score_heavy_holdout(request, plan, candidate, reference)


def test_score_requires_generated_atoms_not_retained_observations():
    request, plan, result, reference = heavy_case()
    request = replace(request, observations=(*request.observations, reference[0]))
    plan = replace(plan, request_digest=request.digest)
    candidate = replace(result.raw, request_digest=request.digest, plan_digest=plan.digest)
    with pytest.raises(HoldoutScoringError, match='overlaps-observation'):
        score_heavy_holdout(request, plan, candidate, reference)


def test_score_is_available_for_complete_geometry_rejections_not_a_pass_filter(tmp_path):
    evidence = {}
    request, authority, _ = controls.observed_input('LYS', 'LYS', 'heavy', tmp_path, evidence)
    plan = compile_peptide_request_plan(request, controls.HEAVY_POLICY)
    result = propose(request, plan)
    reference = controls.heldout_after_proposal('observed', authority, tmp_path)
    candidate = replace(result.raw, additions=tuple(replace(o, position=(o.position[0]+10., o.position[1], o.position[2]))
                                                   for o in result.raw.additions))
    from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate
    assert not audit_candidate(request, plan, candidate).passed
    score = score_heavy_holdout(request, plan, candidate, reference)
    assert score.atom_count == 5 and score.rmsd_angstrom > 5


def test_heldout_extraction_occurs_after_proposal_and_never_reaches_request(tmp_path, monkeypatch):
    evidence = {}
    request, authority, _ = controls.observed_input('LYS', 'LYS', 'heavy', tmp_path, evidence)
    pending = request.ownership.generated
    assert len(pending) == 5 and not pending & {o.identity for o in request.observations}
    state = {'proposal_finished': False, 'reference_extracted': False}
    original_propose = controls.propose
    original_reference = controls.heldout_after_proposal
    def wrapped_propose(r, p):
        assert not state['reference_extracted']
        assert r is request and not any(o.identity in pending for o in r.observations)
        result = original_propose(r, p)
        state['proposal_finished'] = True
        return result
    def wrapped_reference(*args):
        assert state['proposal_finished']
        state['reference_extracted'] = True
        return original_reference(*args)
    monkeypatch.setattr(controls, 'propose', wrapped_propose)
    monkeypatch.setattr(controls, 'heldout_after_proposal', wrapped_reference)
    item = {}
    controls.execute_heavy(request, authority, 'observed', tmp_path, item, evidence)
    assert all(state.values()) and evidence['heavy_scores']['raw']['score'].atom_count == 5


def test_exact_name_score_does_not_hide_symmetric_atom_name_swaps():
    request, plan, result, reference = heavy_case('PHE', 'PHE')
    points = {o.identity.name: o.position for o in result.raw.additions}
    swapped = replace(result.raw, additions=tuple(replace(o, position=points['CD2' if o.identity.name == 'CD1' else 'CD1'])
                                if o.identity.name in {'CD1', 'CD2'} else o for o in result.raw.additions))
    score = score_heavy_holdout(request, plan, swapped, reference)
    assert score.rmsd_angstrom > .1


def test_fixed_state_table_and_distinct_heavy_H_policies():
    assert len(controls.PAIRS) == len({a for a, _ in controls.PAIRS}) == 20
    assert dict(controls.PAIRS)['ASP'] == 'ASH' and dict(controls.PAIRS)['GLU'] == 'GLH'
    assert dict(controls.PAIRS)['HIS'] == 'HIP'
    assert controls.HEAVY_POLICY.bond_tolerance_angstrom == .10
    assert controls.HEAVY_POLICY.angle_tolerance_degrees == 12
    assert controls.GENERATED_POLICY.bond_tolerance_angstrom == .02
    assert controls.GENERATED_POLICY.angle_tolerance_degrees == 5


def test_eighty_frozen_entries_preserve_not_applicable_missing_MET_and_fixed_planes(tmp_path):
    result = controls.run_controls(tmp_path/'controls')
    assert result['declared'] == result['executed_entries'] == 80
    expected = {('synthetic', 'heavy'): {'declared-geometry-passed': 19, 'not-applicable-no-sidechain-heavy-atoms': 1},
                ('synthetic', 'H'): {'declared-geometry-passed': 20},
                ('observed', 'heavy'): {'declared-geometry-passed': 18, 'not-applicable-no-sidechain-heavy-atoms': 1, 'refused': 1},
                ('observed', 'H'): {'declared-geometry-passed': 16, 'unsupported': 3, 'refused': 1}}
    for (arm, operation), counts in expected.items():
        selected = [r for r in result['controls'] if r['arm'] == arm and r['operation'] == operation]
        assert Counter(r['status'] for r in selected) == counts
    no_MET = [r for r in result['controls'] if r['arm'] == 'observed' and r['source_component'] == 'MET']
    assert all(r['status'] == 'refused' and r['reason'] == 'declared interior target not found' for r in no_MET)
    rejected_H = [r['source_component'] for r in result['controls'] if r['arm'] == 'observed' and r['operation'] == 'H'
                  and r['status'] == 'unsupported']
    assert rejected_H == ['ARG', 'ASP', 'TYR']
    assert not result['all20_production_native_quality_or_combined_pipeline_qualified']
    with pytest.raises(FileExistsError):
        controls.run_controls(tmp_path/'controls')
