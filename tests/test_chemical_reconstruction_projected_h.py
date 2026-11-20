"""Projected-parent H is a distinct generated-conditioning role, not a torsion parent."""
from __future__ import annotations

from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_all20_controls import FFBONDED, HEAVY_POLICY, synthetic_input

from dvbfixer.chemical_reconstruction import projected_hydrogens
from dvbfixer.chemical_reconstruction.admission import admit
from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    ContractError,
    Observation,
    OpaqueContextPolicy,
    Origin,
    Ownership,
    PreparedContextRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.heavy_conditioning import (
    HeavyConditionedRequest,
    condition_on_heavy,
)
from dvbfixer.chemical_reconstruction.metric_projection import propose_projection
from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan
from dvbfixer.chemical_reconstruction.planar_port import (
    PlanarPortDeclaration,
    load_planar_port_authority,
    planar_port_constraints,
)
from dvbfixer.chemical_reconstruction.projected_conditioning import (
    ProjectedConditionedRequest,
    condition_on_projected,
)
from dvbfixer.chemical_reconstruction.projected_hydrogen_audit import (
    audit_projected_hydrogen_result,
)
from dvbfixer.chemical_reconstruction.projected_hydrogen_contracts import (
    ProjectedHeavyParent,
    ProjectedHydrogenResult,
)
from dvbfixer.chemical_reconstruction.projection_contracts import append_projection_constraints
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def synthetic_parent(*, fractional=False):
    from dvbfixer.chemical_reconstruction.construction import propose

    full, authority, incoming = synthetic_input('PRO', 'PRO', {})
    generated = frozenset(ExactAtomRef(authority.identity, name) for name in ('CB', 'CG', 'CD'))
    original = replace(full, observations=tuple(o for o in full.observations if o.identity not in generated), ownership=Ownership(generated))
    if fractional:
        outside = Observation(ExactAtomRef(ComponentInstanceRef(1, 'X', '999'), 'Q'), 'C', (50., 50., 50.), Unit.ANGSTROM, .5)
        observations = (*original.observations, outside)
        original = PreparedContextRequest(
            original.source_sha256, original.operation,
            (*original.components, Component(outside.identity.component, ChemistryKind.UNKNOWN)), observations,
            original.targets, original.ownership, original.connections, original.constraints, original.ph,
            Origin('prepared-observed-peptide-v1', 'a'*64), context_policy=OpaqueContextPolicy(
                original.source_sha256, 'a'*64, frozenset({outside.identity}), digest_record(observations)))
    plan = compile_peptide_request_plan(original, HEAVY_POLICY)
    baseline = propose(original, plan)
    port_parameters = load_planar_port_authority(FFBONDED)
    n, ca, cd = (ExactAtomRef(authority.identity, name) for name in ('N', 'CA', 'CD'))
    declaration = PlanarPortDeclaration(incoming, n, ca, cd, ('C', 'N', 'XC', 'CT'), authority.digest, port_parameters.digest)
    request = append_projection_constraints(original, planar_port_constraints(original, authority, declaration, port_parameters))
    projection = propose_projection(original, request, plan, baseline, (ca, n, ExactAtomRef(authority.identity, 'C')))
    assert projection.status == 'projected-geometry-passed'
    parent = ProjectedHeavyParent(original, request, projection, authority, declaration, port_parameters)
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    return parent, parameters


def run(parent, parameters):
    return projected_hydrogens.complete_projected_hydrogens(parent, parameters, backbone_types=None)


def readback(parent, parameters, result):
    return audit_projected_hydrogen_result(parent, parameters, None, result)


def test_distinct_parent_and_conditioning_preserve_all_three_roles():
    parent, parameters = synthetic_parent()
    base = condition_on_projected(parent.request, parent.projection)
    assert isinstance(base, ProjectedConditionedRequest) and not isinstance(base, (PreparedContextRequest, HeavyConditionedRequest))
    assert base.original is parent.request and base.projection is parent.projection
    assert base.observations == parent.original.observations + parent.projection.rounded.additions
    assert base.digest != parent.original.digest and base.source_sha256 != parent.original.source_sha256
    assert not base.ownership.generated and parent.original.ownership.generated
    with pytest.raises(ContractError, match='heavy-checkpoint-binding'):
        condition_on_heavy(parent.request, parent.projection.rounded)
    result = run(parent, parameters)
    assert result.status in {'projected-H-geometry-passed', 'H-rejected'}
    assert not readback(parent, parameters, result)
    assert result.raw is not None and result.rounded is not None
    assert result.raw.generated_heavy is parent.projection.rounded.additions
    assert result.rounded.generated_heavy is parent.projection.rounded.additions
    assert len(result.raw.generated_hydrogens) == 7
    assert {o.identity.name for o in result.raw.generated_hydrogens} == {
        a.name for a in parent.peptide_authority.state.graph.atoms if a.element == 'H'}
    assert isinstance(result.view.request, ProjectedConditionedRequest)
    assert not result.md_ready and not result.native_quality_qualified


@pytest.mark.parametrize('mutation', ['partial', 'nonrounded', 'element', 'unit', 'occupancy', 'strategy'])
def test_generated_checkpoint_corruption_never_becomes_source_observation(mutation):
    parent, _ = synthetic_parent()
    p = parent.projection
    atom = p.rounded.additions[0]
    changes = {'nonrounded': {'position': (atom.position[0]+.00001, *atom.position[1:])},
               'element': {'element': 'N'}, 'unit': {'unit': Unit.NANOMETER}, 'occupancy': {'occupancy': .5}}
    if mutation == 'partial':
        checkpoint = replace(p.rounded, additions=p.rounded.additions[:-1])
    elif mutation == 'strategy':
        from dvbfixer.chemical_reconstruction.proposals import Strategy
        checkpoint = replace(p.rounded, strategy=Strategy.TORSION)
    else:
        checkpoint = replace(p.rounded, additions=(replace(atom, **changes[mutation]), *p.rounded.additions[1:]))
    with pytest.raises(ValueError):
        condition_on_projected(parent.request, replace(p, rounded=checkpoint))


@pytest.mark.parametrize('mutation', ['source', 'port', 'projection-cost', 'projection-status', 'authority', 'parameters', 'types'])
def test_parent_binding_mutations_fail_before_H_construction(mutation, monkeypatch):
    parent, parameters = synthetic_parent()
    if mutation == 'source':
        parent = replace(parent, original=replace(parent.original, observations=tuple(reversed(parent.original.observations))))
    elif mutation == 'port':
        parent = replace(parent, request=replace(parent.request, constraints=parent.request.constraints[:-1]))
    elif mutation == 'projection-cost':
        parent = replace(parent, projection=replace(parent.projection, final_squared_residual=999.))
    elif mutation == 'projection-status':
        parent = replace(parent, projection=replace(parent.projection, status='rejected'))
    elif mutation == 'authority':
        parent = replace(parent, port_declaration=replace(parent.port_declaration, peptide_authority_digest='0'*64))
    elif mutation == 'parameters':
        parent = replace(parent, port_parameters=replace(parent.port_parameters, source='synthetic foreign authority path'))
    else:
        parent = replace(parent, port_declaration=replace(parent.port_declaration, atom_types=('C', 'N', 'CT', 'XC')))
    monkeypatch.setattr(projected_hydrogens, 'propose', lambda *a: pytest.fail('H constructed from invalid parent'))
    result = run(parent, parameters)
    assert result.status == 'parent-refused' and result.view is None and result.raw is None
    assert not readback(parent, parameters, result)


def test_H_failure_has_exactly_one_proposal_and_never_moves_or_projects_heavy(monkeypatch):
    parent, parameters = synthetic_parent()
    before = digest_record(parent)
    ordinary = projected_hydrogens.propose
    calls = []
    def reject(request, plan):
        calls.append(request)
        return replace(ordinary(request, plan), status='rejected', reason='synthetic-H-negative-no-heavy-change')
    monkeypatch.setattr(projected_hydrogens, 'propose', reject)
    result = run(parent, parameters)
    assert len(calls) == 1 and result.status == 'H-rejected' and not result.audit.passed
    assert digest_record(parent) == before
    assert result.raw.generated_heavy == parent.projection.rounded.additions
    assert not readback(parent, parameters, result)


def test_constructor_solver_free_combined_audit_and_strict_DTO_roundtrip(monkeypatch):
    import json

    from scipy import optimize
    from scripts.chemical_reconstruction_ccd_controls import _json

    from dvbfixer.chemical_reconstruction import construction, metric_projection
    from dvbfixer.chemical_reconstruction.research_decode import decode_research_record

    parent, parameters = synthetic_parent()
    result = run(parent, parameters)
    saved_parent = decode_research_record(json.dumps(_json(parent)).encode(), ProjectedHeavyParent)
    saved_result = decode_research_record(json.dumps(_json(result)).encode(), ProjectedHydrogenResult)
    assert digest_record(saved_parent) == digest_record(parent) and digest_record(saved_result) == digest_record(result)
    def forbidden(*args, **kwargs):
        pytest.fail('combined audit called a constructor/solver')
    for module, name in ((construction, 'propose'), (construction, '_internal'), (projected_hydrogens, 'propose'),
                         (projected_hydrogens, 'complete_projected_hydrogens'), (metric_projection, 'propose_projection'),
                         (metric_projection, '_solve'), (optimize, 'least_squares')):
        monkeypatch.setattr(module, name, forbidden)
    assert not readback(saved_parent, parameters, saved_result)


def test_fractional_context_is_proven_and_conditioning_is_not_prepared_source():
    from types import SimpleNamespace

    from dvbfixer.chemical_reconstruction.source_preparation import audit_source_transport

    parent, parameters = synthetic_parent(fractional=True)
    result = run(parent, parameters)
    request = result.view.request
    assert request.allowed_fractional_atoms == parent.original.context_policy.fractional_atoms
    assert not any(f.code == 'partial-occupancy-or-altloc-requires-explicit-selection' for f in admit(request).findings)
    assert audit_source_transport(None, SimpleNamespace(request=request)).findings == ('prepared-context-request-required',)
    atom = next(iter(request.allowed_fractional_atoms))
    with pytest.raises(ValueError):
        replace(request, ownership=Ownership(request.ownership.generated | {atom}))
    with pytest.raises(ValueError):
        replace(request, original=request)
    assert not readback(parent, parameters, result)


@pytest.mark.parametrize('mutation', ['source-digest', 'port-digest', 'parent-digest', 'heavy', 'H', 'rounding', 'CIP', 'pair', 'status', 'plan'])
def test_saved_combined_corruption_is_rejected(mutation):
    parent, parameters = synthetic_parent()
    result = run(parent, parameters)
    if mutation in {'source-digest', 'port-digest', 'parent-digest'}:
        key = {'source-digest': 'original_request_digest', 'port-digest': 'port_request_digest', 'parent-digest': 'parent_digest'}[mutation]
        result = replace(result, **{key: '0'*64})
    elif mutation in {'heavy', 'H'}:
        key = 'generated_heavy' if mutation == 'heavy' else 'generated_hydrogens'
        atoms = getattr(result.raw, key)
        result = replace(result, raw=replace(result.raw, **{key: atoms[:-1]}))
    elif mutation == 'rounding':
        result = replace(result, hydrogen_result=replace(result.hydrogen_result, rounded=result.hydrogen_result.raw))
    elif mutation == 'CIP':
        result = replace(result, audit=replace(result.audit, raw_cip=replace(result.audit.raw_cip, findings=('forged',))))
    elif mutation == 'pair':
        coverage = result.audit.pair_coverage
        assert coverage.hydrogen_heavy > 0
        result = replace(result, audit=replace(result.audit, pair_coverage=replace(coverage, hydrogen_heavy=0)))
    elif mutation == 'status':
        result = replace(result, status='H-refused')
    else:
        view = result.view
        result = replace(result, view=replace(view, plan=replace(view.plan, policy=replace(view.plan.policy, angle_tolerance_degrees=14.))))
    assert readback(parent, parameters, result)


def test_missing_frozen_parent_refuses_before_output_or_H_science(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_projected_h_controls as controls

    monkeypatch.setattr(controls, 'check_pins', lambda: {'entries': [
        {'parent_artifact': str(tmp_path/'missing-parent'), 'parent_sha256': '0'*64}]})
    monkeypatch.setattr(controls, 'complete_projected_hydrogens', lambda *a, **k: pytest.fail('fabricated missing projection parent'))
    output = tmp_path/'new-output'
    with pytest.raises(FileNotFoundError):
        controls.run_controls(output)
    assert not output.exists()


def test_partial_H_failure_keeps_negative_raw_without_complete_pair_claim(monkeypatch):
    from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate

    parent, parameters = synthetic_parent()
    ordinary = projected_hydrogens.propose
    def partial(request, plan):
        result = ordinary(request, plan)
        raw = replace(result.raw, additions=result.raw.additions[:-1])
        return replace(result, status='construction-failed', reason='synthetic-partial-H', raw=raw,
                       raw_audit=audit_candidate(request, plan, raw), rounded=None, rounded_audit=None)
    monkeypatch.setattr(projected_hydrogens, 'propose', partial)
    result = run(parent, parameters)
    assert result.status == 'H-rejected' and result.raw is not None and result.rounded is None
    assert result.audit.pair_coverage is None and 'staged-pair-H-inventory-incomplete' in result.audit.findings
    assert not readback(parent, parameters, result)
