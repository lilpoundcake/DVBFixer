"""Bounded generated-only projection is not H, native, or production qualification."""
from __future__ import annotations

from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_all20_controls import FFBONDED, HEAVY_POLICY, synthetic_input

from dvbfixer.chemical_reconstruction import metric_projection
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    Constraint,
    ConstraintKind,
    Origin,
    Ownership,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan
from dvbfixer.chemical_reconstruction.planar_port import (
    PlanarPortDeclaration,
    load_planar_port_authority,
    planar_port_constraints,
)
from dvbfixer.chemical_reconstruction.projection_audit import (
    audit_projection_result,
    measure_projection,
)
from dvbfixer.chemical_reconstruction.projection_contracts import append_projection_constraints
from dvbfixer.domain.structure_identity import ExactAtomRef


def case(mask=('CB', 'CG', 'CD'), *, transform=None):
    full, authority, incoming = synthetic_input('PRO', 'PRO', {})
    if transform:
        observations = tuple(replace(o, position=transform(o.position)) for o in full.observations)
        pin = digest_record(('synthetic-transformed-projection-v1', observations))
        full = replace(full, observations=observations, source_sha256=pin, origin=Origin('synthetic-transformed-projection-v1', pin))
    generated = frozenset(ExactAtomRef(authority.identity, name) for name in mask)
    original = replace(full, observations=tuple(o for o in full.observations if o.identity not in generated), ownership=Ownership(generated))
    plan = compile_peptide_request_plan(original, HEAVY_POLICY)
    baseline = propose(original, plan)
    parameters = load_planar_port_authority(FFBONDED)
    n, ca, cd = (ExactAtomRef(authority.identity, name) for name in ('N', 'CA', 'CD'))
    declaration = PlanarPortDeclaration(incoming, n, ca, cd, ('C', 'N', 'XC', 'CT'), authority.digest, parameters.digest)
    extra = planar_port_constraints(original, authority, declaration, parameters)
    request = append_projection_constraints(original, extra)
    anchors = ca, n, ExactAtomRef(authority.identity, 'C')
    return original, request, plan, baseline, anchors


def test_real_port_projection_keeps_all_observations_and_whole_ring_mask():
    inputs = case()
    original, request, _, baseline, _ = inputs
    before = digest_record(original.observations)
    result = metric_projection.propose_projection(*inputs)
    assert baseline.status == 'declared-geometry-passed'
    assert result.status == 'projected-geometry-passed'
    assert result.audit.passed and result.solver_success
    assert request.observations is original.observations and digest_record(request.observations) == before
    assert {a.identity for a in result.raw.additions} == original.ownership.generated
    assert result.audit.maximum_raw_movement_angstrom <= 2.
    assert not audit_projection_result(*inputs, result)
    assert not result.md_ready and not result.heavy_H_or_native_quality_qualified
    assert propose(request, result.plan).reason == 'strategy-needs-dedicated-constructor'


@pytest.mark.parametrize('mask', [('CD',), ('CG', 'CD'), ('CB', 'CG', 'CD')])
def test_partial_and_complete_generated_ring_subsets_have_exact_ownership(mask):
    inputs = case(mask)
    result = metric_projection.propose_projection(*inputs)
    assert result.status in {'projected-geometry-passed', 'rejected'}
    assert result.raw is not None and {a.identity.name for a in result.raw.additions} == set(mask)
    assert not audit_projection_result(*inputs, result)


def test_constructor_free_readback_never_imports_or_calls_a_solver(monkeypatch):
    from scipy import optimize

    inputs = case()
    result = metric_projection.propose_projection(*inputs)
    def forbidden(*args, **kwargs):
        pytest.fail('audit called coordinate construction')
    monkeypatch.setattr(optimize, 'least_squares', forbidden)
    monkeypatch.setattr(metric_projection, '_solve', forbidden)
    monkeypatch.setattr(metric_projection, 'propose_projection', forbidden)
    assert not audit_projection_result(*inputs, result)


def test_proper_rigid_frame_equivariance_without_template_or_source_refitting():
    first = metric_projection.propose_projection(*case())
    def transform(p):
        return -p[1]+5., p[0]-2., p[2]+7.
    second = metric_projection.propose_projection(*case(transform=transform))
    assert first.status == second.status == 'projected-geometry-passed'
    for a, b in zip(first.raw.additions, second.raw.additions, strict=True):
        assert b.position == pytest.approx(transform(a.position), abs=1e-5)


def test_fixed_source_constraint_defect_cannot_be_repaired_by_generated_coordinates():
    original, _, plan, baseline, anchors = case()
    constraint = Constraint(ConstraintKind.DISTANCE, anchors[:2], Unit.ANGSTROM, 4., 4.2)
    request = append_projection_constraints(original, (constraint,))
    result = metric_projection.propose_projection(original, request, plan, baseline, anchors)
    assert result.status == 'rejected' and not result.audit.passed
    assert request.observations is original.observations
    assert not audit_projection_result(original, request, plan, baseline, anchors, result)


def test_optional_optimizer_absence_is_retained_refusal_not_fallback(monkeypatch):
    inputs = case()
    def unavailable(*args):
        raise ImportError('synthetic absent SciPy control')
    monkeypatch.setattr(metric_projection, '_solve', unavailable)
    result = metric_projection.propose_projection(*inputs)
    assert result.status == 'solver-failed' and result.raw is None and result.plan is not None
    assert not audit_projection_result(*inputs, result)


def test_nonconvergence_is_not_accepted_and_never_restarts(monkeypatch):
    calls = []
    def no_convergence(function, variables):
        calls.append(variables)
        return metric_projection.ProjectionSolution((0.,)*variables, False, 0, 'synthetic max-evaluation control', 500)
    monkeypatch.setattr(metric_projection, '_solve', no_convergence)
    inputs = case()
    result = metric_projection.propose_projection(*inputs)
    assert len(calls) == 1 and result.status == 'rejected'
    assert 'projection-solver-not-converged' in result.audit.findings
    assert not audit_projection_result(*inputs, result)


@pytest.mark.parametrize('mutation', ['raw', 'rounding', 'cost', 'success', 'nfev', 'audit', 'plan', 'status'])
def test_independent_auditor_rejects_saved_transport_and_metadata_corruption(mutation):
    inputs = case()
    result = metric_projection.propose_projection(*inputs)
    if mutation == 'raw':
        atom = result.raw.additions[0]
        raw = replace(result.raw, additions=(replace(atom, position=(20., 20., 20.)), *result.raw.additions[1:]))
        result = replace(result, raw=raw)
    elif mutation == 'rounding':
        result = replace(result, rounded=result.raw)
    elif mutation == 'cost':
        result = replace(result, final_squared_residual=999.)
    elif mutation == 'success':
        result = replace(result, solver_success=False)
    elif mutation == 'nfev':
        result = replace(result, nfev=True)
    elif mutation == 'audit':
        result = replace(result, audit=replace(result.audit, findings=('forged',)))
    elif mutation == 'plan':
        term = result.plan.terms[0]
        result = replace(result, plan=replace(result.plan, terms=(replace(term, value=term.value+.01), *result.plan.terms[1:])))
    else:
        result = replace(result, status='rejected')
    assert audit_projection_result(*inputs, result)


def test_movement_cap_is_independent_of_solver_cost_or_convergence():
    original, request, plan, baseline, anchors = inputs = case()
    result = metric_projection.propose_projection(*inputs)
    atom = result.raw.additions[0]
    moved = replace(result.raw, additions=(replace(atom, position=(20., 20., 20.)), *result.raw.additions[1:]))
    audit = measure_projection(request, result.plan, moved, moved, True)
    assert not audit.passed and 'projection-per-atom-movement-cap' in audit.findings
    assert audit.maximum_raw_movement_angstrom > 2.


def test_residual_budget_and_unsupported_torsion_fail_closed_before_solving(monkeypatch):
    original, _, plan, baseline, anchors = case()
    monkeypatch.setattr(metric_projection, '_solve', lambda *args: pytest.fail('solver called for refused plan'))
    extra = Constraint(ConstraintKind.DISTANCE, anchors[:2], Unit.ANGSTROM, 1., 2.)
    request = append_projection_constraints(original, (extra,)*512)
    result = metric_projection.propose_projection(original, request, plan, baseline, anchors)
    assert result.status == 'refused' and 'budget' in result.reason
    refs = (*anchors, next(iter(original.ownership.generated)))
    torsion = Constraint(ConstraintKind.TORSION, refs, Unit.DEGREE, -180., 180.)
    request = append_projection_constraints(original, (torsion,))
    result = metric_projection.propose_projection(original, request, plan, baseline, anchors)
    assert result.status == 'refused' and 'unqualified-constraint' in result.reason


@pytest.mark.parametrize('mutation', ['types', 'incoming', 'inventory', 'charge', 'link'])
def test_planar_port_requires_explicit_source_graph_typing_and_actual_link(mutation):
    full, authority, incoming = synthetic_input('PRO', 'PRO', {})
    cd = ExactAtomRef(authority.identity, 'CD')
    original = replace(full, observations=tuple(o for o in full.observations if o.identity != cd), ownership=Ownership(frozenset({cd})))
    parameters = load_planar_port_authority(FFBONDED)
    declaration = PlanarPortDeclaration(incoming, ExactAtomRef(authority.identity, 'N'), ExactAtomRef(authority.identity, 'CA'), cd,
                                       ('C', 'N', 'XC', 'CT'), authority.digest, parameters.digest)
    if mutation == 'types':
        declaration = replace(declaration, atom_types=('C', 'N', 'CT', 'XC'))
    elif mutation == 'incoming':
        declaration = replace(declaration, incoming=cd)
    elif mutation == 'inventory':
        declaration = replace(declaration, peptide_authority_digest='0'*64)
    elif mutation == 'charge':
        # A different declared target state is not silently substituted.
        _, other, _ = synthetic_input('HIS', 'HIP', {})
        authority = other
    else:
        original = replace(original, connections=())
    with pytest.raises(ValueError):
        planar_port_constraints(original, authority, declaration, parameters)


def test_strict_projection_DTO_roundtrip_has_distinct_strategy_and_lineage():
    import json

    from scripts.chemical_reconstruction_ccd_controls import _json

    from dvbfixer.chemical_reconstruction.projection_contracts import ProjectionResult
    from dvbfixer.chemical_reconstruction.research_decode import decode_research_record

    inputs = case()
    result = metric_projection.propose_projection(*inputs)
    decoded = decode_research_record(json.dumps(_json(result)).encode(), ProjectionResult)
    assert digest_record(decoded) == digest_record(result)
    assert not audit_projection_result(*inputs, decoded)


def test_frozen_three_source_control_dispositions_and_no_reference_before_projection(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_projection_controls as controls

    old_projection, old_reference = controls.propose_projection, controls.heldout_after_proposal
    completed = []
    def proposal(*args, **kwargs):
        result = old_projection(*args, **kwargs)
        completed.append(result.status)
        return result
    def reference(*args, **kwargs):
        assert completed
        return old_reference(*args, **kwargs)
    monkeypatch.setattr(controls, 'propose_projection', proposal)
    monkeypatch.setattr(controls, 'heldout_after_proposal', reference)
    result = controls.run_controls(tmp_path/'controls')
    assert result['declared'] == 3 and result['H_proposals'] == 0
    by_id = {row['id']: row for row in result['controls']}
    assert by_id['1UBQ']['status'] == 'projected-geometry-passed'
    assert by_id['1FR2']['reason'] == 'actual-incoming-template-unresolved:HIS'
    assert by_id['1A3N']['reason'] == 'explicit-altloc-selection-required'
    assert len(completed) == 1
    with pytest.raises(FileExistsError):
        controls.run_controls(tmp_path/'controls')


def test_source_binding_variable_budget_and_frozen_solver_configuration():
    original, request, plan, baseline, anchors = inputs = case()
    result = metric_projection.propose_projection(*inputs)
    with pytest.raises(ValueError, match='frozen-plan'):
        replace(result.plan, max_nfev=1000)
    with pytest.raises(ValueError, match='frozen-plan'):
        replace(result.plan, per_atom_movement_cap_angstrom=3.)
    changed = replace(request, observations=tuple(reversed(request.observations)))
    bad = metric_projection.propose_projection(original, changed, plan, baseline, anchors)
    assert bad.status == 'refused' and 'binding' in bad.reason
    generated = frozenset(ExactAtomRef(next(iter(original.targets)), f'FAKE{i}') for i in range(17))
    oversized = replace(original, ownership=Ownership(generated))
    constraint = Constraint(ConstraintKind.DISTANCE, anchors[:2], Unit.ANGSTROM, 1., 2.)
    bad = metric_projection.propose_projection(oversized, append_projection_constraints(oversized, (constraint,)), plan, baseline, anchors)
    assert bad.status == 'refused' and 'bounded-fixed-heavy' in bad.reason
