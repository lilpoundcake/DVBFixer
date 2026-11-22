"""Read-only diagnosis cannot become an implicit ring repair or feasibility claim."""
from __future__ import annotations

import math
from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_all20_controls import HEAVY_POLICY, synthetic_input

from dvbfixer.chemical_reconstruction.candidate_audit import measure_additions
from dvbfixer.chemical_reconstruction.closure_diagnostics import (
    _constraint_value,
    fixed_neighbour_support,
    is_cycle_edge,
    measure_terms,
    necessary_support,
    trace_frames,
)
from dvbfixer.chemical_reconstruction.contracts import Constraint, ConstraintKind, Ownership, Unit
from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan


@pytest.fixture
def case():
    request, authority, _ = synthetic_input('PRO', 'PRO', {})
    additions = tuple(o for o in request.observations if o.identity.component == authority.identity and o.identity.name == 'CD')
    request = replace(request, observations=tuple(o for o in request.observations if o not in additions),
                      ownership=Ownership(frozenset(o.identity for o in additions)))
    return request, compile_peptide_request_plan(request, HEAVY_POLICY), additions


def test_ideal_fixture_terms_match_existing_audit_without_constructing(case):
    request, plan, additions = case
    report = measure_terms(request, plan, additions)
    assert report.audit_findings == measure_additions(request, plan, additions).findings
    assert not report.unexplained_findings
    assert all(row.passed for row in report.terms)
    assert {row.role for row in report.terms} == {'fixed-only', 'generated-involving'}


def test_changed_saved_pose_is_explained_not_repaired(case):
    request, plan, additions = case
    changed = (replace(additions[0], position=tuple(v+2 for v in additions[0].position)),)
    report = measure_terms(request, plan, changed)
    assert report.audit_findings == measure_additions(request, plan, changed).findings
    assert 'authority-bond-geometry' in report.explained_findings
    assert 'authority-angle-geometry' in report.explained_findings
    assert all(row.role == 'generated-involving' for row in report.terms if not row.passed)
    assert changed[0].position == tuple(v+2 for v in additions[0].position)


def test_fixed_only_failure_is_not_attributed_to_missing_atom(case):
    request, _, additions = case
    fixed = tuple(o.identity for o in request.observations if o.identity.component in request.targets)[:2]
    q = Constraint(ConstraintKind.DISTANCE, fixed, Unit.ANGSTROM, lower=100., upper=101.)
    request = replace(request, constraints=(*request.constraints, q))
    # No compilation/admission success is invented for this diagnosis fixture.
    from dvbfixer.chemical_reconstruction.proposals import GeometryPlan, Strategy
    plan = GeometryPlan(request.digest, Strategy.INTERNAL, HEAVY_POLICY)
    report = measure_terms(request, plan, additions)
    row = report.terms[-1]
    assert not row.passed and row.role == 'fixed-only'
    assert 'fixed-constraint-conflict' in report.explained_findings
    assert 'candidate-constraint-conflict' in report.scalar_only_findings


def test_frame_trace_exposes_non_length_enforced_cycle_edge(case):
    request, plan, additions = case
    frames, edges = trace_frames(request, plan, additions)
    assert len(frames) == 1 and frames[0].atom == additions[0].identity
    assert frames[0].anchors == plan.placements[0].anchors
    assert frames[0].anchor_roles == ('fixed-observation',)*3
    assert all(abs(m.delta) < 1e-10 for m in frames[0].metrics)
    assert any(e.is_cycle for e in edges)
    assert all(set(e.atoms) & request.ownership.generated for e in edges)
    assert not hasattr(frames[0], 'candidate')


def test_cycle_classifier_does_not_depend_on_PRO_atom_names(case):
    request, _, additions = case
    a = additions[0].identity
    b, c = (replace(a, name=name) for name in ('X1', 'X2'))
    assert not is_cycle_edge(((a, b), (b, c)), (a, b))
    assert is_cycle_edge(((a, b), (b, c), (c, a)), (a, b))
    with pytest.raises(ValueError):
        is_cycle_edge(((a, b),), (b, c))


def test_necessary_support_is_not_joint_feasibility(case):
    request, plan, _ = case
    rows = fixed_neighbour_support(request, plan.policy)
    assert len(rows) == 1
    support = rows[0][-1]
    assert not isinstance(support, str) and support.distance_shell_possible and support.angle_interval_possible
    assert 'not a jointly feasible pose' in support.scope


def test_distance_shell_rejection_and_distinct_angle_rejection():
    far = necessary_support(5., 1., 1., .1, 90., 5.)
    assert not far.distance_shell_possible and not far.angle_interval_possible
    wrong_angle = necessary_support(1., 1., 1., .001, 170., 1.)
    assert wrong_angle.distance_shell_possible and not wrong_angle.angle_interval_possible


def test_conservative_intervals_contain_known_scalar_witnesses():
    # Scalar law-of-cosines witnesses only, not a generated Cartesian proposal.
    for r1 in (1.4, 1.5, 1.6):
        for r2 in (1.3, 1.4, 1.5):
            for angle in (85., 90., 95.):
                distance = math.sqrt(r1*r1+r2*r2-2*r1*r2*math.cos(math.radians(angle)))
                support = necessary_support(distance, 1.5, 1.4, .1, 90., 5.)
                assert support.distance_shell_possible and support.angle_interval_possible
                assert support.conservative_cosine_interval[0] <= math.cos(math.radians(angle)) <= support.conservative_cosine_interval[1]


@pytest.mark.parametrize('values', [
    (0., 1., 1., .1, 90., 5.), (-1., 1., 1., .1, 90., 5.),
    (1., 1., 1., 1., 90., 5.), (1., 1., 1., .1, 0., 5.),
    (1., 1., 1., .1, 90., 0.), (math.nan, 1., 1., .1, 90., 5.),
    (math.inf, 1., 1., .1, 90., 5.), (1e308, 1., 1., .1, 90., 5.),
])
def test_invalid_or_overflowing_scalar_bounds_refuse(values):
    with pytest.raises(ValueError):
        necessary_support(*values)


@pytest.mark.parametrize('change', ['empty', 'duplicate', 'element', 'unit', 'nonfinite'])
def test_malformed_or_partial_saved_coordinates_refused(case, change):
    request, plan, additions = case
    with pytest.raises(ValueError):
        if change == 'empty':
            changed = ()
        elif change == 'duplicate':
            changed = additions+additions
        elif change == 'nonfinite':
            changed = (replace(additions[0], position=(math.nan, 0., 0.)),)
        else:
            changed = (replace(additions[0], **({'element': 'N'} if change == 'element' else {'unit': Unit.NANOMETER})),)
        measure_terms(request, plan, changed)


def test_wrong_plan_binding_and_duplicate_placements_refuse(case):
    request, plan, additions = case
    with pytest.raises(ValueError):
        measure_terms(request, replace(plan, request_digest='0'*64), additions)
    with pytest.raises(ValueError):
        trace_frames(request, replace(plan, placements=plan.placements*2), additions)


def test_unknown_constraint_kind_is_not_guessed(case):
    _, _, additions = case
    q = Constraint(ConstraintKind.DISTANCE, (additions[0].identity, replace(additions[0].identity, name='X')), Unit.ANGSTROM, lower=1., upper=2.)
    with pytest.raises(ValueError):
        _constraint_value(replace(q, kind='unrecognized'), ((0., 0., 0.), (1., 0., 0.)))


def test_diagnosis_primitives_do_not_call_any_coordinate_constructor(case, monkeypatch):
    from scipy import optimize

    from dvbfixer.chemical_reconstruction import construction
    from dvbfixer.protein_completion import geometry

    def forbidden(*args, **kwargs):
        pytest.fail('read-only diagnosis constructed coordinates')
    for module, name in ((construction, 'propose'), (construction, '_internal'),
                         (geometry, 'extend'), (geometry, 'two_angle_tip'), (optimize, 'least_squares')):
        monkeypatch.setattr(module, name, forbidden)
    request, plan, additions = case
    measure_terms(request, plan, additions)
    trace_frames(request, plan, additions)
    fixed_neighbour_support(request, plan.policy)


def test_missing_or_corrupt_input_refuses_before_new_output(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_closure_diagnosis as controls

    monkeypatch.setattr(controls, 'PROTOCOL', tmp_path/'missing')
    with pytest.raises(FileNotFoundError):
        controls.run_diagnosis(tmp_path/'new-output')
    assert not (tmp_path/'new-output').exists()
    p = tmp_path/'bad'
    p.write_bytes(b'not the frozen source')
    with pytest.raises(ValueError, match='pin'):
        controls._pinned(p, '0'*64)


def test_saved_diagnostic_report_tampering_is_detected(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_closure_diagnosis as controls

    monkeypatch.setattr(controls, 'build_report', lambda: {'schema': 'synthetic-report', 'new_coordinate_proposals': 0})
    (tmp_path/'diagnosis.json').write_text('{"schema":"synthetic-report","new_coordinate_proposals":1}')
    with pytest.raises(ValueError, match='report-mismatch'):
        controls.readback_diagnosis(tmp_path)
