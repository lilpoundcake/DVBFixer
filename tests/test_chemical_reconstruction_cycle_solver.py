"""Frozen R1 numerical unit controls, NOT the tracked-source 30-cell study."""
from __future__ import annotations

import inspect
import math
from dataclasses import replace

import pytest

from dvbfixer.chemical_reconstruction import cycle_solver
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Constraint,
    ConstraintKind,
    Observation,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.cycle_contract import compile_cycle_problem
from dvbfixer.chemical_reconstruction.cycle_initialization import initialize_cycle_problem
from dvbfixer.chemical_reconstruction.cycle_residuals import (
    CycleNumericalError,
    residual_values,
    squared_residual,
)
from dvbfixer.chemical_reconstruction.cycle_solve_plan import (
    SETTINGS,
    CycleSolvePlan,
    _branch_plan,
    compile_cycle_solve_plan,
    point_order,
    record_bytes,
)
from dvbfixer.chemical_reconstruction.cycle_solver import (
    CycleConstructionResult,
    CycleEvaluationTracker,
    CycleExecutionAborted,
    CycleSolution,
    CycleSolverStats,
    solve_cycle_plan,
)
from dvbfixer.chemical_reconstruction.cycle_solver_audit import (
    CycleConstructionAudit,
    audit_cycle_construction,
    measure_cycle_construction,
)
from dvbfixer.chemical_reconstruction.research_decode import decode_research_record
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from tests.test_chemical_reconstruction_cycle_audit import source_fixture


def case(mask=('CD',)):
    problem, prepared, metadata, witness = source_fixture(mask=mask)
    initializer = initialize_cycle_problem(problem)
    return compile_cycle_solve_plan(initializer, prepared, metadata), prepared, witness


def mock_witness_solver(plan, witness, monkeypatch, *, converged=True):
    """Known synthetic coordinate transport only, not a numerical proposal claim."""
    indices = []
    values = {o.identity: o.position for o in witness}
    frame = plan.initializer.system.frame
    def solve(function, variables):
        index = len(indices)
        indices.append(index)
        function(tuple(0. for _ in range(variables)))
        start = plan.initializer.starts[index]
        x = tuple(float(sum((values[o.identity][k]-o.position[k])*axis[k] for k in range(3)))
                  for o in start.additions for axis in frame.axes)
        return CycleSolution(x, CycleSolverStats(converged, 1 if converged else 0, 'synthetic mock', 1, 1))
    monkeypatch.setattr(cycle_solver, '_solve', solve)
    return indices


def test_plan_is_source_authenticated_and_frozen_without_an_old_parent_gate():
    plan, prepared, _ = case(('CB', 'CD'))
    assert len(plan.branches) == len(plan.initial_audits) == 8
    assert len(plan.initializer.problem.coverage.blocks) == 2
    assert not all(a.passed for a in plan.initial_audits)
    assert all(b.ready for b in plan.branches)
    assert plan.settings == SETTINGS and not plan.constructor_approved
    assert list(inspect.signature(solve_cycle_plan).parameters) == ['plan', 'on_branch']
    assert prepared.decode() != plan.metadata.masked_pdb
    assert 'prepared' not in {f for f in inspect.signature(CycleSolvePlan).parameters}


@pytest.mark.parametrize('mask', [('CD',), ('CB', 'CD'), ('CB', 'CG', 'CD')])
def test_real_synthetic_solver_executes_all_eight_and_independent_audit_retains_every_outcome(mask):
    plan, prepared, _ = case(mask)
    seen = []
    result = solve_cycle_plan(plan, on_branch=seen.append)
    assert tuple(seen) == result.branches and len(seen) == 8
    report = measure_cycle_construction(plan, prepared, result)
    assert len(report.branches) == 8 and report.status in {'selected-declared-geometry-only', 'rejected'}
    assert not result.constructor_approved and not result.md_ready
    assert not report.H_parent_approved and not report.md_ready and not report.source_coverage_verified
    for branch in result.branches:
        assert branch.solver_calls <= 1 and branch.optimizer_residual_calls <= 50001
        assert branch.status in {'coordinates-returned', 'numerical-failed'}
        if branch.status == 'coordinates-returned':
            assert branch.external_residual_calls == 2
            assert branch.optimizer_residual_calls >= branch.stats.nfev
            assert branch.last_finite_x == branch.returned_x
    audit_cycle_construction(plan, prepared, result, report)


def test_selection_is_objective_then_movement_then_index_without_native_scores(monkeypatch):
    plan, prepared, witness = case()
    calls = mock_witness_solver(plan, witness, monkeypatch)
    result = solve_cycle_plan(plan)
    assert calls == list(range(8))
    report = measure_cycle_construction(plan, prepared, result)
    eligible = [(squared_residual(b.final_residuals), a.pair.raw_maximum_movement_angstrom, b.index)
                for b, a in zip(result.branches, report.branches) if a.selectable]
    assert eligible and report.selected_index == min(eligible)[2]
    assert all(b.solver_calls == 1 for b in result.branches)


def test_geometry_pass_without_solver_convergence_is_not_selected(monkeypatch):
    plan, prepared, witness = case()
    mock_witness_solver(plan, witness, monkeypatch, converged=False)
    result = solve_cycle_plan(plan)
    report = measure_cycle_construction(plan, prepared, result)
    assert all(a.pair.passed for a in report.branches)
    assert report.selected_index is None and report.status == 'rejected'


def test_callback_budget_counts_blocked_breach_and_keeps_last_finite_snapshot():
    plan, _, _ = case()
    tracker = CycleEvaluationTracker(plan, 0)
    x = (0., 0., 0.)
    tracker.evaluate(x, optimizer=False)
    tracker.optimizer_calls = 50000  # Direct counter-boundary test, no changed solver policy.
    with pytest.raises(CycleNumericalError, match='call-budget'):
        tracker.callback(x)
    assert tracker.optimizer_calls == 50001 and tracker.last_x == x and tracker.external_calls == 1


def test_budget_breach_precedes_nonfinite_input_classification():
    plan, _, _ = case()
    tracker = CycleEvaluationTracker(plan, 0)
    tracker.evaluate((0., 0., 0.), optimizer=False)
    tracker.optimizer_calls = 50000
    with pytest.raises(CycleNumericalError, match='call-budget'):
        tracker.callback((float('nan'), 0., 0.))
    assert tracker.optimizer_calls == 50001 and tracker.last_x == (0., 0., 0.)


def test_missing_scipy_is_unavailable_for_all_eight_with_no_optimization(monkeypatch):
    plan, prepared, _ = case()
    def missing(*args):
        raise ModuleNotFoundError('synthetic missing SciPy', name='scipy')
    monkeypatch.setattr(cycle_solver, '_solve', missing)
    result = solve_cycle_plan(plan)
    assert len(result.branches) == 8 and all(b.status == 'unavailable' and b.solver_calls == 0 for b in result.branches)
    report = measure_cycle_construction(plan, prepared, result)
    assert report.status == 'unavailable' and report.selected_index is None


def test_numerical_failure_retains_partial_not_a_scored_or_accepted_proposal(monkeypatch):
    plan, prepared, _ = case()
    def fail(function, variables):
        function(tuple(0. for _ in range(variables)))
        raise CycleNumericalError('synthetic numeric failure')
    monkeypatch.setattr(cycle_solver, '_solve', fail)
    result = solve_cycle_plan(plan)
    assert all(b.status == 'numerical-failed' and b.last_finite_x is not None and not b.raw for b in result.branches)
    assert measure_cycle_construction(plan, prepared, result).selected_index is None


@pytest.mark.parametrize('kind', ['programming', 'unrelated-dependency'])
def test_unexpected_failure_aborts_with_partial_evidence_not_a_chemical_rejection(kind, monkeypatch):
    plan, prepared, _ = case()
    def fail(function, variables):
        function(tuple(0. for _ in range(variables)))
        if kind == 'unrelated-dependency':
            raise ModuleNotFoundError('unexpected dependency', name='other')
        raise RuntimeError('synthetic bug')
    monkeypatch.setattr(cycle_solver, '_solve', fail)
    with pytest.raises(CycleExecutionAborted) as error:
        solve_cycle_plan(plan)
    result = error.value.partial
    assert len(result.branches) == 1 and result.branches[0].status == 'aborted' and result.aborted_reason
    assert result.branches[0].last_finite_x is not None
    assert measure_cycle_construction(plan, prepared, result).status == 'aborted'


def test_sink_failure_preserves_completed_branch_without_fabricating_unrun_attempts(monkeypatch):
    plan, prepared, witness = case()
    calls = mock_witness_solver(plan, witness, monkeypatch)
    def sink(branch):
        raise OSError('synthetic quota error')
    with pytest.raises(CycleExecutionAborted) as error:
        solve_cycle_plan(plan, on_branch=sink)
    assert calls == [0] and len(error.value.partial.branches) == 1
    report = measure_cycle_construction(plan, prepared, error.value.partial)
    assert report.status == 'aborted' and report.selected_index is None


@pytest.mark.parametrize('value', [float('nan'), float('inf'), 10.])
def test_nonfinite_or_out_of_box_solver_return_is_a_retained_failure(value, monkeypatch):
    plan, prepared, _ = case()
    def invalid(function, variables):
        function(tuple(0. for _ in range(variables)))
        return CycleSolution((value, 0., 0.), CycleSolverStats(True, 1, 'synthetic invalid return', 1, 1))
    monkeypatch.setattr(cycle_solver, '_solve', invalid)
    result = solve_cycle_plan(plan)
    assert all(b.status == 'numerical-failed' and not b.raw for b in result.branches)
    assert measure_cycle_construction(plan, prepared, result).selected_index is None


def test_saved_readback_disables_optimizer_initializer_and_materializer(monkeypatch):
    import numpy as np

    from dvbfixer.chemical_reconstruction import cycle_initialization

    plan, prepared, witness = case()
    mock_witness_solver(plan, witness, monkeypatch)
    result = solve_cycle_plan(plan)
    report = measure_cycle_construction(plan, prepared, result)
    def forbidden(*args, **kwargs):
        pytest.fail('saved readback constructed/optimized coordinates')
    for module, name in ((cycle_solver, '_solve'), (cycle_solver, 'vector_additions'),
                         (cycle_initialization, 'initialize_cycle_problem'), (cycle_initialization, '_direct_solve'),
                         (cycle_initialization, '_make_start'), (np.linalg, 'solve')):
        monkeypatch.setattr(module, name, forbidden)
    decoded_plan = decode_research_record(record_bytes(plan), CycleSolvePlan)
    decoded_result = decode_research_record(record_bytes(result), CycleConstructionResult)
    decoded_report = decode_research_record(record_bytes(report), CycleConstructionAudit)
    audit_cycle_construction(decoded_plan, prepared, decoded_result, decoded_report)


@pytest.mark.parametrize('mutation', ['pairs', 'far-count', 'far-proof', 'settings', 'initial-audit', 'branch-order'])
def test_plan_cache_tampering_fails_exact_binding(mutation):
    plan, _, _ = case()
    with pytest.raises(ValueError):
        if mutation == 'settings':
            replace(plan, settings=replace(plan.settings, max_nfev=601))
        elif mutation == 'initial-audit':
            replace(plan, initial_audits=(replace(plan.initial_audits[0], source_metadata_digest='0'*64), *plan.initial_audits[1:]))
        elif mutation == 'branch-order':
            replace(plan, branches=tuple(reversed(plan.branches)))
        else:
            b = plan.branches[0]
            changed = replace(b, pairs=((0, 1),)) if mutation == 'pairs' else replace(b, omitted_far_count=b.omitted_far_count+1) if mutation == 'far-count' else replace(b, omitted_far_digest='0'*64)
            replace(plan, branches=(changed, *plan.branches[1:]))


@pytest.mark.parametrize('mutation', ['initial', 'final', 'partial', 'calls', 'status', 'raw', 'rounded', 'environment', 'plan', 'selection'])
def test_trace_and_selection_tampering_is_rejected(mutation, monkeypatch):
    plan, prepared, witness = case()
    mock_witness_solver(plan, witness, monkeypatch)
    result = solve_cycle_plan(plan)
    report = measure_cycle_construction(plan, prepared, result)
    with pytest.raises(ValueError):
        if mutation == 'selection':
            bad = replace(report, selected_index=(report.selected_index+1) % 8)
            audit_cycle_construction(plan, prepared, result, bad)
            return
        if mutation == 'environment':
            altered = replace(result, environment=())
        elif mutation == 'plan':
            altered = replace(result, plan_digest='0'*64)
        else:
            b = result.branches[0]
            if mutation == 'initial':
                b = replace(b, initial_residuals=(0.,))
            elif mutation == 'final':
                b = replace(b, final_residuals=(0.,))
            elif mutation == 'partial':
                b = replace(b, last_finite_residuals=(0.,))
            elif mutation == 'calls':
                b = replace(b, optimizer_residual_calls=0)
            elif mutation == 'status':
                b = replace(b, status='unavailable')
            elif mutation == 'raw':
                b = replace(b, raw=(replace(b.raw[0], position=(0., 0., 0.)),))
            else:
                b = replace(b, rounded=(replace(b.rounded[0], position=(0., 0., 0.)),))
            altered = replace(result, branches=(b, *result.branches[1:]))
        audit_cycle_construction(plan, prepared, altered, report)


def test_nonzero_lower_plane_bound_is_refused_before_optimizer():
    plan, prepared, _ = case()
    problem = plan.initializer.problem
    target = next(iter(problem.source.targets))
    # The existing Constraint owner already refuses this before R0/R1. Do not
    # weaken it merely to exercise the numerical compiler's defensive check.
    with pytest.raises(ValueError, match='plane distance tolerance'):
        Constraint(ConstraintKind.PLANARITY, tuple(ExactAtomRef(target, n) for n in ('N', 'CA', 'CG', 'CD')),
                   Unit.ANGSTROM, .001, .05)


def test_pair_budget_refuses_without_a_truncated_numerical_vector():
    plan, _, _ = case(('CB', 'CG', 'CD'))
    problem = plan.initializer.problem
    context = ComponentInstanceRef(1, 'z', '999')
    # 3 generated atoms * 22,000 opaque observations exceed the frozen 65,536
    # scalar ceiling; this is an abstract source-contract budget unit fixture.
    extras = tuple(Observation(ExactAtomRef(context, f'Q{i}'), 'C', (0., 0., 0.), Unit.ANGSTROM) for i in range(22000))
    source = replace(problem.source, components=(*problem.source.components, Component(context, ChemistryKind.UNKNOWN)),
                     observations=(*problem.source.observations, *extras))
    initializer = initialize_cycle_problem(compile_cycle_problem(source, problem.anchors))
    branch = _branch_plan(initializer, 0)
    assert branch.residual_count > 65536 and not branch.ready and not branch.pairs


def test_far_pair_proof_is_complete_and_vector_does_not_change_with_evaluation():
    extra = Observation(ExactAtomRef(ComponentInstanceRef(1, 'z', '999'), 'Q'), 'C', (50., 50., 50.), Unit.ANGSTROM)
    problem, prepared, metadata, _ = source_fixture(extra=(extra,))
    plan = compile_cycle_solve_plan(initialize_cycle_problem(problem), prepared, metadata)
    for branch in plan.branches:
        assert branch.omitted_far_count >= 1 and branch.minimum_omitted_distance_angstrom > 8.7
        order = point_order(plan.initializer)
        assert all(extra.identity not in (order[a], order[b]) for a, b in branch.pairs)
        assert len(residual_values(plan, branch.index, (0., 0., 0.))) == branch.residual_count
        assert len(residual_values(plan, branch.index, (.01, .01, .01))) == branch.residual_count


def test_no_complete_source_or_reference_argument_can_be_passed_to_numeric_entry():
    plan, prepared, witness = case()
    for field, value in (('prepared', prepared), ('reference', witness), ('seed', 7), ('candidate', witness)):
        with pytest.raises(TypeError):
            solve_cycle_plan(plan, **{field: value})


def test_frozen_scipy_options_and_call_count_are_explicit(monkeypatch):
    from types import SimpleNamespace

    import numpy as np
    from scipy import optimize

    captured = []
    def fake(function, x, **kwargs):
        captured.append(kwargs)
        function(np.array(x, dtype=float))
        return SimpleNamespace(x=np.array(x, dtype=float), success=True, status=1, nfev=1, njev=1, message='mock')
    monkeypatch.setattr(optimize, 'least_squares', fake)
    plan, prepared, _ = case()
    result = solve_cycle_plan(plan)
    assert len(captured) == 8
    expected = dict(bounds=(-8/math.sqrt(3), 8/math.sqrt(3)), method='trf', jac='2-point', tr_solver='exact',
                    x_scale=1., diff_step=1e-6, loss='linear', ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=600)
    assert all(digest_record(tuple(sorted(row.items()))) == digest_record(tuple(sorted(expected.items()))) for row in captured)
    measure_cycle_construction(plan, prepared, result)


def test_forged_initial_geometry_cache_is_rejected_by_external_source_readback(monkeypatch):
    plan, prepared, witness = case()
    mock_witness_solver(plan, witness, monkeypatch)
    result = solve_cycle_plan(plan)
    altered = replace(plan, initial_audits=(replace(plan.initial_audits[0], findings=()), *plan.initial_audits[1:]))
    with pytest.raises(ValueError, match='solve-plan-binding'):
        measure_cycle_construction(altered, prepared, replace(result, plan_digest=altered.digest))


def test_real_generated_context_pairs_and_fixed_metrics_are_included():
    plan, _, _ = case(('CB', 'CG', 'CD'))
    problem = plan.initializer.problem
    branch = plan.branches[0]
    values = residual_values(plan, 0, tuple(0. for _ in range(9)))
    assert len(values) == len(problem.coverage.metrics)+problem.coverage.constraint_scalar_count+len(branch.pairs)
    assert any(t.role == 'fixed-only' for t in problem.coverage.metrics)
    assert branch.pairs == tuple(sorted(branch.pairs))


@pytest.mark.parametrize('kind', ['exclusion-upper', 'nanometer-distance'])
def test_constraint_residuals_have_frozen_semantics_not_a_dropped_upper_bound(kind):
    plan, prepared, _ = case()
    problem = plan.initializer.problem
    target = next(iter(problem.source.targets))
    pair = tuple(ExactAtomRef(target, n) for n in ('N', 'CD'))
    q = (Constraint(ConstraintKind.EXCLUSION, pair, Unit.ANGSTROM, .1, .2) if kind == 'exclusion-upper'
         else Constraint(ConstraintKind.DISTANCE, pair, Unit.NANOMETER, .14, .16))
    other = compile_cycle_problem(problem.source, problem.anchors, extra_constraints=(q,))
    from dvbfixer.chemical_reconstruction.cycle_source import bind_cycle_source
    newer = compile_cycle_solve_plan(initialize_cycle_problem(other), prepared, bind_cycle_source(other, prepared))
    from dvbfixer.chemical_reconstruction.cycle_residuals import vector_points
    points = vector_points(newer, 0, (0., 0., 0.))
    distance = math.dist(*(points[a] for a in pair))
    expected = max(0., .1-distance)/.1+max(0., distance-.2)/.2 if kind == 'exclusion-upper' else (distance/10-.15)/.01
    offset = len(other.coverage.metrics)+other.coverage.constraint_scalar_count-1
    assert residual_values(newer, 0, (0., 0., 0.))[offset] == pytest.approx(expected)


def test_budget_breach_does_not_become_a_partial_proposal_or_selection(monkeypatch):
    plan, prepared, _ = case()
    def exhausted(function, variables):
        tracker = function.__self__
        tracker.optimizer_calls = 50000  # Synthetic boundary injection, not a sampled trajectory.
        function(tuple(0. for _ in range(variables)))
        pytest.fail('ceiling failed')
    monkeypatch.setattr(cycle_solver, '_solve', exhausted)
    result = solve_cycle_plan(plan)
    assert all(b.status == 'numerical-failed' and b.optimizer_residual_calls == 50001 and not b.raw for b in result.branches)
    assert measure_cycle_construction(plan, prepared, result).selected_index is None


def test_preexisting_context_D_cannot_be_hidden_by_numeric_convergence(monkeypatch):
    base, _, witness = case()
    problem = base.initializer.problem
    target = next(iter(problem.source.targets))
    atoms = {o.identity.name: o for o in (*problem.source.observations, *witness) if o.identity.component == target}
    context = ComponentInstanceRef(1, 'z', '999')
    extra = tuple(Observation(ExactAtomRef(context, n), atoms[n].element,
        (50.-atoms[n].position[0], 50.+atoms[n].position[1], 50.+atoms[n].position[2]), Unit.ANGSTROM)
        for n in ('N', 'CA', 'C', 'CB'))
    problem, prepared, metadata, witness = source_fixture(extra=extra)
    plan = compile_cycle_solve_plan(initialize_cycle_problem(problem), prepared, metadata)
    mock_witness_solver(plan, witness, monkeypatch)
    result = solve_cycle_plan(plan)
    report = measure_cycle_construction(plan, prepared, result)
    assert all(b.stats.success for b in result.branches)
    assert all(a.pair.raw.global_l_status == 'rejected' for a in report.branches)
    assert report.status == 'rejected' and report.selected_index is None


def test_snapshot_numeric_representation_and_solver_success_flags_are_strict(monkeypatch):
    plan, prepared, witness = case()
    mock_witness_solver(plan, witness, monkeypatch)
    result = solve_cycle_plan(plan)
    with pytest.raises(ValueError):
        replace(result.branches[0].stats, success=1)
    with pytest.raises(ValueError):
        replace(result.branches[0], returned_x=tuple(0 for _ in range(3)))
    with pytest.raises(ValueError):
        replace(plan, settings=replace(plan.settings, x_scale=1))
