"""Synthetic initializer controls only: no tracked-source R1 solver trials."""
from __future__ import annotations

import inspect
import json
import math
from dataclasses import replace

import numpy as np
import pytest
from scripts.chemical_reconstruction_ccd_controls import _json

from dvbfixer.chemical_reconstruction import cycle_initialization as initialization
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Observation,
    OpaqueContextPolicy,
    Origin,
    PreparedContextRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.cycle_contract import compile_cycle_problem
from dvbfixer.chemical_reconstruction.cycle_initialization import (
    CycleInitializerRecord,
    audit_cycle_initializer,
    equation_residual,
    initialize_cycle_problem,
)
from dvbfixer.chemical_reconstruction.proposals import GeometryPlan, Strategy
from dvbfixer.chemical_reconstruction.research_decode import decode_research_record
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from tests.test_chemical_reconstruction_cycle_contract import mapped_graph, source_case


def problem_case(name='PRO', mask=None, state=None):
    source, anchors = source_case(name, mask, state)
    return compile_cycle_problem(source, anchors)


@pytest.mark.parametrize('mask', [('CD',), ('CG',), ('CB',), ('CG', 'CD'), ('CB', 'CD'), ('CB', 'CG'), ('CB', 'CG', 'CD')])
def test_eight_starts_for_every_nonempty_subset_without_a_parent_gate(mask):
    problem = problem_case(mask=mask)
    record = initialize_cycle_problem(problem)
    assert record.status == 'initialized-unvalidated'
    assert record.problem is problem
    assert record.system.variables == problem.coverage.variables
    assert record.system.anchors == problem.anchors
    assert len(record.starts) == 8
    assert record.maximum_equation_residual <= 1e-10
    assert record.coordinates_produced
    assert not record.geometry_approved and not record.solver_ready and not record.md_ready
    assert audit_cycle_initializer(record) == ()
    for index, start in enumerate(record.starts):
        assert start.index == index and start.sign == (1 if index < 4 else -1) and start.quarter_turn == index % 4
        assert start.status == 'coordinates-initialized'
        assert tuple(o.identity for o in start.additions) == problem.coverage.variables
        assert {o.identity for o in start.additions} == problem.source.ownership.generated
        assert start.problem_digest == problem.digest
        assert not hasattr(start, 'strategy') and not isinstance(start, GeometryPlan)
    if set(mask) == {'CB', 'CD'}:
        assert len(problem.coverage.blocks) == 2 and len(record.system.matrix) == 2


@pytest.mark.parametrize('name,state', [('PHE', 'PHE'), ('TYR', 'TYR'), ('TRP', 'TRP'), ('HIS', 'HIP')])
def test_complete_aromatic_fused_and_branched_scaffolds_are_not_PRO_patches(name, state):
    problem = problem_case(name, state=state)
    record = initialize_cycle_problem(problem)
    assert record.status == 'initialized-unvalidated'
    assert all(len(start.additions) == len(problem.coverage.variables) for start in record.starts)
    assert record.system.length_scale_angstrom > 0
    if name == 'TYR':
        assert all('OH' in {o.identity.name for o in start.additions} for start in record.starts)
    if name == 'TRP':
        assert len(record.system.matrix) == 10
    assert not record.geometry_approved


def test_one_variable_harmonic_point_is_the_weighted_observed_boundary_average():
    problem = problem_case(mask=('CD',))
    record = initialize_cycle_problem(problem)
    variable = problem.coverage.variables[0]
    neighbours = []
    observed = {o.identity: o.angstrom for o in problem.source.observations}
    for metric in problem.coverage.metrics:
        if metric.kind == 'bond' and variable in metric.atoms:
            other = next(a for a in metric.atoms if a != variable)
            neighbours.append((observed[other], 1/metric.target**2))
    assert len(neighbours) == 2  # Both graph boundaries, not a selected spanning edge.
    expected = tuple(sum(p[k]*weight for p, weight in neighbours)/sum(w for _, w in neighbours) for k in range(3))
    actual = record.system.frame.source(record.harmonic[0])
    assert actual == pytest.approx(expected, abs=1e-12)


def test_orbits_match_frozen_independent_formula_without_deduplication():
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    system = record.system
    h = record.harmonic[0]
    # n=1: theta=pi, v=(-1,0,1/2). Check four signed rotations independently.
    for index, (sign, k) in enumerate((s, k) for s in (1, -1) for k in range(4)):
        phi = k*math.pi/2
        v = (-1., -.5*math.sin(phi), .5*math.cos(phi))
        expected = system.frame.source(tuple(h[j]+sign*system.length_scale_angstrom/2*v[j] for j in range(3)))
        assert record.starts[index].additions[0].position == pytest.approx(expected, abs=1e-12)
    assert len(record.starts) == 8


def test_frame_is_right_handed_and_preserves_exact_identity_order():
    problem = problem_case(mask=('CB', 'CD'))
    record = initialize_cycle_problem(problem)
    frame = record.system.frame
    assert np.linalg.det(np.array(frame.axes)) == pytest.approx(1., abs=1e-12)
    for observation in problem.source.observations:
        assert frame.source(frame.local(observation.angstrom)) == pytest.approx(observation.angstrom, abs=1e-12)
    assert all(a.component.chain_id == 'a' and a.component.insertion_code == 'B' and a.component.occurrence == 2
               for a in record.system.variables)
    assert any(c.identity.chain_id == 'A' for c in record.problem.source.components)


def test_full_context_changes_provenance_not_the_target_boundary_solve():
    problem = problem_case(mask=('CD',))
    old = initialize_cycle_problem(problem)
    source = problem.source
    observations = tuple(replace(o, position=(70., 70., 70.)) if o.identity.component not in source.targets else o
                         for o in source.observations)
    newer = compile_cycle_problem(replace(source, observations=observations), problem.anchors)
    new = initialize_cycle_problem(newer)
    assert old.problem.digest != new.problem.digest and old.digest != new.digest
    assert old.system.observations_digest != new.system.observations_digest
    assert old.system.matrix == new.system.matrix and old.system.rhs == new.system.rhs
    assert tuple(s.additions for s in old.starts) == tuple(s.additions for s in new.starts)
    assert new.problem.source.observations is observations


def test_fractional_opaque_context_proof_and_all_observations_remain_exact():
    problem = problem_case(mask=('CD',))
    source = problem.source
    outside = ExactAtomRef(ComponentInstanceRef(1, 'q', '999'), 'Q')
    observations = (*source.observations, Observation(outside, 'C', (50., 50., 50.), Unit.ANGSTROM, .5))
    origin = Origin('prepared-observed-peptide-v1', 'b'*64)
    policy = OpaqueContextPolicy(source.source_sha256, 'b'*64, frozenset({outside}), digest_record(observations))
    components = (*source.components, Component(outside.component, ChemistryKind.UNKNOWN))
    prepared = PreparedContextRequest(source.source_sha256, source.operation, components, observations, source.targets,
        source.ownership, source.connections, source.constraints, source.ph, origin, context_policy=policy)
    record = initialize_cycle_problem(compile_cycle_problem(prepared, problem.anchors))
    assert record.problem.source is prepared and prepared.context_policy is policy
    assert prepared.observations is observations
    assert all(o.identity != outside for start in record.starts for o in start.additions)


def test_ideal_coordinate_translation_cannot_be_used_as_a_fitted_start():
    problem = problem_case(mask=('CD',))
    state = next(c.resolved for c in problem.source.components if c.identity in problem.source.targets)
    shifted_graph = replace(state.graph, atoms=tuple(replace(a, ideal_position=tuple(x+20. for x in a.ideal_position))
        for a in state.graph.atoms))
    shifted = compile_cycle_problem(mapped_graph(problem.source, shifted_graph), problem.anchors)
    old, new = initialize_cycle_problem(problem), initialize_cycle_problem(shifted)
    assert old.system.graph_digest != new.system.graph_digest
    assert old.harmonic == new.harmonic
    assert tuple(s.additions for s in old.starts) == tuple(s.additions for s in new.starts)


def test_initializer_accepts_only_the_problem_not_seed_reference_or_external_coordinates():
    assert list(inspect.signature(initialize_cycle_problem).parameters) == ['problem']
    problem = problem_case(mask=('CD',))
    for key in ('seed', 'native', 'reference', 'candidate', 'baseline', 'initial_positions'):
        with pytest.raises(TypeError):
            initialize_cycle_problem(problem, **{key: ()})
    with pytest.raises(ValueError, match='typed-problem'):
        initialize_cycle_problem(problem.source)


def test_no_optimizer_pseudoinverse_template_or_internal_constructor_is_called(monkeypatch):
    from scipy import optimize

    from dvbfixer.chemical_reconstruction import authority_geometry, construction, metric_projection

    problem = problem_case('TRP')
    def forbidden(*args, **kwargs):
        pytest.fail('initializer used an undeclared optimizer/pose/constructor')
    for module, name in ((optimize, 'least_squares'), (np.linalg, 'pinv'), (np.linalg, 'lstsq'),
                         (construction, 'propose'), (construction, '_internal'), (construction, '_template'),
                         (metric_projection, '_solve'), (authority_geometry, 'geometry_view')):
        monkeypatch.setattr(module, name, forbidden)
    assert initialize_cycle_problem(problem).status == 'initialized-unvalidated'


def test_saved_readback_and_strict_transport_do_not_solve_or_construct(monkeypatch):
    record = initialize_cycle_problem(problem_case(mask=('CB', 'CD')))
    payload = json.dumps(_json(record)).encode()
    def forbidden(*args, **kwargs):
        pytest.fail('saved initializer audit constructed coordinates or solved')
    for name in ('_direct_solve', '_make_start', 'initialize_cycle_problem'):
        monkeypatch.setattr(initialization, name, forbidden)
    monkeypatch.setattr(np.linalg, 'solve', forbidden)
    restored = decode_research_record(payload, CycleInitializerRecord)
    assert restored.digest == record.digest and audit_cycle_initializer(restored) == ()
    with pytest.raises(ValueError):
        decode_research_record(payload, GeometryPlan)
    assert {s.value for s in Strategy} == {'graph-internal-authority-v1', 'class-a-rigid-authority-v1',
        'graph-fragment-torsion-grid-v1', 'generated-heavy-metric-projection-v1'}


@pytest.mark.parametrize('mutation', ['system', 'rhs', 'frame', 'scale', 'harmonic', 'equation', 'source', 'order',
                                    'addition', 'element', 'occupancy', 'missing', 'algorithm', 'status'])
def test_cached_system_start_and_status_tampering_is_rejected(mutation):
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    with pytest.raises(ValueError):
        if mutation == 'system':
            replace(record, system=replace(record.system, problem_digest='0'*64))
        elif mutation == 'rhs':
            replace(record, system=replace(record.system, rhs=((0., 0., 0.),)))
        elif mutation == 'frame':
            replace(record, system=replace(record.system, frame=replace(record.system.frame, axes=tuple(reversed(record.system.frame.axes)))))
        elif mutation == 'scale':
            replace(record, system=replace(record.system, length_scale_angstrom=1.))
        elif mutation == 'harmonic':
            replace(record, harmonic=((0., 0., 0.),))
        elif mutation == 'equation':
            replace(record, maximum_equation_residual=1.)
        elif mutation == 'source':
            replace(record, problem=replace(record.problem, source=replace(record.problem.source, source_sha256='0'*64)))
        elif mutation == 'order':
            replace(record, starts=tuple(reversed(record.starts)))
        elif mutation in {'addition', 'element', 'occupancy'}:
            start = record.starts[0]
            addition = start.additions[0]
            change = {'position': (0., 0., 0.)} if mutation == 'addition' else {'element': 'O'} if mutation == 'element' else {'occupancy': 1}
            replace(record, starts=(replace(start, additions=(replace(addition, **change),)), *record.starts[1:]))
        elif mutation == 'missing':
            replace(record, starts=record.starts[:-1])
        elif mutation == 'algorithm':
            replace(record, algorithm='other')
        else:
            replace(record, status='refused')


def test_independent_deterministic_initializer_replay_is_byte_identical_not_a_new_trial():
    problem = problem_case('TRP')
    first = initialize_cycle_problem(problem)
    second = initialize_cycle_problem(problem)
    assert json.dumps(_json(first), sort_keys=True) == json.dumps(_json(second), sort_keys=True)


def test_singular_system_keeps_refusal_and_never_tries_another_algorithm(monkeypatch):
    problem = problem_case(mask=('CD',))
    calls = []
    def singular(*args, **kwargs):
        calls.append(1)
        raise np.linalg.LinAlgError('synthetic singular')
    monkeypatch.setattr(np.linalg, 'solve', singular)
    record = initialize_cycle_problem(problem)
    assert calls == [1] and record.status == 'refused' and record.reason == 'cycle-initializer-singular-system'
    assert record.system is not None and not record.starts and not record.harmonic
    assert not record.coordinates_produced and audit_cycle_initializer(record) == ()


def test_numpy_missing_is_environmental_unavailable_not_chemistry_rejection(monkeypatch):
    def missing(*args):
        raise ModuleNotFoundError('no numpy', name='numpy')
    monkeypatch.setattr(initialization, '_direct_solve', missing)
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    assert record.status == 'unavailable' and not record.starts and record.system is not None
    assert not record.geometry_approved and audit_cycle_initializer(record) == ()


@pytest.mark.parametrize('kind', ['dependency', 'programming'])
def test_unexpected_errors_propagate_without_becoming_scientific_refusals(kind, monkeypatch):
    def error(*args):
        if kind == 'dependency':
            raise ModuleNotFoundError('other missing module', name='unrelated_module')
        raise RuntimeError('synthetic bug')
    monkeypatch.setattr(initialization, '_direct_solve', error)
    with pytest.raises(ModuleNotFoundError if kind == 'dependency' else RuntimeError):
        initialize_cycle_problem(problem_case(mask=('CD',)))


def test_large_equation_residual_retains_finite_failed_scaffold_but_no_start(monkeypatch):
    monkeypatch.setattr(initialization, '_direct_solve', lambda system: ((0., 0., 0.),))
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    assert record.status == 'refused' and record.reason == 'cycle-initializer-equation-residual'
    assert record.harmonic == ((0., 0., 0.),) and record.maximum_equation_residual > 1e-10
    assert not record.starts and audit_cycle_initializer(record) == ()


def test_nonfinite_linear_solution_is_retained_as_failure_not_zero_filled(monkeypatch):
    monkeypatch.setattr(np.linalg, 'solve', lambda *args: np.full((1, 3), float('nan')))
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    assert record.status == 'refused' and record.reason == 'cycle-initializer-nonfinite-solution'
    assert not record.starts and not record.harmonic


def test_nonfinite_observed_frame_refuses_before_a_linear_solve(monkeypatch):
    problem = problem_case(mask=('CD',))
    monkeypatch.setattr(initialization, '_unit', lambda point: (float('nan'), 0., 0.))
    monkeypatch.setattr(np.linalg, 'solve', lambda *args: pytest.fail('invalid frame reached solve'))
    record = initialize_cycle_problem(problem)
    assert record.status == 'refused' and record.system is None and not record.starts


def test_failed_orbit_is_labelled_without_losing_other_starts_or_substituting(monkeypatch):
    original = initialization.offset_position
    def one_failed(system, harmonic, atom_index, branch_index):
        if branch_index == 3:
            return (100001., 0., 0.)
        return original(system, harmonic, atom_index, branch_index)
    monkeypatch.setattr(initialization, 'offset_position', one_failed)
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    assert len(record.starts) == 8 and record.starts[3].status == 'failed' and not record.starts[3].additions
    assert all(s.status == 'coordinates-initialized' for i, s in enumerate(record.starts) if i != 3)
    assert audit_cycle_initializer(record) == () and not record.geometry_approved


def test_equation_readback_refuses_nonfinite_or_wrong_shape():
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    with pytest.raises(ValueError):
        equation_residual(record.system, ())
    with pytest.raises(ValueError):
        equation_residual(record.system, ((float('nan'), 0., 0.),))


def test_strict_transport_rejects_undeclared_native_score_or_solver_fields():
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    payload = _json(record)
    for key in ('native_reference', 'rmsd', 'solver_success'):
        with pytest.raises(ValueError):
            decode_research_record(json.dumps(dict(payload, **{key: 0})).encode(), CycleInitializerRecord)


def test_integer_float_and_bool_numeric_equality_cannot_bypass_bindings():
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    start = record.starts[0]
    with pytest.raises(ValueError):
        replace(start, index=True)
    with pytest.raises(ValueError):
        replace(start, sign=True)
    # Independently verify strict occupancy representation (1 vs 1.0).
    addition = replace(start.additions[0], occupancy=1)
    assert addition == start.additions[0]
    with pytest.raises(ValueError):
        replace(record, starts=(replace(start, additions=(addition,)), *record.starts[1:]))


def test_all_failed_orbits_do_not_claim_coordinates_or_solver_readiness(monkeypatch):
    monkeypatch.setattr(initialization, 'offset_position', lambda *args: (100001., 0., 0.))
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    assert len(record.starts) == 8 and all(s.status == 'failed' and not s.additions for s in record.starts)
    assert not record.coordinates_produced and not record.solver_ready and not record.source_coverage_verified


def test_local_linear_algebra_values_do_not_inherit_source_coordinate_bounds():
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    # RHS has inverse-length units, not source-coordinate units. Its shape DTO must
    # allow finite values >1e5; the parent record still binds it to the actual source.
    system = replace(record.system, rhs=((200000., 0., 0.),))
    assert system.rhs[0][0] == 200000.
    with pytest.raises(ValueError, match='system-source-binding'):
        replace(record, system=system)


def test_failure_status_cannot_hide_successful_additions():
    record = initialize_cycle_problem(problem_case(mask=('CD',)))
    with pytest.raises(ValueError):
        replace(record, status='unavailable', reason='cycle-initializer-numpy-unavailable')
    with pytest.raises(ValueError):
        replace(record.starts[0], status='failed', reason='synthetic refusal')
