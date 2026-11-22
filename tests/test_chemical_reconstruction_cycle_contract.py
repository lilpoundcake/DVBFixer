"""Planning complete cyclic subsets is not numerical construction qualification."""
from __future__ import annotations

import json
from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_all20_controls import (
    load_parent,
    sidechain_heavy_mask,
    synthetic_input,
)
from scripts.chemical_reconstruction_ccd_controls import _json

from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Constraint,
    ConstraintKind,
    Observation,
    OpaqueContextPolicy,
    Operation,
    Origin,
    Ownership,
    PreparedContextRequest,
    ReconstructionRequest,
    Unit,
    canonical_graph_digest,
    digest_record,
)
from dvbfixer.chemical_reconstruction.cycle_contract import (
    CycleConstructionProblem,
    CyclePlanningAdmission,
    compile_cycle_problem,
    constraint_request,
    plan_cycle_construction,
)
from dvbfixer.chemical_reconstruction.proposals import GeometryPlan
from dvbfixer.chemical_reconstruction.research_decode import decode_research_record
from dvbfixer.domain.nonprotein_reconstruction import (
    AuthorityMode,
    AuthorityProvenance,
    ChemicalAtom,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def source_case(name='PRO', mask=None, state=None):
    full, authority, _ = synthetic_input(name, state or name, {})
    wanted = sidechain_heavy_mask(authority) if mask is None else frozenset(mask)
    generated = frozenset(ExactAtomRef(authority.identity, n) for n in wanted)
    source = replace(full, observations=tuple(o for o in full.observations if o.identity not in generated), ownership=Ownership(generated))
    anchors = tuple(ExactAtomRef(authority.identity, n) for n in ('CA', 'N', 'C'))
    return source, anchors


@pytest.mark.parametrize('mask', [('CD',), ('CG',), ('CB',), ('CG', 'CD'), ('CB', 'CD'), ('CB', 'CG'), ('CB', 'CG', 'CD')])
def test_all_nonempty_PRO_subsets_are_complete_problems_not_one_atom_patches(mask):
    source, anchors = source_case(mask=mask)
    problem = compile_cycle_problem(source, anchors)
    assert {a.name for a in problem.coverage.variables} == set(mask)
    assert len(problem.coverage.cycle_edges) == 5
    assert {a for block in problem.coverage.blocks for a in block.atoms} == source.ownership.generated
    bonds = {frozenset(m.atoms) for m in problem.coverage.metrics if m.kind == 'bond'}
    assert all(frozenset(edge) in bonds for edge in problem.coverage.cycle_edges)
    assert all(set(edge) & source.ownership.generated for edge in problem.coverage.generated_cycle_edges)
    if set(mask) == {'CB', 'CD'}:
        assert len(problem.coverage.blocks) == 2  # Fixed CG separates islands; do not solve independently.
    assert not isinstance(problem, GeometryPlan) and not hasattr(problem, 'strategy')
    assert not problem.coordinates_produced and not problem.geometry_approved and not problem.source_coverage_verified
    assert not problem.solver_defined and not problem.md_ready


@pytest.mark.parametrize('name,state,cycle_count', [('PHE', 'PHE', 6), ('TYR', 'TYR', 6), ('TRP', 'TRP', 10), ('HIS', 'HIP', 5)])
def test_aromatic_fused_and_ring_branch_metadata_without_chemical_qualification(name, state, cycle_count):
    source, anchors = source_case(name=name, state=state)
    problem = compile_cycle_problem(source, anchors)
    assert len(problem.coverage.cycle_edges) == cycle_count
    assert set(problem.coverage.variables) == source.ownership.generated
    if name == 'TYR':
        assert any(a.name == 'OH' for a in problem.coverage.variables)
        assert all(a.name != 'OH' for e in problem.coverage.cycle_edges for a in e)
    assert not problem.geometry_approved


def test_isolated_resolved_free_molecule_in_retained_opaque_context():
    authority = load_parent('PRO')
    identity = ComponentInstanceRef(1, 'x', '12', 'A')
    mask = frozenset(ExactAtomRef(identity, n) for n in ('CB', 'CG', 'CD'))
    observations = tuple(Observation(ExactAtomRef(identity, a.name), a.element, a.ideal_position, Unit.ANGSTROM)
                         for a in authority.state.graph.atoms if a.element not in {'H', 'D'} and ExactAtomRef(identity, a.name) not in mask)
    outside = ComponentInstanceRef(1, 'X', '13')
    observations += (Observation(ExactAtomRef(outside, 'Z1'), 'C', (50., 50., 50.), Unit.ANGSTROM),)
    source = ReconstructionRequest('a'*64, Operation.HEAVY,
        (Component(identity, ChemistryKind.ISOLATED, (authority.state,)), Component(outside, ChemistryKind.UNKNOWN)),
        observations, frozenset({identity}), Ownership(mask), constraints=authority.constraints(identity, include_hydrogens=False,
        minimum_volume=.05, planar_tolerance_angstrom=.05))
    problem = compile_cycle_problem(source, tuple(ExactAtomRef(identity, n) for n in ('CA', 'N', 'C')))
    assert problem.coverage.unassessed_context == (outside,) and not problem.source_coverage_verified


def test_source_and_new_provenance_are_exact_and_separate_from_projection():
    source, anchors = source_case(mask=('CD',))
    problem = compile_cycle_problem(source, anchors)
    request = constraint_request(source, problem.extra_constraints)
    assert source is problem.source and request.observations is source.observations
    assert request.components is source.components and request.connections is source.connections
    assert request.ownership is source.ownership and request.ph == source.ph
    assert request.constraints[:len(source.constraints)] == source.constraints
    assert request.digest == problem.coverage.request_digest != source.digest
    assert dict(request.origin.settings)['cycle_construction_source'] == source.digest
    assert not any(k.startswith('metric_projection_') for k, _ in request.origin.settings)
    with pytest.raises(ValueError, match='nested'):
        compile_cycle_problem(request, anchors)


def test_prepared_fractional_context_proof_is_retained_not_owned_or_relabelled():
    source, anchors = source_case(mask=('CD',))
    outside = ComponentInstanceRef(1, 'q', '999')
    observation = Observation(ExactAtomRef(outside, 'Q'), 'C', (50., 50., 50.), Unit.ANGSTROM, .5)
    observations = (*source.observations, observation)
    prepared = PreparedContextRequest(source.source_sha256, source.operation,
        (*source.components, Component(outside, ChemistryKind.UNKNOWN)), observations, source.targets,
        source.ownership, source.connections, source.constraints, source.ph, Origin('prepared-observed-peptide-v1', 'b'*64),
        context_policy=OpaqueContextPolicy(source.source_sha256, 'b'*64, frozenset({observation.identity}), digest_record(observations)))
    problem = compile_cycle_problem(prepared, anchors)
    effective = constraint_request(prepared, ())
    assert isinstance(effective, PreparedContextRequest) and effective.context_policy is prepared.context_policy
    assert effective.observations is prepared.observations and problem.source is prepared
    assert observation.identity not in problem.coverage.variables
    with pytest.raises(ValueError):
        compile_cycle_problem(replace(source, observations=observations,
            components=(*source.components, Component(outside, ChemistryKind.UNKNOWN))), anchors)


def mapped_graph(source, graph):
    target = next(iter(source.targets))
    component = next(c for c in source.components if c.identity == target)
    state = replace(component.resolved, graph=graph, authority=AuthorityProvenance(AuthorityMode.USER_MAPPED,
        'synthetic-cycle-contract-test', canonical_graph_digest(graph), None))
    return replace(source, components=tuple(replace(c, states=(state,), selected_state=None) if c.identity == target else c for c in source.components))


def test_unanchored_generated_ring_fragment_refuses_without_absolute_ideal_pose():
    source, anchors = source_case(name='PHE', mask=('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'))
    state = next(c.resolved for c in source.components if c.identity in source.targets)
    graph = replace(state.graph, bonds=tuple(b for b in state.graph.bonds if {b.atom1, b.atom2} != {'CB', 'CG'}), component_count=2)
    admission = plan_cycle_construction(mapped_graph(source, graph), anchors)
    assert admission.status == 'refused' and 'unanchored' in admission.reason and admission.problem is None


def test_element_token_case_uses_existing_admission_semantics_not_chain_normalization():
    source, anchors = source_case(mask=('CD',))
    state = next(c.resolved for c in source.components if c.identity in source.targets)
    graph = replace(state.graph, atoms=tuple(replace(a, element=a.element.lower()) for a in state.graph.atoms))
    problem = compile_cycle_problem(mapped_graph(source, graph), anchors)
    assert len(problem.coverage.heavy_nodes) == 7  # Lowercase h is not a heavy variable.
    assert all(a.component.chain_id == 'a' for a in problem.coverage.heavy_nodes)
    assert any(c.identity.chain_id == 'A' for c in problem.source.components)


def test_atom_and_component_names_are_not_a_PRO_or_chain_case_rule():
    source, anchors = source_case(mask=('CD',))
    target = next(iter(source.targets))
    state = next(c.resolved for c in source.components if c.identity == target)
    names = {a.name: f'Z{i}' for i, a in enumerate(state.graph.atoms)}
    def ref(a):
        return replace(a, name=names[a.name]) if a.component == target else a
    graph = replace(state.graph, component_id='SYNTHETIC', state_id='synthetic-cycle',
        atoms=tuple(replace(a, name=names[a.name]) for a in state.graph.atoms),
        bonds=tuple(replace(b, atom1=names[b.atom1], atom2=names[b.atom2]) for b in state.graph.bonds))
    renamed = mapped_graph(source, graph)
    renamed = replace(renamed, observations=tuple(replace(o, identity=ref(o.identity)) for o in renamed.observations),
        ownership=Ownership(frozenset(ref(a) for a in source.ownership.generated)),
        constraints=tuple(replace(q, atoms=tuple(ref(a) for a in q.atoms)) for q in source.constraints),
        connections=tuple(replace(c, link=replace(c.link, component_atom=ref(c.link.component_atom), other_atom=ref(c.link.other_atom))) for c in source.connections))
    problem = compile_cycle_problem(renamed, tuple(ref(a) for a in anchors))
    assert len(problem.coverage.cycle_edges) == 5
    assert all(a.component.chain_id == 'a' and a.component.insertion_code == 'B' and a.component.occurrence == 2 for a in problem.coverage.heavy_nodes)
    assert any(c.identity.chain_id == 'A' for c in problem.source.components)


@pytest.mark.parametrize('name', ['LYS', 'THR', 'ALA'])
def test_acyclic_work_is_not_silently_admitted_to_the_cycle_route(name):
    source, anchors = source_case(name=name)
    result = plan_cycle_construction(source, anchors)
    assert result.status == 'refused' and result.reason == 'cycle-no-generated-involving-cycle'


@pytest.mark.parametrize('mutation', ['empty', 'mobile', 'temporary', 'wrong-operation', 'missing-heavy', 'target-H', 'extra-target', 'huge-mask'])
def test_unsupported_inventory_or_ownership_refuses(mutation):
    source, anchors = source_case(mask=('CD',))
    target = next(iter(source.targets))
    if mutation == 'empty':
        source = replace(source, ownership=Ownership())
    elif mutation in {'mobile', 'temporary'}:
        fixed = next(o.identity for o in source.observations if o.identity.component == target)
        source = replace(source, ownership=Ownership(source.ownership.generated, **{mutation: frozenset({fixed if mutation == 'mobile' else replace(fixed, name='TEMP')})}))
    elif mutation == 'wrong-operation':
        source = replace(source, operation=Operation.HYDROGENS)
    elif mutation == 'missing-heavy':
        source = replace(source, observations=tuple(o for o in source.observations if o.identity.name != 'CG'))
    elif mutation == 'target-H':
        source = replace(source, observations=(*source.observations, Observation(ExactAtomRef(target, 'HA'), 'H', (0., 0., 0.), Unit.ANGSTROM)))
    elif mutation == 'extra-target':
        source = replace(source, targets=frozenset(c.identity for c in source.components))
    else:
        source = replace(source, ownership=Ownership(frozenset(ExactAtomRef(target, f'V{i}') for i in range(17))))
    result = plan_cycle_construction(source, anchors)
    assert result.status == 'refused' and result.problem is None and not result.coordinates_produced


@pytest.mark.parametrize('mutation', ['duplicate', 'outside', 'generated', 'collinear'])
def test_observed_frame_is_explicit_and_never_replaced(mutation):
    source, anchors = source_case(mask=('CD',))
    if mutation == 'duplicate':
        anchors = (anchors[0], anchors[0], anchors[2])
    elif mutation == 'outside':
        anchors = (*anchors[:2], next(o.identity for o in source.observations if o.identity.component not in source.targets))
    elif mutation == 'generated':
        anchors = (*anchors[:2], next(iter(source.ownership.generated)))
    else:
        # Duplicate fixed positions also cause an explicit source-constraint refusal,
        # never a search for replacement anchors.
        source = replace(source, observations=tuple(replace(o, position=(0., 0., 0.)) if o.identity in anchors else o for o in source.observations))
    assert plan_cycle_construction(source, anchors).status == 'refused'


def test_nonfinite_frame_arithmetic_refuses_instead_of_defining_a_solver_frame(monkeypatch):
    from dvbfixer.chemical_reconstruction import cycle_contract

    source, anchors = source_case(mask=('CD',))
    monkeypatch.setattr(cycle_contract, '_unit', lambda *a: (float('nan'), 0., 0.))
    result = plan_cycle_construction(source, anchors)
    assert result.status == 'refused' and result.reason == 'cycle-degenerate-observed-frame'


@pytest.mark.parametrize('mutation', ['torsion', 'zero-width', 'fixed-conflict', 'projection-origin'])
def test_constraint_semantics_and_provenance_are_not_silently_ignored(mutation):
    source, anchors = source_case(mask=('CD',))
    missing = next(iter(source.ownership.generated))
    if mutation == 'projection-origin':
        source = replace(source, origin=replace(source.origin, settings=(('metric_projection_parent', source.digest),)))
        extra = ()
    elif mutation == 'torsion':
        extra = (Constraint(ConstraintKind.TORSION, (*anchors, missing), Unit.DEGREE, lower=-180., upper=180.),)
    elif mutation == 'zero-width':
        extra = (Constraint(ConstraintKind.DISTANCE, (anchors[0], missing), Unit.ANGSTROM, lower=1., upper=1.),)
    else:
        extra = (Constraint(ConstraintKind.DISTANCE, anchors[:2], Unit.ANGSTROM, lower=100., upper=101.),)
    assert plan_cycle_construction(source, anchors, extra_constraints=extra).status == 'refused'


def test_generated_external_link_is_not_inferred_or_dropped():
    source, anchors = source_case(mask=('CD',))
    connection = source.connections[0]
    changed = replace(connection, link=replace(connection.link, other_atom=next(iter(source.ownership.generated))))
    source = replace(source, connections=(changed,))
    result = plan_cycle_construction(source, anchors)
    assert result.status == 'refused' and 'generated-external' in result.reason
    assert source.connections == (changed,)


def test_fixed_authority_conflict_does_not_request_observation_motion():
    source, anchors = source_case(mask=('CD',))
    source = replace(source, observations=tuple(replace(o, position=(o.position[0]+10., *o.position[1:]))
        if o.identity.component in source.targets and o.identity.name == 'CG' else o for o in source.observations))
    result = plan_cycle_construction(source, anchors)
    assert result.status == 'refused' and result.reason == 'cycle-fixed-authority-bond-conflict'


def test_target_node_and_total_scalar_budgets_refuse_without_truncation():
    source, anchors = source_case(mask=('CD',))
    problem = compile_cycle_problem(source, anchors)
    with pytest.raises(ValueError, match='budget'):
        replace(problem.coverage, constraint_scalar_count=512)
    state = next(c.resolved for c in source.components if c.identity in source.targets)
    extras = tuple(ChemicalAtom(f'EX{i}', 'C', 0, (100.+i, 50., 50.)) for i in range(60))
    graph = replace(state.graph, atoms=(*state.graph.atoms, *extras), component_count=61)
    source = mapped_graph(source, graph)
    target = next(iter(source.targets))
    source = replace(source, observations=(*source.observations, *(Observation(ExactAtomRef(target, a.name), a.element,
        a.ideal_position, Unit.ANGSTROM) for a in extras)))
    result = plan_cycle_construction(source, anchors)
    assert result.status == 'refused' and result.reason == 'cycle-target-node-budget'


@pytest.mark.parametrize('mutation', ['variables', 'edges', 'cycles', 'metrics', 'boundary', 'roles', 'digest', 'policy', 'source'])
def test_self_consistent_transport_does_not_allow_coverage_tampering(mutation):
    source, anchors = source_case(mask=('CD',))
    problem = compile_cycle_problem(source, anchors)
    with pytest.raises(ValueError):
        if mutation == 'policy':
            replace(problem, policy=replace(problem.policy, bond_tolerance_angstrom=.2))
        elif mutation == 'source':
            replace(problem, source=replace(source, observations=tuple(reversed(source.observations))))
        else:
            coverage = problem.coverage
            if mutation == 'variables':
                coverage = replace(coverage, variables=())
            elif mutation in {'edges', 'cycles', 'metrics'}:
                key = {'edges': 'heavy_edges', 'cycles': 'cycle_edges', 'metrics': 'metrics'}[mutation]
                coverage = replace(coverage, **{key: getattr(coverage, key)[:-1]})
            elif mutation == 'boundary':
                block = coverage.blocks[0]
                coverage = replace(coverage, blocks=(replace(block, boundary_edges=block.boundary_edges[:-1]),))
            elif mutation == 'roles':
                coverage = replace(coverage, constraint_roles=())
            else:
                coverage = replace(coverage, graph_digest='0'*64)
            replace(problem, coverage=coverage)


def test_standalone_DTO_roundtrip_does_not_expand_old_geometry_plan_variants():
    source, anchors = source_case(mask=('CD',))
    result = plan_cycle_construction(source, anchors)
    restored = decode_research_record(json.dumps(_json(result)).encode(), CyclePlanningAdmission)
    assert restored == result and not restored.geometry_approved
    problem = decode_research_record(json.dumps(_json(result.problem)).encode(), CycleConstructionProblem)
    assert problem.digest == result.problem.digest
    with pytest.raises(ValueError):
        decode_research_record(json.dumps(_json(problem)).encode(), GeometryPlan)
    payload = _json(problem)
    payload['native_reference'] = []
    with pytest.raises(ValueError):
        decode_research_record(json.dumps(payload).encode(), CycleConstructionProblem)


def test_numeric_representation_changes_cannot_evade_exact_digest_bindings():
    source, anchors = source_case(mask=('CD',))
    problem = compile_cycle_problem(source, anchors)
    with pytest.raises(ValueError):
        replace(problem, policy=replace(problem.policy, angle_tolerance_degrees=12))
    metrics = tuple(replace(m, tolerance=12) if m.kind == 'angle' else m for m in problem.coverage.metrics)
    assert metrics == problem.coverage.metrics  # Python equality alone is insufficient.
    with pytest.raises(ValueError):
        replace(problem, coverage=replace(problem.coverage, metrics=metrics))
    from scripts.chemical_reconstruction_cycle_contract_controls import _canonical_json
    assert _canonical_json({'planned': True}) != _canonical_json({'planned': 1})
    assert _canonical_json({'count': 1}) != _canonical_json({'count': 1.0})


def test_unsupported_parent_gate_stays_closed_even_when_source_problem_is_planned(monkeypatch):
    from dvbfixer.chemical_reconstruction.construction import proposal_admission
    from dvbfixer.chemical_reconstruction.metric_projection import propose_projection
    from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan
    from dvbfixer.chemical_reconstruction.projection_contracts import (
        POLICY,
        append_projection_constraints,
    )
    from dvbfixer.chemical_reconstruction.proposals import ProposalResult

    source, anchors = source_case(mask=('CD',))
    extra = (Constraint(ConstraintKind.EXCLUSION, (anchors[0], next(iter(source.ownership.generated))), Unit.ANGSTROM, lower=.01, upper=100.),)
    planned = plan_cycle_construction(source, anchors, extra_constraints=extra)
    assert planned.status == 'planned-contract-only'
    plan = compile_peptide_request_plan(source, POLICY)
    rejected = ProposalResult('rejected', proposal_admission(source, plan), plan.digest, 'synthetic-negative-no-geometry-parent')
    old_request = append_projection_constraints(source, extra)
    result = propose_projection(source, old_request, plan, rejected, anchors)
    assert result.status == 'refused' and result.plan is None


def test_compilation_never_calls_initializer_solver_or_coordinate_constructor(monkeypatch):
    from scipy import optimize

    from dvbfixer.chemical_reconstruction import construction, metric_projection
    from dvbfixer.protein_completion import geometry

    def forbidden(*args, **kwargs):
        pytest.fail('cycle contract generated coordinates')
    for module, name in ((construction, 'propose'), (construction, '_internal'), (metric_projection, 'propose_projection'),
                         (geometry, 'extend'), (geometry, 'two_angle_tip'), (optimize, 'least_squares')):
        monkeypatch.setattr(module, name, forbidden)
    source, anchors = source_case(name='TRP')
    problem = compile_cycle_problem(source, anchors)
    assert len(problem.coverage.cycle_edges) == 10 and not problem.solver_defined


def test_missing_protocol_refuses_before_output_or_planning(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_cycle_contract_controls as controls

    monkeypatch.setattr(controls, 'PROTOCOL', tmp_path/'missing')
    with pytest.raises(FileNotFoundError):
        controls.run_controls(tmp_path/'output')
    assert not (tmp_path/'output').exists()


def test_source_refusals_preserve_nine_cell_denominator_without_initializer(tmp_path, monkeypatch):
    from scripts import chemical_reconstruction_cycle_contract_controls as controls

    def refuse(*args):
        raise ValueError('synthetic-source-refusal')
    monkeypatch.setattr(controls, 'prepare_cell', refuse)
    monkeypatch.setattr(controls, 'plan_cycle_construction', lambda *a, **k: pytest.fail('source failure compiled'))
    summary = controls.run_controls(tmp_path/'output')
    assert len(summary['controls']) == 9 and summary['new_coordinate_proposals'] == 0
    assert all(r['status'] == 'source-refused' for r in summary['controls'])
    assert controls.readback_controls(tmp_path/'output') == summary
