"""Explicit PRO/HIP declarations and source-only metrics; no state inference."""
from __future__ import annotations

from dataclasses import replace
from itertools import combinations

import pytest
from scripts.chemical_reconstruction_declared_state_controls import (
    CASES,
    FFBONDED,
    FIXED_POLICY,
    GENERATED_POLICY,
    compile_view,
    declaration_inputs,
    run_controls,
    synthetic_base,
)

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import ConstraintKind, canonical_graph_digest
from dvbfixer.chemical_reconstruction.declared_peptide import (
    compile_declared_hydrogens,
    derive_declared_peptide,
)
from dvbfixer.chemical_reconstruction.hydrogen_completion import audit_hydrogen_cip
from dvbfixer.chemical_reconstruction.joint_hydrogen_geometry import (
    JointMetricGeometryPlan,
    JointSourceMetricGeometryPlan,
)
from dvbfixer.chemical_reconstruction.metric_geometry import (
    GeneratedMetricGeometryPlan,
    SourceMetricGeometryPlan,
    metric_findings,
    metric_key,
)
from dvbfixer.chemical_reconstruction.peptide_state import derive_local_peptide
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils import variants


def case(state='PRO'):
    source, _, leaving = next(row for row in CASES if row[1] == state)
    base, authority, incoming = synthetic_base(source, state, leaving, {})
    view = compile_view(base, authority, incoming, {})
    return base, authority, incoming, view


@pytest.mark.parametrize('state,count,charge', [('PRO', 7, 0), ('HIP', 8, 1)])
def test_declared_states_keep_exact_unchanged_source_chemistry_and_complete_H(state, count, charge):
    base, authority, _, view = case(state)
    graph = authority.state.graph
    source_atoms = {a.name: a for a in authority.parent.state.graph.atoms}
    assert all(a == source_atoms[a.name] for a in graph.atoms)
    assert all(b in authority.parent.state.graph.bonds for b in graph.bonds)
    assert graph.net_charge == charge and authority.inventory.state_name == state
    assert authority.removed_atoms == ({'OXT', 'HXT', 'H'} if state == 'PRO' else {'OXT', 'HXT', 'H2'})
    assert view.request.observations is base.observations and view.request.components == base.components
    assert view.request.connections == base.connections
    assert dict(view.request.origin.settings)['peptide_state_declaration'] == authority.declaration.digest
    result = propose(view.request, view.plan)
    assert result.status == 'declared-geometry-passed'
    assert len(result.raw.additions) == count and result.raw_audit.passed and result.rounded_audit.passed
    assert audit_hydrogen_cip(view, authority, result.raw).passed
    assert audit_hydrogen_cip(view, authority, result.rounded).passed
    assert not result.md_ready and not result.geometry_approved


@pytest.mark.parametrize('state', ['PRO', 'HIP'])
def test_old_derivation_still_refuses(state):
    _, authority, _, _ = case(state)
    with pytest.raises(ValueError):
        derive_local_peptide(authority.parent, authority.witness.source_units[0][1], authority.identity, authority.inventory)


@pytest.mark.parametrize('field', ['target', 'requested_state', 'parent_digest', 'alanine_digest', 'inventory_digest', 'substitutions'])
def test_every_declaration_binding_and_directed_role_is_required(field):
    _, authority, _, _ = case('HIP')
    declaration = authority.declaration
    changes = {'target': ComponentInstanceRef(1, 'A', '82', 'B', occurrence=2), 'requested_state': 'HID',
               'parent_digest': '0'*64, 'alanine_digest': '0'*64, 'inventory_digest': '0'*64,
               'substitutions': tuple(reversed(declaration.substitutions))}
    with pytest.raises(ValueError):
        derive_declared_peptide(authority.parent, authority.witness.source_units[0][1], authority.identity,
                                authority.inventory, replace(declaration, **{field: changes[field]}))


@pytest.mark.parametrize('mutation', ['state', 'witness', 'removed'])
def test_constructed_authority_cannot_be_replaced_with_unbound_projection(mutation):
    _, authority, _, _ = case()
    if mutation == 'state':
        changes = {'state': authority.parent.state}
    elif mutation == 'witness':
        changes = {'witness': replace(authority.witness, extra_bonds=())}
    else:
        changes = {'removed_atoms': frozenset()}
    with pytest.raises(ValueError, match='projection-binding'):
        replace(authority, **changes)


def test_existing_variant_owner_not_a_duplicate_mapping_is_consulted(monkeypatch):
    _, authority, _, _ = case('HIP')
    monkeypatch.setitem(variants.VARIANT_TO_PARENT, 'HIP', 'ALA')
    with pytest.raises(ValueError, match='parent-name-mismatch'):
        derive_declared_peptide(authority.parent, authority.witness.source_units[0][1], authority.identity,
                                authority.inventory, authority.declaration)


@pytest.mark.parametrize('state', ['HID', 'HIE'])
def test_neutral_names_cannot_edit_unchanged_cationic_graph(state):
    parent, ala, inventory, declaration = declaration_inputs('HIS', state, ComponentInstanceRef(1, 'A', '82', 'B'), 'H2')
    before = parent.state.graph
    with pytest.raises(ValueError, match='unchanged-graph-inventory-mismatch'):
        derive_declared_peptide(parent, ala, declaration.target, inventory.canonical, declaration)
    assert parent.state.graph is before and before.net_charge == 1


def test_PRO_leaving_atom_must_exist_and_have_source_leaving_flag():
    parent, ala, inventory, declaration = declaration_inputs('PRO', 'PRO', ComponentInstanceRef(1, 'A', '82'), 'H2')
    with pytest.raises(ValueError, match='removal outside'):
        derive_declared_peptide(parent, ala, declaration.target, inventory.canonical, declaration)
    _, _, _, declared = case()
    authority = declared.request.components[0].resolved
    assert authority is not None
    # HA exists but is not a CCD leaving atom; graph valence may not be guessed.
    parent, ala, inventory, declaration = declaration_inputs('PRO', 'PRO', ComponentInstanceRef(1, 'A', '82'), 'HA')
    with pytest.raises(ValueError, match='not a declared source leaving'):
        derive_declared_peptide(parent, ala, declaration.target, inventory.canonical, declaration)


@pytest.mark.parametrize('state', ['PRO', 'HIP'])
def test_backbone_H_typing_presence_is_explicit(state):
    base, authority, incoming, _ = case(state)
    types = ('N', 'H', 'XC') if state == 'PRO' else None
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    with pytest.raises(ValueError):
        compile_declared_hydrogens(base, authority, incoming, parameters, backbone_types=types,
                                   fixed_policy=FIXED_POLICY, generated_policy=GENERATED_POLICY)


@pytest.mark.parametrize('mutation', ['missing-link', 'unknown-order', 'wrong-kind', 'missing-observation', 'partial', 'generated'])
def test_incoming_carbon_cannot_be_an_unbound_or_mobile_witness(mutation):
    base, authority, incoming, _ = case()
    if mutation == 'missing-link':
        base = replace(base, connections=())
    elif mutation == 'unknown-order':
        base = replace(base, connections=(replace(base.connections[0], link=replace(base.connections[0].link, order=None)),))
    elif mutation == 'wrong-kind':
        from dvbfixer.chemical_reconstruction.contracts import LinkKind
        base = replace(base, connections=(replace(base.connections[0], kind=LinkKind.UNKNOWN,
                           link=replace(base.connections[0].link, order=None)),))
    elif mutation == 'missing-observation':
        base = replace(base, observations=tuple(o for o in base.observations if o.identity != incoming))
    elif mutation == 'partial':
        base = replace(base, observations=tuple(replace(o, occupancy=.5) if o.identity == incoming else o for o in base.observations))
    else:
        base = replace(base, ownership=replace(base.ownership, generated=frozenset({incoming})))
    with pytest.raises(ValueError):
        compile_view(base, authority, incoming, {})


def test_PRO_no_imaginary_backbone_H_and_complete_source_metric_coverage():
    base, authority, incoming, view = case()
    assert isinstance(view.plan, JointSourceMetricGeometryPlan) and isinstance(view.plan, SourceMetricGeometryPlan)
    assert not any(a.name == 'H' for a in authority.state.graph.atoms)
    assert not any(g.kind == 'backbone-amide' for g in view.groups)
    assert not any(a.name == 'H' for a in view.request.ownership.generated)
    assert all(incoming not in m.atoms for m in view.plan.metrics)
    assert len(view.plan.metrics) == 25
    graph = authority.state.graph
    target = authority.identity
    adjacent = {a.name: set() for a in graph.atoms}
    expected = set()
    generated = {a.name for a in view.request.ownership.generated}
    for b in graph.bonds:
        adjacent[b.atom1].add(b.atom2)
        adjacent[b.atom2].add(b.atom1)
        if {b.atom1, b.atom2} & generated:
            expected.add(metric_key((ExactAtomRef(target, b.atom1), ExactAtomRef(target, b.atom2))))
    for centre, neighbours in adjacent.items():
        for a, b in combinations(neighbours, 2):
            if {a, centre, b} & generated:
                expected.add(metric_key(tuple(ExactAtomRef(target, n) for n in (a, centre, b))))
    assert {metric_key(m.atoms) for m in view.plan.metrics} == expected
    assert all(m.source_sha256 == canonical_graph_digest(graph) for m in view.plan.metrics)
    assert view.plan.policy == FIXED_POLICY and view.plan.generated_policy == GENERATED_POLICY
    plane = next(c for c in view.request.constraints if c.kind == ConstraintKind.PLANARITY and incoming in c.atoms)
    assert {a.name for a in plane.atoms} == {'N', 'C', 'CA', 'CD'} and plane.upper == .05
    assert not metric_findings(view.request, view.plan)
    assert view.request.observations == base.observations


@pytest.mark.parametrize('mutation', ['value', 'missing', 'extra-observed', 'pin', 'graph-pin', 'mask'])
def test_source_metric_plan_cannot_be_used_as_arbitrary_override(mutation):
    _, _, _, view = case()
    plan, request = view.plan, view.request
    first = plan.metrics[0]
    if mutation == 'value':
        plan = replace(plan, metrics=(replace(first, value=first.value+.000001), *plan.metrics[1:]))
    elif mutation == 'missing':
        plan = replace(plan, metrics=plan.metrics[1:])
    elif mutation == 'extra-observed':
        from dvbfixer.chemical_reconstruction.contracts import Unit
        a, b = (ExactAtomRef(next(iter(request.targets)), n) for n in ('N', 'CA'))
        plan = replace(plan, metrics=(*plan.metrics, replace(first, atoms=(a, b), value=1.4, unit=Unit.ANGSTROM)))
    elif mutation == 'pin':
        plan = replace(plan, metrics=(replace(first, source_sha256='0'*64), *plan.metrics[1:]))
    elif mutation == 'graph-pin':
        plan = replace(plan, authority_graph_digest='0'*64)
    else:
        request = replace(request, ownership=replace(request.ownership, generated=request.ownership.generated-{next(iter(request.ownership.generated))}))
        plan = replace(plan, request_digest=request.digest)
    assert metric_findings(request, plan)
    result = propose(request, plan)
    assert result.status != 'declared-geometry-passed' and result.raw is None


def test_old_override_guard_and_nonempty_invariant_stay_intact():
    _, _, _, view = case()
    old = JointMetricGeometryPlan(view.request.digest, view.plan.strategy, view.plan.policy, view.plan.placements,
                                 view.plan.metrics, generated_policy=view.plan.generated_policy, joint_groups=view.plan.joint_groups)
    assert 'metric-generated-atom-strategy-unqualified' in metric_findings(view.request, old)
    with pytest.raises(ValueError, match='nonempty'):
        GeneratedMetricGeometryPlan(view.request.digest, view.plan.strategy, view.plan.policy, view.plan.placements, (),
                                    generated_policy=view.plan.generated_policy)


@pytest.mark.parametrize('state', ['PRO', 'HIP'])
def test_declared_geometry_is_proper_rigid_equivariant_and_corruption_rejected(state):
    _, authority, _, view = case(state)
    result = propose(view.request, view.plan)
    def transform(p):
        return (-p[1]+4., p[0]-3., p[2]+7.)
    request = replace(view.request, observations=tuple(replace(o, position=transform(o.position)) for o in view.request.observations))
    other = propose(request, replace(view.plan, request_digest=request.digest))
    assert other.status == 'declared-geometry-passed'
    for a, b in zip(result.raw.additions, other.raw.additions):
        assert a.identity == b.identity and b.position == pytest.approx(transform(a.position), abs=1e-10)
    h = result.raw.additions[0]
    bad = replace(result.raw, additions=(replace(h, position=(h.position[0]+.3, h.position[1], h.position[2])), *result.raw.additions[1:]))
    assert not audit_candidate(view.request, view.plan, bad).passed
    assert audit_hydrogen_cip(view, authority, result.raw).passed


def test_eight_frozen_controls_preserve_neutral_state_and_bad_declaration_refusals(tmp_path):
    summary = run_controls(tmp_path/'controls')
    assert summary['declared'] == summary['executed_entries'] == 8
    assert [r['status'] for r in summary['controls']] == ['declared-geometry-passed']*4+['refused']*4
    assert [r['expected_hydrogens'] for r in summary['controls'][:4]] == [7, 8, 7, 8]
    assert [r.get('target') for r in summary['controls'][2:4]] == ['19', '68']
    with pytest.raises(FileExistsError):
        run_controls(tmp_path/'controls')
