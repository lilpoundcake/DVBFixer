"""Explicit planar-N/thiol groups; no hidden state or source-geometry correction."""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scripts.chemical_reconstruction_extended_h_controls import (
    FIXTURES,
    compile_view,
    run_controls,
    synthetic_base,
)
from scripts.chemical_reconstruction_hydrogen_controls import (
    FFBONDED,
    FIXED_POLICY,
    GENERATED_POLICY,
    backbone_types,
    synthetic_view,
)

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import ConstraintKind, canonical_graph_digest
from dvbfixer.chemical_reconstruction.hydrogen_completion import (
    audit_hydrogen_cip,
    compile_complete_hydrogens,
    compile_extended_hydrogens,
)
from dvbfixer.chemical_reconstruction.measurement import angle_degrees, satisfied, torsion_degrees
from dvbfixer.chemical_reconstruction.proposals import PlacementKind


def base_case(name):
    evidence = {}
    base, authority, incoming = synthetic_base(name, evidence)
    return base, authority, incoming


def compile_case(base, authority, incoming, *, extended=True):
    compiler = compile_extended_hydrogens if extended else compile_complete_hydrogens
    return compiler(base, authority, incoming, load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256),
                    backbone_types=backbone_types(authority.inventory.state_name), fixed_policy=FIXED_POLICY,
                    generated_policy=GENERATED_POLICY)


@pytest.mark.parametrize('name,count', [('ASN', 6), ('GLN', 8), ('ARG', 13), ('CYS', 5), ('TRP', 10), ('GLY', 3)])
def test_explicit_nominal_complete_inventory_rounding_and_CIP(name, count):
    base, authority, incoming = base_case(name)
    evidence = {}
    view = compile_view(base, authority, incoming, evidence)
    result = propose(view.request, view.plan)
    assert result.status == 'declared-geometry-passed', result
    assert len(result.raw.additions) == len(result.rounded.additions) == count
    assert result.raw_audit.passed and result.rounded_audit.passed
    assert audit_hydrogen_cip(view, authority, result.raw).passed
    assert audit_hydrogen_cip(view, authority, result.rounded).passed
    assert view.request.observations is base.observations
    assert all(o.identity in view.request.ownership.generated for o in result.raw.additions)
    assert not result.md_ready and not result.geometry_approved
    assert dict(view.request.origin.settings)['complete_H_policy'] == 'complete-H-groups-v2-explicit-sp2-N-thiol'


@pytest.mark.parametrize('name', ['ASN', 'GLN', 'ARG', 'CYS', 'TRP'])
def test_old_compiler_does_not_silently_gain_new_group_support(name):
    base, authority, incoming = base_case(name)
    with pytest.raises(ValueError, match='hydrogen-group-unqualified'):
        compile_case(base, authority, incoming, extended=False)


@pytest.mark.parametrize('name', ['ALA', 'THR', 'ILE', 'LYS', 'PHE'])
def test_old_five_default_results_stay_equal_and_new_policy_is_explicit(name):
    old, authority, _ = synthetic_view(name)
    # Use the original heavy-only input rather than pretending the full-H request
    # (with its extra constraints/origin) was the historical compiler input.
    from scripts.chemical_reconstruction_amide_controls import synthetic_case

    from dvbfixer.chemical_reconstruction.contracts import (
        Operation,
        Origin,
        Ownership,
        digest_record,
    )
    raw, _, _ = synthetic_case(name)
    obs = tuple(o for o in raw.observations if o.identity.component not in raw.targets or o.element not in {'H', 'D'})
    pin = digest_record(('synthetic-complete-H-heavy-only-v1', name, 121.96, obs))
    base = replace(raw, source_sha256=pin, operation=Operation.HEAVY, ownership=Ownership(frozenset()), observations=obs,
                   origin=Origin('synthetic-complete-H-heavy-only-v1', pin))
    incoming = raw.connections[0].link.component_atom
    assert compile_case(base, authority, incoming, extended=False) == old
    extended = compile_case(base, authority, incoming)
    assert extended.groups == old.groups and extended.plan.placements == old.plan.placements
    assert extended.request.observations == old.request.observations
    assert extended.request.digest != old.request.digest
    assert propose(extended.request, extended.plan).raw.additions == propose(old.request, old.plan).raw.additions


@pytest.mark.parametrize('name', ['ASN', 'GLN', 'ARG', 'TRP'])
def test_nitrogen_needs_exact_existing_SP2_neighbour_plane_authority(name):
    base, authority, incoming = base_case(name)
    authority = replace(authority, parent=replace(authority.parent, planar_groups=()))
    with pytest.raises(ValueError, match='hydrogen-group-unqualified'):
        compile_case(base, authority, incoming)


@pytest.mark.parametrize('name,parent,kind', [('ASN', 'ND2', 'charge'), ('GLN', 'NE2', 'H-bond-order'),
                                           ('CYS', 'SG', 'charge'), ('ASN', 'ND2', 'deuterium')])
def test_wrong_state_or_H_bond_chemistry_is_not_repaired(name, parent, kind):
    base, authority, incoming = base_case(name)
    graph = authority.state.graph
    h = next(b.atom2 if b.atom1 == parent else b.atom1 for b in graph.bonds
             if parent in {b.atom1, b.atom2} and next(a.element for a in graph.atoms
                 if a.name == (b.atom2 if b.atom1 == parent else b.atom1)) == 'H')
    if kind == 'H-bond-order':
        graph = replace(graph, bonds=tuple(replace(b, order=2.) if {b.atom1, b.atom2} == {parent, h} else b for b in graph.bonds))
    else:
        atoms = tuple(replace(a, formal_charge=-1) if kind == 'charge' and a.name == parent else
                      replace(a, element='D') if kind == 'deuterium' and a.name == h else a for a in graph.atoms)
        graph = replace(graph, atoms=atoms, net_charge=sum(a.formal_charge for a in atoms))
    state = replace(authority.state, graph=graph, authority=replace(authority.state.authority, source_digest=canonical_graph_digest(graph)))
    authority = replace(authority, state=state)
    base = replace(base, components=tuple(replace(c, states=(state,)) if c.identity == authority.identity else c for c in base.components))
    with pytest.raises(ValueError):
        compile_case(base, authority, incoming)


@pytest.mark.parametrize('mutation', ['generated', 'mobile', 'missing-heavy', 'measured-H', 'occupancy'])
def test_extended_compiler_retains_fixed_complete_heavy_and_no_measured_H_boundary(mutation):
    base, authority, incoming = base_case('ASN')
    atom = base.observations[0]
    if mutation == 'generated':
        base = replace(base, ownership=replace(base.ownership, generated=frozenset({atom.identity})))
    elif mutation == 'mobile':
        base = replace(base, ownership=replace(base.ownership, mobile=frozenset({atom.identity})))
    elif mutation == 'missing-heavy':
        base = replace(base, observations=base.observations[1:])
    elif mutation == 'measured-H':
        base = replace(base, observations=(replace(atom, element='H'), *base.observations[1:]))
    else:
        base = replace(base, observations=(replace(atom, occupancy=.5), *base.observations[1:]))
    with pytest.raises(ValueError):
        compile_case(base, authority, incoming)


@pytest.mark.parametrize('name', ['ASN', 'GLN', 'ARG'])
def test_named_NH2_uses_deterministic_heavy_chain_and_preserves_H_H_geometry(name):
    base, authority, incoming = base_case(name)
    view = compile_view(base, authority, incoming, {})
    result = propose(view.request, view.plan)
    points = {o.identity: o.position for o in (*view.request.observations, *result.raw.additions)}
    ideals = {a.name: a.ideal_position for a in authority.state.graph.atoms}
    groups = [g for g in view.groups if g.kind == 'planar-NH2']
    assert len(groups) == (2 if name == 'ARG' else 1)
    for g in groups:
        assert all(a in {o.identity for o in base.observations} for a in g.anchors)
        for h in g.hydrogens:
            placement = next(p for p in view.plan.placements if p.atom == h)
            assert placement.kind == PlacementKind.CHAIN and placement.anchors == g.anchors
            assert torsion_degrees(tuple(points[a] for a in (*g.anchors, h))) == pytest.approx(
                torsion_degrees(tuple(ideals[a.name] for a in (*g.anchors, h))), abs=1e-9)
        h1, h2 = g.hydrogens
        assert angle_degrees(points[h1], points[g.parent], points[h2]) == pytest.approx(
            angle_degrees(ideals[h1.name], ideals[g.parent.name], ideals[h2.name]), abs=1e-9)
        assert any(c.kind == ConstraintKind.PLANARITY and set(c.atoms) == {*g.anchors, *g.hydrogens}
                   and c.upper == .005 for c in view.request.constraints)


def test_thiol_has_exact_source_named_phase_and_length():
    base, authority, incoming = base_case('CYS')
    view = compile_view(base, authority, incoming, {})
    result = propose(view.request, view.plan)
    group = next(g for g in view.groups if g.kind == 'thiol')
    h, = group.hydrogens
    points = {o.identity: o.position for o in (*view.request.observations, *result.raw.additions)}
    ideals = {a.name: a.ideal_position for a in authority.state.graph.atoms}
    assert torsion_degrees(tuple(points[a] for a in (*group.anchors, h))) == pytest.approx(
        torsion_degrees(tuple(ideals[a.name] for a in (*group.anchors, h))), abs=1e-9)
    bond = next(b for b in authority.state.graph.bonds if {b.atom1, b.atom2} == {group.parent.name, h.name})
    assert np.linalg.norm(np.asarray(points[h])-np.asarray(points[group.parent])) == pytest.approx(bond.ideal_length, abs=1e-12)


def test_corrupt_NH2_plane_is_rejected_by_measurement_only_audit():
    base, authority, incoming = base_case('ASN')
    view = compile_view(base, authority, incoming, {})
    candidate = propose(view.request, view.plan).raw
    group = next(g for g in view.groups if g.kind == 'planar-NH2')
    points = {o.identity: np.asarray(o.position) for o in base.observations}
    a, b, c = group.anchors
    normal = np.cross(points[b]-points[a], points[c]-points[a])
    normal /= np.linalg.norm(normal)
    bad = replace(candidate, additions=tuple(replace(o, position=tuple(np.asarray(o.position)+.2*normal))
                  if o.identity == group.hydrogens[0] else o for o in candidate.additions))
    audit = audit_candidate(view.request, view.plan, bad)
    assert not audit.passed and 'candidate-constraint-conflict' in audit.findings
    bad_points = {o.identity: o.position for o in (*view.request.observations, *bad.additions)}
    plane = next(c for c in view.request.constraints if c.kind == ConstraintKind.PLANARITY
                 and set(c.atoms) == {*group.anchors, *group.hydrogens})
    assert not satisfied(plane, tuple(bad_points[a] for a in plane.atoms))


@pytest.mark.parametrize('name', ['ASN', 'ARG', 'CYS', 'TRP'])
def test_extended_groups_are_proper_rigid_equivariant(name):
    base, authority, incoming = base_case(name)
    view = compile_view(base, authority, incoming, {})
    result = propose(view.request, view.plan)
    def transform(p):
        return (-p[1]+6., p[0]-3., p[2]+2.)
    request = replace(view.request, observations=tuple(replace(o, position=transform(o.position)) for o in base.observations))
    other = propose(request, replace(view.plan, request_digest=request.digest))
    assert other.status == 'declared-geometry-passed'
    for a, b in zip(result.raw.additions, other.raw.additions):
        assert a.identity == b.identity and b.position == pytest.approx(transform(a.position), abs=1e-10)


def test_sixteen_declared_entries_and_chemical_fixed_input_refusals(tmp_path):
    import json
    output = tmp_path/'controls'
    summary = run_controls(output)
    assert summary['declared'] == summary['executed_entries'] == 16
    statuses = [r['status'] for r in summary['controls']]
    assert statuses == ['declared-geometry-passed']*6+['refused']*2+[
        'declared-geometry-passed', 'declared-geometry-passed', 'unsupported',
        'declared-geometry-passed', 'declared-geometry-passed', 'declared-geometry-passed', 'refused', 'refused']
    arg = json.loads((output/'observed-ARG/evidence.json').read_text())
    assert arg['result']['admission']['status'] == 'invalid-input'
    assert arg['result']['raw'] is None
    assert summary['controls'][6]['reason'] == 'removal outside named source components'
    assert summary['controls'][7]['reason'] == 'explicit-same-name-RTP-state-unavailable:HIS'
    assert len(list(FIXTURES.glob('*.cif'))) == 8
    with pytest.raises(FileExistsError):
        run_controls(output)
