"""Joint orientation tests; known-case comparison is not broad H qualification."""
from __future__ import annotations

from collections import Counter
from dataclasses import replace

import numpy as np
import pytest
from scripts.chemical_reconstruction_hydrogen_controls import synthetic_view
from scripts.chemical_reconstruction_joint_h_controls import ROOT, run_controls, source_view

from dvbfixer.chemical_reconstruction import construction, joint_hydrogen_geometry
from dvbfixer.chemical_reconstruction.authority_geometry import geometry_view
from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import canonical_graph_digest
from dvbfixer.chemical_reconstruction.hydrogen_completion import (
    audit_hydrogen_cip,
    compile_joint_hydrogens,
)
from dvbfixer.chemical_reconstruction.joint_hydrogen_geometry import (
    DirectionGroup,
    fit_joint_directions,
    joint_group_findings,
    joint_group_points,
)
from dvbfixer.chemical_reconstruction.measurement import angle_degrees
from dvbfixer.chemical_reconstruction.metric_geometry import GeneratedMetricGeometryPlan
from dvbfixer.chemical_reconstruction.proposals import PlacementKind


def joint(name="ALA"):
    old, authority, _ = synthetic_view(name)
    return compile_joint_hydrogens(old), old, authority


def test_origin_constrained_proper_rotation_reuses_antipodal_owner(monkeypatch) -> None:
    from dvbfixer.model.diffusion import geometry
    owner = geometry.weighted_kabsch
    seen = []
    def capture(a, b, weights=None):
        seen.append((a.copy(), b.copy(), weights))
        return owner(a, b, weights)
    monkeypatch.setattr(geometry, "weighted_kabsch", capture)
    fit = fit_joint_directions(((2., 0., 0.), (0., 4., 0.)), ((0., 3., 0.), (-5., 0., 0.)))
    assert fit.apply((1., 0., 0.)) == pytest.approx((0., 1., 0.))
    assert fit.apply((0., 1., 0.)) == pytest.approx((-1., 0., 0.))
    rotation = np.asarray(fit.rotation)
    assert np.linalg.det(rotation) == pytest.approx(1.)
    assert rotation.T @ rotation == pytest.approx(np.eye(3))
    assert fit.translation_norm < 1e-10 and fit.directional_chord_rms < 1e-12
    assert len(seen) == 1 and seen[0][2] is None
    for cloud in seen[0][:2]:
        assert cloud.shape == (4, 3) and cloud.mean(axis=0) == pytest.approx(np.zeros(3))
        assert cloud[::2] == pytest.approx(-cloud[1::2])


@pytest.mark.parametrize("source,target", [
    (((0., 0., 0.), (0., 1., 0.)), ((1., 0., 0.), (0., 1., 0.))),
    (((1., 0., 0.), (2., 0., 0.)), ((1., 0., 0.), (0., 1., 0.))),
    (((1., 0., 0.), (0., 1., 0.)), ((1., 0., 0.), (1., 1e-9, 0.))),
    (((float('nan'), 0., 0.), (0., 1., 0.)), ((1., 0., 0.), (0., 1., 0.))),
    (((1., 0., 0.),), ((1., 0., 0.),)),
])
def test_degenerate_or_invalid_frames_refuse_without_an_alternate_pair(source, target) -> None:
    with pytest.raises(ValueError):
        fit_joint_directions(source, target)


def test_mirrored_input_never_returns_an_improper_rotation() -> None:
    fit = fit_joint_directions(((1., 0., 0.), (0., 1., 0.), (0., 0., 1.)),
                               ((-1., 0., 0.), (0., 1., 0.), (0., 0., 1.)))
    assert np.linalg.det(np.asarray(fit.rotation)) == pytest.approx(1.)
    assert fit.directional_chord_rms > .1


@pytest.mark.parametrize("bad", ["reflection", "translation", "nonorthogonal"])
def test_broken_fit_owner_output_cannot_move_or_reflect_observations(bad, monkeypatch) -> None:
    from dvbfixer.model.diffusion import geometry
    rotation = np.eye(3)
    translation = np.zeros(3)
    if bad == "reflection":
        rotation[0, 0] = -1.
    elif bad == "nonorthogonal":
        rotation[0, 0] = 2.
    else:
        translation[0] = .1
    monkeypatch.setattr(geometry, "weighted_kabsch", lambda *a: geometry.RigidTransform(rotation, translation))
    with pytest.raises(ValueError, match="proper origin-constrained"):
        fit_joint_directions(((1., 0., 0.), (0., 1., 0.)), ((1., 0., 0.), (0., 1., 0.)))


@pytest.mark.parametrize("name", ["ALA", "THR", "ILE", "LYS", "PHE"])
def test_joint_nominal_preserves_request_constraints_and_complete_inventory(name) -> None:
    view, old, authority = joint(name)
    assert view.request is old.request and view.groups == old.groups
    assert view.plan.policy == old.plan.policy and view.plan.generated_policy == old.plan.generated_policy
    assert view.plan.digest != old.plan.digest
    result = propose(view.request, view.plan)
    assert result.status == "declared-geometry-passed", result
    assert len(result.raw.additions) == len(view.request.ownership.generated)
    assert result.rounded_audit.passed and audit_hydrogen_cip(view, authority, result.raw).passed
    assert not result.geometry_approved and not result.md_ready


def test_one_fit_per_whole_group_and_independent_measurement_auditor(monkeypatch) -> None:
    view, _, _ = joint("LYS")
    original = construction.joint_group_points
    calls = []
    def capture(request, group):
        calls.append(group.parent)
        return original(request, group)
    monkeypatch.setattr(construction, "joint_group_points", capture)
    result = propose(view.request, view.plan)
    assert Counter(calls) == Counter(g.parent for g in view.plan.joint_groups)
    monkeypatch.setattr(joint_hydrogen_geometry, "fit_joint_directions", lambda *a: pytest.fail("auditor called fitter"))
    assert audit_candidate(view.request, view.plan, result.raw).passed
    h = result.raw.additions[0]
    bad = replace(result.raw, additions=(replace(h, position=(h.position[0]+.3, h.position[1], h.position[2])), *result.raw.additions[1:]))
    assert not audit_candidate(view.request, view.plan, bad).passed


def test_all_three_heavy_neighbours_influence_fit_not_only_old_cone_pair() -> None:
    view, _, _ = joint("ALA")
    group = view.plan.joint_groups[0]
    third = group.heavy_neighbours[2]
    assert third not in next(p.anchors for p in view.plan.placements if p.atom in group.hydrogens)
    before = joint_group_points(view.request, group)
    observations = tuple(replace(o, position=(o.position[0]+.07, o.position[1]-.02, o.position[2]+.03)) if o.identity == third else o
                         for o in view.request.observations)
    after = joint_group_points(replace(view.request, observations=observations), group)
    assert np.linalg.norm(np.asarray(after[group.hydrogens[0]])-np.asarray(before[group.hydrogens[0]])) > 1e-8


@pytest.mark.parametrize("name", ["ILE", "LYS"])
def test_known_observed_angular_deficits_resolve_without_changing_the_old_result(name) -> None:
    original = (ROOT / "tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    old, authority, _, _, _, _ = source_view(original, name, None)
    assert propose(old.request, old.plan).status == "rejected"
    view = compile_joint_hydrogens(old)
    result = propose(view.request, view.plan)
    assert view.request is old.request and result.status == "declared-geometry-passed", result
    assert result.raw_audit.passed and result.rounded_audit.passed
    assert audit_hydrogen_cip(view, authority, result.raw).passed and audit_hydrogen_cip(view, authority, result.rounded).passed
    ideals, bonds = geometry_view(view.request)
    points = {o.identity: o.position for o in view.request.observations} | {o.identity: o.position for o in result.raw.additions}
    for g in view.plan.joint_groups:
        for h in g.hydrogens:
            assert np.linalg.norm(np.asarray(points[h])-np.asarray(points[g.parent])) == pytest.approx(bonds[frozenset((h, g.parent))], abs=1e-12)
        if len(g.hydrogens) == 2:
            a, b = g.hydrogens
            assert angle_degrees(points[a], points[g.parent], points[b]) == pytest.approx(angle_degrees(ideals[a], ideals[g.parent], ideals[b]), abs=1e-10)


def test_proper_rigid_equivariance_with_all_groups() -> None:
    view, _, _ = joint("LYS")
    result = propose(view.request, view.plan)
    def transform(p):
        return (-p[1]+10., p[0]-7., p[2]+3.)
    request = replace(view.request, observations=tuple(replace(o, position=transform(o.position)) for o in view.request.observations))
    other = propose(request, replace(view.plan, request_digest=request.digest))
    assert other.status == "declared-geometry-passed"
    for before, after in zip(result.raw.additions, other.raw.additions):
        assert after.identity == before.identity and after.position == pytest.approx(transform(before.position), abs=1e-10)


@pytest.mark.parametrize("change", ["missing", "duplicate", "wrong-kind"])
def test_group_plan_membership_cannot_be_silently_incomplete_or_replaced(change) -> None:
    view, _, _ = joint("LYS")
    with pytest.raises(ValueError):
        if change == "missing":
            replace(view.plan, joint_groups=view.plan.joint_groups[:-1])
        elif change == "duplicate":
            replace(view.plan, joint_groups=(*view.plan.joint_groups, view.plan.joint_groups[0]))
        else:
            placements = tuple(replace(p, kind=PlacementKind.CONE) if p.kind == PlacementKind.JOINT else p for p in view.plan.placements)
            replace(view.plan, placements=placements)


def test_plain_plan_does_not_silently_fall_back_to_independent_cones() -> None:
    view, _, _ = joint("ALA")
    result = propose(view.request, view.plan)
    plain = GeneratedMetricGeometryPlan(view.request.digest, view.plan.strategy, view.plan.policy, view.plan.placements, view.plan.metrics,
                                        generated_policy=view.plan.generated_policy)
    assert propose(view.request, plain).reason == "joint-group-plan-required"
    report = audit_candidate(view.request, plain, replace(result.raw, plan_digest=plain.digest))
    assert not report.passed and "joint-group-plan-required" in report.findings


@pytest.mark.parametrize("mutation", ["element", "charge", "order", "observation", "occupancy", "ownership", "mobile-anchor", "generated-anchor"])
def test_joint_authority_and_fixed_ownership_fail_closed(mutation) -> None:
    view, _, _ = joint("ALA")
    r = view.request
    group = view.plan.joint_groups[0]
    if mutation in {"element", "charge", "order"}:
        component = next(c for c in r.components if c.identity == group.parent.component)
        state = component.resolved
        graph = state.graph
        if mutation == "order":
            graph = replace(graph, bonds=tuple(replace(b, order=2.) if {b.atom1, b.atom2} == {group.parent.name, group.hydrogens[0].name} else b for b in graph.bonds))
        else:
            changes = {"element": "N"} if mutation == "element" else {"formal_charge": 1}
            atoms = tuple(replace(a, **changes) if a.name == group.parent.name else a for a in graph.atoms)
            graph = replace(graph, atoms=atoms, net_charge=sum(a.formal_charge for a in atoms))
        state = replace(state, graph=graph, authority=replace(state.authority, source_digest=canonical_graph_digest(graph)))
        r = replace(r, components=tuple(replace(c, states=(state,)) if c.identity == component.identity else c for c in r.components))
    elif mutation == "observation":
        r = replace(r, observations=tuple(o for o in r.observations if o.identity != group.heavy_neighbours[2]))
    elif mutation == "occupancy":
        r = replace(r, observations=tuple(replace(o, occupancy=.5) if o.identity == group.parent else o for o in r.observations))
    elif mutation == "ownership":
        r = replace(r, ownership=replace(r.ownership, generated=r.ownership.generated-{group.hydrogens[0]}))
    elif mutation == "mobile-anchor":
        r = replace(r, ownership=replace(r.ownership, mobile=frozenset({group.parent})))
    else:
        r = replace(r, ownership=replace(r.ownership, generated=r.ownership.generated | {group.parent}))
    assert joint_group_findings(r, replace(view.plan, request_digest=r.digest))


def test_group_identity_shape_is_exact() -> None:
    view, _, _ = joint("ALA")
    g = view.plan.joint_groups[0]
    with pytest.raises(ValueError):
        replace(g, hydrogens=())
    with pytest.raises(ValueError):
        replace(g, heavy_neighbours=tuple(reversed(g.heavy_neighbours)))
    with pytest.raises(ValueError):
        DirectionGroup(g.parent, g.heavy_neighbours[:2], (g.hydrogens[0], g.hydrogens[0]))


def test_frozen_twelve_entries_preserve_refusals_and_nonoverwriting_evidence(tmp_path) -> None:
    output = tmp_path / "controls"
    result = run_controls(output)
    assert result["declared"] == result["executed_entries"] == 12
    observed = [r for r in result["controls"] if r["arm"] == "observed"]
    assert [r["status"] for r in observed] == ["declared-geometry-passed"]*3+["unsupported", "declared-geometry-passed", "refused"]
    assert result["controls"][5]["status"] == "rejected"
    assert all(r.get('request_and_constraints_unchanged', True) for r in result['controls'])
    with pytest.raises(FileExistsError):
        run_controls(output)
