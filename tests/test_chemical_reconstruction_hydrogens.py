"""Complete-H research regressions, not measured/pro-R/pro-S H qualification."""
from __future__ import annotations

from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_amide_controls import FFBONDED
from scripts.chemical_reconstruction_ccd_controls import ROOT
from scripts.chemical_reconstruction_hydrogen_controls import (
    FIXED_POLICY,
    GENERATED_POLICY,
    backbone_types,
    run_controls,
    synthetic_view,
)
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import load_ccd

from dvbfixer.chemical_reconstruction import construction
from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Observation,
    Operation,
    Ownership,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.hydrogen_completion import (
    audit_hydrogen_cip,
    compile_complete_hydrogens,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.measurement import _cross, _dot, _sub, _unit
from dvbfixer.chemical_reconstruction.metric_geometry import MetricGeometryPlan, audit_policy_for
from dvbfixer.chemical_reconstruction.observed_peptide import (
    prepare_observed_probe,
    select_interior_target,
)
from dvbfixer.chemical_reconstruction.peptide_state import derive_local_peptide
from dvbfixer.chemical_reconstruction.proposals import AuditPolicy, PlacementKind
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def source_view(name, *, path="tests/fixtures/protein_context/1UBQ.pdb", number=None):
    original = (ROOT / path).read_bytes()
    prepared = prepare_source(original)
    target = select_interior_target(original, name, chain="A", number=number)
    inventory = transport_inventory(raw_inventory(name), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
    authority = derive_local_peptide(load_ccd(name), load_ccd("ALA"), target, inventory.canonical)
    probe = prepare_source_probe(prepared, authority, withheld=frozenset())
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    view = compile_complete_hydrogens(probe.request, authority, probe.real_peptide_partners[0], parameters,
                                       backbone_types=backbone_types(name), fixed_policy=FIXED_POLICY, generated_policy=GENERATED_POLICY)
    return view, authority, prepared, probe


@pytest.mark.parametrize("name,count", [("ALA", 5), ("THR", 7), ("ILE", 11), ("LYS", 13), ("PHE", 9)])
def test_complete_synthetic_inventory_geometry_and_conditional_cip(name, count) -> None:
    view, authority, _ = synthetic_view(name)
    r = view.request
    before = r.observations
    result = propose(r, view.plan)
    assert result.status == "declared-geometry-passed", result
    assert len(r.ownership.generated) == len(result.raw.additions) == len(result.rounded.additions) == count
    assert all(o.element == "H" for o in result.raw.additions)
    assert all(o.element != "H" for o in r.observations if o.identity.component in r.targets)
    assert audit_hydrogen_cip(view, authority, result.raw).passed
    assert audit_hydrogen_cip(view, authority, result.rounded).passed
    assert r.observations == before
    assert not result.geometry_approved and not result.md_ready and not result.source_coverage_verified


@pytest.mark.parametrize("name", ["ALA", "THR", "ILE", "LYS", "PHE"])
def test_groups_use_only_observed_heavy_anchors_and_cover_each_h_once(name) -> None:
    view, _, _ = synthetic_view(name)
    observed = {o.identity for o in view.request.observations if o.element != "H"}
    assert all(set(p.anchors) <= observed for p in view.plan.placements)
    assert {p.atom for p in view.plan.placements} == view.request.ownership.generated
    assert len(view.plan.placements) == len(view.request.ownership.generated)
    assert sum(len(g.hydrogens) for g in view.groups) == len(view.request.ownership.generated)
    assert {g.kind for g in view.groups} <= {"backbone-amide", "aromatic-CH", "tetrahedral-CH", "methylene", "terminal-three-H", "hydroxyl"}
    assert all(p.kind == PlacementKind.PLANAR for p in view.plan.placements if p.atom.name == "H")


def test_distorted_amide_keeps_complete_raw_failure_without_using_cip_as_geometry_pass() -> None:
    view, authority, _ = synthetic_view("ALA", incoming_angle=100.)
    result = propose(view.request, view.plan)
    assert result.status == "rejected" and len(result.raw.additions) == 5
    assert not result.raw_audit.passed and not result.rounded_audit.passed
    assert audit_hydrogen_cip(view, authority, result.raw).passed


def test_complete_groups_are_proper_rigid_frame_equivariant() -> None:
    view, authority, _ = synthetic_view("THR")
    original = propose(view.request, view.plan).raw
    def transform(p):
        return (-p[1]+10., p[0]-7., p[2]+3.)
    observations = tuple(replace(o, position=transform(o.position)) for o in view.request.observations)
    digest = digest_record(("rigid synthetic H regression", observations))
    request = replace(view.request, observations=observations, source_sha256=digest,
                       origin=replace(view.request.origin, request_digest=digest))
    moved = replace(view, request=request, plan=replace(view.plan, request_digest=request.digest))
    result = propose(request, moved.plan)
    assert result.status == "declared-geometry-passed", result
    for before, after in zip(original.additions, result.raw.additions):
        assert after.identity == before.identity and after.position == pytest.approx(transform(before.position), abs=1e-10)
    assert audit_hydrogen_cip(moved, authority, result.raw).passed


def test_missing_H_or_metric_cannot_be_reported_as_complete_and_context_clash_is_retained() -> None:
    view, authority, _ = synthetic_view("ALA")
    candidate = propose(view.request, view.plan).raw
    missing = replace(candidate, additions=candidate.additions[:-1])
    assert not audit_hydrogen_cip(view, authority, missing).passed
    assert not audit_candidate(view.request, view.plan, missing).passed
    incomplete_plan = replace(view.plan, metrics=view.plan.metrics[1:])
    assert propose(view.request, incomplete_plan).status == "invalid-plan"
    p = candidate.additions[0].position
    ref = ExactAtomRef(ComponentInstanceRef(1, "Z", "99"), "C1")
    request = replace(view.request, components=(*view.request.components, Component(ref.component, ChemistryKind.UNKNOWN)),
                       observations=(*view.request.observations, Observation(ref, "C", (p[0]+.2, p[1], p[2]), Unit.ANGSTROM)))
    result = propose(request, replace(view.plan, request_digest=request.digest))
    assert result.status == "rejected" and "generated-context-minimum-separation" in result.raw_audit.findings


def test_generated_metric_bounds_are_stricter_than_fixed_heavy_and_bound_into_digest() -> None:
    view, _, _ = synthetic_view("ALA")
    r, plan = view.request, view.plan
    target = next(iter(r.targets))
    n, h, ca = (ExactAtomRef(target, s) for s in ("N", "H", "CA"))
    assert audit_policy_for(plan, (n, ca), r.ownership.generated) == FIXED_POLICY
    assert audit_policy_for(plan, (n, h), r.ownership.generated) == GENERATED_POLICY
    loose = MetricGeometryPlan(r.digest, plan.strategy, plan.policy, plan.placements, plan.metrics)
    assert loose.digest != plan.digest
    candidate = propose(r, plan).raw
    observed = {o.identity: o.position for o in r.observations}
    additions = tuple(replace(o, position=tuple(observed[n][i]+(o.position[i]-observed[n][i])*1.03 for i in range(3)))
                      if o.identity == h else o for o in candidate.additions)
    bad = replace(candidate, additions=additions)
    assert not audit_candidate(r, plan, bad).passed
    # Negative audit-scope comparison only; not another qualifying proposal.
    assert audit_candidate(r, loose, replace(bad, plan_digest=loose.digest)).passed
    with pytest.raises(ValueError):
        replace(plan, generated_policy=AuditPolicy(.2, 12., .7, .02))


def test_swap_of_methylene_names_is_not_independent_pro_R_pro_S_evidence() -> None:
    view, authority, _ = synthetic_view("LYS")
    candidate = propose(view.request, view.plan).raw
    group = next(g for g in view.groups if g.kind == "methylene")
    a, b = group.hydrogens
    points = {o.identity: o.position for o in candidate.additions}
    swapped = replace(candidate, additions=tuple(replace(o, position=points[b if o.identity == a else a]) if o.identity in {a, b} else o
                                                 for o in candidate.additions))
    assert not audit_candidate(view.request, view.plan, swapped).passed
    assert audit_hydrogen_cip(view, authority, swapped).passed  # Identical H do not define chemical CIP.


def test_cip_and_geometry_have_distinct_negative_controls_without_proposal_constructor(monkeypatch) -> None:
    view, authority, _ = synthetic_view("ALA")
    candidate = propose(view.request, view.plan).raw
    monkeypatch.setattr(construction, "_internal", lambda *a: pytest.fail("audit called proposal constructor"))
    assert audit_candidate(view.request, view.plan, candidate).passed
    assert audit_hydrogen_cip(view, authority, candidate).passed
    target = next(iter(view.request.targets))
    points = {o.identity.name: o.position for o in view.request.observations if o.identity.component == target}
    first, second, third = (points[n] for n in ("N", "C", "CB"))
    normal = _unit(_cross(_sub(second, first), _sub(third, first)))
    additions = []
    for o in candidate.additions:
        if o.identity.name == "HA":
            delta = _dot(_sub(o.position, first), normal)
            o = replace(o, position=tuple(o.position[i]-2*delta*normal[i] for i in range(3)))
        additions.append(o)
    corrupt = replace(candidate, additions=tuple(additions))
    # RDKit retains S from the unchanged heavy orientation despite this bad H.
    # CIP is not a substitute for complete-ligand geometry measurements.
    assert audit_hydrogen_cip(view, authority, corrupt).passed
    assert not audit_candidate(view.request, view.plan, corrupt).passed
    def mirror(p):
        return (-p[0], p[1], p[2])
    observations = tuple(replace(o, position=mirror(o.position)) for o in view.request.observations)
    request = replace(view.request, observations=observations)
    mirrored = replace(view, request=request, plan=replace(view.plan, request_digest=request.digest))
    d_candidate = replace(candidate, request_digest=request.digest, plan_digest=mirrored.plan.digest,
                           additions=tuple(replace(o, position=mirror(o.position)) for o in candidate.additions))
    assert not audit_hydrogen_cip(mirrored, authority, d_candidate).passed


@pytest.mark.parametrize("mutation", ["measured-H", "missing-heavy", "fractional-heavy", "wrong-types"])
def test_incomplete_or_unqualified_input_is_not_silently_repaired(mutation) -> None:
    view, authority, evidence = synthetic_view("ALA")
    base = replace(view.request, operation=Operation.HEAVY, ownership=Ownership(frozenset()),
                   origin=replace(view.request.origin, settings=()))
    target = authority.identity
    if mutation == "measured-H":
        base = replace(base, observations=(*base.observations, Observation(ExactAtomRef(target, "H"), "H", (0., 0., 0.), Unit.ANGSTROM)))
    elif mutation == "missing-heavy":
        base = replace(base, observations=tuple(o for o in base.observations if o.identity != ExactAtomRef(target, "CB")))
    elif mutation == "fractional-heavy":
        base = replace(base, observations=tuple(replace(o, occupancy=.5) if o.identity == ExactAtomRef(target, "CB") else o for o in base.observations))
    types = ("X", "H", "XC") if mutation == "wrong-types" else backbone_types("ALA")
    incoming = next(g.anchors[1] for g in view.groups if g.kind == "backbone-amide")
    with pytest.raises(ValueError):
        compile_complete_hydrogens(base, authority, incoming, evidence["parameters"], backbone_types=types,
                                    fixed_policy=FIXED_POLICY, generated_policy=GENERATED_POLICY)


def test_fixed_observed_D_ca_is_refused_not_reflected() -> None:
    view, _, _ = synthetic_view("ALA")
    r = view.request
    target = next(iter(r.targets))
    points = {o.identity.name: o.position for o in r.observations if o.identity.component == target}
    centre = points["CA"]
    normal = _unit(_cross(_sub(points["N"], centre), _sub(points["C"], centre)))
    delta = _dot(_sub(points["CB"], centre), normal)
    reflected = tuple(points["CB"][i]-2*delta*normal[i] for i in range(3))
    observations = tuple(replace(o, position=reflected) if o.identity == ExactAtomRef(target, "CB") else o for o in r.observations)
    request = replace(r, observations=observations)
    result = propose(request, replace(view.plan, request_digest=request.digest))
    assert result.status == "unsupported" and result.raw is None
    assert "fixed-constraint-conflict" in {f.code for f in result.admission.findings}
    assert request.observations == observations


@pytest.mark.parametrize("name,path", [("ALA", "tests/fixtures/protein_context/1UBQ.pdb"), ("THR", "tests/fixtures/regressions/1FR2.pdb")])
def test_retained_source_complete_H_passes_without_any_observation_change(name, path) -> None:
    view, authority, prepared, probe = source_view(name, path=path)
    assert view.request.observations == probe.request.observations
    assert audit_source_transport(prepared, replace(probe, request=view.request)).passed
    result = propose(view.request, view.plan)
    assert result.status == "declared-geometry-passed", result
    assert audit_hydrogen_cip(view, authority, result.raw).passed and audit_hydrogen_cip(view, authority, result.rounded).passed
    with pytest.raises(ValueError, match="amide-hydrogen-geometry-authority-unqualified"):
        prepare_observed_probe(b"", authority, withheld=frozenset(), operation=Operation.HYDROGENS)


@pytest.mark.parametrize("name", ["ILE", "LYS"])
def test_observed_negative_geometry_results_remain_failures_despite_correct_cip(name) -> None:
    view, authority, _, _ = source_view(name)
    result = propose(view.request, view.plan)
    assert result.status == "rejected" and "authority-angle-geometry" in result.raw_audit.findings
    assert audit_hydrogen_cip(view, authority, result.raw).passed


def test_PHE_fixed_ring_constraint_is_admission_refusal_not_H_construction_failure() -> None:
    view, _, _, _ = source_view("PHE", number="45")
    result = propose(view.request, view.plan)
    assert result.status == "unsupported" and result.raw is None
    assert "fixed-constraint-conflict" in {f.code for f in result.admission.findings}


def test_twelve_frozen_entries_and_nonoverwriting_evidence(tmp_path) -> None:
    output = tmp_path / "controls"
    result = run_controls(output)
    assert result["declared"] == result["executed_entries"] == 12
    synthetic = [r for r in result["controls"] if r["arm"] == "synthetic"]
    observed = [r for r in result["controls"] if r["arm"] == "observed"]
    assert [r["status"] for r in synthetic] == ["declared-geometry-passed"]*5+["rejected"]
    assert [r["status"] for r in observed] == ["declared-geometry-passed", "rejected", "rejected", "unsupported", "declared-geometry-passed", "refused"]
    assert not result["full_polymer_or_general_H_qualification"] and not result["production_enabled"]
    with pytest.raises(FileExistsError):
        run_controls(output)
