"""Synthetic metric/H controls; no observed-source or general H qualification."""
from __future__ import annotations

from dataclasses import replace

import pytest
from scripts.chemical_reconstruction_amide_controls import (
    FFBONDED,
    POLICY,
    run_controls,
    synthetic_case,
)
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import load_ccd

from dvbfixer.chemical_reconstruction import construction
from dvbfixer.chemical_reconstruction.amide_authority import (
    FFBONDED_SHA256,
    _select,
    load_amide_authority,
)
from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    ContractError,
    Observation,
    Operation,
    Unit,
    canonical_graph_digest,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.metric_geometry import GeometryMetric, metric_findings
from dvbfixer.chemical_reconstruction.observed_peptide import prepare_observed_probe
from dvbfixer.chemical_reconstruction.peptide_state import derive_local_peptide
from dvbfixer.chemical_reconstruction.proposals import GeometryPlan, Strategy
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


@pytest.mark.parametrize("name", ["ALA", "THR", "ILE", "LYS", "PHE"])
def test_frozen_synthetic_amide_has_two_independent_passing_audits(name) -> None:
    request, plan, _ = synthetic_case(name)
    before = request.observations
    result = propose(request, plan)
    assert result.status == "declared-geometry-passed", result
    assert result.raw_audit.passed and result.rounded_audit.passed
    assert len(result.raw.additions) == 1 and result.raw.additions[0].element == "H"
    assert request.observations == before
    assert not result.geometry_approved and not result.md_ready and not result.source_coverage_verified


def test_distorted_fixed_angle_retains_raw_and_rounded_rejection() -> None:
    request, plan, _ = synthetic_case("ALA", incoming_angle=100.)
    before = request.observations
    result = propose(request, plan)
    assert result.status == "rejected"
    assert result.raw.additions and result.rounded.additions
    assert not result.raw_audit.passed and not result.rounded_audit.passed
    assert "explicit-metric-geometry" in result.raw_audit.findings
    assert request.observations == before


def test_raw_proposal_is_rigid_frame_equivariant() -> None:
    request, plan, _ = synthetic_case("THR")
    result = propose(request, plan)
    def transform(p):
        return (-p[1]+10., p[0]-7., p[2]+3.)
    moved = replace(request, observations=tuple(replace(o, position=transform(o.position)) for o in request.observations))
    moved_plan = replace(plan, request_digest=moved.digest)
    other = propose(moved, moved_plan)
    assert other.status == "declared-geometry-passed", other
    assert other.raw.additions[0].position == pytest.approx(transform(result.raw.additions[0].position), abs=1e-12)


def test_free_ccd_h_position_is_not_used_for_construction_or_metric_audit() -> None:
    request, plan, _ = synthetic_case("ALA")
    result = propose(request, plan)
    target = next(c for c in request.components if c.resolved is not None)
    state = target.resolved
    graph = replace(state.graph, atoms=tuple(replace(a, ideal_position=(100., -200., 300.)) if a.name == "H" else a
                                            for a in state.graph.atoms))
    state = replace(state, graph=graph, authority=replace(state.authority, source_digest=canonical_graph_digest(graph)))
    altered = replace(request, components=tuple(replace(c, states=(state,)) if c.identity == target.identity else c
                                               for c in request.components))
    other = propose(altered, replace(plan, request_digest=altered.digest))
    assert other.status == "declared-geometry-passed", other
    assert other.raw.additions == result.raw.additions


@pytest.mark.parametrize("missing", [0, 1, 2])
def test_each_metric_is_mandatory_in_constructor_and_independent_auditor(missing) -> None:
    request, plan, _ = synthetic_case("ALA")
    result = propose(request, plan)
    incomplete = replace(plan, metrics=tuple(m for i, m in enumerate(plan.metrics) if i != missing))
    bad = propose(request, incomplete)
    assert bad.status == "invalid-plan"
    candidate = replace(result.raw, plan_digest=incomplete.digest)
    report = audit_candidate(request, incomplete, candidate)
    assert not report.passed and "planar-needs-explicit-length-and-both-angles" in report.findings


def test_duplicate_reversed_metric_and_observed_only_overrides_refuse() -> None:
    request, plan, _ = synthetic_case("ALA")
    with pytest.raises(ContractError, match="duplicate"):
        replace(plan, metrics=(*plan.metrics, replace(plan.metrics[0], atoms=tuple(reversed(plan.metrics[0].atoms)))))
    n, _, ca = plan.placements[0].anchors
    override = GeometryMetric((n, ca), 1.449, Unit.ANGSTROM, "negative control", "1"*64)
    result = propose(request, replace(plan, metrics=(*plan.metrics, override)))
    assert result.status == "invalid-plan" and result.reason == "metric-cannot-override-only-observed-geometry"


@pytest.mark.parametrize("mutation", ["missing-link", "wrong-link", "no-plane", "collinear"])
def test_external_geometry_admission_does_not_guess_or_repair(mutation) -> None:
    request, plan, _ = synthetic_case("ALA")
    if mutation == "missing-link":
        request = replace(request, connections=())
    elif mutation == "wrong-link":
        connection = request.connections[0]
        wrong = replace(connection.link, other_atom=plan.placements[0].anchors[2])
        request = replace(request, connections=(replace(connection, link=wrong),))
    elif mutation == "no-plane":
        request = replace(request, constraints=request.constraints[:-1])
    else:
        centre, incoming, alpha = plan.placements[0].anchors
        points = {o.identity: o.position for o in request.observations}
        point = tuple(2*points[centre][i]-points[alpha][i] for i in range(3))
        request = replace(request, observations=tuple(replace(o, position=point) if o.identity == incoming else o
                                                     for o in request.observations))
    result = propose(request, replace(plan, request_digest=request.digest))
    assert result.status in {"invalid-plan", "construction-failed"}, result
    assert result.status != "declared-geometry-passed"


@pytest.mark.parametrize("mutation", ["position", "length", "angle", "request-binding"])
def test_measurement_only_audit_rejects_corruption(mutation, monkeypatch) -> None:
    request, plan, _ = synthetic_case("ALA")
    candidate = propose(request, plan).raw
    monkeypatch.setattr(construction, "planar_exterior_tip", lambda *a, **k: pytest.fail("auditor called constructor"))
    assert audit_candidate(request, plan, candidate).passed
    if mutation == "position":
        o = candidate.additions[0]
        candidate = replace(candidate, additions=(replace(o, position=(o.position[0], o.position[1], o.position[2]+.2)),))
    elif mutation in {"length", "angle"}:
        index = 0 if mutation == "length" else 1
        metric = plan.metrics[index]
        metrics = tuple(replace(m, value=m.value+(0.1 if index == 0 else 20.)) if i == index else m
                        for i, m in enumerate(plan.metrics))
        plan = replace(plan, metrics=metrics)
        candidate = replace(candidate, plan_digest=plan.digest)
        assert metric.value != plan.metrics[index].value
    else:
        candidate = replace(candidate, request_digest="0"*64)
    assert not audit_candidate(request, plan, candidate).passed


def test_retained_opaque_context_clash_is_not_deleted_or_excluded() -> None:
    request, plan, _ = synthetic_case("ALA")
    p = propose(request, plan).raw.additions[0].position
    ref = ExactAtomRef(ComponentInstanceRef(1, "Z", "99"), "C1")
    request = replace(request, components=(*request.components, Component(ref.component, ChemistryKind.UNKNOWN)),
                       observations=(*request.observations, Observation(ref, "C", (p[0]+.2, p[1], p[2]), Unit.ANGSTROM)))
    result = propose(request, replace(plan, request_digest=request.digest))
    assert result.status == "rejected"
    assert "generated-context-minimum-separation" in result.raw_audit.findings
    assert ref in {o.identity for o in request.observations}


def test_baseline_plan_digest_and_wp2b_h_refusal_stay_unchanged() -> None:
    plan = GeometryPlan("0"*64, Strategy.INTERNAL, POLICY)
    assert plan.digest == "289498f8dae84679ba36d431efdb92e760ee2f3793c4a137b4959343f2562fe1"
    inv = transport_inventory(raw_inventory("ALA"), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
    a = derive_local_peptide(load_ccd("ALA"), load_ccd("ALA"), ComponentInstanceRef(1, "A", "2"), inv.canonical)
    with pytest.raises(ValueError, match="amide-hydrogen-geometry-authority-unqualified"):
        prepare_observed_probe(b"", a, withheld=frozenset(), operation=Operation.HYDROGENS)


@pytest.mark.parametrize("mutation", ["function", "types", "stiffness", "pin", "value"])
def test_unqualified_parameter_record_refuses(mutation) -> None:
    authority = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    with pytest.raises(ValueError):
        if mutation == "function":
            replace(authority, bond_row=(*authority.bond_row[:2], "2", *authority.bond_row[3:]))
        elif mutation == "types":
            replace(authority, bond_row=("X", *authority.bond_row[1:]))
        elif mutation == "stiffness":
            replace(authority, alpha_angle_row=(*authority.alpha_angle_row[:-1], "400.0"))
        elif mutation == "value":
            replace(authority, bond_row=(*authority.bond_row[:-2], "0.102", authority.bond_row[-1]))
        else:
            replace(authority, source_sha256="0"*64)


def test_parameter_loader_checks_pin_before_use_and_requires_unique_rows(tmp_path) -> None:
    path = tmp_path / "bad.itp"
    path.write_text("[ bondtypes ]\nH N 1 .101 363171.2\n")
    with pytest.raises(ValueError, match="pin-mismatch"):
        load_amide_authority(path, expected_sha256=FFBONDED_SHA256)
    with pytest.raises(ValueError, match="pinned"):
        load_amide_authority(path, expected_sha256="0"*64)
    for text in ("[ bondtypes ]\n", "[ bondtypes ]\nH N 1 .101 363171.2\nN H 1 .101 363171.2\n"):
        with pytest.raises(ValueError, match="not-unique"):
            _select(text, "bondtypes", ("H", "N"))


def test_plain_plan_cannot_bypass_metric_requirement_in_auditor() -> None:
    request, plan, _ = synthetic_case("ALA")
    result = propose(request, plan)
    plain = GeometryPlan(request.digest, plan.strategy, plan.policy, plan.placements)
    report = audit_candidate(request, plain, replace(result.raw, plan_digest=plain.digest))
    assert not report.passed and "planar-needs-explicit-metric-plan" in report.findings


def test_planar_anchors_must_be_actual_observations_not_only_link_metadata() -> None:
    request, plan, _ = synthetic_case("ALA")
    incoming = plan.placements[0].anchors[1]
    altered = replace(request, observations=tuple(o for o in request.observations if o.identity != incoming))
    assert "planar-anchor-observations-required" in metric_findings(altered, replace(plan, request_digest=altered.digest))


def test_frozen_six_control_replay_is_not_observed_qualification(tmp_path) -> None:
    result = run_controls(tmp_path / "controls")
    assert result["declared"] == result["executed_entries"] == 6
    assert [r["status"] for r in result["controls"]] == ["declared-geometry-passed"]*5 + ["rejected"]
    assert not result["observed_pdb_reconstruction"] and not result["h_or_complex_qualification"]
    with pytest.raises(FileExistsError):
        run_controls(tmp_path / "controls")
