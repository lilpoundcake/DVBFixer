"""Post-assembly anchor motion is explicit; moved donors are still forbidden."""
from dataclasses import asdict, replace
from types import SimpleNamespace

import pytest

from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    DiffusionContractError,
    ResidueIdentity,
)
from dvbfixer.model.diffusion.gap_portfolio import _coordinates, digest
from dvbfixer.model.diffusion.portfolio_anchors import (
    anchor_repair_scope,
    audit_portfolio_anchor_repair,
    repair_portfolio_anchors,
)
from dvbfixer.model.diffusion.portfolio_contract import (
    PortfolioAnchorArtifacts,
    PortfolioAnchoredEvidence,
    PortfolioEvidence,
    decode_portfolio_evidence,
)
from dvbfixer.model.diffusion.portfolio_evidence import (
    inspect_portfolio_anchor_inputs,
    json_bytes,
    portfolio_resume_key,
    strict_json,
)
from dvbfixer.model.diffusion.portfolio_publication import (
    publish_portfolio_bundle,
    resume_portfolio_bundle,
)
from tests.test_diffusion_anchor_acceptance import shift
from tests.test_diffusion_gap_portfolio import two_gaps
from tests.test_diffusion_portfolio_publication import setup_evidence


def frames():
    request, source, good = two_gaps()
    selected = shift(good, "65", 1.5, atom_name="C")
    return request, source, good, selected, shift(good, "66", 0.04)


def test_repair_restarts_from_assembly_not_failed_local_coordinates(monkeypatch):
    request, source, good, selected, final = frames()
    stage1 = shift(selected, "65", 0.001, atom_name="N")
    indices, residues, atoms, anchors = anchor_repair_scope(request, source, selected, stage1)
    assert indices == (0,)
    assert ResidueIdentity("C", "66") in anchors
    captured = []

    def kernel(text, **kwargs):
        assert text == selected and text != stage1
        assert kwargs["anchor_movement_limits"] == AnchorMovementLimits()
        assert kwargs["generated_atoms"] == atoms
        assert all(a.residue_number == "65" for a in atoms)
        allowed = set(atoms) | {r.atom for r in kwargs["observed_anchor_restraints"]}
        xyz = {a: x for a, x in _coordinates(final)[0].items() if a in allowed}
        assert kwargs["candidate_score"](xyz) == (0, 0, 0)
        return SimpleNamespace(coordinates_angstrom=xyz, platform="Reference", search_diagnostics=())

    monkeypatch.setattr("dvbfixer.model.diffusion.boundary_refinement.refine_generated_region", kernel)
    result = repair_portfolio_anchors(request, source, selected, stage1, refinement_seed=7,
        retain_trial=lambda *args: captured.append(args))
    assert result.coordinate_text == final
    assert result.report["audit"]["summary"]["passed"]
    assert result.report["input_sha256"] == digest(selected)
    assert result.report["stage1_sha256"] == digest(stage1)
    assert not result.report["new_sampling"] and captured
    # A second clean generated gap remains byte-identical.
    assert {a: line for a, line in _coordinates(final)[1].items() if a.residue_number == "69"} == {
        a: line for a, line in _coordinates(selected)[1].items() if a.residue_number == "69"}


@pytest.mark.parametrize("mutation", ["clean-gap", "unselected-observation", "header", "metadata", "stage1-observation"])
def test_masks_and_noncoordinate_records_are_independently_enforced(mutation):
    request, source, good, selected, final = frames()
    stage1 = selected
    if mutation == "clean-gap":
        final = shift(final, "69", 0.001, atom_name="CA")
    elif mutation == "unselected-observation":
        final = shift(final, "67", 0.001, atom_name="CA")
    elif mutation == "header":
        final = "REMARK changed\n" + final
    elif mutation == "metadata":
        final = final.replace("  1.00", "  0.90", 1)
    else:
        stage1 = shift(stage1, "66", 0.001, atom_name="CA")
    with pytest.raises(DiffusionContractError):
        audit_portfolio_anchor_repair(request, source, selected, stage1, final)


def test_anchor_limits_remain_original_frame_per_anchor_not_aligned():
    request, source, good, selected, final = frames()
    candidate = shift(good, "66", 0.6)
    audit = audit_portfolio_anchor_repair(request, source, selected, selected, candidate)
    assert not audit.summary.passed
    assert "anchor-repair-per-anchor-movement-limit" in audit.summary.hard_gate_failures
    assert audit.movement_limits == AnchorMovementLimits()


@pytest.mark.parametrize("kwargs", [{"refinement_seed": True}, {"max_iterations": 0}, {"restart_count": 9}])
def test_execution_budgets_fail_before_kernel(kwargs):
    request, source, good, selected, final = frames()
    with pytest.raises(DiffusionContractError, match="budget"):
        repair_portfolio_anchors(request, source, selected, selected, **({"refinement_seed": 7} | kwargs))


def test_clean_assembly_or_successful_generated_only_attempt_never_authorizes_anchors():
    request, source, good, selected, final = frames()
    for raw, first in ((good, selected), (selected, good)):
        with pytest.raises(DiffusionContractError, match="unnecessary"):
            anchor_repair_scope(request, source, raw, first)


def test_preexisting_source_findings_cannot_authorize_observed_motion(monkeypatch):
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    request, source, good, selected, final = frames()
    findings = assessment_adapter(request, source)(selected).issues
    monkeypatch.setattr("dvbfixer.model.diffusion.portfolio_anchors._findings", lambda *args: findings)
    with pytest.raises(DiffusionContractError, match="justified anchors"):
        anchor_repair_scope(request, source, selected, selected)


def test_global_d_ca_is_never_accepted_by_anchor_motion():
    import numpy as np

    from dvbfixer.ffutils.geometry import ChiralityError
    from dvbfixer.model.diffusion.contract import AtomIdentity
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates

    request, source, good, selected, final = frames()
    xyz = _coordinates(final)[0]
    n, ca, c, cb = (xyz[AtomIdentity("C", "65", "", name)] for name in ("N", "CA", "C", "CB"))
    normal = np.cross(n - ca, c - ca)
    reflected = cb - 2 * np.dot(cb - ca, normal) / np.dot(normal, normal) * normal
    text = rewrite_generated_coordinates(final, {AtomIdentity("C", "65", "", "CB"): reflected})
    with pytest.raises(ChiralityError):
        audit_portfolio_anchor_repair(request, source, selected, selected, text)


def setup_anchored(root, monkeypatch):
    base, request = setup_evidence(root, repair=True)
    selected = (root / "selected.pdb").read_text()
    good = (root / "candidate.pdb").read_text()
    final = shift(good, "66", 0.04)
    # Retained failed generated-only attempt, not a synthetic success receipt.
    (root / "candidate.pdb").write_text(selected)
    operation = strict_json((root / "repair.json").read_bytes())
    operation["output_sha256"] = digest(selected)
    (root / "repair.json").write_bytes(json_bytes(operation))
    base = replace(base, candidate=replace(base.candidate, sha256=digest(selected)),
        repair_report=replace(base.repair_report, sha256=digest(json_bytes(operation).decode())))
    inputs = inspect_portfolio_anchor_inputs(base, root)
    assert not inputs.final_assessment.geometry_passed

    def kernel(text, **kwargs):
        allowed = set(kwargs["generated_atoms"]) | {r.atom for r in kwargs["observed_anchor_restraints"]}
        xyz = {a: x for a, x in _coordinates(final)[0].items() if a in allowed}
        kwargs["candidate_score"](xyz)
        return SimpleNamespace(coordinates_angstrom=xyz, platform="Reference", search_diagnostics=())

    monkeypatch.setattr("dvbfixer.model.diffusion.boundary_refinement.refine_generated_region", kernel)
    result = repair_portfolio_anchors(request, (root / "source.pdb").read_text(), selected, selected, refinement_seed=7)

    def write(name, data):
        (root / name).write_bytes(data)
        return ArtifactReference(name, digest(data.decode()))

    report = write("anchor-report.json", json_bytes(result.report))
    candidate = write("anchored.pdb", result.coordinate_text.encode())
    sources = {name: write(name + ".py", b"# unit control, not sampling evidence\n") for name in
        ("adapter_source", "kernel_source", "conditioning_source", "acceptance_source", "policy_source")}
    values = base.to_dict()
    values.update(schema_version=2, candidate=asdict(candidate), anchor=asdict(
        PortfolioAnchorArtifacts(base.candidate, report, **sources)))
    return PortfolioAnchoredEvidence.from_dict(values), base, request


def test_separate_schema_publication_dat_and_readonly_resume(tmp_path, monkeypatch):
    root = tmp_path / "run"
    evidence, base, request = setup_anchored(root, monkeypatch)
    bundle = tmp_path / "bundle"
    with pytest.raises(DiffusionContractError):
        PortfolioEvidence.from_json(evidence.to_json())
    assert isinstance(decode_portfolio_evidence(evidence.to_dict()), PortfolioAnchoredEvidence)
    assert portfolio_resume_key(base) != portfolio_resume_key(evidence)
    publish_portfolio_bundle(evidence, root, bundle)
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in bundle.iterdir()}

    def forbidden(*args, **kwargs):
        pytest.fail("resume must not run a refinement")

    monkeypatch.setattr("dvbfixer.model.diffusion.boundary_refinement.refine_generated_region", forbidden)
    published = resume_portfolio_bundle(bundle, evidence)
    assert published.audit.anchor_audit["summary"]["passed"]
    assert published.dat.added_keys() == {(a.chain, a.residue_number, a.insertion_code, a.atom_name)
                                         for a in request.generated_atoms}
    assert not published.dat.added_keys() & {(a.chain, a.residue_number, a.insertion_code, a.atom_name)
                                           for a in request.fixed_atoms}
    assert before == {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in bundle.iterdir()}
    with pytest.raises(DiffusionContractError):
        inspect_portfolio_anchor_inputs(evidence, root)


@pytest.mark.parametrize("mutation", ["mask", "limits", "restart", "receipt", "stage1", "donor"])
def test_rehashed_anchor_claims_never_authorize_expanded_mobility(tmp_path, monkeypatch, mutation):
    root = tmp_path / "run"
    evidence, base, request = setup_anchored(root, monkeypatch)
    operation = strict_json((root / "anchor-report.json").read_bytes())
    if mutation == "mask":
        operation["mobile_generated_atoms"].extend([asdict(a) for a in request.generated_atoms if a.residue_number == "69"])
    elif mutation == "limits":
        operation["movement_limits"]["ca_displacement_angstrom_max"] = 1.0
    elif mutation == "restart":
        operation["input_sha256"] = evidence.anchor.stage1.sha256 + "0"
    elif mutation == "receipt":
        operation["audit"]["movements"][0]["ca_displacement_angstrom"] = 0
    elif mutation == "stage1":
        text = shift((root / evidence.anchor.stage1.path).read_text(), "66", 0.001, atom_name="CA")
        (root / evidence.anchor.stage1.path).write_text(text)
        evidence = replace(evidence, anchor=replace(evidence.anchor, stage1=replace(evidence.anchor.stage1, sha256=digest(text))))
    else:
        evidence = replace(evidence, parents=(replace(evidence.parents[0], seed=999), *evidence.parents[1:]))
    data = json_bytes(operation)
    (root / "anchor-report.json").write_bytes(data)
    evidence = replace(evidence, anchor=replace(evidence.anchor, report=replace(evidence.anchor.report, sha256=digest(data.decode()))))
    with pytest.raises((DiffusionContractError, ValueError)):
        publish_portfolio_bundle(evidence, root, tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()
