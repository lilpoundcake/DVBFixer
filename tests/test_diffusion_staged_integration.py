"""Versioned runner -> independent validation -> publication/resume regressions."""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from types import SimpleNamespace

import pytest

pytest.importorskip("openmm")
pytest.importorskip("pdbfixer")

from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.boundary_refinement import BoundaryNumericalError
from dvbfixer.model.diffusion.contract import ArtifactReference, SamplerTrace
from dvbfixer.model.diffusion.pipeline import publish_diffusion_bundle
from dvbfixer.model.diffusion.publication import PublicationError, read_published_bundle
from dvbfixer.model.diffusion.refinement_evidence import (
    StagedRefinementEvidence,
    bind_evidence,
    evidence_reference,
)
from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy
from dvbfixer.model.diffusion.runner_refinement import refine_runner_result
from dvbfixer.model.diffusion.validate import (
    DiffusionValidationError,
    build_validated_result,
    validate_runner_result,
)
from tests.test_diffusion_anchor_acceptance import shift
from tests.test_diffusion_validate import (
    _candidate,
    _request,
    _runner_result,
    _source_and_candidate,
    _trace_bytes,
)


def setup_run(tmp_path, *, raw_bad=False, relax=True):
    source, good, fixed, generated = _source_and_candidate()
    request = _request(source.encode(), fixed, generated)
    request = replace(request, backend_options=RefinementPolicy(relax_anchors=relax).backend_options())
    (tmp_path / request.normalized_pdb.path).parent.mkdir(parents=True, exist_ok=True)
    (tmp_path / request.normalized_pdb.path).write_text(source)
    request_path = tmp_path / "request.json"
    request_path.write_text(request.to_json())
    raw = shift(good, "69", 1.5, atom_name="C") if raw_bad else good
    (tmp_path / "raw.pdb").write_text(raw)
    (tmp_path / "sampler-trace.json").write_bytes(_trace_bytes(fixed, generated))
    result = _runner_result((_candidate("raw", "raw.pdb", raw.encode(), generated, score=0.1),))
    return request_path, request, result, raw, good


def install_minimizer(monkeypatch, good, calls, *, stage2=False, numerical=False):
    def minimize(text, **kwargs):
        calls.append((text, kwargs))
        if stage2 and len(calls) == 1:
            if numerical:
                raise BoundaryNumericalError("No finite all-L minimum")
            selected = _structure(text, "raw.pdb")
        else:
            selected = _structure(shift(good, "70", 0.04) if stage2 else good, "selected.pdb")
        atoms = set(kwargs["generated_atoms"]) | {r.atom for r in kwargs["observed_anchor_restraints"]}
        return SimpleNamespace(platform="Reference", coordinates_angstrom={
            atom: selected.atoms[atom].coordinates_angstrom for atom in atoms})
    monkeypatch.setattr("dvbfixer.model.diffusion.runner_refinement.refine_generated_region", minimize)


@pytest.mark.parametrize("stage", ["raw", "generated-only", "generated-and-anchors", "numerical"])
def test_runner_public_validation_and_resume_preserve_staged_custody(tmp_path, monkeypatch, stage):
    request_path, request, result, raw, good = setup_run(tmp_path, raw_bad=stage != "raw")
    calls = []
    install_minimizer(monkeypatch, good, calls, stage2=stage in {"generated-and-anchors", "numerical"},
                      numerical=stage == "numerical")
    refined = refine_runner_result(request_path, result, tmp_path / "refined", platform_name="Reference")
    expected = "generated-and-anchors" if stage == "numerical" else stage
    assert len(calls) == {"raw": 0, "generated-only": 1, "generated-and-anchors": 2}[expected]
    assert all(text == raw for text, _ in calls)
    if len(calls) == 2:
        assert calls[1][1]["observed_anchor_restraints"]
    validated = build_validated_result(request, refined, workspace=tmp_path)
    assert validated.validation_summaries[0].passed, validated.validation_summaries[0].hard_gate_failures
    assert validated.candidates[0].generated_atoms == request.generated_atoms
    assert any("Observed anchor C:70" in w for w in validated.validation_summaries[0].warnings) == (len(calls) == 2)
    bundle = tmp_path / "bundle"
    publish_diffusion_bundle(request, validated, workspace=tmp_path, destination_bundle=bundle)
    published = read_published_bundle(bundle, request=request, source_root=tmp_path, work_parent=tmp_path)
    assert published.selected().pdb_bytes == (tmp_path / refined.candidates[0].coordinate_artifact.path).read_bytes()
    trace_path = next(bundle.glob("*.sampler-trace.json"))
    reference = evidence_reference(SamplerTrace.from_json(trace_path.read_text()))
    evidence = StagedRefinementEvidence.from_json((bundle / reference.path).read_text())
    assert evidence.selected_stage == expected
    assert (bundle / evidence.report.path).is_file()
    assert evidence.raw.path.startswith("raw-staged.")
    (bundle / evidence.raw.path).write_text("REMARK tampered\n")
    with pytest.raises(PublicationError, match="SHA-256|digest"):
        read_published_bundle(bundle, request=request, source_root=tmp_path, work_parent=tmp_path)


def test_optout_never_attempts_anchor_motion(tmp_path, monkeypatch):
    request_path, request, result, _, good = setup_run(tmp_path, raw_bad=True, relax=False)
    calls = []
    install_minimizer(monkeypatch, good, calls, stage2=True)
    refined = refine_runner_result(request_path, result, tmp_path / "refined", platform_name="Reference")
    assert len(calls) == 1
    validated = build_validated_result(request, refined, workspace=tmp_path)
    assert not validated.validation_summaries[0].passed


def test_missing_or_misbound_evidence_cannot_waive_strict_preservation(tmp_path, monkeypatch):
    request_path, request, result, _, good = setup_run(tmp_path)
    with pytest.raises(DiffusionValidationError, match="required"):
        build_validated_result(request, result, workspace=tmp_path)
    install_minimizer(monkeypatch, good, [])
    refined = refine_runner_result(request_path, result, tmp_path / "refined", platform_name="Reference")
    trace_path = tmp_path / refined.candidates[0].sampler_trace_artifact.path
    reference = evidence_reference(SamplerTrace.from_json(trace_path.read_text()))
    path = tmp_path / reference.path
    record = json.loads(path.read_text())
    record["request_sha256"] = hashlib.sha256(b"wrong request").hexdigest()
    path.write_text(json.dumps(record))
    with pytest.raises(DiffusionValidationError, match="SHA-256"):
        build_validated_result(request, refined, workspace=tmp_path)


def test_archive_replay_cli_is_explicitly_not_neural_qualification(tmp_path, monkeypatch):
    from scripts.replay_staged_refinement import main

    request_path, request, _, _, _ = setup_run(tmp_path)
    request = replace(request, seeds=(request.seeds[0],), candidate_count=1)
    request_path.write_text(request.to_json())
    output = tmp_path / "replay"
    monkeypatch.setattr("sys.argv", ["replay", "--request", str(request_path),
        "--source", str(tmp_path / request.normalized_pdb.path),
        "--raw", str(tmp_path / "raw.pdb"), "--output", str(output)])
    with pytest.raises(SystemExit) as exited:
        main()
    assert exited.value.code == 0
    report = json.loads((output / "report.json").read_text())
    assert report["status"] == "success" and report["sampler_qualification"] is False


def test_model_zbs_homology_expose_independent_default_enabled_switches():
    from dvbfixer.homology import parse_args as homology_args
    from dvbfixer.model.cli import parse_args as model_args
    from dvbfixer.zbs import parse_args as zbs_args

    options = ["--no-relax-anchors", "--no-anchor-bfactor-weighting"]
    for parse, base in (
        (model_args, ["input.pdb", "--backend", "diffusion"]),
        (zbs_args, ["input.pdb", "--model-backend", "diffusion"]),
        (homology_args, ["target.fasta", "--template-plan", "plan.json", "--backend", "diffusion"]),
    ):
        default = parse(base)
        assert not default.no_relax_anchors and not default.no_anchor_bfactor_weighting
        disabled = parse([*base, *options])
        assert disabled.no_relax_anchors and disabled.no_anchor_bfactor_weighting
    with pytest.raises(SystemExit):
        model_args(["input.pdb", "--backend", "diffusion", "--diffusion-request", "request.json", *options])


def test_seed_scoped_runner_evidence_is_verified_and_rebound_for_publication(tmp_path, monkeypatch):
    request_path, parent, result, _, good = setup_run(tmp_path)
    single = replace(parent, seeds=(result.candidates[0].seed,), candidate_count=1)
    request_path.write_text(single.to_json())
    install_minimizer(monkeypatch, good, [])
    refined = refine_runner_result(request_path, result, tmp_path / "refined", platform_name="Reference")
    validated = build_validated_result(parent, refined, workspace=tmp_path)
    assert validated.validation_summaries[0].passed
    bundle = tmp_path / "bundle"
    publish_diffusion_bundle(parent, validated, workspace=tmp_path, destination_bundle=bundle)
    assert read_published_bundle(bundle, request=parent, source_root=tmp_path, work_parent=tmp_path).status == "success"
    published = StagedRefinementEvidence.from_json(next(bundle.glob("*.refinement-evidence.json")).read_text())
    assert published.request_sha256 == hashlib.sha256(parent.to_json().encode()).hexdigest()


def test_staged_runner_cannot_relabel_observed_anchors_as_generated(tmp_path, monkeypatch):
    request_path, request, result, _, good = setup_run(tmp_path)
    install_minimizer(monkeypatch, good, [])
    refined = refine_runner_result(request_path, result, tmp_path / "refined", platform_name="Reference")
    candidate = refined.candidates[0]
    forged = replace(refined, candidates=(replace(candidate,
        generated_atoms=(*candidate.generated_atoms, request.fixed_atoms[0])),))
    with pytest.raises(DiffusionValidationError, match="generated mask"):
        build_validated_result(request, forged, workspace=tmp_path)


@pytest.mark.parametrize("raw_bad", [False, True])
def test_sampler_raw_diagnostics_do_not_require_future_refinement_evidence(tmp_path, raw_bad):
    _, request, result, _, _ = setup_run(tmp_path, raw_bad=raw_bad)
    diagnostic = validate_runner_result(request, result, workspace=tmp_path, phase="raw")[0]
    assert diagnostic.summary.passed is not raw_bad
    if raw_bad:
        assert "junction-peptide-connectivity" in diagnostic.summary.hard_gate_failures
    # The diagnostic-only phase does not change the default publication gate.
    with pytest.raises(DiffusionValidationError, match="required"):
        build_validated_result(request, result, workspace=tmp_path)


def test_raw_diagnostics_never_waive_observed_motion_or_accept_refined_traces(tmp_path, monkeypatch):
    _, request, result, _, good = setup_run(tmp_path)
    moved = shift(good, "70", 0.04)
    (tmp_path / "raw.pdb").write_text(moved)
    candidate = result.candidates[0]
    moved_result = replace(result, candidates=(replace(candidate,
        coordinate_artifact=replace(candidate.coordinate_artifact, sha256=hashlib.sha256(moved.encode()).hexdigest())),))
    diagnostic = validate_runner_result(request, moved_result, workspace=tmp_path, phase="raw")[0]
    assert not diagnostic.summary.passed
    assert "fixed-heavy-atom-max-displacement" in diagnostic.summary.hard_gate_failures
    request_path, request, result, _, good = setup_run(tmp_path / "second")
    install_minimizer(monkeypatch, good, [])
    refined = refine_runner_result(request_path, result, request_path.parent / "refined", platform_name="Reference")
    with pytest.raises(DiffusionValidationError, match="post-refinement"):
        validate_runner_result(request, refined, workspace=request_path.parent, phase="raw")


def test_rehashed_runner_claim_cannot_authorize_an_unselected_observed_atom(tmp_path, monkeypatch):
    request_path, request, result, _, good = setup_run(tmp_path, raw_bad=True)
    install_minimizer(monkeypatch, good, [], stage2=True)
    refined = refine_runner_result(request_path, result, tmp_path / "refined", platform_name="Reference")
    candidate = refined.candidates[0]
    selected_path = tmp_path / candidate.coordinate_artifact.path
    selected_path.write_text(shift(selected_path.read_text(), "64", 0.001, atom_name="CA"))
    selected = ArtifactReference(candidate.coordinate_artifact.path,
        hashlib.sha256(selected_path.read_bytes()).hexdigest())
    trace_path = tmp_path / candidate.sampler_trace_artifact.path
    trace = SamplerTrace.from_json(trace_path.read_text())
    reference = evidence_reference(trace)
    path = tmp_path / reference.path
    record = replace(StagedRefinementEvidence.from_json(path.read_text()), selected=selected)
    path.write_text(record.to_json())
    trace = bind_evidence(trace, ArtifactReference(reference.path, hashlib.sha256(path.read_bytes()).hexdigest()))
    trace_path.write_text(trace.to_json())
    forged = replace(refined, candidates=(replace(candidate, coordinate_artifact=selected,
        sampler_trace_artifact=ArtifactReference(candidate.sampler_trace_artifact.path,
            hashlib.sha256(trace_path.read_bytes()).hexdigest())),))
    validated = build_validated_result(request, forged, workspace=tmp_path)
    assert not validated.validation_summaries[0].passed
    assert "anchor-repair-unselected-observed-record-changed" in validated.validation_summaries[0].hard_gate_failures
