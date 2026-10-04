"""Tests for the versioned backend-neutral diffusion request/result contract."""

from __future__ import annotations

import json
from dataclasses import replace

import pytest

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    BackendOption,
    BackendProvenance,
    DiffusionCandidate,
    DiffusionContractError,
    DiffusionRequest,
    DiffusionResult,
    DiffusionStatus,
    GapRegion,
    Metric,
    ResidueIdentity,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResourceMetrics,
    RunnerResult,
    SamplerTrace,
    SamplingAblationMode,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
    ValidationSummary,
)
from dvbfixer.model.diffusion.trace import build_sampler_trace, validate_sampler_trace_context

DIGEST = "a" * 64


def _request() -> DiffusionRequest:
    residues = tuple(ResidueIdentity("D", str(number)) for number in range(10, 17))
    generated_residues = residues[1:6]
    fixed_atoms = (
        AtomIdentity("D", "10", "", "CA"),
        AtomIdentity("D", "16", "", "CA"),
    )
    generated_atoms = tuple(
        AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, "CA")
        for residue in generated_residues
    )
    return DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("input/normalized.pdb", DIGEST),
        target_sequences=(TargetSequence("D", "FSGSKSG"),),
        sequence_placements=(
            SequencePlacement(
                chain="D",
                target_length=7,
                observed_target_indices=(0, 6),
                observed_residues=(residues[0], residues[-1]),
            ),
        ),
        gaps=(
            GapRegion(
                chain="D",
                target_interval=TargetInterval(1, 6),
                left_anchor=residues[0],
                right_anchor=residues[-1],
                generated_residues=generated_residues,
                movable_junction_residues=(residues[0], *generated_residues, residues[-1]),
            ),
        ),
        fixed_atoms=fixed_atoms,
        generated_atoms=generated_atoms,
        retained_explicit_links=(),
        candidate_count=2,
        seeds=(7, 11),
        backend_options=(BackendOption("precision", "float32"),),
    )


def _success_result() -> DiffusionResult:
    candidate = DiffusionCandidate(
        candidate_id="candidate-0001",
        seed=7,
        coordinate_artifact=ArtifactReference("candidates/candidate-0001.pdb", DIGEST),
        sampler_trace_artifact=ArtifactReference("candidates/candidate-0001.trace.json", DIGEST),
        generated_atoms=(AtomIdentity("D", "11", "", "CA"),),
        generated_residues=(ResidueIdentity("D", "11"),),
        raw_backend_score=0.75,
        score_provenance="fake-runner:test-score-v1",
    )
    return DiffusionResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=(candidate,),
        validation_summaries=(
            ValidationSummary(
                passed=True,
                metrics=(Metric("fixed-heavy-atom-rmsd", 0.0, "angstrom"),),
            ),
        ),
        runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="fake",
            runner_protocol_version=2,
            engine_repository="https://example.invalid/fake",
            engine_revision="test-revision",
        ),
    )


def test_request_round_trip_is_strict_and_deterministic() -> None:
    request = _request()

    encoded = request.to_json()
    decoded = DiffusionRequest.from_json(encoded)

    assert decoded == request
    assert decoded.to_json() == encoded
    assert json.loads(encoded)["schema_version"] == DIFFUSION_SCHEMA_VERSION
    assert decoded.target_sequences[0].chain == "D"


def test_runner_result_round_trip_has_no_external_validation_claims() -> None:
    request = _request()
    candidate = RunnerCandidate(
        candidate_id="candidate-0001",
        seed=7,
        coordinate_artifact=ArtifactReference("candidates/candidate-0001.pdb", DIGEST),
        sampler_trace_artifact=ArtifactReference("candidates/candidate-0001.trace.json", DIGEST),
        generated_atoms=request.generated_atoms,
        generated_residues=request.gaps[0].generated_residues,
        raw_backend_score=0.75,
        score_provenance="fake-runner:test-score-v1",
    )
    result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=(candidate,),
        runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="fake",
            runner_protocol_version=2,
            engine_repository="https://example.invalid/fake",
            engine_revision="test-revision",
        ),
    )

    decoded = RunnerResult.from_json(result.to_json())

    assert decoded == result
    assert "validation_summaries" not in decoded.to_dict()

    trace_free = result.to_dict()
    del trace_free["candidates"][0]["sampler_trace_artifact"]
    with pytest.raises(DiffusionContractError, match="sampler_trace_artifact"):
        RunnerResult.from_dict(trace_free)


def test_result_round_trip_preserves_status_and_nested_types() -> None:
    result = _success_result()

    decoded = DiffusionResult.from_json(result.to_json())

    assert decoded == result
    assert decoded.status is DiffusionStatus.SUCCESS
    assert decoded.candidates[0].seed == 7
    assert isinstance(decoded.candidates[0].generated_atoms[0], AtomIdentity)
    assert isinstance(decoded.validation_summaries[0].metrics[0], Metric)


def test_backend_provenance_round_trip_preserves_optional_runtime_evidence() -> None:
    provenance = BackendProvenance(
        backend="test-engine",
        runner_protocol_version=2,
        engine_repository="https://example.invalid/engine",
        engine_revision="pinned-revision",
        source_license="BSD-3-Clause",
        checkpoint_sha256="b" * 64,
        checkpoint_license="upstream-review-required",
        container_digest="sha256:" + "c" * 64,
        environment_hash="d" * 64,
        environment_identity="lockfile:environment.lock",
        device="cuda:0",
        precision="float32",
        framework="torch",
        framework_version="2.5.1",
        cuda_version="12.4",
        driver_version="550.54",
        deterministic_algorithms=False,
        deterministic_flags=("CUBLAS_WORKSPACE_CONFIG=:4096:8",),
        known_nondeterministic_operations=("scatter_add",),
    )

    result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.FAILED,
        candidates=(),
        runner_diagnostics=RunnerDiagnostics(exit_code=1, timed_out=False),
        backend_provenance=provenance,
        message="failed",
    )
    round_trip = RunnerResult.from_json(result.to_json())

    assert round_trip.backend_provenance == provenance
    assert round_trip.backend_provenance.deterministic_algorithms is False


def test_schema_v4_resource_metrics_and_sampler_trace_are_strict() -> None:
    metrics = RunnerResourceMetrics(
        wall_time_seconds=12.5,
        model_load_seconds=3.0,
        peak_ram_bytes=1024,
        peak_vram_bytes=2048,
    )
    result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.FAILED,
        candidates=(),
        runner_diagnostics=RunnerDiagnostics(exit_code=1, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="fake",
            runner_protocol_version=4,
            engine_repository="builtin://fake",
            engine_revision="test",
        ),
        resource_metrics=metrics,
        message="failed",
    )

    assert RunnerResult.from_json(result.to_json()).resource_metrics == metrics
    with pytest.raises(DiffusionContractError, match="finite and non-negative"):
        RunnerResourceMetrics(wall_time_seconds=float("nan"))
    with pytest.raises(DiffusionContractError, match="non-negative"):
        RunnerResourceMetrics(peak_vram_bytes=-1)

    fixed = AtomIdentity("D", "10", "", "CA")
    generated = AtomIdentity("D", "11", "", "CA")
    trace = build_sampler_trace(
        profile="test",
        engine_repository="builtin://fake",
        engine_revision="test",
        patch_identity="none",
        atom_order=(fixed, generated),
        fixed_atoms=(fixed,),
        device="cpu",
        fallback_disabled=True,
        denoising_update_count=1,
        projection_errors_angstrom=(0.0,),
        final_fixed_coordinate_restoration=False,
    )
    assert trace.ablation_mode is SamplingAblationMode.REINJECTION
    assert trace.represented_fixed_atoms == (fixed,)
    assert trace.sampler_evidence_complete
    assert trace.fixed_tolerance_angstrom == 0.01
    assert trace.max_projection_error_angstrom == 0.0
    assert SamplerTrace.from_json(trace.to_json()) == trace

    with pytest.raises(DiffusionContractError, match="exceeds fixed tolerance"):
        build_sampler_trace(
            profile="test",
            engine_repository="builtin://fake",
            engine_revision="test",
            patch_identity="none",
            atom_order=(fixed, generated),
            fixed_atoms=(fixed,),
            device="cpu",
            fallback_disabled=True,
            denoising_update_count=1,
            projection_errors_angstrom=(0.02,),
            final_fixed_coordinate_restoration=False,
        )
    relaxed = build_sampler_trace(
        profile="test",
        engine_repository="builtin://fake",
        engine_revision="test",
        patch_identity="none",
        atom_order=(fixed, generated),
        fixed_atoms=(fixed,),
        device="cpu",
        fallback_disabled=True,
        denoising_update_count=1,
        projection_errors_angstrom=(0.02,),
        final_fixed_coordinate_restoration=False,
        fixed_tolerance_angstrom=0.03,
    )
    assert relaxed.fixed_tolerance_angstrom == 0.03

    with pytest.raises(DiffusionContractError, match="cannot claim"):
        replace(trace, ablation_mode=SamplingAblationMode.TEMPLATE_ONLY)

    unknown = trace.to_dict()
    unknown["future_field"] = True
    with pytest.raises(DiffusionContractError, match="unknown fields"):
        SamplerTrace.from_dict(unknown)

    missing = trace.to_dict()
    del missing["device"]
    with pytest.raises(DiffusionContractError, match="missing fields"):
        SamplerTrace.from_dict(missing)


def test_sampler_trace_enforces_production_profile_semantics() -> None:
    fixed = AtomIdentity("D", "10", "", "CA")
    generated = AtomIdentity("D", "11", "", "CA")
    common = {
        "engine_repository": "builtin://test",
        "engine_revision": "test",
        "patch_identity": "none",
        "atom_order": (fixed, generated),
        "fixed_atoms": (fixed,),
        "fallback_disabled": True,
    }

    with pytest.raises(DiffusionContractError, match="per-step reinjection"):
        build_sampler_trace(
            profile="protpardelle-1c-mps",
            device="mps",
            denoising_update_count=1,
            projection_errors_angstrom=(0.0,),
            final_fixed_coordinate_restoration=True,
            **common,
        )
    with pytest.raises(DiffusionContractError, match="final fixed-coordinate"):
        build_sampler_trace(
            profile="protpardelle-1c-mps",
            device="mps",
            denoising_update_count=500,
            final_fixed_coordinate_restoration=False,
            **common,
        )
    with pytest.raises(DiffusionContractError, match="reinjection capabilities"):
        build_sampler_trace(
            profile="protenix-v1-cuda",
            device="cuda:0",
            denoising_update_count=1,
            final_fixed_coordinate_restoration=False,
            **common,
        )
    apple = build_sampler_trace(
        profile="protpardelle-1c-mps",
        device="mps",
        denoising_update_count=500,
        final_fixed_coordinate_restoration=True,
        **common,
    )
    assert apple.denoising_update_count == 500
    assert apple.callback_update_count == 0
    assert apple.steps == ()
    assert apple.max_projection_error_angstrom is None
    passthrough_fixed = AtomIdentity("D", "12", "", "OXT")
    apple_with_passthrough = build_sampler_trace(
        profile="protpardelle-1c-mps",
        engine_repository="builtin://test",
        engine_revision="test",
        patch_identity="none",
        atom_order=(fixed, generated),
        fixed_atoms=(fixed, passthrough_fixed),
        represented_fixed_atoms=(fixed,),
        device="mps",
        fallback_disabled=True,
        denoising_update_count=500,
        final_fixed_coordinate_restoration=True,
    )
    assert apple_with_passthrough.represented_fixed_atoms == (fixed,)
    with pytest.raises(DiffusionContractError, match="exactly 500"):
        build_sampler_trace(
            profile="protpardelle-1c-mps",
            device="mps",
            denoising_update_count=499,
            final_fixed_coordinate_restoration=True,
            **common,
        )

    protenix = build_sampler_trace(
        profile="protenix-v1-cuda",
        device="cuda:0",
        denoising_update_count=200,
        projection_errors_angstrom=(0.0,) * 200,
        final_fixed_coordinate_restoration=False,
        **common,
    )
    assert protenix.denoising_update_count == 200
    assert protenix.callback_update_count == 200
    assert len(protenix.steps) == 200
    with pytest.raises(DiffusionContractError, match="every denoising update"):
        build_sampler_trace(
            profile="protenix-v1-cuda",
            device="cuda:0",
            denoising_update_count=199,
            projection_errors_angstrom=(0.0,) * 199,
            final_fixed_coordinate_restoration=False,
            **common,
        )

    with pytest.raises(DiffusionContractError, match="every fixed atom"):
        build_sampler_trace(
            profile="protenix-v1-cuda",
            device="cuda:0",
            denoising_update_count=200,
            projection_errors_angstrom=(0.0,) * 200,
            represented_fixed_atoms=(),
            final_fixed_coordinate_restoration=False,
            **common,
        )
    incomplete = build_sampler_trace(
        profile="unknown-source-legacy-refinement",
        device="unknown-input-sampler-device",
        denoising_update_count=None,
        represented_fixed_atoms=(),
        final_fixed_coordinate_restoration=False,
        sampler_evidence_complete=False,
        **common,
    )
    assert incomplete.callback_update_count is None
    assert incomplete.denoising_update_count is None
    assert not incomplete.sampler_evidence_complete
    with pytest.raises(DiffusionContractError, match="atom identity digest"):
        replace(
            protenix,
            steps=(
                replace(protenix.steps[0], atom_identity_sha256="0" * 64),
                *protenix.steps[1:],
            ),
        )


def test_sampler_trace_context_requires_provenance_device_match() -> None:
    request = _request()
    trace = build_sampler_trace(
        profile="test",
        engine_repository="builtin://test",
        engine_revision="revision",
        patch_identity="none",
        atom_order=(*request.fixed_atoms, *request.generated_atoms),
        fixed_atoms=request.fixed_atoms,
        device="cpu",
        fallback_disabled=True,
        denoising_update_count=0,
        final_fixed_coordinate_restoration=True,
    )
    provenance = BackendProvenance(
        backend="test",
        runner_protocol_version=4,
        engine_repository="builtin://test",
        engine_revision="revision",
        device="mps",
    )

    with pytest.raises(DiffusionContractError, match="device"):
        validate_sampler_trace_context(trace, request, provenance)


def test_sampler_trace_context_requires_exact_generated_sampler_atoms() -> None:
    request = _request()
    provenance = BackendProvenance(
        backend="test",
        runner_protocol_version=4,
        engine_repository="builtin://test",
        engine_revision="revision",
        device="cpu",
    )
    for atom_order in (
        request.fixed_atoms,
        (*request.fixed_atoms, *request.generated_atoms, AtomIdentity("D", "99", "", "CA")),
    ):
        trace = build_sampler_trace(
            profile="test",
            engine_repository="builtin://test",
            engine_revision="revision",
            patch_identity="none",
            atom_order=atom_order,
            fixed_atoms=request.fixed_atoms,
            device="cpu",
            fallback_disabled=True,
            denoising_update_count=0,
            final_fixed_coordinate_restoration=True,
        )
        with pytest.raises(DiffusionContractError, match="requested generated atoms"):
            validate_sampler_trace_context(trace, request, provenance)


def test_backend_provenance_and_candidates_reject_invalid_seed_metadata() -> None:
    with pytest.raises(DiffusionContractError, match="candidate seed"):
        RunnerCandidate(
            candidate_id="candidate",
            seed=-1,
            coordinate_artifact=ArtifactReference("candidate.pdb", DIGEST),
            sampler_trace_artifact=ArtifactReference("candidate.trace.json", DIGEST),
            generated_atoms=(),
            generated_residues=(),
            raw_backend_score=None,
            score_provenance="test",
        )

    with pytest.raises(DiffusionContractError, match="deterministic_flags"):
        BackendProvenance(
            backend="fake",
            runner_protocol_version=2,
            engine_repository="builtin://fake",
            engine_revision="test",
            deterministic_flags=("flag", "flag"),
        )


def test_case_distinct_chains_and_insertion_codes_remain_distinct() -> None:
    upper = AtomIdentity("D", "82", "", "CA")
    lower = AtomIdentity("d", "82", "", "CA")
    inserted = AtomIdentity("D", "82", "A", "CA")

    assert len({upper, lower, inserted}) == 3
    assert upper.chain == "D"
    assert lower.chain == "d"
    assert inserted.insertion_code == "A"


def test_unknown_schema_version_is_rejected() -> None:
    raw = _request().to_dict()
    raw["schema_version"] = DIFFUSION_SCHEMA_VERSION + 1

    with pytest.raises(DiffusionContractError, match="unsupported schema_version"):
        DiffusionRequest.from_dict(raw)


def test_unknown_and_missing_fields_are_rejected() -> None:
    unknown = _request().to_dict()
    unknown["future_field"] = True
    with pytest.raises(DiffusionContractError, match="unknown fields"):
        DiffusionRequest.from_dict(unknown)

    missing = _request().to_dict()
    del missing["normalized_pdb"]
    with pytest.raises(DiffusionContractError, match="missing fields"):
        DiffusionRequest.from_dict(missing)


def test_nested_unknown_field_is_rejected() -> None:
    raw = _request().to_dict()
    raw["gaps"][0]["silent_policy"] = "guess"

    with pytest.raises(DiffusionContractError, match="GapRegion contains unknown fields"):
        DiffusionRequest.from_dict(raw)


def test_request_rejects_overlapping_masks_and_seed_mismatch() -> None:
    request = _request()
    common = AtomIdentity("D", "10", "", "CA")

    with pytest.raises(DiffusionContractError, match="must be disjoint"):
        DiffusionRequest(
            schema_version=DIFFUSION_SCHEMA_VERSION,
            normalized_pdb=request.normalized_pdb,
            target_sequences=request.target_sequences,
            sequence_placements=request.sequence_placements,
            gaps=request.gaps,
            fixed_atoms=(common,),
            generated_atoms=(common,),
            retained_explicit_links=(),
            candidate_count=1,
            seeds=(7,),
        )

    with pytest.raises(DiffusionContractError, match="seeds length"):
        DiffusionRequest(
            schema_version=DIFFUSION_SCHEMA_VERSION,
            normalized_pdb=request.normalized_pdb,
            target_sequences=request.target_sequences,
            sequence_placements=request.sequence_placements,
            gaps=request.gaps,
            fixed_atoms=request.fixed_atoms,
            generated_atoms=request.generated_atoms,
            retained_explicit_links=(),
            candidate_count=2,
            seeds=(7,),
        )


def test_artifact_paths_must_be_relative_and_contained() -> None:
    with pytest.raises(DiffusionContractError, match="relative and contained"):
        ArtifactReference("../escape.pdb", DIGEST)
    with pytest.raises(DiffusionContractError, match="relative and contained"):
        ArtifactReference("/absolute/output.pdb", DIGEST)
    with pytest.raises(DiffusionContractError, match="relative and contained"):
        ArtifactReference(".", DIGEST)
    with pytest.raises(DiffusionContractError, match="sha256"):
        ArtifactReference("candidate.pdb", "not-a-digest")


def test_result_status_controls_candidate_publication() -> None:
    provenance = BackendProvenance(
        backend="fake",
        runner_protocol_version=2,
        engine_repository="https://example.invalid/fake",
        engine_revision="test-revision",
    )
    diagnostics = RunnerDiagnostics(exit_code=2, timed_out=False, stderr="unsupported")

    unsupported = DiffusionResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.UNSUPPORTED,
        candidates=(),
        validation_summaries=(),
        runner_diagnostics=diagnostics,
        backend_provenance=provenance,
        message="terminal gaps are unsupported",
    )
    assert DiffusionResult.from_json(unsupported.to_json()) == unsupported

    with pytest.raises(DiffusionContractError, match="cannot publish candidates"):
        DiffusionResult(
            schema_version=DIFFUSION_SCHEMA_VERSION,
            status=DiffusionStatus.FAILED,
            candidates=_success_result().candidates,
            validation_summaries=(),
            runner_diagnostics=diagnostics,
            backend_provenance=provenance,
            message="runner failed",
        )
