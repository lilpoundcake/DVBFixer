"""Tests for the versioned backend-neutral diffusion request/result contract."""

from __future__ import annotations

import json

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
    RunnerResult,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
    ValidationSummary,
)

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


def test_backend_provenance_and_candidates_reject_invalid_seed_metadata() -> None:
    with pytest.raises(DiffusionContractError, match="candidate seed"):
        RunnerCandidate(
            candidate_id="candidate",
            seed=-1,
            coordinate_artifact=ArtifactReference("candidate.pdb", DIGEST),
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
