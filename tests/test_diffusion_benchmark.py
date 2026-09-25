"""Tests for CPU-only repeatability and withheld-coordinate metrics."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pytest

from dvbfixer.model.diffusion.benchmark import (
    BackendCaseResult,
    BenchmarkOutcome,
    DiffusionBenchmarkError,
    ModellerComparison,
    StratumRun,
    assess_same_seed_repeatability,
    candidate_quality,
    compare_paired_backends,
    load_pdb_coordinates,
    pairwise_gap_backbone_rmsds,
    summarize_backend_cases,
    summarize_ensemble,
    summarize_strata,
    wilson_score_interval,
)
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.runner import (
    prepare_diffusion_workspace,
    run_prepared_diffusion_runner,
)


def _coordinates(offset: float = 0.0) -> dict[AtomIdentity, np.ndarray]:
    return {
        AtomIdentity("D", "10", "", "CA"): np.asarray((0.0, 0.0, 0.0)),
        AtomIdentity("D", "11", "", "N"): np.asarray((1.0 + offset, 0.0, 0.0)),
        AtomIdentity("D", "11", "", "CA"): np.asarray((2.0 + offset, 0.0, 0.0)),
        AtomIdentity("D", "11", "", "C"): np.asarray((3.0 + offset, 0.0, 0.0)),
        AtomIdentity("D", "11", "", "O"): np.asarray((4.0 + offset, 0.0, 0.0)),
        AtomIdentity("D", "11", "", "CB"): np.asarray((2.0 + offset, 1.0, 0.0)),
        AtomIdentity("D", "12", "A", "N"): np.asarray((5.0 + offset, 0.0, 0.0)),
        AtomIdentity("D", "12", "A", "CA"): np.asarray((6.0 + offset, 0.0, 0.0)),
        AtomIdentity("D", "12", "A", "C"): np.asarray((7.0 + offset, 0.0, 0.0)),
        AtomIdentity("D", "12", "A", "O"): np.asarray((8.0 + offset, 0.0, 0.0)),
        AtomIdentity("D", "13", "", "CA"): np.asarray((9.0, 0.0, 0.0)),
    }


def _pdb_bytes() -> bytes:
    return (
        b"ATOM      1  CA  ALA D  10       0.000   0.000   0.000  1.00 20.00           C  \n"
        b"ATOM      2  CA  GLY D  11       1.000   0.000   0.000  1.00 20.00           C  \n"
        b"ATOM      3  CA  SER D  12A      2.000   0.000   0.000  1.00 20.00           C  \n"
        b"ATOM      4  CA  ALA D  13       3.000   0.000   0.000  1.00 20.00           C  \n"
        b"END\n"
    )


def _request(input_bytes: bytes) -> DiffusionRequest:
    residues = (
        ResidueIdentity("D", "10"),
        ResidueIdentity("D", "11"),
        ResidueIdentity("D", "12", "A"),
        ResidueIdentity("D", "13"),
    )
    return DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference(
            "input/normalized.pdb",
            hashlib.sha256(input_bytes).hexdigest(),
        ),
        target_sequences=(TargetSequence("D", "AGSA"),),
        sequence_placements=(
            SequencePlacement(
                chain="D",
                target_length=4,
                observed_target_indices=(0, 3),
                observed_residues=(residues[0], residues[-1]),
            ),
        ),
        gaps=(
            GapRegion(
                chain="D",
                target_interval=TargetInterval(1, 3),
                left_anchor=residues[0],
                right_anchor=residues[-1],
                generated_residues=residues[1:3],
                movable_junction_residues=residues,
            ),
        ),
        fixed_atoms=(
            AtomIdentity("D", "10", "", "CA"),
            AtomIdentity("D", "13", "", "CA"),
        ),
        generated_atoms=(
            AtomIdentity("D", "11", "", "CA"),
            AtomIdentity("D", "12", "A", "CA"),
        ),
        retained_explicit_links=(),
        candidate_count=2,
        seeds=(7, 11),
    )


def test_repeatability_classifies_same_seed_by_stable_identity() -> None:
    reference = _coordinates()
    repeated = {key: value.copy() for key, value in reference.items()}
    repeated[AtomIdentity("D", "11", "", "CA")] += np.asarray((0.005, 0.0, 0.0))

    deterministic = assess_same_seed_repeatability(reference, repeated)

    assert deterministic.classification == "deterministic"
    assert deterministic.coordinate_rmsd_angstrom < 0.01

    repeated[AtomIdentity("D", "11", "", "CA")] += np.asarray((0.1, 0.0, 0.0))
    nondeterministic = assess_same_seed_repeatability(reference, repeated)
    assert nondeterministic.classification == "nondeterministic"


def test_repeatability_rejects_identity_mismatch() -> None:
    reference = _coordinates()
    candidate = dict(reference)
    candidate.pop(AtomIdentity("D", "12", "A", "CA"))

    with pytest.raises(DiffusionBenchmarkError, match="identical requested atoms"):
        assess_same_seed_repeatability(reference, candidate)

    candidate = dict(reference)
    candidate[AtomIdentity("d", "99", "", "CA")] = np.asarray((0.0, 0.0, 0.0))
    with pytest.raises(DiffusionBenchmarkError, match="identical requested atoms"):
        assess_same_seed_repeatability(reference, candidate)


def test_candidate_quality_preserves_insertion_code_and_metric_scopes() -> None:
    reference = _coordinates()
    candidate = {key: value.copy() for key, value in reference.items()}
    candidate[AtomIdentity("D", "11", "", "CB")] += np.asarray((1.0, 0.0, 0.0))
    candidate[AtomIdentity("D", "10", "", "CA")] += np.asarray((0.02, 0.0, 0.0))
    generated = {
        ResidueIdentity("D", "11"),
        ResidueIdentity("D", "12", "A"),
    }

    quality = candidate_quality(
        "candidate-1",
        reference,
        candidate,
        generated_residues=generated,
        fixed_atoms={AtomIdentity("D", "10", "", "CA")},
        anchor_residues={
            ResidueIdentity("D", "10"),
            ResidueIdentity("D", "13"),
        },
        closure_passed=True,
        validation_passed=True,
    )

    assert quality.gap_backbone_rmsd_angstrom == 0.0
    assert quality.gap_all_heavy_rmsd_angstrom > 0.0
    assert quality.fixed_heavy_rmsd_angstrom == pytest.approx(0.02)
    assert quality.anchor_heavy_rmsd_angstrom > 0.0


def test_ensemble_reports_top1_oracle_diversity_and_enrichment() -> None:
    reference = _coordinates()
    generated = {
        ResidueIdentity("D", "11"),
        ResidueIdentity("D", "12", "A"),
    }
    first = _coordinates(0.2)
    second = _coordinates(0.1)
    first_quality = candidate_quality(
        "first",
        reference,
        first,
        generated_residues=generated,
        fixed_atoms={AtomIdentity("D", "10", "", "CA")},
        anchor_residues={ResidueIdentity("D", "10"), ResidueIdentity("D", "13")},
        closure_passed=False,
        validation_passed=False,
    )
    second_quality = candidate_quality(
        "second",
        reference,
        second,
        generated_residues=generated,
        fixed_atoms={AtomIdentity("D", "10", "", "CA")},
        anchor_residues={ResidueIdentity("D", "10"), ResidueIdentity("D", "13")},
        closure_passed=True,
        validation_passed=True,
    )
    diversity = pairwise_gap_backbone_rmsds((first, second), generated)

    summary = summarize_ensemble(
        (first_quality, second_quality),
        ("first", "second"),
        diversity,
    )

    assert summary.candidate_count == 2
    assert summary.passing_count == 1
    assert summary.closure_pass_rate == 0.5
    assert summary.validation_pass_rate == 0.5
    assert summary.top1_gap_backbone_rmsd_angstrom > summary.oracle_gap_backbone_rmsd_angstrom
    assert summary.ranking_enrichment_angstrom > 0.0
    assert summary.mean_pairwise_gap_backbone_rmsd_angstrom == pytest.approx(0.1)


def test_fake_runner_repeatability_is_measured_across_independent_workspaces(
    tmp_path: Path,
) -> None:
    input_bytes = _pdb_bytes()
    source_root = tmp_path / "source"
    input_path = source_root / "input" / "normalized.pdb"
    input_path.parent.mkdir(parents=True)
    input_path.write_bytes(input_bytes)
    request = _request(input_bytes)

    first = prepare_diffusion_workspace(
        request,
        source_root=source_root,
        workspace=tmp_path / "first-workspace",
    )
    first_result = run_prepared_diffusion_runner(
        first,
        (sys.executable, "-m", "dvbfixer.model.diffusion.fake_runner"),
    )
    second = prepare_diffusion_workspace(
        request,
        source_root=source_root,
        workspace=tmp_path / "second-workspace",
    )
    second_result = run_prepared_diffusion_runner(
        second,
        (sys.executable, "-m", "dvbfixer.model.diffusion.fake_runner"),
    )

    first_coordinates = load_pdb_coordinates(
        first.workspace,
        first_result.candidates[0].coordinate_artifact,
    )
    second_coordinates = load_pdb_coordinates(
        second.workspace,
        second_result.candidates[0].coordinate_artifact,
    )
    assessment = assess_same_seed_repeatability(
        first_coordinates,
        second_coordinates,
        atom_identities=set(request.generated_atoms),
    )

    assert first_result.candidates[0].seed == second_result.candidates[0].seed == 7
    assert assessment.classification == "deterministic"
    assert assessment.coordinate_rmsd_angstrom == 0.0
    assert assessment.maximum_displacement_angstrom == 0.0
    assert assessment.atom_count == 2


def test_stratum_summary_reports_outcomes_resources_timeouts_and_crashes() -> None:
    summaries = summarize_strata(
        (
            StratumRun(
                "short-gap",
                BenchmarkOutcome.SUCCESS,
                wall_time_seconds=4.0,
                model_load_seconds=1.0,
                peak_ram_bytes=100,
                peak_vram_bytes=200,
            ),
            StratumRun(
                "short-gap",
                BenchmarkOutcome.FAILED,
                wall_time_seconds=6.0,
                model_load_seconds=3.0,
                peak_ram_bytes=150,
                peak_vram_bytes=250,
                external_process_timed_out=True,
            ),
            StratumRun("short-gap", BenchmarkOutcome.UNSUPPORTED),
            StratumRun(
                "long-gap",
                BenchmarkOutcome.FAILED,
                external_process_crashed=True,
            ),
        )
    )

    long_gap, short_gap = summaries
    assert long_gap.stratum == "long-gap"
    assert long_gap.failure_rate == 1.0
    assert long_gap.crash_rate == 1.0
    assert long_gap.median_wall_time_seconds is None
    assert short_gap.run_count == 3
    assert short_gap.success_rate == pytest.approx(1 / 3)
    assert short_gap.unsupported_rate == pytest.approx(1 / 3)
    assert short_gap.failure_rate == pytest.approx(1 / 3)
    assert short_gap.timeout_rate == pytest.approx(1 / 3)
    assert short_gap.crash_rate == 0.0
    assert short_gap.median_wall_time_seconds == 5.0
    assert short_gap.median_model_load_seconds == 2.0
    assert short_gap.peak_ram_bytes == 150
    assert short_gap.peak_vram_bytes == 250


def test_stratum_run_rejects_inconsistent_failure_evidence() -> None:
    with pytest.raises(ValueError, match="failed outcome"):
        StratumRun(
            "short-gap",
            BenchmarkOutcome.SUCCESS,
            external_process_timed_out=True,
        )
    with pytest.raises(ValueError, match="both"):
        StratumRun(
            "short-gap",
            BenchmarkOutcome.FAILED,
            external_process_timed_out=True,
            external_process_crashed=True,
        )


def test_modeller_comparison_applies_predeclared_acceptance_gates() -> None:
    passing = ModellerComparison(
        experimental_median_gap_backbone_rmsd_angstrom=1.20,
        modeller_median_gap_backbone_rmsd_angstrom=1.00,
        experimental_junction_pass_rate=0.95,
        modeller_junction_pass_rate=0.90,
        experimental_fixed_heavy_rmsd_angstrom=0.005,
        modeller_fixed_heavy_rmsd_angstrom=0.02,
    )
    failing = ModellerComparison(
        experimental_median_gap_backbone_rmsd_angstrom=1.30,
        modeller_median_gap_backbone_rmsd_angstrom=1.00,
        experimental_junction_pass_rate=0.80,
        modeller_junction_pass_rate=0.90,
        experimental_fixed_heavy_rmsd_angstrom=0.02,
        modeller_fixed_heavy_rmsd_angstrom=0.02,
    )

    assert passing.gap_backbone_rmsd_delta_angstrom == pytest.approx(0.20)
    assert passing.passed
    assert failing.gap_backbone_rmsd_delta_angstrom == pytest.approx(0.30)
    assert not failing.passes_rmsd_gate
    assert not failing.passes_junction_gate
    assert not failing.passes_fixed_adherence_gate
    assert not failing.passed


def test_backend_summary_counts_failures_and_reports_wilson_interval() -> None:
    runs = (
        BackendCaseResult("candidate", "case-1", "structure-a", True, True, True, True, 1.0),
        BackendCaseResult("candidate", "case-2", "structure-a", False, True, True, True),
        BackendCaseResult("candidate", "case-3", "structure-b", True, True, True, True, 2.0),
    )

    summary = summarize_backend_cases(runs)
    interval = wilson_score_interval(50, 100)

    assert summary.case_count == 3
    assert summary.passing_count == 2
    assert summary.validation_pass_rate == pytest.approx(2 / 3)
    assert summary.median_valid_gap_backbone_rmsd_angstrom == 1.5
    assert summary.wall_time_observed_count == 0
    assert summary.peak_vram_observed_count == 0
    assert not summary.resource_reporting_complete
    assert summary.eligible
    assert interval.lower == pytest.approx(0.4038, abs=1e-4)
    assert interval.upper == pytest.approx(0.5962, abs=1e-4)


def test_paired_backend_comparison_bootstraps_independence_groups() -> None:
    candidate = (
        BackendCaseResult("candidate", "a-short", "structure-a", True, True, True, True, 1.0),
        BackendCaseResult("candidate", "a-long", "structure-a", True, True, True, True, 1.2),
        BackendCaseResult("candidate", "b-short", "structure-b", True, True, True, True, 0.8),
        BackendCaseResult("candidate", "c-short", "structure-c", True, True, True, True, 1.1),
    )
    baseline = (
        BackendCaseResult("baseline", "a-short", "structure-a", True, True, True, True, 2.0),
        BackendCaseResult("baseline", "a-long", "structure-a", True, True, True, True, 2.2),
        BackendCaseResult("baseline", "b-short", "structure-b", True, True, True, True, 1.8),
        BackendCaseResult("baseline", "c-short", "structure-c", True, True, True, True, 2.1),
    )

    comparison = compare_paired_backends(
        candidate,
        baseline,
        minimum_independence_groups=3,
        minimum_both_valid_independence_groups=3,
        require_complete_resources=False,
        bootstrap_samples=500,
        random_seed=7,
    )

    assert comparison.paired_case_count == 4
    assert comparison.both_passed_count == 4
    assert comparison.pass_rate_difference == 0.0
    assert comparison.pass_rate_difference_interval.lower == 0.0
    assert comparison.sample_size_sufficient
    assert comparison.noninferiority_demonstrated
    assert comparison.median_paired_rmsd_difference_angstrom == pytest.approx(-1.0)
    assert comparison.median_paired_rmsd_difference_interval is not None
    assert comparison.median_paired_rmsd_difference_interval.upper < 0
    assert comparison.decision == "candidate-preferred"


def test_paired_backend_comparison_fails_closed_on_fixed_coordinate_drift() -> None:
    candidate = (
        BackendCaseResult("candidate", "case-1", "structure-a", True, False, True, True, 0.5),
    )
    baseline = (
        BackendCaseResult("baseline", "case-1", "structure-a", True, True, True, True, 1.0),
    )

    comparison = compare_paired_backends(
        candidate,
        baseline,
        bootstrap_samples=20,
    )

    assert not comparison.candidate.eligible
    assert comparison.candidate.ineligibility_reasons == ("fixed-coordinates-not-exact",)
    assert not comparison.noninferiority_demonstrated
    assert comparison.decision == "candidate-ineligible"


def test_paired_backend_comparison_does_not_select_from_pilot_sample() -> None:
    candidate = (
        BackendCaseResult("candidate", "case-1", "structure-a", True, True, True, True, 0.5),
    )
    baseline = (
        BackendCaseResult("baseline", "case-1", "structure-a", True, True, True, True, 1.0),
    )

    comparison = compare_paired_backends(
        candidate,
        baseline,
        bootstrap_samples=20,
    )

    assert comparison.minimum_independence_groups == 180
    assert comparison.independence_group_count == 1
    assert not comparison.sample_size_sufficient
    assert not comparison.noninferiority_demonstrated
    assert comparison.decision == "insufficient-sample"


def test_backend_summary_rejects_unresolved_training_leakage() -> None:
    summary = summarize_backend_cases(
        (
            BackendCaseResult(
                "candidate",
                "case-1",
                "structure-a",
                True,
                True,
                True,
                False,
                0.5,
            ),
        )
    )

    assert not summary.eligible
    assert summary.ineligibility_reasons == ("unresolved-training-leakage",)


def test_paired_backend_comparison_requires_enough_both_valid_groups() -> None:
    candidate = (
        BackendCaseResult("candidate", "case-a", "group-a", True, True, True, True, 0.5),
        BackendCaseResult("candidate", "case-b", "group-b", False, True, True, True),
        BackendCaseResult("candidate", "case-c", "group-c", False, True, True, True),
    )
    baseline = (
        BackendCaseResult("baseline", "case-a", "group-a", True, True, True, True, 1.0),
        BackendCaseResult("baseline", "case-b", "group-b", False, True, True, True),
        BackendCaseResult("baseline", "case-c", "group-c", False, True, True, True),
    )

    comparison = compare_paired_backends(
        candidate,
        baseline,
        minimum_independence_groups=3,
        minimum_both_valid_independence_groups=2,
        require_complete_resources=False,
        bootstrap_samples=50,
    )

    assert comparison.sample_size_sufficient
    assert comparison.noninferiority_demonstrated
    assert comparison.both_valid_independence_group_count == 1
    assert not comparison.rmsd_sample_size_sufficient
    assert comparison.decision == "noninferior-rmsd-insufficient-sample"


def test_paired_backend_comparison_requires_complete_resource_reporting() -> None:
    candidate = (BackendCaseResult("candidate", "case-a", "group-a", True, True, True, True, 0.5),)
    baseline = (BackendCaseResult("baseline", "case-a", "group-a", True, True, True, True, 1.0),)

    comparison = compare_paired_backends(
        candidate,
        baseline,
        minimum_independence_groups=1,
        maximum_independence_groups=1,
        minimum_both_valid_independence_groups=1,
        bootstrap_samples=20,
    )

    assert comparison.noninferiority_demonstrated
    assert not comparison.candidate.resource_reporting_complete
    assert comparison.decision == "resource-reporting-incomplete"


def test_paired_backend_comparison_enforces_preregistered_maximum() -> None:
    candidate = (
        BackendCaseResult("candidate", "case-a", "group-a", True, True, True, True, 0.5),
        BackendCaseResult("candidate", "case-b", "group-b", True, True, True, True, 0.6),
    )
    baseline = (
        BackendCaseResult("baseline", "case-a", "group-a", True, True, True, True, 1.0),
        BackendCaseResult("baseline", "case-b", "group-b", True, True, True, True, 1.1),
    )

    comparison = compare_paired_backends(
        candidate,
        baseline,
        minimum_independence_groups=1,
        maximum_independence_groups=1,
        minimum_both_valid_independence_groups=1,
        require_complete_resources=False,
        bootstrap_samples=20,
    )

    assert not comparison.sample_size_sufficient
    assert comparison.maximum_independence_groups == 1
    assert comparison.decision == "sample-exceeds-preregistered-maximum"
