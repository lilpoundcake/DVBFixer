"""CPU-testable repeatability and withheld-coordinate diffusion metrics."""

from __future__ import annotations

import hashlib
import math
import os
import stat
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import numpy as np

from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    ResidueIdentity,
)

REPEATABILITY_RMSD_MAX_ANGSTROM = 0.01
_BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})


class DiffusionBenchmarkError(RuntimeError):
    """Raised when benchmark structures cannot be compared safely."""


@dataclass(frozen=True, slots=True)
class RepeatabilityAssessment:
    classification: str
    coordinate_rmsd_angstrom: float
    maximum_displacement_angstrom: float
    atom_count: int
    threshold_angstrom: float = REPEATABILITY_RMSD_MAX_ANGSTROM


@dataclass(frozen=True, slots=True)
class CandidateQuality:
    candidate_id: str
    gap_backbone_rmsd_angstrom: float
    gap_all_heavy_rmsd_angstrom: float
    fixed_heavy_rmsd_angstrom: float
    anchor_heavy_rmsd_angstrom: float
    closure_passed: bool
    validation_passed: bool


@dataclass(frozen=True, slots=True)
class EnsembleSummary:
    candidate_count: int
    passing_count: int
    closure_pass_rate: float
    validation_pass_rate: float
    top1_gap_backbone_rmsd_angstrom: float
    oracle_gap_backbone_rmsd_angstrom: float
    mean_pairwise_gap_backbone_rmsd_angstrom: float
    ranking_enrichment_angstrom: float


class BenchmarkOutcome(StrEnum):
    SUCCESS = "success"
    UNSUPPORTED = "unsupported"
    FAILED = "failed"


@dataclass(frozen=True, slots=True)
class StratumRun:
    stratum: str
    outcome: BenchmarkOutcome
    wall_time_seconds: float | None = None
    model_load_seconds: float | None = None
    peak_ram_bytes: int | None = None
    peak_vram_bytes: int | None = None
    external_process_timed_out: bool = False
    external_process_crashed: bool = False

    def __post_init__(self) -> None:
        if not self.stratum:
            raise ValueError("benchmark stratum must not be empty")
        for value in (self.wall_time_seconds, self.model_load_seconds):
            if value is not None and (not math.isfinite(value) or value < 0):
                raise ValueError("benchmark timings must be finite and non-negative")
        for value in (self.peak_ram_bytes, self.peak_vram_bytes):
            if value is not None and value < 0:
                raise ValueError("benchmark memory values must be non-negative")
        if self.external_process_timed_out and self.external_process_crashed:
            raise ValueError("one benchmark run cannot be both timed out and crashed")
        if self.outcome is not BenchmarkOutcome.FAILED and (
            self.external_process_timed_out or self.external_process_crashed
        ):
            raise ValueError("timeout or crash evidence requires a failed outcome")


@dataclass(frozen=True, slots=True)
class StratumSummary:
    stratum: str
    run_count: int
    success_count: int
    unsupported_count: int
    failure_count: int
    success_rate: float
    unsupported_rate: float
    failure_rate: float
    timeout_rate: float
    crash_rate: float
    median_wall_time_seconds: float | None
    median_model_load_seconds: float | None
    peak_ram_bytes: int | None
    peak_vram_bytes: int | None


@dataclass(frozen=True, slots=True)
class ModellerComparison:
    experimental_median_gap_backbone_rmsd_angstrom: float
    modeller_median_gap_backbone_rmsd_angstrom: float
    experimental_junction_pass_rate: float
    modeller_junction_pass_rate: float
    experimental_fixed_heavy_rmsd_angstrom: float
    modeller_fixed_heavy_rmsd_angstrom: float
    rmsd_delta_limit_angstrom: float = 0.25

    def __post_init__(self) -> None:
        values = (
            self.experimental_median_gap_backbone_rmsd_angstrom,
            self.modeller_median_gap_backbone_rmsd_angstrom,
            self.experimental_fixed_heavy_rmsd_angstrom,
            self.modeller_fixed_heavy_rmsd_angstrom,
            self.rmsd_delta_limit_angstrom,
        )
        if any(not math.isfinite(value) or value < 0 for value in values):
            raise ValueError("MODELLER comparison values must be finite and non-negative")
        for rate in (
            self.experimental_junction_pass_rate,
            self.modeller_junction_pass_rate,
        ):
            if not math.isfinite(rate) or not 0 <= rate <= 1:
                raise ValueError("MODELLER comparison rates must be within [0, 1]")

    @property
    def gap_backbone_rmsd_delta_angstrom(self) -> float:
        return (
            self.experimental_median_gap_backbone_rmsd_angstrom
            - self.modeller_median_gap_backbone_rmsd_angstrom
        )

    @property
    def passes_rmsd_gate(self) -> bool:
        return self.gap_backbone_rmsd_delta_angstrom <= self.rmsd_delta_limit_angstrom

    @property
    def passes_junction_gate(self) -> bool:
        return self.experimental_junction_pass_rate >= self.modeller_junction_pass_rate

    @property
    def passes_fixed_adherence_gate(self) -> bool:
        return (
            self.experimental_fixed_heavy_rmsd_angstrom
            < self.modeller_fixed_heavy_rmsd_angstrom
        )

    @property
    def passed(self) -> bool:
        return (
            self.passes_rmsd_gate
            and self.passes_junction_gate
            and self.passes_fixed_adherence_gate
        )


def assess_same_seed_repeatability(
    first: dict[AtomIdentity, np.ndarray],
    second: dict[AtomIdentity, np.ndarray],
    *,
    atom_identities: set[AtomIdentity] | None = None,
    threshold_angstrom: float = REPEATABILITY_RMSD_MAX_ANGSTROM,
) -> RepeatabilityAssessment:
    """Classify two same-seed coordinate sets using stable atom identities."""
    if not math.isfinite(threshold_angstrom) or threshold_angstrom <= 0:
        raise ValueError("repeatability threshold must be finite and positive")
    if atom_identities is None:
        identities = set(first)
        if set(second) != identities:
            raise DiffusionBenchmarkError(
                "repeatability structures do not contain identical requested atoms"
            )
    else:
        identities = set(atom_identities)
        if set(first) & identities != identities or set(second) & identities != identities:
            raise DiffusionBenchmarkError(
                "repeatability structures do not contain identical requested atoms"
            )
    if not identities:
        raise DiffusionBenchmarkError("repeatability atom identity set is empty")
    displacements = _displacements(first, second, identities)
    rmsd = float(np.sqrt(np.mean(displacements * displacements)))
    maximum = float(np.max(displacements))
    return RepeatabilityAssessment(
        classification="deterministic" if rmsd <= threshold_angstrom else "nondeterministic",
        coordinate_rmsd_angstrom=rmsd,
        maximum_displacement_angstrom=maximum,
        atom_count=len(identities),
        threshold_angstrom=threshold_angstrom,
    )


def candidate_quality(
    candidate_id: str,
    reference: dict[AtomIdentity, np.ndarray],
    candidate: dict[AtomIdentity, np.ndarray],
    *,
    generated_residues: set[ResidueIdentity],
    fixed_atoms: set[AtomIdentity],
    anchor_residues: set[ResidueIdentity],
    closure_passed: bool,
    validation_passed: bool,
) -> CandidateQuality:
    """Compute honest coordinate metrics against withheld reference atoms."""
    generated_heavy = {
        identity
        for identity in reference
        if _residue(identity) in generated_residues
        and not _is_hydrogen(identity)
    }
    generated_backbone = {
        identity
        for identity in generated_heavy
        if identity.atom_name in _BACKBONE_ATOMS
    }
    fixed_heavy = {
        identity
        for identity in fixed_atoms
        if identity in reference and not _is_hydrogen(identity)
    }
    anchor_heavy = {
        identity
        for identity in reference
        if _residue(identity) in anchor_residues
        and not _is_hydrogen(identity)
    }
    return CandidateQuality(
        candidate_id=candidate_id,
        gap_backbone_rmsd_angstrom=_rmsd(
            reference,
            candidate,
            generated_backbone,
            "gap backbone",
        ),
        gap_all_heavy_rmsd_angstrom=_rmsd(
            reference,
            candidate,
            generated_heavy,
            "gap all-heavy",
        ),
        fixed_heavy_rmsd_angstrom=_rmsd(
            reference,
            candidate,
            fixed_heavy,
            "fixed heavy",
        ),
        anchor_heavy_rmsd_angstrom=_rmsd(
            reference,
            candidate,
            anchor_heavy,
            "anchor heavy",
        ),
        closure_passed=closure_passed,
        validation_passed=validation_passed,
    )


def summarize_ensemble(
    qualities: tuple[CandidateQuality, ...],
    ranked_candidate_ids: tuple[str, ...],
    pairwise_gap_backbone_rmsds: tuple[float, ...],
) -> EnsembleSummary:
    """Summarize top-1/oracle quality, rates, diversity, and ranking enrichment."""
    if not qualities:
        raise DiffusionBenchmarkError("ensemble contains no candidates")
    by_id = {quality.candidate_id: quality for quality in qualities}
    if len(by_id) != len(qualities):
        raise DiffusionBenchmarkError("ensemble candidate IDs are not unique")
    if set(ranked_candidate_ids) != set(by_id):
        raise DiffusionBenchmarkError(
            "ranked candidate IDs do not match ensemble candidates"
        )
    if any(not math.isfinite(value) or value < 0 for value in pairwise_gap_backbone_rmsds):
        raise DiffusionBenchmarkError("pairwise RMSDs must be finite and non-negative")
    top1 = by_id[ranked_candidate_ids[0]].gap_backbone_rmsd_angstrom
    oracle = min(
        quality.gap_backbone_rmsd_angstrom
        for quality in qualities
    )
    return EnsembleSummary(
        candidate_count=len(qualities),
        passing_count=sum(quality.validation_passed for quality in qualities),
        closure_pass_rate=(
            sum(quality.closure_passed for quality in qualities)
            / len(qualities)
        ),
        validation_pass_rate=(
            sum(quality.validation_passed for quality in qualities)
            / len(qualities)
        ),
        top1_gap_backbone_rmsd_angstrom=top1,
        oracle_gap_backbone_rmsd_angstrom=oracle,
        mean_pairwise_gap_backbone_rmsd_angstrom=(
            float(np.mean(pairwise_gap_backbone_rmsds))
            if pairwise_gap_backbone_rmsds
            else 0.0
        ),
        ranking_enrichment_angstrom=top1 - oracle,
    )


def pairwise_gap_backbone_rmsds(
    candidates: tuple[dict[AtomIdentity, np.ndarray], ...],
    generated_residues: set[ResidueIdentity],
) -> tuple[float, ...]:
    """Return pairwise generated-backbone RMSDs for ensemble diversity."""
    if not candidates:
        return ()
    identities = {
        identity
        for identity in candidates[0]
        if _residue(identity) in generated_residues
        and identity.atom_name in _BACKBONE_ATOMS
        and not _is_hydrogen(identity)
    }
    values: list[float] = []
    for index, first in enumerate(candidates):
        for second in candidates[index + 1 :]:
            values.append(_rmsd(first, second, identities, "pairwise gap backbone"))
    return tuple(values)


def summarize_strata(
    runs: tuple[StratumRun, ...],
) -> tuple[StratumSummary, ...]:
    """Summarize outcomes and observed resources without fabricating missing data."""
    if not runs:
        raise DiffusionBenchmarkError("benchmark run set is empty")
    grouped: dict[str, list[StratumRun]] = {}
    for run in runs:
        grouped.setdefault(run.stratum, []).append(run)
    summaries: list[StratumSummary] = []
    for stratum in sorted(grouped):
        stratum_runs = grouped[stratum]
        run_count = len(stratum_runs)
        success_count = sum(
            run.outcome is BenchmarkOutcome.SUCCESS
            for run in stratum_runs
        )
        unsupported_count = sum(
            run.outcome is BenchmarkOutcome.UNSUPPORTED
            for run in stratum_runs
        )
        failure_count = sum(
            run.outcome is BenchmarkOutcome.FAILED
            for run in stratum_runs
        )
        wall_times = [
            run.wall_time_seconds
            for run in stratum_runs
            if run.wall_time_seconds is not None
        ]
        model_load_times = [
            run.model_load_seconds
            for run in stratum_runs
            if run.model_load_seconds is not None
        ]
        peak_ram = [
            run.peak_ram_bytes
            for run in stratum_runs
            if run.peak_ram_bytes is not None
        ]
        peak_vram = [
            run.peak_vram_bytes
            for run in stratum_runs
            if run.peak_vram_bytes is not None
        ]
        summaries.append(
            StratumSummary(
                stratum=stratum,
                run_count=run_count,
                success_count=success_count,
                unsupported_count=unsupported_count,
                failure_count=failure_count,
                success_rate=success_count / run_count,
                unsupported_rate=unsupported_count / run_count,
                failure_rate=failure_count / run_count,
                timeout_rate=(
                    sum(run.external_process_timed_out for run in stratum_runs)
                    / run_count
                ),
                crash_rate=(
                    sum(run.external_process_crashed for run in stratum_runs)
                    / run_count
                ),
                median_wall_time_seconds=(
                    float(np.median(wall_times)) if wall_times else None
                ),
                median_model_load_seconds=(
                    float(np.median(model_load_times))
                    if model_load_times
                    else None
                ),
                peak_ram_bytes=max(peak_ram) if peak_ram else None,
                peak_vram_bytes=max(peak_vram) if peak_vram else None,
            )
        )
    return tuple(summaries)


def load_pdb_coordinates(
    root: Path,
    artifact: ArtifactReference,
) -> dict[AtomIdentity, np.ndarray]:
    """Read stable-identity coordinates from a verified contained PDB artifact."""
    data = _read_contained(root, artifact)
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise DiffusionBenchmarkError("benchmark PDB is not valid UTF-8") from exc
    coordinates: dict[AtomIdentity, np.ndarray] = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 54:
            raise DiffusionBenchmarkError(
                f"truncated benchmark PDB coordinate on line {line_number}"
            )
        identity = AtomIdentity(
            line[21],
            line[22:26].strip(),
            line[26].strip(),
            line[12:16].strip(),
        )
        if identity in coordinates:
            raise DiffusionBenchmarkError("duplicate benchmark atom identity")
        try:
            point = np.asarray(
                (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ),
                dtype=np.float64,
            )
        except ValueError as exc:
            raise DiffusionBenchmarkError(
                f"malformed benchmark PDB coordinate on line {line_number}"
            ) from exc
        if not np.all(np.isfinite(point)):
            raise DiffusionBenchmarkError("non-finite benchmark PDB coordinate")
        coordinates[identity] = point
    if not coordinates:
        raise DiffusionBenchmarkError("benchmark PDB contains no coordinates")
    return coordinates


def _rmsd(
    reference: dict[AtomIdentity, np.ndarray],
    candidate: dict[AtomIdentity, np.ndarray],
    identities: set[AtomIdentity],
    label: str,
) -> float:
    if not identities:
        raise DiffusionBenchmarkError(f"{label} atom identity set is empty")
    if set(reference) & identities != identities or set(candidate) & identities != identities:
        raise DiffusionBenchmarkError(f"{label} atom identity mismatch")
    displacements = _displacements(reference, candidate, identities)
    return float(np.sqrt(np.mean(displacements * displacements)))


def _displacements(
    first: dict[AtomIdentity, np.ndarray],
    second: dict[AtomIdentity, np.ndarray],
    identities: set[AtomIdentity],
) -> np.ndarray:
    ordered = sorted(identities)
    return np.asarray(
        [
            np.linalg.norm(second[identity] - first[identity])
            for identity in ordered
        ],
        dtype=np.float64,
    )


def _residue(identity: AtomIdentity) -> ResidueIdentity:
    return ResidueIdentity(
        identity.chain,
        identity.residue_number,
        identity.insertion_code,
    )


def _is_hydrogen(identity: AtomIdentity) -> bool:
    return identity.atom_name.lstrip("0123456789").upper().startswith("H")


def _read_contained(root: Path, artifact: ArtifactReference) -> bytes:
    if root.is_symlink() or not root.is_dir():
        raise DiffusionBenchmarkError("benchmark root must be a regular directory")
    resolved_root = root.resolve(strict=True)
    path = resolved_root.joinpath(*Path(artifact.path).parts)
    current = resolved_root
    for part in Path(artifact.path).parts:
        current = current / part
        if current.is_symlink():
            raise DiffusionBenchmarkError("benchmark artifact path contains a symlink")
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise DiffusionBenchmarkError("benchmark artifact does not exist") from exc
    if not resolved.is_relative_to(resolved_root):
        raise DiffusionBenchmarkError("benchmark artifact escapes its root")
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        file_stat = os.fstat(descriptor)
        if not stat.S_ISREG(file_stat.st_mode) or file_stat.st_nlink != 1:
            raise DiffusionBenchmarkError(
                "benchmark artifact is not a private regular file"
            )
        digest = hashlib.sha256()
        data = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            digest.update(chunk)
            data.extend(chunk)
    finally:
        os.close(descriptor)
    if digest.hexdigest() != artifact.sha256.lower():
        raise DiffusionBenchmarkError("benchmark artifact SHA-256 mismatch")
    return bytes(data)
