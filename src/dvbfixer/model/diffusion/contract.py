"""Versioned request/result contract for experimental diffusion runners."""

from __future__ import annotations

import json
from dataclasses import MISSING, dataclass, fields, is_dataclass
from enum import StrEnum
from pathlib import Path
from types import UnionType
from typing import Any, ClassVar, TypeVar, Union, get_args, get_origin, get_type_hints

DIFFUSION_SCHEMA_VERSION = 1


class DiffusionContractError(ValueError):
    """Raised when a request or result violates the diffusion protocol schema."""


class DiffusionStatus(StrEnum):
    SUCCESS = "success"
    UNSUPPORTED = "unsupported"
    FAILED = "failed"


def _required_text(value: str, field_name: str) -> None:
    if not value:
        raise DiffusionContractError(f"{field_name} must not be empty")


def _single_character(value: str, field_name: str, *, allow_blank: bool = False) -> None:
    if len(value) != 1 or (not allow_blank and value == " "):
        raise DiffusionContractError(f"{field_name} must be exactly one non-blank character")


@dataclass(frozen=True, slots=True, order=True)
class ResidueIdentity:
    chain: str
    residue_number: str
    insertion_code: str = ""

    def __post_init__(self) -> None:
        _single_character(self.chain, "chain")
        _required_text(self.residue_number, "residue_number")
        if len(self.insertion_code) > 1:
            raise DiffusionContractError("insertion_code must contain at most one character")


@dataclass(frozen=True, slots=True, order=True)
class AtomIdentity:
    chain: str
    residue_number: str
    insertion_code: str
    atom_name: str

    def __post_init__(self) -> None:
        ResidueIdentity(self.chain, self.residue_number, self.insertion_code)
        _required_text(self.atom_name, "atom_name")
        if len(self.atom_name) > 4:
            raise DiffusionContractError("atom_name must contain at most four characters")


@dataclass(frozen=True, slots=True)
class TargetInterval:
    start: int
    stop: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise DiffusionContractError("target interval start must be non-negative")
        if self.stop <= self.start:
            raise DiffusionContractError("target interval stop must be greater than start")


@dataclass(frozen=True, slots=True)
class TargetSequence:
    chain: str
    sequence: str

    def __post_init__(self) -> None:
        _single_character(self.chain, "target sequence chain")
        _required_text(self.sequence, "target sequence")
        if not self.sequence.isalpha() or self.sequence != self.sequence.upper():
            raise DiffusionContractError("target sequence must contain uppercase alphabetic symbols")


@dataclass(frozen=True, slots=True)
class SequencePlacement:
    chain: str
    target_length: int
    observed_target_indices: tuple[int, ...]
    observed_residues: tuple[ResidueIdentity, ...]

    def __post_init__(self) -> None:
        _single_character(self.chain, "sequence placement chain")
        if self.target_length <= 0:
            raise DiffusionContractError("target_length must be positive")
        if len(self.observed_target_indices) != len(self.observed_residues):
            raise DiffusionContractError(
                "observed_target_indices and observed_residues must have equal length"
            )
        if tuple(sorted(self.observed_target_indices)) != self.observed_target_indices:
            raise DiffusionContractError("observed_target_indices must be strictly ordered")
        if len(set(self.observed_target_indices)) != len(self.observed_target_indices):
            raise DiffusionContractError("observed_target_indices must not contain duplicates")
        if any(index < 0 or index >= self.target_length for index in self.observed_target_indices):
            raise DiffusionContractError("observed target index is outside target_length")
        if any(residue.chain != self.chain for residue in self.observed_residues):
            raise DiffusionContractError("observed residue chain does not match placement chain")
        if len(set(self.observed_residues)) != len(self.observed_residues):
            raise DiffusionContractError("observed_residues must not contain duplicates")


@dataclass(frozen=True, slots=True)
class GapRegion:
    chain: str
    target_interval: TargetInterval
    left_anchor: ResidueIdentity
    right_anchor: ResidueIdentity
    generated_residues: tuple[ResidueIdentity, ...]
    movable_junction_residues: tuple[ResidueIdentity, ...]

    def __post_init__(self) -> None:
        _single_character(self.chain, "gap chain")
        if self.left_anchor.chain != self.chain or self.right_anchor.chain != self.chain:
            raise DiffusionContractError("gap anchors must belong to the target chain")
        if self.left_anchor == self.right_anchor:
            raise DiffusionContractError("left and right anchors must be distinct")
        expected_count = self.target_interval.stop - self.target_interval.start
        if len(self.generated_residues) != expected_count:
            raise DiffusionContractError(
                "generated_residues length must match target_interval length"
            )
        if not self.generated_residues:
            raise DiffusionContractError("gap must contain at least one generated residue")
        if any(residue.chain != self.chain for residue in self.generated_residues):
            raise DiffusionContractError("generated residue chain does not match gap chain")
        if len(set(self.generated_residues)) != len(self.generated_residues):
            raise DiffusionContractError("generated_residues must not contain duplicates")
        if any(residue.chain != self.chain for residue in self.movable_junction_residues):
            raise DiffusionContractError(
                "movable junction residue chain does not match gap chain"
            )
        if len(set(self.movable_junction_residues)) != len(self.movable_junction_residues):
            raise DiffusionContractError(
                "movable_junction_residues must not contain duplicates"
            )
        required_movable = {
            self.left_anchor,
            self.right_anchor,
            *self.generated_residues,
        }
        if not required_movable <= set(self.movable_junction_residues):
            raise DiffusionContractError(
                "movable_junction_residues must contain both anchors and all generated residues"
            )


@dataclass(frozen=True, slots=True)
class ExplicitLink:
    atom1: AtomIdentity
    atom2: AtomIdentity
    source: str

    def __post_init__(self) -> None:
        if self.atom1 == self.atom2:
            raise DiffusionContractError("explicit link endpoints must be distinct")
        _required_text(self.source, "explicit link source")


@dataclass(frozen=True, slots=True)
class BackendOption:
    name: str
    value: str

    def __post_init__(self) -> None:
        _required_text(self.name, "backend option name")


@dataclass(frozen=True, slots=True)
class ArtifactReference:
    path: str
    sha256: str

    def __post_init__(self) -> None:
        _required_text(self.path, "artifact path")
        artifact_path = Path(self.path)
        if (
            artifact_path.is_absolute()
            or not artifact_path.parts
            or ".." in artifact_path.parts
            or artifact_path.name in {"", ".", ".."}
        ):
            raise DiffusionContractError("artifact path must be relative and contained")
        if len(self.sha256) != 64:
            raise DiffusionContractError("artifact sha256 must contain 64 hexadecimal characters")
        try:
            int(self.sha256, 16)
        except ValueError as exc:
            raise DiffusionContractError("artifact sha256 must be hexadecimal") from exc


@dataclass(frozen=True, slots=True)
class Metric:
    name: str
    value: float
    unit: str = ""

    def __post_init__(self) -> None:
        _required_text(self.name, "metric name")


@dataclass(frozen=True, slots=True)
class ValidationSummary:
    passed: bool
    hard_gate_failures: tuple[str, ...] = ()
    metrics: tuple[Metric, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.passed and self.hard_gate_failures:
            raise DiffusionContractError("passing validation cannot contain hard gate failures")
        if not self.passed and not self.hard_gate_failures:
            raise DiffusionContractError("failed validation must name a hard gate failure")
        names = [metric.name for metric in self.metrics]
        if len(names) != len(set(names)):
            raise DiffusionContractError("validation metric names must be unique")


@dataclass(frozen=True, slots=True)
class RunnerDiagnostics:
    exit_code: int | None
    timed_out: bool
    stdout: str = ""
    stderr: str = ""


@dataclass(frozen=True, slots=True)
class BackendProvenance:
    backend: str
    runner_protocol_version: int
    engine_repository: str
    engine_revision: str
    checkpoint_sha256: str = ""
    environment_hash: str = ""

    def __post_init__(self) -> None:
        _required_text(self.backend, "backend")
        if self.runner_protocol_version <= 0:
            raise DiffusionContractError("runner_protocol_version must be positive")
        _required_text(self.engine_repository, "engine_repository")
        _required_text(self.engine_revision, "engine_revision")
        if self.checkpoint_sha256:
            ArtifactReference("checkpoint", self.checkpoint_sha256)


@dataclass(frozen=True, slots=True)
class DiffusionCandidate:
    candidate_id: str
    coordinate_artifact: ArtifactReference
    generated_atoms: tuple[AtomIdentity, ...]
    generated_residues: tuple[ResidueIdentity, ...]
    raw_backend_score: float | None
    score_provenance: str
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _required_text(self.candidate_id, "candidate_id")
        _required_text(self.score_provenance, "score_provenance")
        if len(set(self.generated_atoms)) != len(self.generated_atoms):
            raise DiffusionContractError("generated_atoms must not contain duplicates")
        if len(set(self.generated_residues)) != len(self.generated_residues):
            raise DiffusionContractError("generated_residues must not contain duplicates")


@dataclass(frozen=True, slots=True)
class RunnerCandidate:
    candidate_id: str
    coordinate_artifact: ArtifactReference
    generated_atoms: tuple[AtomIdentity, ...]
    generated_residues: tuple[ResidueIdentity, ...]
    raw_backend_score: float | None
    score_provenance: str
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _required_text(self.candidate_id, "candidate_id")
        _required_text(self.score_provenance, "score_provenance")
        if len(set(self.generated_atoms)) != len(self.generated_atoms):
            raise DiffusionContractError("generated_atoms must not contain duplicates")
        if len(set(self.generated_residues)) != len(self.generated_residues):
            raise DiffusionContractError("generated_residues must not contain duplicates")


T = TypeVar("T", bound="_JsonContract")


@dataclass(frozen=True, slots=True)
class _JsonContract:
    schema_version: int

    _schema_version: ClassVar[int] = DIFFUSION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != self._schema_version:
            raise DiffusionContractError(
                f"unsupported schema_version {self.schema_version}; "
                f"expected {self._schema_version}"
            )

    def to_dict(self) -> dict[str, Any]:
        return _encode(self)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, sort_keys=True, separators=(",", ":"))

    @classmethod
    def from_dict(cls: type[T], raw: dict[str, Any]) -> T:
        if not isinstance(raw, dict):
            raise DiffusionContractError(f"{cls.__name__} must be a JSON object")
        return _decode_dataclass(cls, raw)

    @classmethod
    def from_json(cls: type[T], text: str) -> T:
        try:
            raw = json.loads(text)
        except json.JSONDecodeError as exc:
            raise DiffusionContractError(f"invalid JSON: {exc.msg}") from exc
        if not isinstance(raw, dict):
            raise DiffusionContractError(f"{cls.__name__} must be a JSON object")
        return cls.from_dict(raw)


@dataclass(frozen=True, slots=True)
class DiffusionRequest(_JsonContract):
    normalized_pdb: ArtifactReference
    target_sequences: tuple[TargetSequence, ...]
    sequence_placements: tuple[SequencePlacement, ...]
    gaps: tuple[GapRegion, ...]
    fixed_atoms: tuple[AtomIdentity, ...]
    generated_atoms: tuple[AtomIdentity, ...]
    retained_explicit_links: tuple[ExplicitLink, ...]
    candidate_count: int
    seeds: tuple[int, ...]
    backend_options: tuple[BackendOption, ...] = ()

    def __post_init__(self) -> None:
        super(DiffusionRequest, self).__post_init__()
        if not self.target_sequences:
            raise DiffusionContractError("target_sequences must not be empty")
        if not self.sequence_placements:
            raise DiffusionContractError("sequence_placements must not be empty")
        if not self.gaps:
            raise DiffusionContractError("gaps must not be empty")
        sequence_chains = [sequence.chain for sequence in self.target_sequences]
        placement_chains = [placement.chain for placement in self.sequence_placements]
        if len(sequence_chains) != len(set(sequence_chains)):
            raise DiffusionContractError("target sequence chains must be unique")
        if len(placement_chains) != len(set(placement_chains)):
            raise DiffusionContractError("sequence placement chains must be unique")
        if set(sequence_chains) != set(placement_chains):
            raise DiffusionContractError(
                "target_sequences and sequence_placements must cover the same chains"
            )
        lengths = {sequence.chain: len(sequence.sequence) for sequence in self.target_sequences}
        if any(lengths[placement.chain] != placement.target_length for placement in self.sequence_placements):
            raise DiffusionContractError("sequence placement target_length does not match sequence")
        if any(gap.chain not in lengths for gap in self.gaps):
            raise DiffusionContractError("gap references a chain without a target sequence")
        if self.candidate_count <= 0:
            raise DiffusionContractError("candidate_count must be positive")
        if len(self.seeds) != self.candidate_count:
            raise DiffusionContractError("seeds length must match candidate_count")
        if len(set(self.seeds)) != len(self.seeds):
            raise DiffusionContractError("seeds must be unique")
        if any(seed < 0 for seed in self.seeds):
            raise DiffusionContractError("seeds must be non-negative")
        if len(set(self.fixed_atoms)) != len(self.fixed_atoms):
            raise DiffusionContractError("fixed_atoms must not contain duplicates")
        if len(set(self.generated_atoms)) != len(self.generated_atoms):
            raise DiffusionContractError("generated_atoms must not contain duplicates")
        overlap = set(self.fixed_atoms) & set(self.generated_atoms)
        if overlap:
            raise DiffusionContractError("fixed_atoms and generated_atoms must be disjoint")
        option_names = [option.name for option in self.backend_options]
        if len(option_names) != len(set(option_names)):
            raise DiffusionContractError("backend option names must be unique")


@dataclass(frozen=True, slots=True)
class RunnerResult(_JsonContract):
    status: DiffusionStatus
    candidates: tuple[RunnerCandidate, ...]
    runner_diagnostics: RunnerDiagnostics
    backend_provenance: BackendProvenance
    message: str = ""

    def __post_init__(self) -> None:
        super(RunnerResult, self).__post_init__()
        if not isinstance(self.status, DiffusionStatus):
            try:
                object.__setattr__(self, "status", DiffusionStatus(self.status))
            except ValueError as exc:
                raise DiffusionContractError(f"unknown diffusion status: {self.status}") from exc
        if self.status is DiffusionStatus.SUCCESS and not self.candidates:
            raise DiffusionContractError("successful runner result must contain candidates")
        if self.status is not DiffusionStatus.SUCCESS and self.candidates:
            raise DiffusionContractError(
                "unsupported or failed runner result cannot contain candidates"
            )
        if self.status is DiffusionStatus.UNSUPPORTED and not self.message:
            raise DiffusionContractError("unsupported runner result must include a message")


@dataclass(frozen=True, slots=True)
class DiffusionResult(_JsonContract):
    status: DiffusionStatus
    candidates: tuple[DiffusionCandidate, ...]
    validation_summaries: tuple[ValidationSummary, ...]
    runner_diagnostics: RunnerDiagnostics
    backend_provenance: BackendProvenance
    message: str = ""

    def __post_init__(self) -> None:
        super(DiffusionResult, self).__post_init__()
        if not isinstance(self.status, DiffusionStatus):
            try:
                object.__setattr__(self, "status", DiffusionStatus(self.status))
            except ValueError as exc:
                raise DiffusionContractError(f"unknown diffusion status: {self.status}") from exc
        if self.status is DiffusionStatus.SUCCESS:
            if not self.candidates:
                raise DiffusionContractError("successful result must contain candidates")
            if len(self.candidates) != len(self.validation_summaries):
                raise DiffusionContractError(
                    "successful result must contain one validation summary per candidate"
                )
            if not all(summary.passed for summary in self.validation_summaries):
                raise DiffusionContractError("successful result cannot contain failed validation")
        elif self.candidates:
            raise DiffusionContractError("unsupported or failed result cannot publish candidates")
        if self.status is DiffusionStatus.UNSUPPORTED and not self.message:
            raise DiffusionContractError("unsupported result must include a message")


def _encode(value: Any) -> Any:
    if isinstance(value, StrEnum):
        return value.value
    if is_dataclass(value):
        return {field.name: _encode(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, tuple):
        return [_encode(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"unsupported contract value: {type(value).__name__}")


def _decode_dataclass(cls: type[Any], raw: dict[str, Any]) -> Any:
    expected = {field.name for field in fields(cls)}
    unknown = set(raw) - expected
    missing = {
        field.name
        for field in fields(cls)
        if field.name not in raw
        and field.default is MISSING
        and field.default_factory is MISSING
    }
    if unknown:
        raise DiffusionContractError(
            f"{cls.__name__} contains unknown fields: {', '.join(sorted(unknown))}"
        )
    if missing:
        raise DiffusionContractError(
            f"{cls.__name__} is missing fields: {', '.join(sorted(missing))}"
        )
    hints = get_type_hints(cls)
    values: dict[str, Any] = {}
    for field in fields(cls):
        if field.name in raw:
            values[field.name] = _decode_value(hints[field.name], raw[field.name], field.name)
    try:
        return cls(**values)
    except (TypeError, ValueError) as exc:
        if isinstance(exc, DiffusionContractError):
            raise
        raise DiffusionContractError(f"invalid {cls.__name__}: {exc}") from exc


def _decode_value(annotation: Any, value: Any, field_name: str) -> Any:
    origin = get_origin(annotation)
    args = get_args(annotation)

    if origin is tuple:
        if not isinstance(value, list):
            raise DiffusionContractError(f"{field_name} must be a JSON array")
        item_type = args[0]
        return tuple(_decode_value(item_type, item, field_name) for item in value)
    if origin in {Union, UnionType}:
        if value is None and type(None) in args:
            return None
        non_none = [arg for arg in args if arg is not type(None)]
        if len(non_none) == 1:
            return _decode_value(non_none[0], value, field_name)
    if isinstance(annotation, type) and issubclass(annotation, StrEnum):
        if not isinstance(value, str):
            raise DiffusionContractError(f"{field_name} must be a string")
        try:
            return annotation(value)
        except ValueError as exc:
            raise DiffusionContractError(f"{field_name} has unknown value {value!r}") from exc
    if isinstance(annotation, type) and is_dataclass(annotation):
        if not isinstance(value, dict):
            raise DiffusionContractError(f"{field_name} must be a JSON object")
        return _decode_dataclass(annotation, value)
    if annotation is bool:
        if type(value) is not bool:
            raise DiffusionContractError(f"{field_name} must be a boolean")
        return value
    if annotation is int:
        if type(value) is not int:
            raise DiffusionContractError(f"{field_name} must be an integer")
        return value
    if annotation is float:
        if type(value) not in {int, float}:
            raise DiffusionContractError(f"{field_name} must be a number")
        return float(value)
    if annotation is str:
        if not isinstance(value, str):
            raise DiffusionContractError(f"{field_name} must be a string")
        return value
    raise DiffusionContractError(f"unsupported schema annotation for {field_name}")
