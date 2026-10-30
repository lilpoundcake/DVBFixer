"""Separate, deny-by-default wire contract for experimental atom completion.

This contract does not confer scientific admission or enable any existing runner.
Schema-4 gap requests are neither rewritten nor reinterpreted.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from dvbfixer.model.diffusion.atom_completion import (
    AtomCompletionMasks,
    CompletionKind,
    plan_atom_completion,
)
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    DiffusionContractError,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetSequence,
    TemplateOwnership,
    _JsonContract,
)
from dvbfixer.model.diffusion.protein_atoms import (
    _CANONICAL_ONE_TO_THREE,
    canonical_heavy_atom_identities,
)

ATOM_COMPLETION_SCHEMA_VERSION = 1
MAX_TARGET_RESIDUES = 100_000
MAX_CANDIDATES = 64


@dataclass(frozen=True, slots=True)
class AtomCompletionRequest(_JsonContract):
    normalized_pdb: ArtifactReference
    target_sequences: tuple[TargetSequence, ...]
    sequence_placements: tuple[SequencePlacement, ...]
    gaps: tuple[GapRegion, ...]
    masks: AtomCompletionMasks
    seeds: tuple[int, ...]
    profile: str
    template_ownership: tuple[TemplateOwnership, ...] = ()

    _schema_version: ClassVar[int] = ATOM_COMPLETION_SCHEMA_VERSION

    @property
    def generated_atoms(self) -> tuple[AtomIdentity, ...]:
        return self.masks.gap_generated_atoms + self.masks.missing_atoms

    def __post_init__(self) -> None:
        super(AtomCompletionRequest, self).__post_init__()
        if type(self.schema_version) is not int:
            raise DiffusionContractError("schema_version must be an integer")
        if not self.profile.strip():
            raise DiffusionContractError("atom completion requires an explicit runner profile")
        if self.normalized_pdb.sha256.lower() != self.masks.source_sha256:
            raise DiffusionContractError("atom-completion source digest differs from mask binding")
        if not self.masks.partial_residues:
            raise DiffusionContractError("atom completion requires at least one partial residue")
        if not 1 <= len(self.seeds) <= MAX_CANDIDATES or len(set(self.seeds)) != len(self.seeds):
            raise DiffusionContractError("atom completion requires 1-64 unique seeds")
        if any(type(seed) is not int or not 0 <= seed < 2**32 for seed in self.seeds):
            raise DiffusionContractError("atom-completion seeds must be unsigned 32-bit integers")
        targets = {item.chain: item.sequence for item in self.target_sequences}
        placements = {item.chain: item for item in self.sequence_placements}
        if (not targets or len(targets) != len(self.target_sequences)
                or len(placements) != len(self.sequence_placements) or targets.keys() != placements.keys()):
            raise DiffusionContractError("unique targets and placements must cover the same chains")
        if sum(map(len, targets.values())) > MAX_TARGET_RESIDUES:
            raise DiffusionContractError("atom-completion target exceeds contract resource bound")
        if any(letter not in _CANONICAL_ONE_TO_THREE for seq in targets.values() for letter in seq):
            raise DiffusionContractError("atom completion requires canonical target sequences")
        observed: dict[tuple[str, int], ResidueIdentity] = {}
        for chain, placement in placements.items():
            if placement.target_length != len(targets[chain]):
                raise DiffusionContractError("placement length differs from target sequence")
            observed.update(((chain, index), residue) for index, residue in zip(
                placement.observed_target_indices, placement.observed_residues, strict=True))
        fixed_residues = {_residue(atom) for atom in self.masks.fixed_atoms}
        if not set(observed.values()) <= fixed_residues:
            raise DiffusionContractError("observed placements must be represented by fixed atoms")
        for partial in self.masks.partial_residues:
            key = partial.residue.chain, partial.target_index
            if observed.get(key) != partial.residue:
                raise DiffusionContractError("partial residue must match its observed target placement")
            if _CANONICAL_ONE_TO_THREE[targets[key[0]][key[1]]] != partial.residue_name:
                raise DiffusionContractError("partial residue name differs from target sequence")
        generated_positions: dict[tuple[str, int], ResidueIdentity] = {}
        generated_residues: set[ResidueIdentity] = set()
        gap_atoms: set[AtomIdentity] = set()
        for gap in self.gaps:
            sequence = targets.get(gap.chain)
            if sequence is None or gap.target_interval.stop > len(sequence):
                raise DiffusionContractError("gap exceeds target sequence")
            start, stop = gap.target_interval.start, gap.target_interval.stop
            if ((gap.left_anchor is None) != (start == 0)
                    or (gap.right_anchor is None) != (stop == len(sequence))):
                raise DiffusionContractError("terminal gap anchors must match target sequence boundaries")
            if (gap.left_anchor != observed.get((gap.chain, start - 1))
                    or gap.right_anchor != observed.get((gap.chain, stop))):
                raise DiffusionContractError("gap anchors must match adjacent observed target positions")
            for index, residue in zip(range(start, stop), gap.generated_residues, strict=True):
                key = gap.chain, index
                if key in observed or key in generated_positions or residue in fixed_residues or residue in generated_residues:
                    raise DiffusionContractError("whole-gap ownership overlaps observed or generated ownership")
                generated_positions[key] = residue
                generated_residues.add(residue)
                gap_atoms.update(canonical_heavy_atom_identities(residue, _CANONICAL_ONE_TO_THREE[sequence[index]]))
        if gap_atoms != set(self.masks.gap_generated_atoms):
            raise DiffusionContractError("gap mask must contain exactly every canonical whole-gap atom")
        if len(observed) + len(generated_positions) != sum(map(len, targets.values())):
            raise DiffusionContractError("placements and gaps must cover every target position exactly")
        ownership_keys = [(item.chain, item.target_index) for item in self.template_ownership]
        if len(set(ownership_keys)) != len(ownership_keys) or any(
            observed.get((item.chain, item.target_index)) != item.residue for item in self.template_ownership
        ):
            raise DiffusionContractError("template ownership must match unique observed target placements")


def audit_completion_request_source(request: AtomCompletionRequest, source: bytes) -> None:
    """Recompute the partial masks from digest-bound source bytes before execution."""
    actual = plan_atom_completion(
        source, source_sha256=request.normalized_pdb.sha256,
        target_sequences=request.target_sequences, sequence_placements=request.sequence_placements,
        fixed_atoms=request.masks.fixed_atoms, gap_generated_atoms=request.masks.gap_generated_atoms,
    )
    if actual != request.masks:
        raise DiffusionContractError("atom-completion masks differ from independent source inventory")


@dataclass(frozen=True, slots=True)
class AtomCompletionCapabilities(_JsonContract):
    profile: str
    sidechain: bool = False
    backbone: bool = False
    whole_gaps: bool = False
    mixed_gaps: bool = False
    max_target_residues: int = 0
    max_generated_atoms: int = 0
    max_candidates: int = 0

    _schema_version: ClassVar[int] = ATOM_COMPLETION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        super(AtomCompletionCapabilities, self).__post_init__()
        if type(self.schema_version) is not int:
            raise DiffusionContractError("schema_version must be an integer")
        if not self.profile.strip():
            raise DiffusionContractError("atom-completion capability requires a runner profile")
        for name in ("sidechain", "backbone", "whole_gaps", "mixed_gaps"):
            if type(getattr(self, name)) is not bool:
                raise DiffusionContractError("atom-completion features must be explicit booleans")
        for name, ceiling in (("max_target_residues", MAX_TARGET_RESIDUES),
                              ("max_generated_atoms", MAX_TARGET_RESIDUES * 14),
                              ("max_candidates", MAX_CANDIDATES)):
            value = getattr(self, name)
            if type(value) is not int or not 0 <= value <= ceiling:
                raise DiffusionContractError(f"invalid atom-completion resource bound: {name}")


def require_completion_capability(request: AtomCompletionRequest, capability: AtomCompletionCapabilities) -> None:
    """Wire support only: independent source/scope/geometry admission is still required."""
    if request.profile != capability.profile:
        raise DiffusionContractError("atom-completion runner profile mismatch")
    for partial in request.masks.partial_residues:
        supported = capability.backbone if partial.kind is CompletionKind.BACKBONE else capability.sidechain
        if not supported:
            raise DiffusionContractError(f"runner does not support partial {partial.kind.value} completion")
    if request.gaps and not (capability.whole_gaps and capability.mixed_gaps):
        raise DiffusionContractError("runner does not support mixed partial-residue and whole-gap completion")
    if (sum(len(item.sequence) for item in request.target_sequences) > capability.max_target_residues
            or len(request.generated_atoms) > capability.max_generated_atoms
            or len(request.seeds) > capability.max_candidates):
        raise DiffusionContractError("atom-completion request exceeds runner resource bounds")


def _residue(atom: AtomIdentity) -> ResidueIdentity:
    return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)
