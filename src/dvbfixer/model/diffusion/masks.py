"""Sequence placement and deterministic masks for diffusion gap requests."""

from __future__ import annotations

from dataclasses import dataclass

from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
)
from dvbfixer.sequence_alignment import align_observed_to_reference


class DiffusionMaskError(ValueError):
    """Raised when a sequence placement or generated region is unsupported."""


@dataclass(frozen=True, slots=True)
class ObservedResidue:
    identity: ResidueIdentity
    one_letter_code: str
    atoms: tuple[str, ...]

    def __post_init__(self) -> None:
        if len(self.one_letter_code) != 1 or not self.one_letter_code.isalpha():
            raise DiffusionMaskError("one_letter_code must be one alphabetic character")
        if not self.atoms:
            raise DiffusionMaskError("observed residue must contain at least one atom")
        if any(not atom or len(atom) > 4 for atom in self.atoms):
            raise DiffusionMaskError("observed atom names must contain one to four characters")
        if len(set(self.atoms)) != len(self.atoms):
            raise DiffusionMaskError("observed residue atom names must be unique")


@dataclass(frozen=True, slots=True)
class DiffusionMasks:
    sequence_placement: SequencePlacement
    gaps: tuple[GapRegion, ...]
    fixed_atoms: tuple[AtomIdentity, ...]
    generated_residues: tuple[ResidueIdentity, ...]


def build_sequence_placement(
    chain: str,
    observed_residues: tuple[ObservedResidue, ...],
    target_sequence: str,
) -> SequencePlacement:
    """Place observed residues with the repository's shared affine alignment."""
    if not observed_residues:
        raise DiffusionMaskError("observed_residues must not be empty")
    if any(residue.identity.chain != chain for residue in observed_residues):
        raise DiffusionMaskError("observed residue chain does not match target chain")

    observed_sequence = "".join(residue.one_letter_code.upper() for residue in observed_residues)
    alignment = align_observed_to_reference(observed_sequence, target_sequence)
    if alignment is None:
        raise DiffusionMaskError("observed sequence cannot be placed in target sequence")
    if alignment.ambiguous:
        raise DiffusionMaskError("observed sequence placement is ambiguous")

    return SequencePlacement(
        chain=chain,
        target_length=len(target_sequence),
        observed_target_indices=alignment.positions,
        observed_residues=tuple(residue.identity for residue in observed_residues),
    )


def build_diffusion_masks(
    chain: str,
    observed_residues: tuple[ObservedResidue, ...],
    target_sequence: str,
    *,
    minimum_gap_length: int = 3,
    maximum_gap_length: int = 12,
    movable_flank_residues: int = 1,
) -> DiffusionMasks:
    """Build fixed atoms and two-anchor internal gaps from one placement.

    Generated residue numbering uses target-sequence ordinals. The scientific
    pipeline may later apply its authoritative residue allocator before writing
    PDB output, but mask construction remains deterministic and insertion-code
    aware for every observed residue.
    """
    if minimum_gap_length <= 0 or maximum_gap_length < minimum_gap_length:
        raise DiffusionMaskError("invalid supported gap-length range")
    if movable_flank_residues < 0:
        raise DiffusionMaskError("movable_flank_residues must be non-negative")

    placement = build_sequence_placement(chain, observed_residues, target_sequence)
    by_target = dict(zip(placement.observed_target_indices, observed_residues))

    if placement.observed_target_indices[0] > 0:
        raise DiffusionMaskError("N-terminal one-anchor gaps are unsupported")
    if placement.observed_target_indices[-1] < len(target_sequence) - 1:
        raise DiffusionMaskError("C-terminal one-anchor gaps are unsupported")

    gaps: list[GapRegion] = []
    generated_residues: list[ResidueIdentity] = []
    positions = placement.observed_target_indices
    for left_index, right_index in zip(positions, positions[1:]):
        if right_index == left_index + 1:
            continue
        start = left_index + 1
        stop = right_index
        gap_length = stop - start
        if not minimum_gap_length <= gap_length <= maximum_gap_length:
            raise DiffusionMaskError(
                f"gap length {gap_length} is outside supported range "
                f"{minimum_gap_length}-{maximum_gap_length}"
            )

        left_anchor = by_target[left_index].identity
        right_anchor = by_target[right_index].identity
        generated = tuple(
            ResidueIdentity(chain, str(target_index + 1))
            for target_index in range(start, stop)
        )
        generated_residues.extend(generated)

        movable_observed: list[ResidueIdentity] = []
        left_window_start = max(
            placement.observed_target_indices[0],
            left_index - movable_flank_residues + 1,
        )
        for target_index in range(left_window_start, left_index + 1):
            observed = by_target.get(target_index)
            if observed is not None:
                movable_observed.append(observed.identity)
        right_window_stop = min(
            placement.observed_target_indices[-1],
            right_index + movable_flank_residues - 1,
        )
        for target_index in range(right_index, right_window_stop + 1):
            observed = by_target.get(target_index)
            if observed is not None:
                movable_observed.append(observed.identity)
        movable = tuple(dict.fromkeys([*movable_observed, *generated]))
        gaps.append(
            GapRegion(
                chain=chain,
                target_interval=TargetInterval(start, stop),
                left_anchor=left_anchor,
                right_anchor=right_anchor,
                generated_residues=generated,
                movable_junction_residues=movable,
            )
        )

    if not gaps:
        raise DiffusionMaskError("target sequence contains no supported internal gaps")

    fixed_atoms = tuple(
        AtomIdentity(
            residue.identity.chain,
            residue.identity.residue_number,
            residue.identity.insertion_code,
            atom_name,
        )
        for residue in observed_residues
        for atom_name in residue.atoms
    )
    return DiffusionMasks(
        sequence_placement=placement,
        gaps=tuple(gaps),
        fixed_atoms=fixed_atoms,
        generated_residues=tuple(generated_residues),
    )
