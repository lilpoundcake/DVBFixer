"""Source-bound atom masks for future diffusion completion of partial residues.

Planning only: this module never rebuilds coordinates, changes existing requests,
or grants a backend capability. Whole-residue gap ownership stays separate.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from enum import StrEnum

from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    DiffusionRequest,
    ResidueIdentity,
    SequencePlacement,
    TargetSequence,
)
from dvbfixer.model.diffusion.protein_atoms import (
    _CANONICAL_ONE_TO_THREE,
    canonical_heavy_atom_identities,
)

REVISION = "partial-residue-atom-masks-v1"
BACKBONE = frozenset({"N", "CA", "C", "O"})


class AtomCompletionError(ValueError):
    """An atom-completion inventory cannot be bound unambiguously to its source."""


class CompletionKind(StrEnum):
    SIDECHAIN = "sidechain"
    BACKBONE = "backbone"


@dataclass(frozen=True, slots=True)
class PartialResidueMask:
    residue: ResidueIdentity
    target_index: int
    residue_name: str
    observed_atoms: tuple[AtomIdentity, ...]
    missing_atoms: tuple[AtomIdentity, ...]

    def __post_init__(self) -> None:
        if self.target_index < 0:
            raise AtomCompletionError("partial residue target index must be non-negative")
        expected = set(canonical_heavy_atom_identities(self.residue, self.residue_name))
        observed, missing = set(self.observed_atoms), set(self.missing_atoms)
        if not observed or not missing:
            raise AtomCompletionError("partial residues require both observed and missing canonical atoms")
        if len(observed) != len(self.observed_atoms) or len(missing) != len(self.missing_atoms):
            raise AtomCompletionError("partial residue atom masks contain duplicate identities")
        if observed & missing or (observed | missing) != expected:
            raise AtomCompletionError("partial residue masks must partition its canonical heavy atoms exactly")

    @property
    def kind(self) -> CompletionKind:
        return (CompletionKind.BACKBONE if any(atom.atom_name in BACKBONE for atom in self.missing_atoms)
                else CompletionKind.SIDECHAIN)


@dataclass(frozen=True, slots=True)
class AtomCompletionMasks:
    source_sha256: str
    fixed_atoms: tuple[AtomIdentity, ...]
    gap_generated_atoms: tuple[AtomIdentity, ...]
    partial_residues: tuple[PartialResidueMask, ...]
    revision: str = REVISION

    def __post_init__(self) -> None:
        if self.revision != REVISION:
            raise AtomCompletionError("unknown partial-residue mask revision")
        if len(self.source_sha256) != 64 or any(c not in "0123456789abcdef" for c in self.source_sha256):
            raise AtomCompletionError("partial-residue masks require a lowercase SHA-256 source binding")
        fixed, gaps, missing = set(self.fixed_atoms), set(self.gap_generated_atoms), set(self.missing_atoms)
        if len(fixed) != len(self.fixed_atoms) or len(gaps) != len(self.gap_generated_atoms) or len(missing) != len(self.missing_atoms):
            raise AtomCompletionError("atom ownership masks contain duplicate identities")
        if fixed & gaps or fixed & missing or gaps & missing:
            raise AtomCompletionError("observed, whole-gap and partial-residue atom ownership must be disjoint")
        residues = [item.residue for item in self.partial_residues]
        if len(set(residues)) != len(residues):
            raise AtomCompletionError("partial residue identities must be unique")
        positions = [(item.residue.chain, item.target_index) for item in self.partial_residues]
        if len(set(positions)) != len(positions):
            raise AtomCompletionError("partial residue target positions must be unique")
        gap_residues = {ResidueIdentity(a.chain, a.residue_number, a.insertion_code) for a in gaps}
        if gap_residues & set(residues):
            raise AtomCompletionError("partial residues must not also be owned by whole-residue gaps")
        if any(not set(item.observed_atoms) <= fixed for item in self.partial_residues):
            raise AtomCompletionError("every observed atom in a partial residue must remain fixed")

    @property
    def missing_atoms(self) -> tuple[AtomIdentity, ...]:
        return tuple(atom for residue in self.partial_residues for atom in residue.missing_atoms)


def plan_atom_completion(
    source: bytes,
    *,
    source_sha256: str,
    target_sequences: tuple[TargetSequence, ...],
    sequence_placements: tuple[SequencePlacement, ...],
    fixed_atoms: tuple[AtomIdentity, ...],
    gap_generated_atoms: tuple[AtomIdentity, ...] = (),
) -> AtomCompletionMasks:
    """Inventory only absent canonical atoms; preserve every observed source atom.

    Gapless inputs can use the same planning service before a runner contract is
    built. Backbone defects are classified, not silently accepted for sampling.
    OXT is not a canonical side-chain atom and is never invented by this plan.
    HETATM observations remain fixed; planning does not qualify heterogen context.
    """
    from dvbfixer.model.diffusion.scope import DiffusionScopeError, _parse_coordinate_records

    digest = hashlib.sha256(source).hexdigest()
    if digest != source_sha256.lower():
        raise AtomCompletionError("atom-completion source digest mismatch")
    if not target_sequences or not sequence_placements:
        raise AtomCompletionError("atom-completion inventory requires observed protein targets")
    try:
        records, names, altlocs = _parse_coordinate_records(source.decode("utf-8").splitlines())
    except UnicodeDecodeError as exc:
        raise AtomCompletionError("atom-completion source must be UTF-8 PDB text") from exc
    except DiffusionScopeError as exc:
        raise AtomCompletionError(f"invalid atom-completion source: {exc}") from exc
    if altlocs:
        raise AtomCompletionError("select alternate locations before atom-completion planning")
    heavy = {record.identity for record in records if record.element != "H"}
    if heavy != set(fixed_atoms):
        raise AtomCompletionError("fixed atom mask must cover every observed source heavy atom exactly")
    targets = {target.chain: target for target in target_sequences}
    if len(targets) != len(target_sequences) or len({p.chain for p in sequence_placements}) != len(sequence_placements):
        raise AtomCompletionError("target and placement chains must be unique")
    if targets.keys() != {p.chain for p in sequence_placements}:
        raise AtomCompletionError("targets and placements must cover the same case-sensitive chains")
    protein_residues = {ResidueIdentity(r.identity.chain, r.identity.residue_number, r.identity.insertion_code)
                        for r in records if r.record_name == "ATOM"}
    if not protein_residues:
        raise AtomCompletionError("atom-completion inventory requires observed protein atoms")
    placed = [residue for placement in sequence_placements for residue in placement.observed_residues]
    if len(set(placed)) != len(placed) or set(placed) != protein_residues:
        raise AtomCompletionError("placements must cover each observed protein residue exactly once")
    by_residue: dict[ResidueIdentity, set[AtomIdentity]] = {}
    for record in records:
        if record.element != "H":
            residue = ResidueIdentity(record.identity.chain, record.identity.residue_number, record.identity.insertion_code)
            by_residue.setdefault(residue, set()).add(record.identity)
    partials = []
    for placement in sequence_placements:
        target = targets[placement.chain]
        if placement.target_length != len(target.sequence):
            raise AtomCompletionError("placement length differs from target sequence")
        for index, residue in zip(placement.observed_target_indices, placement.observed_residues, strict=True):
            expected_name = _CANONICAL_ONE_TO_THREE.get(target.sequence[index])
            if expected_name != names[residue]:
                raise AtomCompletionError("atom completion cannot reinterpret a substitution or noncanonical residue")
            canonical = canonical_heavy_atom_identities(residue, expected_name)
            observed = by_residue[residue]
            extras = observed - set(canonical)
            if any(atom.atom_name != "OXT" or index != len(target.sequence) - 1 for atom in extras):
                raise AtomCompletionError("observed residue contains unsupported heavy atoms or nonterminal OXT")
            missing = tuple(atom for atom in canonical if atom not in observed)
            if missing:
                partials.append(PartialResidueMask(residue, index, expected_name,
                    tuple(atom for atom in canonical if atom in observed), missing))
    return AtomCompletionMasks(digest, fixed_atoms, gap_generated_atoms, tuple(partials))


def plan_request_atom_completion(request: DiffusionRequest, source: bytes) -> AtomCompletionMasks:
    """Inspect an existing gap request without changing its frozen atom ownership."""
    return plan_atom_completion(source, source_sha256=request.normalized_pdb.sha256,
        target_sequences=request.target_sequences, sequence_placements=request.sequence_placements,
        fixed_atoms=request.fixed_atoms, gap_generated_atoms=request.generated_atoms)


def audit_atom_completion_masks(plan: AtomCompletionMasks, request: DiffusionRequest, source: bytes) -> None:
    """Recompute masks from original bytes, never trust a caller's missing list."""
    if plan != plan_request_atom_completion(request, source):
        raise AtomCompletionError("atom-completion masks differ from independent source inventory")
