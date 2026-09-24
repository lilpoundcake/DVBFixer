"""Fail-closed mapping from the Boltz-2 atom axis to DVBFixer identities."""

from __future__ import annotations

from typing import Any

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest, ResidueIdentity

ONE_TO_THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}


def target_residue_map(request: DiffusionRequest) -> dict[int, ResidueIdentity]:
    """Map every zero-based target position to its exact requested identity."""
    if len(request.target_sequences) != 1 or len(request.sequence_placements) != 1:
        raise ValueError("Boltz checkpoint smoke requires exactly one target chain")
    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    if placement.chain != target.chain or placement.target_length != len(target.sequence):
        raise ValueError("target sequence and sequence placement are inconsistent")

    mapping = dict(zip(placement.observed_target_indices, placement.observed_residues))
    for gap in request.gaps:
        if gap.chain != target.chain:
            raise ValueError("all gaps must belong to the single target chain")
        for target_index, residue in zip(
            range(gap.target_interval.start, gap.target_interval.stop),
            gap.generated_residues,
        ):
            if target_index in mapping:
                raise ValueError(f"target residue index {target_index} is mapped more than once")
            mapping[target_index] = residue

    expected = set(range(len(target.sequence)))
    if set(mapping) != expected:
        raise ValueError("request does not map every target residue exactly once")
    return mapping


def model_atom_axis(
    tokenized: Any,
    request: DiffusionRequest,
    *,
    chain_name_by_asym_id: dict[int, str],
) -> tuple[tuple[AtomIdentity, ...], tuple[str, ...], tuple[int, ...]]:
    """Build the exact unpadded model atom axis from Boltz token traversal."""
    target = request.target_sequences[0]
    residue_map = target_residue_map(request)
    if set(chain_name_by_asym_id.values()) != {target.chain}:
        raise ValueError("processed Boltz chains do not match the requested target chain")

    identities: list[AtomIdentity] = []
    residue_names: list[str] = []
    token_indices: list[int] = []
    target_index = 0
    expected_atom_start = 0
    seen_tokens: set[int] = set()
    for token_position, token in enumerate(tokenized.tokens):
        token_index = int(token["token_idx"])
        if token_index != token_position or token_index in seen_tokens:
            raise ValueError("Boltz token indices are not unique and contiguous")
        seen_tokens.add(token_index)

        asym_id = int(token["asym_id"])
        if chain_name_by_asym_id.get(asym_id) != target.chain:
            raise ValueError("Boltz token belongs to an unexpected chain")
        if target_index >= len(target.sequence):
            raise ValueError("Boltz token axis contains too many target residues")
        if int(token["res_idx"]) != target_index:
            raise ValueError("Boltz residue indices do not match target sequence order")

        expected_residue_name = ONE_TO_THREE.get(target.sequence[target_index])
        residue_name = str(token["res_name"])
        if expected_residue_name is None or residue_name != expected_residue_name:
            raise ValueError(
                f"Boltz residue {target_index + 1} is {residue_name}, "
                f"expected {expected_residue_name}"
            )
        atom_start = int(token["atom_idx"])
        atom_count = int(token["atom_num"])
        if atom_start != expected_atom_start:
            raise ValueError("Boltz token atom spans are not contiguous and ordered")
        atoms = tokenized.structure.atoms[atom_start : atom_start + atom_count]
        if len(atoms) != atom_count or atom_count <= 0:
            raise ValueError("Boltz token has an invalid atom span")

        residue = residue_map[target_index]
        atom_names = [str(atom["name"]).strip() for atom in atoms]
        if any(not name for name in atom_names) or len(set(atom_names)) != len(atom_names):
            raise ValueError("Boltz residue has blank or duplicate atom names")
        for atom_name in atom_names:
            identities.append(
                AtomIdentity(
                    residue.chain,
                    residue.residue_number,
                    residue.insertion_code,
                    atom_name,
                )
            )
            residue_names.append(residue_name)
            token_indices.append(token_index)
        expected_atom_start += atom_count
        target_index += 1

    if target_index != len(target.sequence):
        raise ValueError("Boltz token axis omits target residues")
    axis = tuple(identities)
    if len(set(axis)) != len(axis):
        raise ValueError("Boltz atom axis maps to duplicate DVBFixer identities")
    required = set(request.fixed_atoms) | set(request.generated_atoms)
    missing = required - set(axis)
    if missing:
        raise ValueError(f"Boltz atom axis omits {len(missing)} requested atoms")
    return axis, tuple(residue_names), tuple(token_indices)


def validate_feature_axis(
    atom_axis: tuple[AtomIdentity, ...],
    token_indices: tuple[int, ...],
    features: dict[str, Any],
) -> tuple[int, ...]:
    """Validate Boltz feature ordering and return atomic numbers for real atoms."""
    pad_mask = _numpy(features["atom_pad_mask"])
    if pad_mask.ndim != 1 or not np.isin(pad_mask, (0, 1)).all():
        raise ValueError("atom_pad_mask must be a one-dimensional binary mask")
    real_indices = np.flatnonzero(pad_mask.astype(bool))
    if not np.array_equal(real_indices, np.arange(len(atom_axis))):
        raise ValueError("Boltz real atoms must be an unpadded prefix matching the mapped axis")

    atom_to_token = _numpy(features["atom_to_token"])
    if atom_to_token.ndim != 2 or atom_to_token.shape[0] != len(pad_mask):
        raise ValueError("atom_to_token has an invalid shape")
    feature_tokens = atom_to_token[: len(atom_axis)].argmax(axis=-1)
    if not np.array_equal(feature_tokens, np.asarray(token_indices)):
        raise ValueError("atom_to_token disagrees with token atom spans")

    encoded_names = _numpy(features["ref_atom_name_chars"])
    if encoded_names.ndim != 3 or encoded_names.shape[:2] != (len(pad_mask), 4):
        raise ValueError("ref_atom_name_chars has an invalid shape")
    decoded_names = tuple(_decode_atom_name(row) for row in encoded_names[: len(atom_axis)])
    expected_names = tuple(identity.atom_name for identity in atom_axis)
    if decoded_names != expected_names:
        raise ValueError("ref_atom_name_chars disagrees with the mapped atom identities")

    elements = _numpy(features["ref_element"])
    if elements.ndim != 2 or elements.shape[0] != len(pad_mask):
        raise ValueError("ref_element has an invalid shape")
    atomic_numbers = tuple(int(value) for value in elements[: len(atom_axis)].argmax(axis=-1))
    if any(number <= 0 for number in atomic_numbers):
        raise ValueError("Boltz real atom has an invalid atomic number")
    return atomic_numbers


def _numpy(value: Any) -> np.ndarray:
    if hasattr(value, "detach"):
        value = value.detach().cpu().numpy()
    return np.asarray(value)


def _decode_atom_name(encoded: np.ndarray) -> str:
    codes = encoded.argmax(axis=-1)
    try:
        return "".join(chr(int(code) + 32) for code in codes if code).strip()
    except (TypeError, ValueError) as exc:
        raise ValueError("invalid encoded Boltz atom name") from exc
