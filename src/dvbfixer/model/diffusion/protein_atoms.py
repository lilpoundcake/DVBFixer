"""Canonical protein heavy-atom identities, independent of scientific runtimes.

Scope admission, atom-mask planning and wire contracts share this one catalog.
"""

from __future__ import annotations

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity

ATOM37_ORDER = (
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD", "CD1", "CD2",
    "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3", "NE", "NE1", "NE2",
    "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT",
)

_CANONICAL_ONE_TO_THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}
CANONICAL_HEAVY_ATOMS: dict[str, tuple[str, ...]] = {
    "ALA": ("N", "CA", "CB", "C", "O"),
    "ARG": ("N", "CA", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "C", "O"),
    "ASN": ("N", "CA", "CB", "CG", "OD1", "ND2", "C", "O"),
    "ASP": ("N", "CA", "CB", "CG", "OD1", "OD2", "C", "O"),
    "CYS": ("N", "CA", "CB", "SG", "C", "O"),
    "GLN": ("N", "CA", "CB", "CG", "CD", "OE1", "NE2", "C", "O"),
    "GLU": ("N", "CA", "CB", "CG", "CD", "OE1", "OE2", "C", "O"),
    "GLY": ("N", "CA", "C", "O"),
    "HIS": ("N", "CA", "CB", "CG", "ND1", "CE1", "NE2", "CD2", "C", "O"),
    "ILE": ("N", "CA", "CB", "CG2", "CG1", "CD1", "C", "O"),
    "LEU": ("N", "CA", "CB", "CG", "CD1", "CD2", "C", "O"),
    "LYS": ("N", "CA", "CB", "CG", "CD", "CE", "NZ", "C", "O"),
    "MET": ("N", "CA", "CB", "CG", "SD", "CE", "C", "O"),
    "PHE": ("N", "CA", "CB", "CG", "CD1", "CE1", "CZ", "CE2", "CD2", "C", "O"),
    "PRO": ("N", "CD", "CG", "CB", "CA", "C", "O"),
    "SER": ("N", "CA", "CB", "OG", "C", "O"),
    "THR": ("N", "CA", "CB", "CG2", "OG1", "C", "O"),
    "TRP": ("N", "CA", "CB", "CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2", "C", "O"),
    "TYR": ("N", "CA", "CB", "CG", "CD1", "CE1", "CZ", "OH", "CE2", "CD2", "C", "O"),
    "VAL": ("N", "CA", "CB", "CG1", "CG2", "C", "O"),
}


def canonical_heavy_atom_identities(residue: ResidueIdentity, residue_name: str) -> tuple[AtomIdentity, ...]:
    """Return the complete canonical heavy-atom identity set in template order."""
    try:
        atom_names = CANONICAL_HEAVY_ATOMS[residue_name]
    except KeyError as exc:
        raise ValueError(f"unknown canonical residue name: {residue_name}") from exc
    return tuple(AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, name) for name in atom_names)
