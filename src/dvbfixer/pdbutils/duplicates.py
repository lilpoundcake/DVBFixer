"""Detect accidentally concatenated coordinate-identical protein chains."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

_PROTEIN = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM",
    "LYN", "MSE", "NLN", "OLS", "OLT",
}


def duplicate_protein_chain_coordinates(
    pdb_path: str | Path, *, tolerance: float = 0.001,
) -> list[tuple[str, str, int]]:
    """Return chain pairs whose protein atoms occupy identical coordinates.

    Atom names and residue names must match in file order. Residue numbers may
    differ because concatenation/export bugs sometimes renumber the duplicate.
    The default tolerance corresponds to PDB coordinate precision.
    """
    chains: dict[str, list[tuple[str, str, float, float, float]]] = defaultdict(list)
    with Path(pdb_path).open() as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                continue
            resname = line[17:20].strip()
            if resname not in _PROTEIN:
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            chains[line[21]].append((resname, line[12:16].strip(), *xyz))

    duplicates: list[tuple[str, str, int]] = []
    chain_ids = list(chains)
    for index, left_id in enumerate(chain_ids):
        left = chains[left_id]
        if not left:
            continue
        for right_id in chain_ids[index + 1:]:
            right = chains[right_id]
            if len(left) != len(right):
                continue
            if all(
                la[:2] == ra[:2]
                and abs(la[2] - ra[2]) <= tolerance
                and abs(la[3] - ra[3]) <= tolerance
                and abs(la[4] - ra[4]) <= tolerance
                for la, ra in zip(left, right)
            ):
                duplicates.append((left_id or "_", right_id or "_", len(left)))
    return duplicates
