"""Shared constants and helpers for ``dvbfixer diagnose``."""

from __future__ import annotations

# Waters. Kept separate from ions (which participate in real
# coordination chemistry and shouldn't be blindly skipped).
WATER_RESIDUES: frozenset[str] = frozenset({
    "HOH", "WAT", "SOL", "H2O",
    "TIP3", "TIP4", "TIP5",
    "SPC", "SPCE",
    "DOD",  # D2O
})


def is_water(resname: str) -> bool:
    return resname.strip().upper() in WATER_RESIDUES


# Non-standard amino acids the chemistry checks should recognise as
# proteinogenic. Adds selenocysteine (SEC), pyrrolysine (PYL),
# selenomethionine (MSE), and common phospho-residues.
NON_STANDARD_AAS: frozenset[str] = frozenset({
    "SEC", "PYL", "MSE",
    "SEP", "TPO", "PTR",   # phosphoserine / threonine / tyrosine
    "CSO", "CSD", "CME",   # oxidised cysteines
    "HYP",                  # hydroxyproline (already handled for cis-PRO)
})


# Elements that can act as H-bond acceptors AND donors (when they
# carry an H, they are donors; without an H they can still hydrogen
# bond as acceptors).
HBOND_HEAVY_ELEMENTS: frozenset[str] = frozenset({"N", "O", "S", "F"})
