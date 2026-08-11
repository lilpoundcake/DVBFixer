"""Canonical residue-name classifications shared across DVBfixer workflows."""

from __future__ import annotations

PDB_SUGAR_RESNAMES = frozenset({
    "NAG", "NDG", "BGL", "BMA", "MAN", "FUC", "FUL", "GAL", "BGC",
    "GLC", "SIA", "NGA", "A2G", "AFU", "AMA", "BGA", "AGA", "XYS",
    "XYP", "RIB", "RIP", "ARA", "GCU", "IDS", "RAM", "NAN",
})

CHARMM_SUGAR_RESNAMES = frozenset({
    "BGLC", "AGLC", "BMAN", "AMAN", "BGAL", "AGAL", "BGLCNA",
    "AGLCNA", "BGALNA", "AGALNA", "AFUC", "BFUC", "ANE5", "BNE5",
    "ANE5AC", "BNE5AC", "BXYL", "AXYL", "ARIB", "BRIB", "AGLCA",
    "BGLCA", "AIDO", "BIDO", "AIDOA", "BIDOA", "ARHA", "BRHA",
})

SIALIC_RESNAMES = frozenset({"SIA", "ANE5", "BNE5", "ANE5AC", "BNE5AC"})
GLYCAM_PROTEIN_RESNAMES = frozenset({"NLN", "OLS", "OLT"})
GLYCAM_CAP_RESNAMES = frozenset({"ROH", "OME", "TBT", "CMET"})

_GLYCAM_LINKAGE_CHARS = frozenset("0123456789VWUZXYTSRQPvwuzxytsr")
_GLYCAM_ANOMER_CHARS = frozenset({"A", "B"})


def is_glycam_sugar(resname: str) -> bool:
    """Return whether *resname* is a GLYCAM sugar or reducing-end cap."""
    return (
        resname in GLYCAM_CAP_RESNAMES
        or (
            len(resname) == 3
            and resname[0] in _GLYCAM_LINKAGE_CHARS
            and resname[1].isalpha()
            and resname[2] in _GLYCAM_ANOMER_CHARS
        )
    )


def is_sugar_resname(resname: str) -> bool:
    """Return whether *resname* uses a supported PDB, CHARMM, or GLYCAM name."""
    return (
        resname in PDB_SUGAR_RESNAMES
        or resname in CHARMM_SUGAR_RESNAMES
        or is_glycam_sugar(resname)
    )


def is_pdb_or_glycam_sugar(resname: str) -> bool:
    """Return whether a residue needs the AMBER/GLYCAM carbohydrate path."""
    return resname in PDB_SUGAR_RESNAMES or is_glycam_sugar(resname)
