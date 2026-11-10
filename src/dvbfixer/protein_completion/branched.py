"""Illustrative paired VAL/THR branch geometry, not empirical rotamer statistics.

THR CA/OG1/CG2 signed volume at CB is positive for natural (2S,3R)
threonine. VAL has no CB stereocentre (two identical methyl groups); the
negative CA/CG1/CG2 volume records PDB stereospecific atom-name convention.
"""
from __future__ import annotations

from dataclasses import dataclass

from dvbfixer.protein_completion.contracts import Point
from dvbfixer.protein_completion.geometry import cross, dot, sub

BRANCHED_REVISION = "canonical-branched-sidechain-staggered-v1"
TETRAHEDRAL_ANGLE = 109.47122063449069


@dataclass(frozen=True)
class PairedBranch:
    first: str
    first_length: float
    second: str
    second_length: float
    second_chi_offset: float
    expected_named_volume_sign: int


BRANCHES = {
    "VAL": PairedBranch("CG1", 1.53, "CG2", 1.53, 120., -1),
    "THR": PairedBranch("OG1", 1.43, "CG2", 1.53, -120., 1),
}


def named_beta_volume(ca: Point, cb: Point, first: Point, second: Point) -> float:
    """Independent measurement; does not reconstruct or reflect observations."""
    return dot(cross(sub(ca, cb), sub(first, cb)), sub(second, cb))
