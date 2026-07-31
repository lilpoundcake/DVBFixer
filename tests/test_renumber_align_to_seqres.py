"""Regression test for the 2026-07-31 audit's `renumber.py` fix.

`align_to_seqres`'s unbounded forward search let a single point
mutation (ATOM residue name differs from wild-type SEQRES at that one
position) overshoot indefinitely looking for a same-name match — and
since `j` was never advanced on a miss, every residue after the
mutation got silently misplaced to the sequence tail instead of its
true position.
"""
from __future__ import annotations

from dvbfixer.renumber import align_to_seqres


def test_align_to_seqres_tolerates_point_mutation() -> None:
    seqres = ["ALA", "GLY", "SER", "THR", "VAL"]
    atom_residues = [
        (1, " ", "ALA"),
        (2, " ", "GLY"),
        (3, " ", "ALA"),  # point mutation: wild-type SEQRES is SER
        (4, " ", "THR"),
        (5, " ", "VAL"),
    ]

    mapping = align_to_seqres(atom_residues, seqres)

    assert mapping[(1, " ")] == 1
    assert mapping[(2, " ")] == 2
    assert mapping[(3, " ")] == 3, (
        "point-mutated residue must map to its true position, not fall "
        "through to the sequence tail"
    )
    assert mapping[(4, " ")] == 4, (
        "residues AFTER the mutation must not be shifted just because "
        "the mutation itself needed a substitution"
    )
    assert mapping[(5, " ")] == 5


def test_align_to_seqres_still_resyncs_genuine_gap() -> None:
    """A real ATOM/SEQRES gap (missing residues, not a substitution)
    must still resync via the forward search."""
    seqres = ["ALA", "GLY", "SER", "THR", "VAL", "LEU"]
    atom_residues = [
        (1, " ", "ALA"),
        (2, " ", "GLY"),
        # SER (3), THR (4) missing from ATOM records — genuine gap
        (5, " ", "VAL"),
        (6, " ", "LEU"),
    ]

    mapping = align_to_seqres(atom_residues, seqres)

    assert mapping[(1, " ")] == 1
    assert mapping[(2, " ")] == 2
    assert mapping[(5, " ")] == 5
    assert mapping[(6, " ")] == 6


def test_align_to_seqres_multiple_point_mutations() -> None:
    """Several separated point mutations must each resolve independently,
    not compound into a single runaway misalignment."""
    seqres = ["ALA", "GLY", "SER", "THR", "VAL", "LEU", "PRO"]
    atom_residues = [
        (1, " ", "ALA"),
        (2, " ", "GLY"),
        (3, " ", "CYS"),  # mutation: wild-type SER
        (4, " ", "THR"),
        (5, " ", "VAL"),
        (6, " ", "TRP"),  # mutation: wild-type LEU
        (7, " ", "PRO"),
    ]

    mapping = align_to_seqres(atom_residues, seqres)

    assert [mapping[(i, " ")] for i in range(1, 8)] == [1, 2, 3, 4, 5, 6, 7]
