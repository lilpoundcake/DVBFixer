"""Tests for deterministic diffusion sequence placement and masks."""

from __future__ import annotations

import pytest

from dvbfixer.model.diffusion.contract import ResidueIdentity
from dvbfixer.model.diffusion.masks import (
    DiffusionMaskError,
    ObservedResidue,
    build_diffusion_masks,
    build_sequence_placement,
)


def _residue(chain: str, number: str, code: str, *, icode: str = "") -> ObservedResidue:
    return ObservedResidue(
        identity=ResidueIdentity(chain, number, icode),
        one_letter_code=code,
        atoms=("N", "CA", "C", "O"),
    )


def test_build_masks_reuses_shared_internal_gap_placement() -> None:
    observed = (
        _residue("D", "10", "A"),
        _residue("D", "11", "C"),
        _residue("D", "20", "D"),
        _residue("D", "21", "E"),
    )

    masks = build_diffusion_masks("D", observed, "ACGGGSDE")

    assert masks.sequence_placement.observed_target_indices == (0, 1, 6, 7)
    assert len(masks.gaps) == 1
    gap = masks.gaps[0]
    assert (gap.target_interval.start, gap.target_interval.stop) == (2, 6)
    assert gap.left_anchor == ResidueIdentity("D", "11")
    assert gap.right_anchor == ResidueIdentity("D", "20")
    assert gap.generated_residues == tuple(
        ResidueIdentity("D", str(number)) for number in range(3, 7)
    )
    assert gap.movable_junction_residues == (
        ResidueIdentity("D", "11"),
        ResidueIdentity("D", "20"),
        *gap.generated_residues,
    )
    assert len(masks.fixed_atoms) == 16


def test_movable_flank_window_is_bounded_and_deterministic() -> None:
    observed = tuple(
        _residue("D", str(number), code)
        for number, code in zip((10, 11, 12, 20, 21, 22), "ABCDEF")
    )

    masks = build_diffusion_masks(
        "D",
        observed,
        "ABCGGGDEF",
        movable_flank_residues=2,
    )

    assert masks.gaps[0].movable_junction_residues == (
        ResidueIdentity("D", "11"),
        ResidueIdentity("D", "12"),
        ResidueIdentity("D", "20"),
        ResidueIdentity("D", "21"),
        ResidueIdentity("D", "4"),
        ResidueIdentity("D", "5"),
        ResidueIdentity("D", "6"),
    )


def test_insertion_codes_and_case_distinct_chains_are_preserved() -> None:
    upper = (
        _residue("D", "81", "A"),
        _residue("D", "82", "C", icode="A"),
        _residue("D", "90", "D"),
        _residue("D", "91", "E"),
    )
    lower = tuple(
        _residue("d", str(number), code)
        for number, code in zip((1, 2, 10, 11), "ACDE")
    )

    upper_masks = build_diffusion_masks("D", upper, "ACGGGDE")
    lower_masks = build_diffusion_masks("d", lower, "ACGGGDE")

    assert upper_masks.sequence_placement.observed_residues[1].insertion_code == "A"
    assert upper_masks.fixed_atoms[4].chain == "D"
    assert lower_masks.fixed_atoms[4].chain == "d"
    assert upper_masks.fixed_atoms[4] != lower_masks.fixed_atoms[4]


def test_ambiguous_placement_is_rejected() -> None:
    observed = (
        _residue("A", "1", "A"),
        _residue("A", "2", "C"),
    )

    with pytest.raises(DiffusionMaskError, match="ambiguous"):
        build_sequence_placement("A", observed, "ACGGGAC")


def test_terminal_one_anchor_gaps_are_rejected() -> None:
    n_terminal_missing = (
        _residue("A", "4", "D"),
        _residue("A", "5", "E"),
        _residue("A", "6", "F"),
    )
    with pytest.raises(DiffusionMaskError, match="N-terminal"):
        build_diffusion_masks("A", n_terminal_missing, "GGGDEF")

    c_terminal_missing = (
        _residue("A", "1", "A"),
        _residue("A", "2", "C"),
        _residue("A", "3", "D"),
    )
    with pytest.raises(DiffusionMaskError, match="C-terminal"):
        build_diffusion_masks("A", c_terminal_missing, "ACDGGG")


def test_gap_length_outside_configured_slice_is_rejected() -> None:
    observed = tuple(
        _residue("A", str(number), code)
        for number, code in zip((1, 2, 3, 20, 21, 22), "MKTWQH")
    )

    with pytest.raises(DiffusionMaskError, match="outside supported range"):
        build_diffusion_masks(
            "A",
            observed,
            "MKT" + "A" * 10 + "WQH",
            maximum_gap_length=9,
        )
