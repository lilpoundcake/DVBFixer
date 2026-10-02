from __future__ import annotations

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates


def _atom_line(serial: int, residue_number: int, atom_name: str, x: float) -> str:
    return (
        f"ATOM  {serial:5d}  {atom_name:<3} ALA A{residue_number:4d}    "
        f"{x:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00  0.00           C\n"
    )


def test_rewrite_generated_coordinates_preserves_non_coordinate_text() -> None:
    fixed = _atom_line(1, 1, "CA", 1.0)
    generated = _atom_line(2, 2, "CA", 4.0)
    text = "HEADER preserved\n" + fixed + generated + "END\n"

    rewritten = rewrite_generated_coordinates(
        text,
        {AtomIdentity("A", "2", "", "CA"): np.asarray([7.0, 8.0, 9.0])},
    )

    original_lines = text.splitlines(keepends=True)
    rewritten_lines = rewritten.splitlines(keepends=True)
    assert rewritten_lines[:2] == original_lines[:2]
    assert rewritten_lines[2][:30] == original_lines[2][:30]
    assert rewritten_lines[2][30:54] == "   7.000   8.000   9.000"
    assert rewritten_lines[2][54:] == original_lines[2][54:]
    assert rewritten_lines[3] == original_lines[3]


def test_rewrite_generated_coordinates_rejects_missing_atoms() -> None:
    coordinates = {AtomIdentity("A", "2", "", "CA"): np.asarray([7.0, 8.0, 9.0])}
    with pytest.raises(ValueError, match="omits 1 generated atoms"):
        rewrite_generated_coordinates("END\n", coordinates)
