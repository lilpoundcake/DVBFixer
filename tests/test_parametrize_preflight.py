from __future__ import annotations

from pathlib import Path

import pytest

from dvbfixer.parametrize import _validate_single_small_molecule_input


def _atom(serial: int, residue: int, element: str = "C") -> str:
    return (
        f"HETATM{serial:5d}  {element:<2}  LIG A{residue:4d}    "
        f"   0.000   0.000   0.000  1.00  0.00          {element:>2}\n"
    )


def test_parametrize_rejects_multi_residue_pdb(tmp_path: Path) -> None:
    path = tmp_path / "complex.pdb"
    path.write_text(_atom(1, 1) + _atom(2, 2))
    with pytest.raises(ValueError, match="one isolated small-molecule"):
        _validate_single_small_molecule_input(path)


def test_parametrize_rejects_metal_pdb(tmp_path: Path) -> None:
    path = tmp_path / "metal.pdb"
    path.write_text(_atom(1, 1, "FE"))
    with pytest.raises(ValueError, match="metal-containing"):
        _validate_single_small_molecule_input(path)
