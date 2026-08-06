"""Biological-assembly extraction tests for ``dvbfixer split``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dvbfixer.biological_assembly import (
    AssemblyError,
    parse_biological_assemblies,
    render_biological_assembly,
)
from dvbfixer.split_chains import CHAIN_IDS, main


def _chains(path: Path) -> set[str]:
    return {
        line[21] for line in path.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    }


def _synthetic_pdb() -> list[str]:
    return [
        "HEADER    BIOMT TEST\n",
        "REMARK 350 BIOMOLECULE: 1\n",
        "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A\n",
        "REMARK 350   BIOMT1   1  0.000000 -1.000000  0.000000       10.00000\n",
        "REMARK 350   BIOMT2   1  1.000000  0.000000  0.000000       20.00000\n",
        "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000       30.00000\n",
        "ATOM      1  CA  ALA A  10       1.000   2.000   3.000  1.00 20.00           C  \n",
        "ANISOU    1  CA  ALA A  10     1000   2000   3000    100    200    300       C  \n",
        "END\n",
    ]


def test_8xj0_all_assemblies(
    biological_assembly_pdb: Path, tmp_workdir: Path,
) -> None:
    prefix = tmp_workdir / "fab.pdb"
    main([str(biological_assembly_pdb), "--assembly", "all", "-o", str(prefix)])
    expected = {"1": {"A", "B"}, "2": {"C", "D"}, "3": {"E", "F"}, "4": {"G", "H"}}
    for identifier, chains in expected.items():
        output = tmp_workdir / f"fab_assembly_{identifier}.pdb"
        assert output.is_file()
        assert _chains(output) == chains
        text = output.read_text()
        assert "REMARK 999 DVBfixer BIOLOGICAL ASSEMBLY" in text
        assert sum(line.startswith("SSBOND") for line in text.splitlines()) == 6


def test_8xj0_one_assembly_honours_exact_output(
    biological_assembly_pdb: Path, tmp_workdir: Path,
) -> None:
    output = tmp_workdir / "chosen.pdb"
    main([str(biological_assembly_pdb), "--assembly", "2", "-o", str(output)])
    assert _chains(output) == {"C", "D"}


def test_biomt_rotation_translation_and_anisou() -> None:
    lines = _synthetic_pdb()
    assembly = parse_biological_assemblies(lines)["1"]
    rendered = render_biological_assembly(lines, assembly, CHAIN_IDS)
    atom = next(line for line in rendered if line.startswith("ATOM"))
    xyz = np.array([float(atom[30:38]), float(atom[38:46]), float(atom[46:54])])
    assert xyz == pytest.approx([8.0, 21.0, 33.0])
    anisou = next(line for line in rendered if line.startswith("ANISOU"))
    assert [int(anisou[start:start + 7]) for start in (28, 35, 42)] == [2000, 1000, 3000]
    assert atom[22:26] == "  10"  # assembly mode preserves deposited numbering

    renumbered = render_biological_assembly(lines, assembly, CHAIN_IDS, renumber=True)
    renumbered_atom = next(line for line in renumbered if line.startswith("ATOM"))
    assert renumbered_atom[22:26] == "   1"


def test_missing_and_incomplete_assembly_records_fail() -> None:
    with pytest.raises(AssemblyError, match="no biological assemblies"):
        parse_biological_assemblies(["ATOM      1  CA  ALA A   1\n"])
    with pytest.raises(AssemblyError, match="incomplete"):
        parse_biological_assemblies([
            "REMARK 350 BIOMOLECULE: 1\n",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A\n",
            "REMARK 350   BIOMT1   1  1.000000 0.000000 0.000000 0.00000\n",
        ])
