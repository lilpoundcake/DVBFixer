"""Unit tests for the LYN HZ3 → HZ1 GROMACS-compat rename helpers.

Coverage:
- ``rename_lyn_hz_for_gromacs(topology)`` — topology-level helper
- ``rename_lyn_hz_for_gromacs_in_pdb_text(path)`` — text-level helper
- Idempotence of both.
- Non-LYN LYS is untouched.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("openmm", reason="needs OpenMM Topology")

from openmm.app import Element, Topology  # noqa: E402

from dvbfixer.ffutils.variants import (  # noqa: E402
    rename_lyn_hz_for_gromacs,
    rename_lyn_hz_for_gromacs_in_pdb_text,
)


def _lyn_topology(atom_names: list[str]) -> Topology:
    top = Topology()
    chain = top.addChain("A")
    r = top.addResidue("LYN", chain)
    for name in atom_names:
        sym = name[0]  # first char is element
        top.addAtom(name, Element.getBySymbol(sym), r)
    return top


def test_topology_helper_renames_hz3_to_hz1() -> None:
    top = _lyn_topology(["N", "CA", "CB", "CG", "CD", "CE", "NZ", "HZ2", "HZ3"])
    n = rename_lyn_hz_for_gromacs(top)
    assert n == 1
    names = [a.name for a in top.atoms()]
    assert "HZ1" in names
    assert "HZ3" not in names
    assert "HZ2" in names   # untouched


def test_topology_helper_idempotent() -> None:
    top = _lyn_topology(["N", "CA", "NZ", "HZ2", "HZ3"])
    n1 = rename_lyn_hz_for_gromacs(top)
    n2 = rename_lyn_hz_for_gromacs(top)
    assert n1 == 1
    assert n2 == 0
    names = [a.name for a in top.atoms()]
    assert "HZ1" in names
    assert "HZ3" not in names


def test_topology_helper_lys_untouched() -> None:
    """A charged LYS with HZ1+HZ2+HZ3 must NOT be modified."""
    top = Topology()
    chain = top.addChain("A")
    r = top.addResidue("LYS", chain)
    for name in ["NZ", "HZ1", "HZ2", "HZ3"]:
        top.addAtom(name, Element.getBySymbol(name[0]), r)
    n = rename_lyn_hz_for_gromacs(top)
    assert n == 0
    names = [a.name for a in top.atoms()]
    assert names == ["NZ", "HZ1", "HZ2", "HZ3"]


_LYN_PDB = """\
ATOM      1  N   LYN A 100      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  LYN A 100      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  NZ  LYN A 100      15.000  10.000  10.000  1.00  0.00           N
ATOM      4  HZ2 LYN A 100      15.500  10.500  10.500  1.00  0.00           H
ATOM      5  HZ3 LYN A 100      15.500   9.500  10.500  1.00  0.00           H
ATOM      6  N   LYS A 101      20.000  10.000  10.000  1.00  0.00           N
ATOM      7  HZ3 LYS A 101      25.500   9.500  10.500  1.00  0.00           H
END
"""


def test_text_helper_renames_only_lyn_hz3(tmp_path: Path) -> None:
    """Text-level helper renames LYN HZ3 → HZ1 but leaves LYS HZ3 alone."""
    pdb = tmp_path / "lyn_test.pdb"
    pdb.write_text(_LYN_PDB)
    n = rename_lyn_hz_for_gromacs_in_pdb_text(pdb)
    assert n == 1
    content = pdb.read_text()
    # LYN got renamed.
    assert " HZ1 LYN A 100" in content
    assert " HZ3 LYN A 100" not in content
    # LYS HZ3 must be preserved.
    assert " HZ3 LYS A 101" in content


def test_text_helper_idempotent(tmp_path: Path) -> None:
    pdb = tmp_path / "lyn_test.pdb"
    pdb.write_text(_LYN_PDB)
    rename_lyn_hz_for_gromacs_in_pdb_text(pdb)
    n2 = rename_lyn_hz_for_gromacs_in_pdb_text(pdb)
    assert n2 == 0
