"""Tests for `dvbfixer renumber`.

Regression coverage for a bare/minimal ``TER\\n`` record (no serial,
resname, or padding — just 4 characters) silently losing its trailing
newline during renumbering. `line[11:]` on such a short line returns
an empty string (Python doesn't error on an out-of-range slice), so
the constructed output line had no newline at all and ran directly
into the next physical line (e.g. ``TER    4748HETATM 4749  C1  ...``).
Every downstream parser's line-start check then silently dropped that
atom (and everything else on the merged line) — in the reported case,
the exact anomeric carbon (C1) that forms a real glycosidic bond,
making the bond undetectable by any tool downstream no matter how
good its CONECT/distance logic was.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def test_renumber_preserves_bare_ter_line_atom_count(tmp_workdir: Path) -> None:
    """A synthetic minimal repro: a bare `TER\\n` immediately followed
    by a HETATM line must not merge the two into one physical line."""
    input_pdb = tmp_workdir / "bare_ter.pdb"
    input_pdb.write_text(
        "ATOM      1  N   ALA A   1      11.104  13.207  10.454  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1      11.804  12.500   9.400  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1      13.300  12.500   9.400  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1      13.900  13.500   9.400  1.00  0.00           O\n"
        "TER\n"
        "HETATM    5  C1  NAG D   1      20.000  20.000  20.000  1.00  0.00           C\n"
        "HETATM    6  C2  NAG D   1      21.000  20.000  20.000  1.00  0.00           C\n"
        "END\n"
    )
    output = tmp_workdir / "out.pdb"
    proc = subprocess.run(
        ["dvbfixer", "renumber", str(input_pdb), "-o", str(output)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"renumber failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    lines = output.read_text().splitlines()
    atom_het_lines = [line for line in lines
                      if line.startswith(("ATOM  ", "HETATM"))]
    # 4 protein atoms + 2 HETATM atoms = 6 total. If the TER/HETATM
    # merge bug regresses, C1 (or the whole merged line) goes missing
    # and this count drops.
    assert len(atom_het_lines) == 6, (
        f"expected 6 ATOM/HETATM lines, got {len(atom_het_lines)}: "
        f"{atom_het_lines}"
    )
    atom_names = {line[12:16].strip() for line in atom_het_lines
                  if line.startswith("HETATM")}
    assert atom_names == {"C1", "C2"}, (
        f"expected both NAG atoms present, got {atom_names}"
    )
    # A regression merges the bare TER directly into the next line
    # ("TER    4748HETATM ..." — no space, no newline between them).
    for line in lines:
        assert "TERHETATM" not in line and "TERATOM" not in line, (
            f"TER line merged into the next record: {line!r}"
        )


@pytest.mark.slow
def test_renumber_preserves_all_atoms_on_real_glycoprotein(
    glycoprot_underannotated_conect_pdb: Path, tmp_workdir: Path,
) -> None:
    """This fixture's own raw PDB has a bare TER immediately before a
    sugar chain's HETATM block — the exact real-world trigger for the
    bug above. Total atom count must be preserved through renumber."""
    output = tmp_workdir / "renum.pdb"
    proc = subprocess.run(
        ["dvbfixer", "renumber", str(glycoprot_underannotated_conect_pdb),
         "-o", str(output)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"renumber failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    n_in = sum(1 for line in glycoprot_underannotated_conect_pdb.read_text()
               .splitlines() if line.startswith(("ATOM  ", "HETATM")))
    n_out = sum(1 for line in output.read_text().splitlines()
                if line.startswith(("ATOM  ", "HETATM")))
    assert n_out == n_in, (
        f"expected {n_in} ATOM/HETATM lines preserved, got {n_out} "
        f"after renumber"
    )
