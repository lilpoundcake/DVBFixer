from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("openmm")


_UNCAPPED = """\
ATOM      1  N   ALA A   2       3.555   3.970   0.000  1.00  0.00           N
ATOM      2  H   ALA A   2       2.733   4.556   0.000  1.00  0.00           H
ATOM      3  CA  ALA A   2       4.853   4.614   0.000  1.00  0.00           C
ATOM      4  HA  ALA A   2       5.408   4.316   0.890  1.00  0.00           H
ATOM      5  CB  ALA A   2       5.661   4.221  -1.232  1.00  0.00           C
ATOM      6  HB1 ALA A   2       5.123   4.521  -2.131  1.00  0.00           H
ATOM      7  HB2 ALA A   2       6.630   4.719  -1.206  1.00  0.00           H
ATOM      8  HB3 ALA A   2       5.809   3.141  -1.241  1.00  0.00           H
ATOM      9  C   ALA A   2       4.713   6.129   0.000  1.00  0.00           C
ATOM     10  O   ALA A   2       3.601   6.653   0.000  1.00  0.00           O
ATOM     11  N   GLY A   3       5.846   6.835   0.000  1.00  0.00           N
ATOM     12  H   GLY A   3       6.737   6.359   0.000  1.00  0.00           H
ATOM     13  CA  GLY A   3       5.846   8.284   0.000  1.00  0.00           C
ATOM     14  HA2 GLY A   3       5.332   8.648  -0.890  1.00  0.00           H
ATOM     15  HA3 GLY A   3       5.332   8.648   0.890  1.00  0.00           H
ATOM     16  C   GLY A   3       7.273   8.814   0.000  1.00  0.00           C
ATOM     17  O   GLY A   3       8.226   8.039   0.000  1.00  0.00           O
ATOM     18  OXT GLY A   3       7.390  10.040   0.000  1.00  0.00           O
TER
END
"""


def _input(tmp_path: Path) -> Path:
    path = tmp_path / "peptide.pdb"
    path.write_text(_UNCAPPED)
    return path


def test_cap_builder_adds_both_caps_and_is_idempotent(tmp_path: Path) -> None:
    from dvbfixer.terminal_caps import add_terminal_caps_to_pdb

    capped = add_terminal_caps_to_pdb(_input(tmp_path))
    text = capped.read_text()
    assert text.count(" ACE ") == 3
    assert text.count(" NME ") == 2
    assert " OXT GLY " not in text
    assert "CONECT" in text
    assert add_terminal_caps_to_pdb(capped) == capped


def test_missing_cap_chain_warns_and_continues(tmp_path: Path, capsys) -> None:
    from dvbfixer.terminal_caps import add_terminal_caps_to_pdb

    output = add_terminal_caps_to_pdb(_input(tmp_path), chain_ids=["Z"])
    assert output.exists()
    assert " ACE " not in output.read_text()
    assert "not found among protein chains: Z" in capsys.readouterr().out


@pytest.mark.parametrize(
    "xmls",
    [
        ["amber19/protein.ff19SB.xml", "amber19/tip3p.xml"],
        ["charmm36.xml", "charmm36/water.xml"],
    ],
    ids=["amber", "charmm"],
)
def test_generated_caps_match_openmm_forcefields(tmp_path: Path, xmls: list[str]) -> None:
    from openmm.app import Modeller, PDBFile

    from dvbfixer.ffutils import create_forcefield_with_openff
    from dvbfixer.terminal_caps import add_terminal_caps_to_pdb

    pdb = PDBFile(str(add_terminal_caps_to_pdb(_input(tmp_path))))
    modeller = Modeller(pdb.topology, pdb.positions)
    ff = create_forcefield_with_openff(xmls, modeller.topology)
    modeller.addHydrogens(ff)
    system = ff.createSystem(modeller.topology)
    assert system.getNumParticles() == len(list(modeller.topology.atoms()))
    assert [r.name for r in modeller.topology.residues()] == ["ACE", "ALA", "GLY", "NME"]


@pytest.mark.parametrize("ff", ["amber", "charmm"])
def test_top_keeps_capped_residues(tmp_path: Path, ff: str) -> None:
    from openmm.app import Modeller, PDBFile

    from dvbfixer.ffutils import create_forcefield_with_openff
    from dvbfixer.terminal_caps import add_terminal_caps_to_pdb
    from dvbfixer.top import main as top_main

    capped = add_terminal_caps_to_pdb(_input(tmp_path))
    pdb = PDBFile(str(capped))
    modeller = Modeller(pdb.topology, pdb.positions)
    xmls = (["charmm36.xml", "charmm36/water.xml"] if ff == "charmm"
            else ["amber19/protein.ff19SB.xml", "amber19/tip3p.xml"])
    forcefield = create_forcefield_with_openff(xmls, modeller.topology)
    modeller.addHydrogens(forcefield)
    hydrated = tmp_path / f"{ff}.pdb"
    with hydrated.open("w") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)

    out_top = tmp_path / f"{ff}.top"
    out_pdb = tmp_path / f"{ff}-top.pdb"
    top_main([str(hydrated), "--ff", ff, "-o", str(out_top), "--pdb", str(out_pdb)])
    residue_names = {line[17:20].strip() for line in out_pdb.read_text().splitlines()
                     if line.startswith("ATOM")}
    assert {"ACE", "ALA", "GLY", "NME"} <= residue_names
    assert len([line for line in out_pdb.read_text().splitlines()
                if line.startswith("ATOM")]) == 29


@pytest.mark.parametrize(
    ("ff", "ace_atoms", "nme_atoms"),
    [
        ("amber", {"CH3", "C", "O", "HH31", "HH32", "HH33"},
         {"N", "H", "CH3", "HH31", "HH32", "HH33"}),
        ("charmm", {"CH3", "C", "O", "HH31", "HH32", "HH33"},
         {"N", "HN", "CH3", "HH31", "HH32", "HH33"}),
    ],
)
def test_zbs_strip_heterogens_preserves_gromacs_named_caps(
    tmp_path: Path, ff: str, ace_atoms: set[str], nme_atoms: set[str]
) -> None:
    """Regression for the real 8B01 command: PDBFixer's stock
    removeHeterogens() must not delete ACE/NME immediately after generation."""
    from dvbfixer.zbs import main as zbs_main

    source = _input(tmp_path)
    output = tmp_path / f"zbs-capped-{ff}.pdb"
    zbs_main([
        str(source), "-o", str(output),
        "--skip-renumber", "--skip-model",
        "--cap-termini", "--strip-heterogens", "--no-solvent",
        "--no-propka", "--no-protassign", "--no-postflight",
        "--no-align-to-input", "--max-iter", "1", "--platform", "Reference",
        "--atom-naming", "gromacs", "--ff", ff,
    ])
    atoms_by_residue: dict[str, set[str]] = {}
    for line in output.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms_by_residue.setdefault(line[17:20].strip(), set()).add(
                line[12:16].strip()
            )
    assert atoms_by_residue["ACE"] == ace_atoms
    assert atoms_by_residue["NME"] == nme_atoms
