"""Trailing waters must not hide an uncapped protein OXT or create false caps."""
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest
from openmm.app import Topology, element
from pdbfixer import PDBFixer

from dvbfixer.diagnose.structural import check_missing_terminals


def fixture(*, chain_id="a", oxt=False, cap=False, cyclic=False):
    top = Topology()
    chain = top.addChain(chain_id)
    protein = top.addResidue("ASP", chain, "136", insertionCode="A")
    atoms = {name: top.addAtom(name, element.oxygen if name.startswith("O") else element.nitrogen if name == "N" else element.carbon, protein)
             for name in ("N", "CA", "C", "O", *(("OXT",) if oxt else ()))}
    if cap:
        other = top.addResidue("NME", chain, "137")
        nitrogen = top.addAtom("N", element.nitrogen, other)
        top.addBond(atoms["C"], nitrogen)
    water = top.addResidue("HOH", chain, "200")
    top.addAtom("O", element.oxygen, water)
    if cyclic:
        other = top.addResidue("ALA", top.addChain("B"), "1")
        nitrogen = top.addAtom("N", element.nitrogen, other)
        top.addBond(atoms["C"], nitrogen)
    return SimpleNamespace(topology=top, missingTerminals={}, missingResidues={}), protein


def test_trailing_water_reports_protein_OXT_without_mutating_fixer():
    fixer, residue = fixture()
    result = check_missing_terminals(fixer)
    assert len(result) == 1
    assert (result[0].chain, result[0].resid, result[0].category) == ("a", "136A", "missing_terminals")
    assert "OXT" in result[0].message and not fixer.missingTerminals
    fixer.missingTerminals = {residue: ["OXT"]}
    assert len(check_missing_terminals(fixer)) == 1


@pytest.mark.parametrize("case", ["present", "cap", "cyclic", "future-gap", "unresolved-seqres", "water-only"])
def test_water_fallback_does_not_invent_free_carboxylates(case):
    fixer, residue = fixture(oxt=case == "present", cap=case == "cap", cyclic=case == "cyclic")
    if case == "future-gap":
        fixer.missingResidues = {(residue.chain.index, 1): ["GLY"]}
    elif case == "unresolved-seqres":
        fixer.sequences = [SimpleNamespace(chainId="a", residues=["ASP", "GLY"])]
    elif case == "water-only":
        top = Topology()
        water = top.addResidue("HOH", top.addChain("a"), "1")
        top.addAtom("O", element.oxygen, water)
        fixer.topology = top
    assert not check_missing_terminals(fixer)


def test_actual_PDBFixer_water_chain_false_negative_is_supplemented():
    source = Path(__file__).parent / "fixtures/8cz8/8cz8_a_u.pdb"
    lines = [line for line in source.read_text().splitlines(keepends=True)
             if line.startswith("ATOM  ") and line[21] == "C" and 6 <= int(line[22:26]) <= 8
             and line[76:78].strip() != "H"]
    text = "".join(lines) + "HETATM  900  O   HOH C 101      50.000  50.000  50.000  1.00 20.00           O\nEND\n"
    fixer = PDBFixer(pdbfile=StringIO(text))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    original = dict(fixer.missingTerminals)
    findings = check_missing_terminals(fixer)
    assert any(f.chain == "C" and f.resid == "8" and "OXT" in f.message for f in findings)
    assert fixer.missingTerminals == original


def test_seqres_chain_case_is_not_collapsed_and_known_N_gap_does_not_hide_C_end():
    fixer, residue = fixture(chain_id="a")
    fixer.sequences = [SimpleNamespace(chainId="A", residues=["ASP", "GLY"])]
    assert len(check_missing_terminals(fixer)) == 1
    fixer.sequences = [SimpleNamespace(chainId="a", residues=["ASP", "GLY"])]
    fixer.missingResidues = {(residue.chain.index, 0): ["GLY"]}
    assert len(check_missing_terminals(fixer)) == 1
