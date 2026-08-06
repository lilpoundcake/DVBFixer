"""Ligand H/valence correctness for `dvbfixer prepare`'s heterogen
H-addition (`prepare/glycan.py`).

PDB files carry no bond-order information, so RDKit's/OpenBabel's
proximity-based bond perception + naive valence-filling silently
over-protonates resonance-delocalized ionizable groups (carboxylate,
sulfonate/sulfate, phosphate) and misses genuine alkenes whose bond
length looks single-bond-like at typical crystallographic resolution.
Regression coverage for both classes of bug, surfaced by
`tests/fixtures/1VCU.pdb` (DAN + EPE).
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from tests.conftest import count_h_on_atom

# RCSB Chemical Component Dictionary connectivity/stereochemistry with the
# physiological ionic states used by prepare: DAN carboxylate and zwitterionic
# HEPES (EPE; hydroxyethyl-substituted piperazine N protonated).
DAN_SMILES = "CC(=O)N[C@@H]1[C@H](C=C(O[C@H]1[C@@H]([C@@H](CO)O)O)C(=O)[O-])O"
EPE_SMILES = "C1CN(CC[NH+]1CCO)CCS(=O)(=O)[O-]"


def _run_prepare(input_pdb: Path, workdir: Path, *extra: str) -> Path:
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / "input.pdb"
    shutil.copy(input_pdb, local)
    output = workdir / "prep.pdb"
    proc = subprocess.run(
        ["dvbfixer", "prepare", str(local), "-o", str(output), *extra],
        capture_output=True, text=True, cwd=str(workdir),
    )
    assert proc.returncode == 0, (
        f"prepare failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    assert output.exists()
    return output


@pytest.mark.slow
def test_prepare_dan_alkene_and_ionizable_groups(
    protein_ligand_gaps_pdb: Path, tmp_workdir: Path,
) -> None:
    """DAN (2,3-didehydro sialic acid analog, DANA) has a real ring
    alkene at C2=C3 that pure-distance bond perception can't recover —
    both RDKit and OpenBabel default to sp3 and over-saturate with H
    unless corrected. Its carboxylate (C1/O1A/O1B) must also come out
    fully ionized (0 H on both O), not the pre-fix inconsistent single
    HO1B."""
    out = _run_prepare(protein_ligand_gaps_pdb, tmp_workdir, "--ff", "amber")
    assert count_h_on_atom(out, "DAN", "C2") == 0, (
        "DAN C2 is sp2 (ring alkene C2=C3, 3 heavy substituents) — "
        "should carry zero H"
    )
    assert count_h_on_atom(out, "DAN", "C3") == 1, (
        "DAN C3 is sp2 (ring alkene C2=C3, 2 heavy substituents + 1 H) "
        "— should carry exactly one H, not two"
    )
    assert count_h_on_atom(out, "DAN", "O1A") == 0
    assert count_h_on_atom(out, "DAN", "O1B") == 0


@pytest.mark.slow
def test_prepare_epe_sulfonate_fully_ionized(
    protein_ligand_gaps_pdb: Path, tmp_workdir: Path,
) -> None:
    """EPE (HEPES buffer)'s sulfonate (-SO3⁻) is fully ionized at
    physiological pH (pKa ≈ -2) — at most one S-OH exists even in the
    neutral acid form, and zero in the ionized form. Every copy in the
    output must show zero H on all three sulfonate oxygens."""
    out = _run_prepare(protein_ligand_gaps_pdb, tmp_workdir, "--ff", "amber")
    for atom_name in ("O1S", "O2S", "O3S"):
        assert count_h_on_atom(out, "EPE", atom_name) == 0, (
            f"EPE {atom_name} (sulfonate) should never be protonated"
        )


@pytest.mark.slow
def test_prepare_1vcu_with_optional_smiles(
    protein_ligand_gaps_pdb: Path, tmp_workdir: Path,
) -> None:
    """The optional SMILES path maps DAN and both EPE instances while the
    existing no-SMILES tests above continue to cover automatic preparation."""
    out = _run_prepare(
        protein_ligand_gaps_pdb, tmp_workdir, "--ff", "amber",
        "--smiles", f"DAN={DAN_SMILES}",
        "--smiles", f"EPE={EPE_SMILES}",
    )

    assert count_h_on_atom(out, "DAN", "C2") == 0
    assert count_h_on_atom(out, "DAN", "C3") == 1
    assert count_h_on_atom(out, "DAN", "O1A") == 0
    assert count_h_on_atom(out, "DAN", "O1B") == 0
    for atom_name in ("O1S", "O2S", "O3S"):
        assert count_h_on_atom(out, "EPE", atom_name) == 0
    assert count_h_on_atom(out, "EPE", "N1") == 0
    assert count_h_on_atom(out, "EPE", "N4") == 2  # one H on each EPE copy

    record = json.loads(out.with_suffix(".dat").read_text())
    mapped_h = [
        atom for atom in record["added_atoms"]
        if atom["resname"] in {"DAN", "EPE"} and atom["element"] == "H"
    ]
    assert mapped_h, "SMILES-generated ligand H atoms must be recorded in .dat"


def test_ionizable_detector_generalizes_beyond_dan_epe() -> None:
    """The carboxylate/sulfonate/sulfate/phosphate detector is
    connectivity-based, not a per-ligand hardcoded list — verify it
    correctly flags a synthetic acetate ion's carboxylate oxygens
    (a completely different molecule from DAN/EPE) and does NOT
    false-positive on a plain alcohol or a simple amide carbonyl."""
    pytest.importorskip("rdkit", reason="needs RDKit")
    from rdkit import Chem

    from dvbfixer.ffutils.ligand_valence import (
        find_ionizable_terminal_oxygens_rdkit,
    )

    # Acetate ion: CH3-COO(-) — both carboxylate O's should be flagged.
    acetate = Chem.MolFromSmiles("CC(=O)[O-]")
    acetate = Chem.AddHs(acetate)
    ionizable = find_ionizable_terminal_oxygens_rdkit(acetate)
    o_indices = {a.GetIdx() for a in acetate.GetAtoms() if a.GetSymbol() == "O"}
    assert ionizable == o_indices, (
        "both carboxylate oxygens of a synthetic acetate ion should be "
        "detected as ionizable, proving the detector generalizes beyond "
        "the DAN/EPE fixture"
    )

    # Ethanol: CH3-CH2-OH — plain alcohol O has only ONE terminal O on
    # its carbon, must NOT be flagged.
    ethanol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert find_ionizable_terminal_oxygens_rdkit(ethanol) == set()

    # Acetamide: CH3-C(=O)-NH2 — carbonyl O has only one terminal O on
    # its carbon (the amide N isn't a terminal O), must NOT be flagged.
    acetamide = Chem.AddHs(Chem.MolFromSmiles("CC(=O)N"))
    assert find_ionizable_terminal_oxygens_rdkit(acetamide) == set()
