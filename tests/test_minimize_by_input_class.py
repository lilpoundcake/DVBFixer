"""Tests for `dvbfixer minimize` across input classes.

Minimize is the guarantor of the chirality invariant (zero D-Cα in
output) via its post-phase-2 unconditional force-reflect. Every
minimize output must:
- Complete without exception.
- Have zero D-Cα in output.
- Preserve SS bonds when present in input.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from tests.conftest import count_d_ca_residues, count_ss_bonds


def _run_prepare_then_minimize(
    input_pdb: Path, workdir: Path, *minimize_extra: str,
) -> Path:
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / "input.pdb"
    shutil.copy(input_pdb, local)
    prep_out = workdir / "prep.pdb"
    p1 = subprocess.run(
        ["dvbfixer", "prepare", str(local), "-o", str(prep_out)],
        capture_output=True, text=True, cwd=str(workdir),
    )
    assert p1.returncode == 0, (
        f"prepare failed: {p1.stdout}\n{p1.stderr}"
    )
    mini_out = workdir / "mini.pdb"
    p2 = subprocess.run(
        ["dvbfixer", "minimize", str(prep_out), "-o", str(mini_out),
         "--no-solvent", *minimize_extra],
        capture_output=True, text=True, cwd=str(workdir),
    )
    assert p2.returncode == 0, (
        f"minimize failed: {p2.stdout}\n{p2.stderr}"
    )
    return mini_out


@pytest.mark.slow
def test_minimize_pure_protein_zero_d_ca(
    pure_protein_small: Path, tmp_workdir: Path,
) -> None:
    out = _run_prepare_then_minimize(pure_protein_small, tmp_workdir)
    assert count_d_ca_residues(out) == 0


@pytest.mark.slow
def test_minimize_ss_bonded_preserves_ss(
    ss_bonded_antibody: Path, tmp_workdir: Path,
) -> None:
    """The SG-SG exception in _drop_spurious_inter_aa_bonds must
    keep disulfides intact through minimize."""
    out = _run_prepare_then_minimize(ss_bonded_antibody, tmp_workdir)
    assert count_d_ca_residues(out) == 0
    assert count_ss_bonds(out) >= 6, (
        "Expected ≥6 SG-SG bonds preserved through minimize"
    )
