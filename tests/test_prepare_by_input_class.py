"""Tests for `dvbfixer prepare` across input classes.

Every prepare output must:
- Complete without exception on every input class.
- Emit a valid .dat file (DatRecord.load succeeds).
- Have zero D-Cα (chirality invariant — enforced downstream in
  minimize, but prepare shouldn't produce D-Cα either).

The tleap-reduce opt-in backend is verified separately on pure-
protein inputs (it rejects non-canonical residues by design).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
from tests.conftest import count_d_ca_residues


def _run_prepare(input_pdb: Path, workdir: Path, *extra: str) -> Path:
    import shutil
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
def test_prepare_pure_protein_default_backend(
    pure_protein_small: Path, tmp_workdir: Path,
) -> None:
    """Default backend (legacy in 0.7.5+) on the smallest pure protein."""
    out = _run_prepare(pure_protein_small, tmp_workdir)
    assert count_d_ca_residues(out) == 0
    # .dat should exist next to output.
    dat = out.with_suffix(".dat")
    assert dat.exists(), f".dat missing: {dat}"
    # Load via DatRecord to verify schema.
    from dvbfixer.ffutils.dat import DatRecord
    rec = DatRecord.load(dat)
    assert rec is not None


@pytest.mark.slow
def test_prepare_pure_protein_tleap_reduce_backend(
    pure_protein_small: Path, tmp_workdir: Path,
) -> None:
    """Opt-in tleap-reduce backend should still work on pure protein."""
    out = _run_prepare(pure_protein_small, tmp_workdir,
                        "--backend", "tleap-reduce")
    assert count_d_ca_residues(out) == 0


@pytest.mark.slow
def test_prepare_ss_bonded_antibody(
    ss_bonded_antibody: Path, tmp_workdir: Path,
) -> None:
    """Antibody with 12 SS bonds — legacy prep should preserve CYX."""
    out = _run_prepare(ss_bonded_antibody, tmp_workdir)
    assert count_d_ca_residues(out) == 0
    # Check for CYX residues (or CYS if legacy prep uses SSBOND
    # records instead of residue renaming).
    from openmm.app import PDBFile
    pdb = PDBFile(str(out))
    cys_family = sum(1 for r in pdb.topology.residues()
                     if r.name in ("CYX", "CYS"))
    assert cys_family >= 12, (
        f"Expected ≥12 CYS/CYX residues, got {cys_family}"
    )
