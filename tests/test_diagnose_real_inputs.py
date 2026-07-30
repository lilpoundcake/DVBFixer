"""Tests for `dvbfixer diagnose` on real PDB inputs.

diagnose is a report-only tool — it should never crash on any
input, and its severity distribution should look sane on the
curated shit/ set (which by construction has known issues)."""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def _run_diagnose(input_pdb: Path) -> tuple[int, str]:
    """Run diagnose; return (exit_code, stdout+stderr)."""
    proc = subprocess.run(
        ["dvbfixer", "diagnose", str(input_pdb)],
        capture_output=True, text=True,
    )
    return proc.returncode, proc.stdout + proc.stderr


def test_diagnose_pure_protein_completes(pure_protein_small: Path) -> None:
    """Smallest protein — diagnose should run cleanly."""
    code, output = _run_diagnose(pure_protein_small)
    assert code in (0, 1), (
        f"diagnose returned {code}; expected 0 (clean) or 1 (issues found)"
    )
    # Should produce some report text.
    assert output.strip(), "diagnose produced no output"


def test_diagnose_ss_bonded_antibody(ss_bonded_antibody: Path) -> None:
    """1DQJ has 12 SS bonds — diagnose should surface them without crash."""
    code, output = _run_diagnose(ss_bonded_antibody)
    assert code in (0, 1)
    # Common diagnose sections should appear.
    assert any(kw in output.lower()
               for kw in ("residue", "atom", "check", "bond")), (
        f"diagnose report unrecognisable:\n{output[:500]}"
    )


@pytest.mark.slow
def test_diagnose_glycoprotein_completes(glycoprot_stress: Path) -> None:
    """Trastuzumab — diagnose should not crash on glycan chains."""
    code, output = _run_diagnose(glycoprot_stress)
    assert code in (0, 1)
    assert output.strip()
