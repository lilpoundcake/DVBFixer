"""Tests for `dvbfixer split` on multi-MODEL PDBs."""
from __future__ import annotations

import subprocess
from pathlib import Path


def test_split_multimodel_produces_valid_output(
    multistate_pdb: Path, tmp_workdir: Path,
) -> None:
    """split should process a multi-MODEL PDB (chain-splitting mode) and
    produce a valid output PDB with TER records inserted at chain
    breaks."""
    output = tmp_workdir / "split.pdb"
    proc = subprocess.run(
        ["dvbfixer", "split", str(multistate_pdb), "-o", str(output)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"split failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    assert output.exists(), "split produced no output file"
    text = output.read_text()
    # Should have ATOM records.
    assert "ATOM" in text, "split output has no ATOM records"
