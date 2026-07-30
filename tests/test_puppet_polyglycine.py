"""Fast unit tests for `dvbfixer puppet` — strip PDB to polyglycine
backbone."""
from __future__ import annotations

import subprocess
from pathlib import Path


def test_puppet_produces_polyglycine_backbone(
    pure_protein_small: Path, tmp_workdir: Path,
) -> None:
    """puppet should reduce every residue to N/CA/C/O + GLY name."""
    output = tmp_workdir / "puppet.pdb"
    proc = subprocess.run(
        ["dvbfixer", "puppet", str(pure_protein_small), "-o", str(output)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"puppet failed: stdout={proc.stdout} stderr={proc.stderr}"
    )
    assert output.exists()

    # Read atoms in output; check every residue is GLY and has only
    # backbone atoms (N, CA, C, O).
    from openmm.app import PDBFile
    pdb = PDBFile(str(output))
    # OXT is preserved on C-terminal residues as part of the backbone
    # (needed for capped-end topology recognition downstream).
    _BB = {"N", "CA", "C", "O", "OXT"}
    non_gly = []
    non_bb = []
    for res in pdb.topology.residues():
        if res.name != "GLY":
            non_gly.append((res.chain.id, res.id, res.name))
        for atom in res.atoms():
            if atom.element.symbol == "H":
                continue  # H atoms may or may not be present
            if atom.name not in _BB:
                non_bb.append((res.chain.id, res.id, atom.name))

    assert not non_gly, f"non-GLY residues in puppet output: {non_gly[:5]}"
    assert not non_bb, (
        f"non-backbone heavy atoms in puppet output: {non_bb[:5]}"
    )
