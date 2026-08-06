"""Baseline end-to-end coverage for `dvbfixer top` / `TopologyBuilder`.

`TopologyBuilder` (top/pipeline.py, ~1175 lines) had ZERO test coverage
before this file — `grep -rln "TopologyBuilder\\|build_glycan_chain\\|
build_glycolipid_chain\\|build_chain\\b" tests/` returned nothing. This
is deliberately added BEFORE moving the class to its own module
(planned follow-up), so the move can be verified against a real
before/after behavioural baseline instead of "trust the diff."

Covers the three real `TopologyBuilder` entry points exercised by
`main()`:
  - `build_chain` — plain protein (`tests/fixtures/8cz8/8cz8_a_u.pdb`).
  - `build_glycan_chain` — the tracked CHARMM glycosylated antibody.
  - The "unrecognized HETATM chain" path (`tests/fixtures/lipid/7x35_r_u.pdb`, PLM)
    — this ligand does NOT reach `build_chain`/`build_glycolipid_chain`
    at all today; `top/pipeline.py`'s chain classifier drops any chain
    with no FF/GLYCAM/ceramide-recognized residue before topology
    building starts, now with an explicit WARNING (previously silent —
    found and fixed alongside this test) rather than a lipid-specific
    builder path. This test locks in that (documented, user-approved)
    current behaviour: protein chains still build correctly, the
    unrecognized chain is dropped with a WARNING, nothing crashes.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("openmm", reason="top's CONECT inference needs OpenMM-adjacent deps")


def _run_top(input_pdb: Path, out_dir: Path, ff: str = "amber") -> tuple[str, str]:
    from dvbfixer.top import main as top_main
    top_out = out_dir / "topol.top"
    pdb_out = out_dir / "topol.pdb"
    top_main([str(input_pdb), "-o", str(top_out), "--pdb", str(pdb_out), "--ff", ff])
    return top_out.read_text(), pdb_out.read_text()


def test_top_pure_protein_baseline(pure_protein_small: Path, tmp_workdir: Path) -> None:
    """`build_chain` on a plain protein: sane atom count, no glycan/lipid
    noise, molecules section lists at least one Protein_chain."""
    top_text, pdb_text = _run_top(pure_protein_small, tmp_workdir)

    assert "[ molecules ]" in top_text
    mol_section = top_text.split("[ molecules ]")[1]
    assert "Protein_chain" in mol_section

    atom_lines = [ln for ln in pdb_text.splitlines() if ln.startswith(("ATOM  ", "HETATM"))]
    assert len(atom_lines) > 100, "suspiciously few atoms in topology-matched PDB output"


def test_top_glycoprotein_baseline(charmm_glycan_pdb: Path, tmp_workdir: Path) -> None:
    """`build_glycan_chain` on a properly-bonded, pre-processed
    CHARMM-GUI glycosylated antibody (`conf.pdb` — CHARMM-GUI's own
    output already carries the glycosidic bonds `detect_glycan_links`
    needs; a raw, un-prepped RCSB-style PDB does not, and that's a
    separate, pre-existing input-readiness gap, not something this
    baseline test is meant to cover): at least one Glycan_
    moleculetype gets built alongside the protein chains."""
    top_text, _ = _run_top(charmm_glycan_pdb, tmp_workdir, ff="charmm")

    assert "[ molecules ]" in top_text
    mol_section = top_text.split("[ molecules ]")[1]
    assert "Protein_chain" in mol_section
    assert "Glycan" in top_text, (
        "expected at least one Glycan_* moleculetype for the tracked "
        "glycosylated antibody fixture"
    )


def test_top_unrecognized_hetatm_chain_warns_and_drops(
    lipid_dir: Path, tmp_workdir: Path, capsys: pytest.CaptureFixture[str],
) -> None:
    """The lipid fixture's PLM ligand has no FF/GLYCAM/ceramide match.
    Today `top` drops the whole chain rather than building a topology
    for it — this must be a loud WARNING (not silence), and protein
    chains must still build correctly around it."""
    src = lipid_dir / "7x35_r_u.pdb"
    assert src.is_file(), f"tracked fixture missing: {src}"

    top_text, _ = _run_top(src, tmp_workdir)
    captured = capsys.readouterr()

    assert "PLM" in captured.err or "PLM" in captured.out, (
        "expected the unrecognized-chain WARNING to name PLM"
    )
    assert "WARNING" in captured.err or "WARNING" in captured.out

    assert "[ molecules ]" in top_text
    mol_section = top_text.split("[ molecules ]")[1]
    assert "Protein_chain" in mol_section, (
        "protein chains must still build a topology even though the "
        "unrecognized PLM chain was dropped"
    )
    assert "PLM" not in top_text, (
        "PLM has no FF entry — it must not appear as a moleculetype "
        "name in the generated topology"
    )
