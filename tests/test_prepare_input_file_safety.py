"""Regression test for a data-loss bug found while writing test coverage
for the 2026-07-31 audit's icode-collapse fix.

`prepare.pipeline.run_pdbfixer` compared `preprocessed_path`/`canon_path`
(always plain `str`, per `_preprocess_glycoprotein_input`/
`_canonicalize_conect_records`'s no-rewrite return paths) against
`input_path` (a `pathlib.Path` object, as passed by `main()`) using `!=`.
`some_str != some_path_object` is ALWAYS `True` in Python regardless of
whether they name the same file — `Path.__eq__` requires matching types.

Confirmed on `main` (pre-existing, not introduced by this branch): a
completely clean PDB (no glycoprotein HETATM/TER fixups needed, no
malformed CONECT records) run through `dvbfixer prepare
--no-infer-conect` unconditionally executed
`Path(canon_path).unlink(missing_ok=True)` in the "no rewrite happened"
case — and since no rewrite happened, `canon_path` names the EXACT SAME
FILE as `input_path`, silently deleting the user's own input file from
disk. `--infer-conect` (the default) usually protects against this by
routing through an disposable intermediate CONECT-materialised temp
copy first, but `--no-infer-conect` on an already-clean input has no
such protection.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("openmm", reason="prepare integration needs OpenMM")
pytest.importorskip("pdbfixer", reason="prepare integration needs PDBFixer")


@pytest.mark.slow
def test_prepare_no_infer_conect_does_not_delete_input_file(tmp_workdir: Path) -> None:
    """`--no-infer-conect` on a clean, complete PDB (no HETATM/TER/CONECT
    fixups needed) must never delete the user's own input file."""
    from dvbfixer.prepare.pipeline import main as prepare_main

    # A minimal, fully-standard single-residue fragment — deliberately
    # "clean" so both _preprocess_glycoprotein_input and
    # _canonicalize_conect_records are no-ops (the exact condition that
    # triggered the bug).
    pdb_text = (
        "ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1      12.560  13.207   2.100  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1      13.090  14.620   2.100  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1      12.350  15.610   2.100  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1      13.090  12.430   3.300  1.00  0.00           C\n"
        "TER\n"
        "END\n"
    )
    input_pdb = tmp_workdir / "clean_input.pdb"
    input_pdb.write_text(pdb_text)
    original_size = input_pdb.stat().st_size
    out = tmp_workdir / "prep.pdb"

    prepare_main([
        str(input_pdb), "-o", str(out), "--no-infer-conect",
        "--no-propka", "--no-protassign",
    ])

    assert input_pdb.exists(), (
        "prepare deleted the user's own input file — the exact data-loss "
        "bug this test guards against"
    )
    assert input_pdb.stat().st_size == original_size, (
        "input file exists but was truncated/overwritten"
    )
    assert out.exists()
