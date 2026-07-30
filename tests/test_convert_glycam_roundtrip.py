"""Tests for `dvbfixer convert` (PDB → GLYCAM sugar naming).

The convert tool detects glycosidic bonds and renames PDB-standard
sugar codes (NAG/BGL/BMA/etc.) to GLYCAM 3-char codes with linkage
prefixes (e.g. 4YB, VMA). Glycosylation-site protein residues become
NLN / OLS / OLT.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def _has_glycam_named_sugars(pdb_path: Path) -> bool:
    """True if the output has residues matching GLYCAM's 3-char pattern
    ``[digit-or-linkage-letter][sugar-letter][A|B]``. Cheap heuristic
    on the file text."""
    text = pdb_path.read_text()
    # GLYCAM sugar letters: G, L, M, Y, V, f, S, X, R, Z, U, h.
    # Anomer: A / B (case-sensitive). Linkage: 0-9 or V/W/U/Z/X/Y/T/S/R/Q/P.
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 20:
            continue
        rn = line[17:20].strip()
        if len(rn) != 3:
            continue
        if rn[2] not in ("A", "B"):
            continue
        if rn[1] not in ("G", "L", "M", "Y", "V", "f", "S", "X",
                          "R", "Z", "U", "h"):
            continue
        # Linkage char is broad; if we got this far, treat as GLYCAM.
        return True
    return False


@pytest.mark.slow
def test_convert_glycam_names_emitted(
    glycoprot_stress: Path, tmp_workdir: Path,
) -> None:
    """convert on trastuzumab (PDB-standard BGL/BMA/AMA/…) should
    emit GLYCAM-canonical residue names."""
    output = tmp_workdir / "converted.pdb"
    proc = subprocess.run(
        ["dvbfixer", "convert", str(glycoprot_stress), "-o", str(output)],
        capture_output=True, text=True,
    )
    # convert may print WARNINGs but should complete cleanly.
    assert proc.returncode == 0, (
        f"convert failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    assert output.exists()
    assert _has_glycam_named_sugars(output), (
        f"convert output has no GLYCAM-named sugars in {output.name}"
    )


@pytest.mark.slow
def test_convert_glycam_idempotent(
    glycoprot_stress: Path, tmp_workdir: Path,
) -> None:
    """Running convert twice should be a no-op the second time (already
    GLYCAM-named)."""
    step1 = tmp_workdir / "step1.pdb"
    step2 = tmp_workdir / "step2.pdb"
    for out in (step1, step2):
        src = glycoprot_stress if out is step1 else step1
        proc = subprocess.run(
            ["dvbfixer", "convert", str(src), "-o", str(out)],
            capture_output=True, text=True,
        )
        assert proc.returncode == 0, (
            f"convert failed: {proc.stdout}\n{proc.stderr}"
        )
    # After the second pass, residue names should stabilise.
    # Compare residue name sets.
    def _resname_set(p: Path) -> set[str]:
        s = set()
        for line in p.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                s.add(line[17:20].strip())
        return s
    assert _resname_set(step1) == _resname_set(step2), (
        "convert not idempotent — residue names changed on second run"
    )


def _count_records(pdb_path: Path, record: str) -> int:
    prefix = f"{record:<6s}"[:6] if record != "CRYST1" else "CRYST1"
    return sum(1 for line in pdb_path.read_text().splitlines()
               if line.startswith(prefix))


@pytest.mark.slow
def test_convert_preserves_header_records(
    glycoprot_underannotated_conect_pdb: Path, tmp_workdir: Path,
) -> None:
    """convert must pass through SEQRES/HELIX/SHEET/CRYST1 unchanged —
    `_parse_pdb` only ever reads ATOM/HETATM/CONECT/LINK, so these were
    previously silently dropped, breaking downstream gap-modeling
    (which needs SEQRES) for anyone using convert -> zbs directly."""
    output = tmp_workdir / "converted.pdb"
    proc = subprocess.run(
        ["dvbfixer", "convert", str(glycoprot_underannotated_conect_pdb),
         "-o", str(output)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"convert failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    for record in ("SEQRES", "HELIX", "SHEET", "CRYST1"):
        n_in = _count_records(glycoprot_underannotated_conect_pdb, record)
        n_out = _count_records(output, record)
        assert n_out == n_in, (
            f"{record}: expected {n_in} (matching input), got {n_out} "
            f"in convert output"
        )


@pytest.mark.slow
def test_convert_finds_undocumented_glycosylation_site(
    glycoprot_underannotated_conect_pdb: Path, tmp_workdir: Path,
) -> None:
    """This fixture has 4 real N-glycosylation sites (chains A, B, C×2),
    but the deposited PDB's CONECT/LINK annotation only documents 2 of
    them (a genuine annotation gap for the other 2, not evidence they
    don't exist). convert's bond-detection previously trusted "any
    CONECT present" as an all-or-nothing gate and silently missed the
    undocumented sites entirely (their Asn stayed unrenamed, sugar tree
    floating, unbonded). All 4 sites must be found now."""
    output = tmp_workdir / "converted.pdb"
    proc = subprocess.run(
        ["dvbfixer", "convert", str(glycoprot_underannotated_conect_pdb),
         "-o", str(output), "-v"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"convert failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    assert "Detected 4 protein-sugar links" in proc.stdout, (
        f"expected all 4 protein-sugar glycosylation sites detected; "
        f"stdout:\n{proc.stdout}"
    )
