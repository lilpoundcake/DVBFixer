"""Regression test for the 2026-07-31 audit's icode-collapse fix.

Before the fix, `prepare.pipeline.run_pdbfixer` collapsed the
`(chain, resid, icode)` variant map down to `(chain, resid)` before
writing the final PDB / `.dat` file — so two residues sharing a resSeq
via insertion code (e.g. a Kabat CDR-loop `H:82`/`H:82A`) silently
overwrote each other's AMBER protonation-variant name. This is exactly
the shape of input this project's primary use case (antibody PDBs)
produces.

The variant names here come from the raw-text capture path (an input
residue already named HIE/HID/HIP/... gets its variant captured
directly, independent of PROPKA/Reduce), so this test is fast and
deterministic — no PROPKA/Reduce subprocess needed.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("openmm", reason="prepare integration needs OpenMM")
pytest.importorskip("pdbfixer", reason="prepare integration needs PDBFixer")


def _atom_line(serial: int, name: str, resname: str, chain: str,
                resseq: int, icode: str, x: float, y: float, z: float) -> str:
    name_field = f" {name:<3s}" if len(name) < 4 else f"{name:<4s}"
    return (
        f"ATOM  {serial:5d} {name_field} {resname:>3s} {chain}{resseq:4d}{icode:1s}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]:>1s}\n"
    )


# Two HIS-family variant residues sharing chain H / resSeq 354, differing
# only by insertion code — HIE at icode "" (no icode), HIP at icode "A".
# Heavy-atom coordinates are a real HIS ring (test/8cz8/8cz8_t_u.pdb chain
# E resid 354), duplicated with the second residue's copy shifted +10 A
# in x so PDBFixer's coincident-atom stripping doesn't touch either.
_HIS_RING_ATOMS = [
    ("N", 31.022, 37.898, -93.853),
    ("CA", 31.183, 38.191, -95.301),
    ("C", 30.098, 37.460, -96.111),
    ("O", 29.647, 38.025, -97.128),
    ("CB", 32.597, 37.833, -95.766),
    ("CG", 32.832, 38.111, -97.210),
    ("ND1", 32.320, 39.235, -97.838),
    ("CD2", 33.507, 37.422, -98.157),
    ("CE1", 32.673, 39.224, -99.106),
    ("NE2", 33.402, 38.124, -99.327),
]


def _build_icode_sibling_fragment() -> str:
    lines = []
    serial = 1
    for name, x, y, z in _HIS_RING_ATOMS:
        lines.append(_atom_line(serial, name, "HIE", "H", 354, " ", x, y, z))
        serial += 1
    for name, x, y, z in _HIS_RING_ATOMS:
        lines.append(_atom_line(serial, name, "HIP", "H", 354, "A", x + 10.0, y, z))
        serial += 1
    lines.append("TER\n")
    lines.append("END\n")
    return "".join(lines)


@pytest.mark.slow
def test_prepare_preserves_distinct_icode_sibling_variants(tmp_workdir: Path) -> None:
    """H:354 (HIE) and H:354A (HIP) must each keep their OWN variant name
    in the output — not collapse to whichever one iteration happened to
    process last."""
    from dvbfixer.prepare.pipeline import main as prepare_main

    input_pdb = tmp_workdir / "icode_siblings.pdb"
    input_pdb.write_text(_build_icode_sibling_fragment())
    out = tmp_workdir / "prep.pdb"

    prepare_main([str(input_pdb), "-o", str(out), "--no-propka", "--no-protassign"])

    assert out.exists()
    text = out.read_text()

    def _resname_at(icode: str) -> str | None:
        for ln in text.splitlines():
            if not ln.startswith(("ATOM  ", "HETATM")):
                continue
            if ln[21] != "H" or ln[22:26].strip() != "354":
                continue
            if ln[26] == (icode or " "):
                return ln[17:20].strip()
        return None

    resname_plain = _resname_at(" ")
    resname_a = _resname_at("A")
    assert resname_plain == "HIE", (
        f"H:354 (no icode) should still be HIE, got {resname_plain!r}"
    )
    assert resname_a == "HIP", (
        f"H:354A should still be HIP, got {resname_a!r} — icode collapse "
        f"would make this collide with H:354's HIE instead"
    )

    dat_path = out.with_suffix(".dat")
    if dat_path.exists():
        import json
        dat = json.loads(dat_path.read_text())
        vo = dat.get("variant_overrides") or {}
        assert vo.get("H:354:") == "HIE"
        assert vo.get("H:354:A") == "HIP"
