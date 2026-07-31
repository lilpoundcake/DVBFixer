"""Tests for `dvbfixer.model.renumber` — residue-number mapping.

Regression for a resSeq-collision bug in `_interpolate_gaps`'s internal
gap "not enough room" branch: when an internal gap's two flanking
residues are numerically adjacent in the INPUT's own numbering (the
depositor never reserved resSeqs for a genuinely-absent loop/linker —
e.g. an scFv construct where the disordered (GGGS)x4 linker between VL
and VH has no electron density, so the author numbered VH's first
residue immediately after VL's last one), the old code numbered the
gap BACKWARD from its right flank (`right - gap_len + k`), colliding
with resSeqs already assigned to `left` and everything before it.
Confirmed on a real structure (`test/8cz8/8cz8_a_u.pdb`, chain C, a
scFv: VL 1-111, missing 16-residue linker, VH 128-230 in the true
numbering) — the old code produced a chain with 16 duplicate
`(chain, resid)` keys and every residue after the linker off by a
constant -16.

Fixed by placing the gap sequentially from `left + 1` (matching the
"enough room" branch) and shifting every already-placed resSeq from
the gap onward forward by the deficit, mirroring the N-terminal
branch's existing shift-the-rest-of-the-chain pattern. The identical
bug existed a second time in `build_resnum_mapping`'s own hand-rolled
fallback gap-loop (align2d-mask path) — replaced with a direct call to
the fixed `_interpolate_gaps` instead of carrying two copies of the
same logic.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from dvbfixer.model.renumber import _interpolate_gaps, build_resnum_mapping

_SCFV_C_FASTA = (
    "DIALTQPASVSGSPGQSITISCTGTSSDIGGYNSVSWYQQHPGKAPKLMIYGVNNRPSGVSNRFSGSKSG"
    "NTASLTISGLQAEDEADYYCSSYDIESATPVFGGGTKLTVLGGGSGGGSGGGSGGGSQVELVQSGAEVK"
    "KPGESLKISCKGSGYSFTSYWIGWVRQAPGKGLEWMGIIDPGDSRTRYSPSFQGQVTISADKSISTAYL"
    "QWSSLKASDTAMYYCARGQLYGGTYMDGWGQGTLVTVSS"
)


def test_build_resnum_mapping_scfv_linker_gap(eightcz8_dir: Path) -> None:
    """The real scFv-linker-gap case: no duplicate resids, linker fills
    112-127, VH resumes at 128 (not the old, colliding 96-111/112...)."""
    src = eightcz8_dir / "8cz8_a_u.pdb"
    if not src.exists():
        pytest.skip(f"fixture missing: {src}")

    lines = src.read_text().splitlines(keepends=True)
    per_chain_masks = [[True] * len(_SCFV_C_FASTA)]
    mapping = build_resnum_mapping(
        per_chain_masks, ["C"], ["C"], lines,
        protein_seq_map={"C": _SCFV_C_FASTA},
    )

    resids = [mapping[("C", n)] for n in range(1, len(_SCFV_C_FASTA) + 1)]
    resseqs = [r[0] for r in resids]

    assert len(resseqs) == len(set(resseqs)), (
        f"duplicate resSeq(s) in mapping: {sorted({r for r in resseqs if resseqs.count(r) > 1})}"
    )
    assert resseqs[0] == 1
    assert resseqs[110] == 111  # last VL residue, unchanged
    assert resseqs[111] == 112  # first linker residue — the old bug put this at 96
    assert resseqs[126] == 127  # last linker residue
    assert resseqs[127] == 128  # first VH residue — the old bug put this at 112
    assert resseqs[-1] == len(_SCFV_C_FASTA)  # final (modeled) C-terminal residue


def test_interpolate_gaps_internal_gap_insufficient_room() -> None:
    """Synthetic, PDB-independent case for the exact collision branch.

    5 known positions (resSeq 1-5), a 3-residue gap with only 0 resSeqs
    of room before the next known residue (resSeq 6), then 2 more known
    positions continuing from resSeq 7-8 in the INPUT's numbering (i.e.
    the author never reserved room for the 3-residue gap).
    """
    full_resids = [
        (1, ' '), (2, ' '), (3, ' '), (4, ' '), (5, ' '),
        None, None, None,
        (6, ' '), (7, ' '),
    ]
    _interpolate_gaps(full_resids, len(full_resids))

    resseqs = [r[0] for r in full_resids]
    assert resseqs == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], resseqs
    assert len(resseqs) == len(set(resseqs))


def test_interpolate_gaps_internal_gap_enough_room_unchanged() -> None:
    """The "enough room" branch must still behave exactly as before —
    this fix only changes the collision branch."""
    full_resids = [(1, ' '), None, None, (10, ' ')]
    _interpolate_gaps(full_resids, len(full_resids))
    assert [r[0] for r in full_resids] == [1, 2, 3, 10]
