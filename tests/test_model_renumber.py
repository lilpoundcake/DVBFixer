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


def _synthetic_pdb_lines(n_protein: int, ligand_resseq: int) -> list[str]:
    """Chain A: `n_protein` ALA residues (resSeq 1..n_protein, one CA atom
    each) + one HETATM ligand ("LIG") at `ligand_resseq` with atoms named
    CA/CB — deliberately colliding with real protein backbone atom names,
    mirroring test/lipid/7x35_r_u.pdb's PLM (palmitic acid uses IUPAC-ish
    names C1..C9,CA,CB,CC... for its alkyl chain)."""
    lines = []
    serial = 1
    for resseq in range(1, n_protein + 1):
        lines.append(
            f"ATOM  {serial:5d}  CA  ALA A{resseq:4d}    "
            f"{float(resseq):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
        )
        serial += 1
    # Same column layout as the ATOM lines above (record field is 6 chars
    # either way — "ATOM  " / "HETATM" — everything after lines up
    # identically): 2-char atom name, 3-char resname, chain, 4-digit resseq.
    for name in ("CA", "CB"):
        lines.append(
            f"HETATM{serial:5d}  {name:<2s}  LIG A{ligand_resseq:4d}    "
            f"{100.0:8.3f}{100.0:8.3f}{100.0:8.3f}  1.00  0.00           C\n"
        )
        serial += 1
    return lines


def test_build_resnum_mapping_hetatm_collision_with_gap_fill() -> None:
    """Regression for the real test/lipid/7x35_r_u.pdb bug: a HETATM
    ligand's ORIGINAL resSeq (as assigned by the earlier, FASTA-blind
    standalone `renumber.py` step — here, deliberately one past the
    ATOM-only residue count) collides with a protein resSeq that only
    exists once the true, FASTA-complete sequence's C-terminal gap gets
    filled in. The ligand must be renumbered out of the way, not
    silently left at a resSeq a real (gap-filled) protein residue now
    also occupies.
    """
    n_protein = 10
    full_seq_len = 13  # 3 more residues than the ATOM records have — a gap
    ligand_resseq = 11  # what upstream renumber.py naively assigned: n_protein + 1

    lines = _synthetic_pdb_lines(n_protein, ligand_resseq)
    # mask length = protein positions + 1 trailing HETATM slot (Modeller's
    # own target alignment includes non-protein residues as template
    # positions too — see CLAUDE.md's "model.py non-protein" note).
    per_chain_masks = [[True] * (full_seq_len + 1)]
    mapping = build_resnum_mapping(
        per_chain_masks, ["A"], ["A"], lines,
        protein_seq_map={"A": "A" * full_seq_len},
    )

    resids = [mapping[("A", n)] for n in range(1, full_seq_len + 2)]
    resseqs = [r[0] for r in resids]

    assert resseqs[:n_protein] == list(range(1, n_protein + 1))
    # Gap-filled protein residues 11, 12, 13 must NOT be stolen by the ligand.
    assert resseqs[n_protein:full_seq_len] == [11, 12, 13]
    # The ligand (last mapping entry, placed at the trailing HETATM slot)
    # must have been moved off resSeq 11 — it must not equal ANY protein
    # resSeq now in use.
    ligand_final_resseq = resseqs[-1]
    assert ligand_final_resseq not in set(resseqs[:full_seq_len]), (
        f"ligand resSeq {ligand_final_resseq} collides with a protein "
        f"residue: {resseqs[:full_seq_len]}"
    )
    assert len(resseqs) == len(set(resseqs)), (
        f"duplicate resSeq(s) in mapping: {resseqs}"
    )
