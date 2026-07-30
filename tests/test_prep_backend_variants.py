"""Unit tests for PROPKA-driven variant assignment in prep_backend.

The tleap/reduce subprocess wrappers are covered by end-to-end fixtures
elsewhere; these tests focus on the pure-Python variant-decision path
that lives entirely in `dvbfixer.prep_backend`:

* HIS tautomer inference from Reduce's placed H atoms
* `assign_amber_variants` file rewrite + terminal skip
* `_patch_variant_hydrogens` H atom additions/deletions for each variant

They exercise the fix for the bug where PROPKA-driven ASH/GLH/LYN/HIP
were silently dropped: the propka_dict collision (chain, resnum) →
(restype, pka) was overwriting side-chain pKas with N+/C- pKas, and
_classify_variant never had a HIS branch.
"""
from __future__ import annotations

from pathlib import Path

import pytest

ASH_MODEL_PDB = """\
ATOM      1  N   ASP A   5      12.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  ASP A   5      13.000  10.000  10.000  1.00  0.00           C
ATOM      3  C   ASP A   5      14.000  10.000  10.000  1.00  0.00           C
ATOM      4  O   ASP A   5      15.000  10.000  10.000  1.00  0.00           O
ATOM      5  CB  ASP A   5      13.500  11.500  10.000  1.00  0.00           C
ATOM      6  CG  ASP A   5      14.500  12.500  10.000  1.00  0.00           C
ATOM      7  OD1 ASP A   5      14.500  13.700  10.000  1.00  0.00           O
ATOM      8  OD2 ASP A   5      15.700  12.100  10.000  1.00  0.00           O
"""


HIP_MODEL_PDB = """\
ATOM      1  N   HIS A  10      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  HIS A  10      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  C   HIS A  10      12.000  10.000  10.000  1.00  0.00           C
ATOM      4  O   HIS A  10      13.000  10.000  10.000  1.00  0.00           O
ATOM      5  CB  HIS A  10      11.500  11.500  10.000  1.00  0.00           C
ATOM      6  CG  HIS A  10      12.500  12.500  10.000  1.00  0.00           C
ATOM      7  ND1 HIS A  10      13.700  12.100  10.000  1.00  0.00           N
ATOM      8  CD2 HIS A  10      12.500  13.900  10.000  1.00  0.00           C
ATOM      9  CE1 HIS A  10      14.500  13.100  10.000  1.00  0.00           C
ATOM     10  NE2 HIS A  10      13.700  14.200  10.000  1.00  0.00           N
ATOM     11  HE2 HIS A  10      13.700  15.200  10.000  1.00  0.00           H
"""


LYN_MODEL_PDB = """\
ATOM      1  N   LYS A  20      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  LYS A  20      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  C   LYS A  20      12.000  10.000  10.000  1.00  0.00           C
ATOM      4  O   LYS A  20      13.000  10.000  10.000  1.00  0.00           O
ATOM      5  CB  LYS A  20      11.500  11.500  10.000  1.00  0.00           C
ATOM      6  CG  LYS A  20      12.500  12.500  10.000  1.00  0.00           C
ATOM      7  CD  LYS A  20      13.500  13.500  10.000  1.00  0.00           C
ATOM      8  CE  LYS A  20      14.500  14.500  10.000  1.00  0.00           C
ATOM      9  NZ  LYS A  20      15.500  15.500  10.000  1.00  0.00           N
ATOM     10  HZ1 LYS A  20      16.000  16.000  10.000  1.00  0.00           H
ATOM     11  HZ2 LYS A  20      16.000  15.000  10.000  1.00  0.00           H
ATOM     12  HZ3 LYS A  20      15.000  16.000  10.000  1.00  0.00           H
"""


# ---------------------------------------------------------------------------
# _infer_his_tautomers_from_atoms
# ---------------------------------------------------------------------------


def test_infer_his_tautomer_hie_from_he2(tmp_path: Path) -> None:
    """HE2 placed → HIE."""
    from dvbfixer.prep_backend import _infer_his_tautomers_from_atoms

    p = tmp_path / "hie.pdb"
    p.write_text(HIP_MODEL_PDB)  # already has HE2 only
    result = _infer_his_tautomers_from_atoms(p)
    assert result == {("A", 10, ""): "HIE"}


def test_infer_his_tautomer_hid_from_hd1(tmp_path: Path) -> None:
    """HD1 placed → HID."""
    from dvbfixer.prep_backend import _infer_his_tautomers_from_atoms

    # Replace HE2 line with an HD1 line to simulate Reduce's other choice.
    pdb = HIP_MODEL_PDB.replace(
        "ATOM     11  HE2 HIS A  10      13.700  15.200  10.000  1.00  0.00           H\n",
        "ATOM     11  HD1 HIS A  10      13.700  11.100  10.000  1.00  0.00           H\n",
    )
    p = tmp_path / "hid.pdb"
    p.write_text(pdb)
    result = _infer_his_tautomers_from_atoms(p)
    assert result == {("A", 10, ""): "HID"}


def test_infer_his_tautomer_hip_when_both_present(tmp_path: Path) -> None:
    """Both HD1 and HE2 → HIP."""
    from dvbfixer.prep_backend import _infer_his_tautomers_from_atoms

    pdb = HIP_MODEL_PDB + (
        "ATOM     12  HD1 HIS A  10      13.700  11.100  10.000  1.00  0.00           H\n"
    )
    p = tmp_path / "hip.pdb"
    p.write_text(pdb)
    result = _infer_his_tautomers_from_atoms(p)
    assert result == {("A", 10, ""): "HIP"}


def test_infer_his_deprotonated_omitted(tmp_path: Path) -> None:
    """No HD1 and no HE2 → key omitted from map."""
    from dvbfixer.prep_backend import _infer_his_tautomers_from_atoms

    # HE2 line removed.
    pdb = HIP_MODEL_PDB.replace(
        "ATOM     11  HE2 HIS A  10      13.700  15.200  10.000  1.00  0.00           H\n",
        "",
    )
    p = tmp_path / "his_deprot.pdb"
    p.write_text(pdb)
    result = _infer_his_tautomers_from_atoms(p)
    assert result == {}


def test_infer_his_icode_is_key_component(tmp_path: Path) -> None:
    """Two HIS at (A, 10, '') and (A, 10, 'A') are distinct."""
    from dvbfixer.prep_backend import _infer_his_tautomers_from_atoms

    pdb = HIP_MODEL_PDB + (
        "ATOM     12  N   HIS A  10A     20.000  10.000  10.000  1.00  0.00           N\n"
        "ATOM     13  CA  HIS A  10A     21.000  10.000  10.000  1.00  0.00           C\n"
        "ATOM     14  HD1 HIS A  10A     22.000  10.000  10.000  1.00  0.00           H\n"
    )
    p = tmp_path / "his_icode.pdb"
    p.write_text(pdb)
    result = _infer_his_tautomers_from_atoms(p)
    assert result == {("A", 10, ""): "HIE", ("A", 10, "A"): "HID"}


# ---------------------------------------------------------------------------
# assign_amber_variants
# ---------------------------------------------------------------------------


def test_assign_amber_variants_writes_variant(tmp_path: Path) -> None:
    """Standard rename path — ASP → ASH is written to the file.

    ASP must be internal (flanking residues on both sides) or the
    terminal-skip logic would revert the rename.
    """
    from dvbfixer.prep_backend import assign_amber_variants

    p = tmp_path / "asp.pdb"
    p.write_text(PROT_MID_CHAIN_PDB)
    applied = assign_amber_variants(
        p, {("A", 2, ""): "ASH"}, verbose=False,
    )
    assert applied == {("A", 2, ""): "ASH"}
    text = p.read_text()
    assert " ASH A   2 " in text
    assert " ASP A   2 " not in text


def test_assign_amber_variants_terminal_skip_for_ash(
    tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    """ASH at a chain terminus is silently reverted (ff14SB has no NASH)."""
    from dvbfixer.prep_backend import assign_amber_variants

    p = tmp_path / "asp_term.pdb"
    # Only-residue chain — the ASP is both N- and C-terminal.
    p.write_text(ASH_MODEL_PDB)
    applied = assign_amber_variants(
        p, {("A", 5, ""): "ASH"}, verbose=False,
    )
    assert applied == {}
    assert " ASP A   5 " in p.read_text()
    out = capsys.readouterr().out
    assert "skipped 1 terminal residue variant rename" in out


def test_assign_amber_variants_terminal_hip_kept(tmp_path: Path) -> None:
    """HIP has NRES/CRES templates in ff14SB — terminal HIP is kept."""
    from dvbfixer.prep_backend import assign_amber_variants

    p = tmp_path / "his_term.pdb"
    p.write_text(HIP_MODEL_PDB)
    applied = assign_amber_variants(
        p, {("A", 10, ""): "HIP"}, verbose=False,
    )
    assert applied == {("A", 10, ""): "HIP"}
    assert " HIP A  10 " in p.read_text()


def test_assign_amber_variants_empty_map(tmp_path: Path) -> None:
    """Empty map is a no-op (returns empty dict, file unchanged)."""
    from dvbfixer.prep_backend import assign_amber_variants

    p = tmp_path / "unchanged.pdb"
    original = ASH_MODEL_PDB
    p.write_text(original)
    applied = assign_amber_variants(p, {}, verbose=False)
    assert applied == {}
    assert p.read_text() == original


# ---------------------------------------------------------------------------
# _patch_variant_hydrogens
# ---------------------------------------------------------------------------


def test_patch_variant_hydrogens_adds_ash_hd2(tmp_path: Path) -> None:
    """ASH renamed residue gains HD2 opposite OD1."""
    from dvbfixer.prep_backend import _patch_variant_hydrogens

    p = tmp_path / "ash.pdb"
    # First rewrite ASP → ASH so the patch function sees the right resname.
    p.write_text(ASH_MODEL_PDB.replace(" ASP ", " ASH "))
    added, dropped = _patch_variant_hydrogens(
        p, {("A", 5, ""): "ASH"}, verbose=False,
    )
    assert added == 1
    assert dropped == 0
    text = p.read_text()
    assert " HD2" in text


def test_patch_variant_hydrogens_adds_hip_hd1(tmp_path: Path) -> None:
    """HIP with only HE2 placed by Reduce gains HD1."""
    from dvbfixer.prep_backend import _patch_variant_hydrogens

    p = tmp_path / "hip.pdb"
    p.write_text(HIP_MODEL_PDB.replace(" HIS ", " HIP "))
    added, dropped = _patch_variant_hydrogens(
        p, {("A", 10, ""): "HIP"}, verbose=False,
    )
    # HE2 already present in fixture, so exactly one atom (HD1) added.
    assert added == 1
    text = p.read_text()
    assert " HD1" in text
    # HE2 remains too.
    assert text.count(" HE2") == 1


def test_patch_variant_hydrogens_drops_lyn_hz1(tmp_path: Path) -> None:
    """LYN drops HZ1 (ff14SB LYN uses HZ2 + HZ3 only)."""
    from dvbfixer.prep_backend import _patch_variant_hydrogens

    p = tmp_path / "lyn.pdb"
    p.write_text(LYN_MODEL_PDB.replace(" LYS ", " LYN "))
    added, dropped = _patch_variant_hydrogens(
        p, {("A", 20, ""): "LYN"}, verbose=False,
    )
    assert added == 0
    assert dropped == 1
    text = p.read_text()
    assert " HZ1" not in text
    assert " HZ2" in text
    assert " HZ3" in text


# ---------------------------------------------------------------------------
# End-to-end variant-decision integration with mocked PROPKA + subprocesses
# ---------------------------------------------------------------------------


def _fake_pka_results(entries):
    """Build the list shape `get_pka_results` returns from tuples like
    (chain, resnum, icode, restype, pka)."""
    return [
        {"chain": c, "resnum": r, "icode": i, "restype": rt, "pka": pka,
         "model_pka": 0.0}
        for (c, r, i, rt, pka) in entries
    ]


def _install_stub_pipeline(monkeypatch, pka_entries, reduce_atoms):
    """Wire up `run_prep`'s subprocess dependencies with in-memory stubs.

    * `run_tleap` / `run_reduce` become file copies (no real subprocess)
    * `run_propka` / `get_pka_results` return the given fake pKa list
    * `_infer_his_tautomers_from_atoms` returns the given map

    `reduce_atoms` is the value `_infer_his_tautomers_from_atoms` should
    return.
    """
    import dvbfixer.prep_backend as pb

    def _fake_tleap(input_pdb, output_pdb, ff=None, extra_leaprc=None,
                   verbose=False):
        Path(output_pdb).write_text(Path(input_pdb).read_text())

    def _fake_reduce(input_pdb, output_pdb, build=True, nuclear=True,
                    verbose=False):
        Path(output_pdb).write_text(Path(input_pdb).read_text())

    def _fake_run_propka(path):
        return object()

    def _fake_get_pka_results(mc):
        return _fake_pka_results(pka_entries)

    def _fake_his(path):
        return dict(reduce_atoms)

    monkeypatch.setattr(pb, "run_tleap", _fake_tleap)
    monkeypatch.setattr(pb, "run_reduce", _fake_reduce)
    monkeypatch.setattr(pb, "_infer_his_tautomers_from_atoms", _fake_his)
    monkeypatch.setattr(
        "dvbfixer.protonate.run_propka", _fake_run_propka, raising=True,
    )
    monkeypatch.setattr(
        "dvbfixer.protonate.get_pka_results", _fake_get_pka_results,
        raising=True,
    )


PROT_MID_CHAIN_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.500   1.500   0.000  1.00  0.00           C
ATOM      6  N   ASP A   2      12.000  10.000  10.000  1.00  0.00           N
ATOM      7  CA  ASP A   2      13.000  10.000  10.000  1.00  0.00           C
ATOM      8  C   ASP A   2      14.000  10.000  10.000  1.00  0.00           C
ATOM      9  O   ASP A   2      15.000  10.000  10.000  1.00  0.00           O
ATOM     10  CB  ASP A   2      13.500  11.500  10.000  1.00  0.00           C
ATOM     11  CG  ASP A   2      14.500  12.500  10.000  1.00  0.00           C
ATOM     12  OD1 ASP A   2      14.500  13.700  10.000  1.00  0.00           O
ATOM     13  OD2 ASP A   2      15.700  12.100  10.000  1.00  0.00           O
ATOM     14  N   HIS A   3      20.000  10.000  10.000  1.00  0.00           N
ATOM     15  CA  HIS A   3      21.000  10.000  10.000  1.00  0.00           C
ATOM     16  C   HIS A   3      22.000  10.000  10.000  1.00  0.00           C
ATOM     17  O   HIS A   3      23.000  10.000  10.000  1.00  0.00           O
ATOM     18  CB  HIS A   3      21.500  11.500  10.000  1.00  0.00           C
ATOM     19  CG  HIS A   3      22.500  12.500  10.000  1.00  0.00           C
ATOM     20  ND1 HIS A   3      23.700  12.100  10.000  1.00  0.00           N
ATOM     21  CD2 HIS A   3      22.500  13.900  10.000  1.00  0.00           C
ATOM     22  CE1 HIS A   3      24.500  13.100  10.000  1.00  0.00           C
ATOM     23  NE2 HIS A   3      23.700  14.200  10.000  1.00  0.00           N
ATOM     24  HE2 HIS A   3      23.700  15.200  10.000  1.00  0.00           H
ATOM     25  N   ALA A   4      30.000   0.000   0.000  1.00  0.00           N
ATOM     26  CA  ALA A   4      31.000   0.000   0.000  1.00  0.00           C
ATOM     27  C   ALA A   4      32.000   0.000   0.000  1.00  0.00           C
ATOM     28  O   ALA A   4      33.000   0.000   0.000  1.00  0.00           O
ATOM     29  CB  ALA A   4      31.500   1.500   0.000  1.00  0.00           C
TER      30      ALA A   4
"""


def test_run_prep_ash_promoted_from_shifted_pka(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PROPKA gives ASP a pKa 8.5 (>7.0) → residue renamed to ASH.

    Regression: pre-fix run_prep would only rename mid-chain ASP if the
    propka_dict entry survived collision with any N+/C- group on the
    same (chain, resnum), and the H patch was never triggered.
    """
    from dvbfixer.prep_backend import run_prep

    _install_stub_pipeline(
        monkeypatch,
        pka_entries=[
            ("A", 2, "", "ASP", 8.5),  # side chain — should promote to ASH
            ("A", 3, "", "HIS", 5.0),  # normal HIS at pH 7 → neutral
            ("A", 1, "", "N+", 8.0),   # N-terminus — was clobbering ASP before
        ],
        reduce_atoms={("A", 3, ""): "HIE"},
    )

    src = tmp_path / "in.pdb"
    src.write_text(PROT_MID_CHAIN_PDB)
    out = tmp_path / "out.pdb"
    result = run_prep(src, out, ph=7.0, assign_variants=True, verbose=False)

    assert result["renames"][("A", 2, "")] == "ASH"
    assert result["renames"][("A", 3, "")] == "HIE"
    text = out.read_text()
    assert " ASH A   2 " in text
    assert " HIE A   3 " in text


def test_run_prep_hip_from_shifted_his_pka(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PROPKA gives HIS a pKa 8.0 → HIP (charged) even though Reduce only
    placed HE2. Prior behavior: HIS always neutral (HIE/HID) because
    Reduce -build essentially never places both HD1 and HE2."""
    from dvbfixer.prep_backend import run_prep

    _install_stub_pipeline(
        monkeypatch,
        pka_entries=[("A", 3, "", "HIS", 8.0)],
        reduce_atoms={("A", 3, ""): "HIE"},  # Reduce says neutral HIE
    )

    src = tmp_path / "in.pdb"
    src.write_text(PROT_MID_CHAIN_PDB)
    out = tmp_path / "out.pdb"
    result = run_prep(src, out, ph=7.0, assign_variants=True, verbose=False)

    assert result["renames"][("A", 3, "")] == "HIP"
    text = out.read_text()
    assert " HIP A   3 " in text
    # Patcher should have added HD1 (HE2 was already in the input).
    assert " HD1" in text


def test_run_prep_reduce_picks_tautomer_for_neutral_his(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PROPKA says HIS is neutral at pH 7 — Reduce's HID/HIE call wins."""
    from dvbfixer.prep_backend import run_prep

    _install_stub_pipeline(
        monkeypatch,
        pka_entries=[("A", 3, "", "HIS", 5.5)],
        reduce_atoms={("A", 3, ""): "HID"},
    )

    src = tmp_path / "in.pdb"
    src.write_text(PROT_MID_CHAIN_PDB)
    out = tmp_path / "out.pdb"
    result = run_prep(src, out, ph=7.0, assign_variants=True, verbose=False)

    assert result["renames"][("A", 3, "")] == "HID"


def test_run_prep_ss_bonds_force_cyx(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An explicit SS pair → CYX regardless of PROPKA CYM pKa."""
    from dvbfixer.prep_backend import run_prep

    cys_pdb = """\
ATOM      1  N   CYS A   5       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  CYS A   5       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   CYS A   5       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   CYS A   5       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  CB  CYS A   5       1.500   1.500   0.000  1.00  0.00           C
ATOM      6  SG  CYS A   5       2.500   2.500   0.000  1.00  0.00           S
ATOM      7  N   ALA A   6       5.000   0.000   0.000  1.00  0.00           N
ATOM      8  CA  ALA A   6       6.000   0.000   0.000  1.00  0.00           C
ATOM      9  C   ALA A   6       7.000   0.000   0.000  1.00  0.00           C
ATOM     10  O   ALA A   6       8.000   0.000   0.000  1.00  0.00           O
ATOM     11  N   CYS A   7       9.000   0.000   0.000  1.00  0.00           N
ATOM     12  CA  CYS A   7      10.000   0.000   0.000  1.00  0.00           C
ATOM     13  C   CYS A   7      11.000   0.000   0.000  1.00  0.00           C
ATOM     14  O   CYS A   7      12.000   0.000   0.000  1.00  0.00           O
ATOM     15  CB  CYS A   7      10.500   1.500   0.000  1.00  0.00           C
ATOM     16  SG  CYS A   7      11.500   2.500   0.000  1.00  0.00           S
TER      17      CYS A   7
"""
    _install_stub_pipeline(
        monkeypatch,
        pka_entries=[
            ("A", 5, "", "CYS", 6.0),  # would be CYM at pH 7 without SS
            ("A", 7, "", "CYS", 5.5),
        ],
        reduce_atoms={},
    )

    src = tmp_path / "in.pdb"
    src.write_text(cys_pdb)
    out = tmp_path / "out.pdb"
    result = run_prep(
        src, out, ph=7.0, assign_variants=True,
        ss_pairs={("A", 5), ("A", 7)}, verbose=False,
    )
    assert result["renames"][("A", 5, "")] == "CYX"
    assert result["renames"][("A", 7, "")] == "CYX"
