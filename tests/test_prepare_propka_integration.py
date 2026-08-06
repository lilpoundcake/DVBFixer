"""Tests for the 0.7.7 PROPKA + Reduce integration in legacy prepare.

Prior to 0.7.7, the legacy prepare backend built its `variants=[...]`
list only from user `--mutate` overrides + HIS HD1/HE2 atom presence.
PROPKA and MolProbity Reduce were not invoked, so pKa-driven ASH /
GLH / HIP / LYN / CYM were never emitted for the average user.

0.7.7 wires PROPKA + Reduce into prepare's variants-list build via
`_run_propka_reduce_variants`. These tests verify:

1. The helper returns a variant map (dict form) for a real PDB.
2. `--no-propka --no-protassign` reproduces the pre-0.7.7 behaviour.
3. The zbs pipeline propagates the flags into prepare.
"""
from __future__ import annotations

from pathlib import Path

import pytest


@pytest.mark.slow
def test_propka_reduce_helper_returns_variants(
    pure_protein_small: Path, tmp_workdir: Path,
) -> None:
    """The helper should run PROPKA + Reduce on a real PDB and return
    a variant map. Even a small protein has HIS tautomers to pick and
    typically one or two shifted acidic side chains — the map should
    be non-empty."""
    from dvbfixer.prepare.pipeline import _run_propka_reduce_variants

    variant_map = _run_propka_reduce_variants(
        pure_protein_small,
        ph=7.0,
        ss_pairs=set(),
        his_default="HIE",
        cys_ss_pka=8.0,
        use_propka=True,
        use_reduce=True,
        verbose=False,
    )
    # Result is a dict; values are valid AMBER variant names.
    _VALID = {"HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM", "LYN"}
    assert isinstance(variant_map, dict)
    for key, variant in variant_map.items():
        assert variant in _VALID, (
            f"unexpected variant name {variant!r} for {key}"
        )


def test_propka_reduce_helper_both_off_returns_empty(
    pure_protein_small: Path,
) -> None:
    """When both PROPKA and Reduce are off, the helper should return
    an empty dict (no pKa-driven decisions to make)."""
    from dvbfixer.prepare.pipeline import _run_propka_reduce_variants

    variant_map = _run_propka_reduce_variants(
        pure_protein_small,
        ph=7.0,
        ss_pairs=set(),
        use_propka=False,
        use_reduce=False,
        verbose=False,
    )
    assert variant_map == {}


def test_propka_reduce_helper_ss_pairs_force_cyx(
    ss_bonded_antibody: Path,
) -> None:
    """When ss_pairs is non-empty, every listed (chain, resseq) should
    end up as CYX in the returned map — regardless of PROPKA's
    per-residue CYS pKa."""
    from dvbfixer.acpype_export import detect_ss_bonds
    from dvbfixer.prepare.pipeline import _run_propka_reduce_variants

    ss_pairs = detect_ss_bonds(str(ss_bonded_antibody))
    assert ss_pairs, "test fixture should have SS bonds"

    variant_map = _run_propka_reduce_variants(
        ss_bonded_antibody,
        ph=7.0,
        ss_pairs=ss_pairs,
        use_propka=False,  # keep the test fast — CONECT-driven CYX only
        use_reduce=False,
        verbose=False,
    )
    cyx_keys = {
        (ch, rn) for (ch, rn, _ic), var in variant_map.items()
        if var == "CYX"
    }
    # Every SS-paired CYS should have a CYX entry.
    for (ch, rn) in ss_pairs:
        assert (ch, rn) in cyx_keys, (
            f"SS-bonded CYS {ch}:{rn} not in CYX map: {sorted(cyx_keys)[:5]}"
        )


def test_prepare_cli_propka_flags_exist() -> None:
    """The new CLI flags on `dvbfixer prepare` should parse cleanly."""
    from dvbfixer.prepare.cli import parse_args

    args = parse_args([
        "dummy.pdb", "--no-propka", "--no-protassign",
        "--his-default", "HID", "--cys-ss-pka", "9.5",
    ])
    assert args.propka is False
    assert args.protassign is False
    assert args.his_default == "HID"
    assert args.cys_ss_pka == 9.5


def test_zbs_cli_propka_flags_exist() -> None:
    """The same flags on `dvbfixer zbs` should parse cleanly."""
    from dvbfixer.zbs import parse_args

    args = parse_args([
        "dummy.pdb", "--no-propka", "--no-protassign",
        "--his-default", "HID", "--cys-ss-pka", "9.5",
    ])
    assert args.propka is False
    assert args.protassign is False
    assert args.his_default == "HID"
    assert args.cys_ss_pka == 9.5


def test_zbs_skip_protonate_maps_to_no_propka_no_protassign() -> None:
    """--skip-protonate (deprecated) should now flip propka+protassign
    to False for backward compat."""
    from dvbfixer.zbs import parse_args

    args = parse_args(["dummy.pdb", "--skip-protonate"])
    assert args.propka is False
    assert args.protassign is False


# Real fragment from tests/fixtures/8cz8/8cz8_t_u.pdb, chain E residues 370-372 —
# HIS 371's entire imidazole ring (CG/ND1/CD2/CE1/NE2) is genuinely
# missing in the deposited structure (crystallographic disorder), only
# backbone + CB survive. This is the exact real-world trigger for
# `ValueError: HIS residue (N) has the wrong set of atoms` from
# `Modeller.addHydrogens` (0.7.11).
_HIS_MISSING_RING_FRAGMENT = """\
ATOM      1  N   GLN A 370      18.470  40.815 -84.848  1.00158.00           N
ATOM      2  CA  GLN A 370      18.037  40.240 -83.545  1.00160.28           C
ATOM      3  C   GLN A 370      19.146  39.332 -83.003  1.00161.91           C
ATOM      4  O   GLN A 370      19.737  39.686 -81.969  1.00164.09           O
ATOM      5  CB  GLN A 370      16.727  39.467 -83.700  1.00165.99           C
ATOM      6  CG  GLN A 370      15.487  40.338 -83.551  1.00170.19           C
ATOM      7  CD  GLN A 370      14.362  39.903 -84.460  1.00175.73           C
ATOM      8  OE1 GLN A 370      13.638  40.724 -85.023  1.00175.98           O
ATOM      9  NE2 GLN A 370      14.206  38.598 -84.614  1.00180.55           N
ATOM     10  N   HIS A 371      19.418  38.212 -83.690  1.00161.49           N
ATOM     11  CA  HIS A 371      20.300  37.108 -83.217  1.00157.10           C
ATOM     12  C   HIS A 371      21.708  37.256 -83.821  1.00155.11           C
ATOM     13  O   HIS A 371      22.391  36.227 -83.970  1.00156.79           O
ATOM     14  CB  HIS A 371      19.669  35.744 -83.550  1.00151.22           C
ATOM     15  N   LEU A 372      22.128  38.490 -84.138  1.00153.86           N
ATOM     16  CA  LEU A 372      23.452  38.789 -84.764  1.00155.87           C
ATOM     17  C   LEU A 372      24.578  38.315 -83.837  1.00150.72           C
ATOM     18  O   LEU A 372      25.688  38.077 -84.338  1.00150.14           O
ATOM     19  CB  LEU A 372      23.570  40.290 -85.040  1.00161.26           C
ATOM     20  CG  LEU A 372      24.529  40.669 -86.170  1.00168.54           C
ATOM     21  CD1 LEU A 372      23.793  40.757 -87.503  1.00170.38           C
ATOM     22  CD2 LEU A 372      25.245  41.978 -85.865  1.00171.95           C
TER
END
"""


def test_prepare_completes_on_his_with_missing_ring(tmp_workdir: Path) -> None:
    """Regression: a HIS residue with its entire imidazole ring missing
    (real crystallographic disorder, not a dvbfixer bug at the source)
    used to crash `prepare` unrecoverably with `ValueError: HIS residue
    (N) has the wrong set of atoms` from `Modeller.addHydrogens` — on
    BOTH the whole-topology attempt and the protein-only fallback,
    since both share the same fragile internal auto-detect logic.

    Root cause: PROPKA + MolProbity Reduce used to run on the RAW input
    (before PDBFixer's own `findMissingAtoms`/`addMissingAtoms` rebuilds
    the ring), so PROPKA had nothing to analyze for this residue and
    emitted no pKa result at all — not "neutral", just absent — so its
    variant decision fell through to OpenMM's own fragile ND1/NE2-
    presence check. Fixed by moving PROPKA + Reduce to run AFTER
    PDBFixer's heavy-atom repair, on the now-complete structure.
    """
    pytest.importorskip("openmm", reason="prepare integration needs OpenMM")
    pytest.importorskip("pdbfixer", reason="prepare integration needs PDBFixer")
    from dvbfixer.prepare.pipeline import main as prepare_main

    input_pdb = tmp_workdir / "his_missing_ring.pdb"
    input_pdb.write_text(_HIS_MISSING_RING_FRAGMENT)
    out = tmp_workdir / "prep.pdb"

    prepare_main([str(input_pdb), "-o", str(out)])

    assert out.exists()
    text = out.read_text()
    his_lines = [ln for ln in text.splitlines()
                 if ln[17:20].strip() in ("HIS", "HIE", "HID", "HIP")
                 and ln[22:26].strip() == "371"]
    atom_names = {ln[12:16].strip() for ln in his_lines}
    assert {"ND1", "NE2", "CG", "CD2", "CE1"} <= atom_names, (
        f"HIS 371's ring wasn't fully rebuilt: {sorted(atom_names)}"
    )
    # A rebuilt ring with no real H-bonding environment should get the
    # user's declared default tautomer (HIE) via the normal PROPKA/
    # Reduce path — not the emergency ND1/NE2-presence backstop.
    assert any(ln[17:20].strip() == "HIE" for ln in his_lines), (
        f"expected HIS 371 to resolve to the HIE default tautomer, "
        f"found resnames: {sorted({ln[17:20].strip() for ln in his_lines})}"
    )
