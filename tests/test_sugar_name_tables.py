"""Regression test for the 2026-07-31 audit's sugar-name-table fixes.

Three independently-maintained sugar-name sets had drifted apart —
`ffutils._PDB_SUGAR_NAMES`, `pdbutils.inference`'s own (now-removed)
copies, and `top.ff_data.PDB_TO_CARB` — missing different subsets of
names from each other, including `BNE5AC` (beta-anomer sialic acid,
a real, distinct CHARMM36 `carb.rtp` residue) from more than one at
once. `top.ff_data.PDB_TO_CARB` in particular was missing it entirely,
silently dropping a BNE5AC residue's glycosidic bonds from GROMACS
topology output with zero warning.
"""
from __future__ import annotations

from dvbfixer.ffutils import is_pdb_sugar_resname
from dvbfixer.top.ff_data import PDB_TO_CARB


def test_bne5ac_recognized_as_sugar() -> None:
    assert is_pdb_sugar_resname("BNE5AC")
    assert is_pdb_sugar_resname("ANE5AC")  # sibling alpha form, sanity check


def test_bne5ac_present_in_pdb_to_carb() -> None:
    assert "BNE5AC" in PDB_TO_CARB
    assert PDB_TO_CARB["BNE5AC"] == "BNE5AC"
    # BNE5 (short alias, mirroring the existing ANE5 -> ANE5AC alias)
    assert PDB_TO_CARB.get("BNE5") == "BNE5AC"


def test_charmm_gui_sugar_names_recognized() -> None:
    """Names present in `top.ff_data.PDB_TO_CARB`'s CHARMM-GUI-native
    block must all be recognized as sugars by the canonical helper —
    these previously lived in a second, independently-drifted copy in
    `pdbutils.inference` that was missing some of them."""
    for name in ("BGLC", "BGAL", "AFUC", "AMAN", "BMAN", "BGLCNA",
                 "BGALNA", "AGALNA", "ANE5AC", "BNE5AC"):
        assert is_pdb_sugar_resname(name), f"{name} not recognized as a sugar"


def test_uncommon_pdb_sugar_names_recognized() -> None:
    """NGA/A2G/XYP/etc. — present in the broader PDB name set that used
    to live only in pdbutils.inference, absent from ffutils' own
    (narrower) copy before consolidation."""
    for name in ("NGA", "A2G", "BGA", "AGA", "XYP", "XYS", "RIB", "RIP",
                 "ARA", "NAN"):
        assert is_pdb_sugar_resname(name), f"{name} not recognized as a sugar"
