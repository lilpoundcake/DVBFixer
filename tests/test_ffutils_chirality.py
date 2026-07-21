"""Unit tests for `dvbfixer.ffutils.geometry.fix_ca_chirality`.

Regression coverage for the 0.6.2 bug where PDBFixer's
``addMissingAtoms`` rebuilds branched-Cβ sidechains (VAL / ILE / THR)
with D stereochemistry. The fix reflects CB through the CA-N-C plane.
"""

from __future__ import annotations

import pytest

pytest.importorskip("openmm", reason="chirality fix needs OpenMM")

from openmm.app import Element, Topology  # noqa: E402
from openmm.unit import Quantity, nanometer  # noqa: E402

from dvbfixer.ffutils.geometry import fix_ca_chirality  # noqa: E402


def _make_residue(
    resname: str,
    n_pos: tuple[float, float, float],
    ca_pos: tuple[float, float, float],
    c_pos: tuple[float, float, float],
    cb_pos: tuple[float, float, float] | None,
) -> tuple[Topology, Quantity, dict[str, int]]:
    """Build a topology holding a single residue with N/CA/C(+ optional CB).

    Positions are given in ångström; converted to nm internally.
    Returns (topology, positions, name→atom_index dict).
    """
    top = Topology()
    chain = top.addChain("A")
    r = top.addResidue(resname, chain)
    index_of: dict[str, int] = {}
    order: list[tuple[str, str, tuple[float, float, float]]] = [
        ("N", "N", n_pos),
        ("CA", "C", ca_pos),
        ("C", "C", c_pos),
    ]
    if cb_pos is not None:
        order.append(("CB", "C", cb_pos))
    for name, sym, _ in order:
        top.addAtom(name, Element.getBySymbol(sym), r)
    for i, (name, _, _) in enumerate(order):
        index_of[name] = i
    positions_nm = [tuple(x / 10.0 for x in p) for _, _, p in order]
    return top, Quantity(positions_nm, nanometer), index_of


def _triple(pos: Quantity, ix: dict[str, int]) -> float:
    """(N - CA) × (C - CA) · (CB - CA) in nm³."""
    def _get(idx: int) -> tuple[float, float, float]:
        p = pos[idx].value_in_unit(nanometer)
        return float(p[0]), float(p[1]), float(p[2])
    ca = _get(ix["CA"])
    vn = tuple(a - b for a, b in zip(_get(ix["N"]), ca))
    vc = tuple(a - b for a, b in zip(_get(ix["C"]), ca))
    vcb = tuple(a - b for a, b in zip(_get(ix["CB"]), ca))
    cross = (
        vn[1] * vc[2] - vn[2] * vc[1],
        vn[2] * vc[0] - vn[0] * vc[2],
        vn[0] * vc[1] - vn[1] * vc[0],
    )
    return cross[0] * vcb[0] + cross[1] * vcb[1] + cross[2] * vcb[2]


def test_l_valine_untouched() -> None:
    """Canonical L geometry: fix_ca_chirality returns 0, CB unchanged."""
    top, pos, ix = _make_residue(
        "VAL",
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        cb_pos=(-0.86, -0.87, -0.87),
    )
    triple_before = _triple(pos, ix)
    assert triple_before > 0, "test fixture is not L geometry"
    original_cb = tuple(pos[ix["CB"]].value_in_unit(nanometer))

    repairs = fix_ca_chirality(top, pos)

    assert repairs == 0
    after_cb = tuple(pos[ix["CB"]].value_in_unit(nanometer))
    assert after_cb == original_cb


def test_d_valine_flipped_to_l() -> None:
    """D-CB → L via mirror through CA-N-C plane."""
    top, pos, ix = _make_residue(
        "VAL",
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        # Mirror CB in z — flips triple sign.
        cb_pos=(-0.86, -0.87, 0.87),
    )
    triple_before = _triple(pos, ix)
    assert triple_before < 0, "test fixture is not D geometry"
    original_ca = tuple(pos[ix["CA"]].value_in_unit(nanometer))
    original_n = tuple(pos[ix["N"]].value_in_unit(nanometer))
    original_c = tuple(pos[ix["C"]].value_in_unit(nanometer))

    repairs = fix_ca_chirality(top, pos)

    assert repairs == 1
    triple_after = _triple(pos, ix)
    assert triple_after > 0
    # Backbone atoms untouched.
    assert tuple(pos[ix["CA"]].value_in_unit(nanometer)) == original_ca
    assert tuple(pos[ix["N"]].value_in_unit(nanometer)) == original_n
    assert tuple(pos[ix["C"]].value_in_unit(nanometer)) == original_c


def test_glycine_skipped() -> None:
    """GLY has no CB — return 0 without error."""
    top, pos, _ = _make_residue(
        "GLY",
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        cb_pos=None,
    )
    repairs = fix_ca_chirality(top, pos)
    assert repairs == 0


def test_residue_missing_cb_skipped() -> None:
    """Non-GLY residue that happens to have no CB — should skip
    silently (defensive, since PDBFixer may not have finished yet)."""
    top, pos, _ = _make_residue(
        "VAL",
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        cb_pos=None,
    )
    repairs = fix_ca_chirality(top, pos)
    assert repairs == 0
