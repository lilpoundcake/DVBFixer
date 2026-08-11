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

from dvbfixer.ffutils.geometry import build_ca_chirality_force, fix_ca_chirality  # noqa: E402


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


def test_signed_volume_guard_is_zero_for_l_and_penalizes_d() -> None:
    from openmm import Context, System, VerletIntegrator
    from openmm.unit import kilojoule_per_mole, picosecond

    top, pos, ix = _make_residue(
        "VAL", (-0.87, 1.21, 0.00), (0.00, 0.00, 0.00),
        (1.44, 0.00, -0.20), (-0.86, -0.87, -0.87),
    )
    force, count = build_ca_chirality_force(top, pos)
    assert count == 1
    system = System()
    for _ in range(4):
        system.addParticle(12.0)
    system.addForce(force)
    context = Context(system, VerletIntegrator(0.001 * picosecond))
    context.setPositions(pos)
    e_l = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    d_pos = Quantity(list(pos.value_in_unit(nanometer)), nanometer)
    cb = d_pos[ix["CB"]].value_in_unit(nanometer)
    d_pos[ix["CB"]] = Quantity((cb[0], cb[1], -cb[2]), nanometer)
    context.setPositions(d_pos)
    state = context.getState(getEnergy=True, getForces=True)
    e_d = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    assert e_l == pytest.approx(0.0, abs=1e-10)
    assert e_d > 0.0
    assert all(value == value for force_vec in state.getForces() for value in force_vec)


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


def test_d_serine_reflects_whole_sidechain() -> None:
    """A D-configured SER with N/CA/C/O/CB/OG must have BOTH CB AND OG
    reflected (whole sidechain), not just CB. The CB-OG bond length
    must survive the reflection."""
    import math
    top, pos, ix = _make_residue(
        "SER",
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        cb_pos=(-0.86, -0.87, 0.87),   # D config (z flipped)
    )
    # Add O and OG so the test exercises the whole-sidechain path.
    top_chain = next(iter(top.chains()))
    res = next(iter(top_chain.residues()))
    top.addAtom("O", Element.getBySymbol("O"), res)
    top.addAtom("OG", Element.getBySymbol("O"), res)

    # Rebuild positions with O (backbone) and OG (sidechain).
    o_pos = (1.6, 1.2, -0.3)
    og_pos = (-1.5, -1.4, 1.5)   # 1.42 Å from D-CB
    # Reconstruct nm-scaled positions from existing (nm) + new (Å→nm).
    from openmm.unit import Quantity as _Q
    existing_nm = [tuple(pos[i].value_in_unit(nanometer))
                   for i in range(4)]  # N, CA, C, CB
    new_positions_nm = existing_nm + [
        tuple(x / 10.0 for x in o_pos),
        tuple(x / 10.0 for x in og_pos),
    ]
    pos = _Q(new_positions_nm, nanometer)
    ix["O"] = 4
    ix["OG"] = 5

    # Distance between CB and OG BEFORE the flip.
    def _get(idx):
        p = pos[idx].value_in_unit(nanometer)
        return float(p[0]), float(p[1]), float(p[2])
    cb_before = _get(ix["CB"])
    og_before = _get(ix["OG"])
    d_before = math.sqrt(sum((a - b) ** 2 for a, b in zip(cb_before, og_before)))

    triple_before = _triple(pos, ix)
    assert triple_before < 0, "test fixture is not D"

    fix_ca_chirality(top, pos)

    triple_after = _triple(pos, ix)
    assert triple_after > 0, "chirality was not flipped"

    # CB-OG bond length preserved (both were reflected together).
    cb_after = _get(ix["CB"])
    og_after = _get(ix["OG"])
    d_after = math.sqrt(sum((a - b) ** 2 for a, b in zip(cb_after, og_after)))
    assert abs(d_after - d_before) < 1e-6, (
        f"CB-OG bond length changed: {d_before:.4f} → {d_after:.4f}"
    )

    # Backbone O unchanged.
    assert _get(ix["O"]) == tuple(x / 10.0 for x in o_pos)


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


# ---------------------------------------------------------------------------
# 0.7.1 defense-in-depth: ChiralityError, restraint, iterative loop
# ---------------------------------------------------------------------------


def _perturb_cb_to_target_triple(
    target_triple_nm3: float,
) -> tuple[Topology, Quantity, dict[str, int]]:
    """Build a VAL fixture whose CB gives ``triple ≈ target_triple_nm3``.

    Uses the standard N/CA/C backbone and picks CB on the negative z
    side (D configuration) with |z| tuned so the triple lands at the
    requested (typically small-negative) value. Useful for probing the
    new -1e-6 deadband boundary.
    """
    # Backbone identical to the other tests.
    n_pos = (-0.87, 1.21, 0.00)
    ca_pos = (0.00, 0.00, 0.00)
    c_pos = (1.44, 0.00, -0.20)
    # Compute cross = (N-CA) × (C-CA) in nm³. Backbone in Å here, so
    # divide by 10 when we go to nm at the end.
    vn = tuple(a - b for a, b in zip(n_pos, ca_pos))
    vc = tuple(a - b for a, b in zip(c_pos, ca_pos))
    cross = (
        vn[1] * vc[2] - vn[2] * vc[1],
        vn[2] * vc[0] - vn[0] * vc[2],
        vn[0] * vc[1] - vn[1] * vc[0],
    )
    # Want (cross · CB) in nm³ = target. cross has units of Å²; CB in
    # Å too, so cross·CB in Å³ → * 1e-3 for nm³. Pick CB = (0,0,z).
    # target_nm3 = (cross_z * z) * 1e-3 → z = target_nm3 * 1e3 / cross_z.
    if abs(cross[2]) < 1e-12:
        raise ValueError("degenerate backbone in fixture")
    z = target_triple_nm3 * 1e3 / cross[2]
    cb_pos = (0.0, 0.0, z)
    return _make_residue("VAL", n_pos=n_pos, ca_pos=ca_pos, c_pos=c_pos,
                         cb_pos=cb_pos)


def test_marginal_d_caught_by_tightened_deadband() -> None:
    """A CB perturbation that gives triple = -5e-5 nm³ (just past the
    old -1e-4 threshold but well past the new -1e-6 one) must now be
    detected and reflected."""
    top, pos, ix = _perturb_cb_to_target_triple(-5e-5)
    t_before = _triple(pos, ix)
    assert -1e-4 < t_before < -1e-6, (
        f"fixture triple {t_before} outside the marginal band"
    )

    repairs = fix_ca_chirality(top, pos)
    assert repairs == 1, "marginal D was not caught by the tightened deadband"
    assert _triple(pos, ix) > 0


def test_assert_all_l_raises_on_d() -> None:
    """assert_all_l raises ChiralityError with the offender listed."""
    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l

    top, pos, _ = _make_residue(
        "VAL",
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        cb_pos=(-0.86, -0.87, 0.87),   # D
    )
    with pytest.raises(ChiralityError) as excinfo:
        assert_all_l(top, pos)
    assert len(excinfo.value.residues) == 1
    ch, rid, name, _tri = excinfo.value.residues[0]
    assert ch == "A"
    assert name == "VAL"


def test_assert_all_l_passes_on_l() -> None:
    """No exception on canonical L geometry."""
    from dvbfixer.ffutils.geometry import assert_all_l

    top, pos, _ = _make_residue(
        "VAL",
        n_pos=(-0.87, 1.21, 0.00),
        ca_pos=(0.00, 0.00, 0.00),
        c_pos=(1.44, 0.00, -0.20),
        cb_pos=(-0.86, -0.87, -0.87),
    )
    assert_all_l(top, pos)  # must not raise


def test_find_d_residues_lists_all() -> None:
    """find_d_residues walks the whole topology, not just the first hit."""
    from dvbfixer.ffutils.geometry import find_d_residues

    top = Topology()
    chain = top.addChain("A")
    positions_nm: list[tuple[float, float, float]] = []
    # 3 residues: L, D, L.
    for i, cb_z in enumerate((-0.087, 0.087, -0.087)):
        r = top.addResidue("VAL", chain)
        for name, sym, p in (
            ("N", "N", (-0.087, 0.121, 0.000)),
            ("CA", "C", (0.000, 0.000, 0.000)),
            ("C", "C", (0.144, 0.000, -0.020)),
            ("CB", "C", (-0.086, -0.087, cb_z)),
        ):
            top.addAtom(name, Element.getBySymbol(sym), r)
            # Shift each residue along x so atoms don't overlap
            # (positional shift affects only presentation; triple is
            # translation-invariant so the chirality signal survives).
            positions_nm.append((p[0] + 2.0 * i, p[1], p[2]))
        r.id = str(i + 1)
    pos = Quantity(positions_nm, nanometer)

    offenders = find_d_residues(top, pos)
    assert len(offenders) == 1
    assert offenders[0][1] == "2"


# The three CustomTorsionForce-related tests
# (test_improper_chirality_restraint_counts_torsions,
#  test_improper_restraint_skips_degenerate_residues,
#  test_chirality_restraint_minimize_pulls_d_to_l)
# were removed alongside `improper_chirality_restraint` itself. Field
# consensus is reflect-then-verify (no runtime harmonic force), and the
# existing reflection tests above cover the shipping primitive.
