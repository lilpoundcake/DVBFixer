"""Unit tests for `dvbfixer.diagnose.steric` — bond-exclusion behaviour
and clash-threshold calibration.

Regression coverage for the 0.5.1 bug where `topology.bonds()` was
empty for HETATMs / non-standard residues and their real covalent
bonds were reported as ERROR clashes.
"""

from __future__ import annotations

import pytest

pytest.importorskip("openmm", reason="steric checks need OpenMM")
pytest.importorskip("scipy", reason="steric checks need scipy")

from openmm.app import Element, Topology  # noqa: E402
from openmm.unit import Quantity, nanometer  # noqa: E402

from dvbfixer.diagnose.report import Severity  # noqa: E402
from dvbfixer.diagnose.steric import clashes_python  # noqa: E402


def _make_topology(
    atom_specs: list[tuple[str, str]],
    positions_a: list[tuple[float, float, float]],
    resname: str = "LIG",
    chain_id: str = "X",
    add_bonds: bool = False,
) -> tuple[Topology, Quantity]:
    """Build a bare topology with one residue holding ``atom_specs``.

    ``atom_specs`` is a list of (atom_name, element_symbol). If
    ``add_bonds`` is False (the default), NO topology bonds are added
    — this mirrors the HETATM / non-standard-residue case where
    OpenMM's PDBFile parser can't infer bonds. Positions are given in
    ångström and internally converted to nanometer.
    """
    top = Topology()
    chain = top.addChain(chain_id)
    res = top.addResidue(resname, chain)
    atoms = []
    for name, sym in atom_specs:
        elem = Element.getBySymbol(sym)
        atoms.append(top.addAtom(name, elem, res))
    if add_bonds:
        for a, b in zip(atoms, atoms[1:]):
            top.addBond(a, b)
    positions_nm = [tuple(x / 10.0 for x in p) for p in positions_a]
    return top, Quantity(positions_nm, nanometer)


def test_bonded_hetatm_cc_not_reported_as_clash() -> None:
    """A two-atom HETATM with C-C at ~1.53 Å (real bond) must be
    recognised as bonded via distance inference and NOT reported as a
    clash, even though ``topology.bonds()`` is empty."""
    top, pos = _make_topology(
        [("C1", "C"), ("C2", "C")],
        [(0.0, 0.0, 0.0), (1.53, 0.0, 0.0)],
        add_bonds=False,
    )
    findings = clashes_python(top, pos)
    assert findings == [], (
        f"Expected zero clash findings (bond inference should exclude "
        f"C-C at 1.53 Å), got: {[str(f) for f in findings]}"
    )


def test_bonded_hetatm_ch_not_reported_as_clash() -> None:
    """C-H at 1.09 Å (typical X-H bond) is recognised as bonded."""
    top, pos = _make_topology(
        [("C1", "C"), ("H1", "H")],
        [(0.0, 0.0, 0.0), (1.09, 0.0, 0.0)],
        add_bonds=False,
    )
    findings = clashes_python(top, pos)
    assert findings == []


def test_bonded_hetatm_ss_disulfide_not_reported() -> None:
    """Disulfide S-S at 2.05 Å (typical bond length) is recognised as
    bonded despite exceeding the heavy-heavy default cutoff."""
    top, pos = _make_topology(
        [("SG", "S"), ("SG", "S")],
        [(0.0, 0.0, 0.0), (2.05, 0.0, 0.0)],
        add_bonds=False,
    )
    findings = clashes_python(top, pos)
    assert findings == []


def test_true_cc_clash_still_error() -> None:
    """Two carbon atoms at 2.0 Å with NO covalent relationship must
    still be reported as an ERROR clash (overlap = 3.40 - 2.0 = 1.4 Å,
    well above the 0.5 Å ERROR threshold). Constructed as two
    single-atom residues so the geometry-inference cutoff (1.90 Å
    heavy-heavy) doesn't accidentally treat them as bonded."""
    top = Topology()
    chain = top.addChain("X")
    r1 = top.addResidue("LG1", chain)
    r2 = top.addResidue("LG2", chain)
    top.addAtom("C1", Element.getBySymbol("C"), r1)
    top.addAtom("C1", Element.getBySymbol("C"), r2)
    # 2.00 Å apart — above the 1.90 Å heavy-heavy bond cutoff, so
    # geometric inference correctly treats them as non-bonded.
    positions_nm = [(0.0, 0.0, 0.0), (0.200, 0.0, 0.0)]
    pos = Quantity(positions_nm, nanometer)
    findings = clashes_python(top, pos)
    assert len(findings) == 1
    assert findings[0].severity == Severity.ERROR
    assert findings[0].category == "clash"


def _make_two_residue_topology(
    r1_atoms: list[tuple[str, str]],
    r2_atoms: list[tuple[str, str]],
    positions_a: list[tuple[float, float, float]],
    r1_bonds: list[tuple[int, int]] | None = None,
    r2_bonds: list[tuple[int, int]] | None = None,
    r1_name: str = "AAA",
    r2_name: str = "BBB",
) -> tuple[Topology, Quantity]:
    """Two-residue topology with explicit intra-residue bonds. Indices
    in ``r1_bonds`` / ``r2_bonds`` are 0-based within each residue's
    ``atoms`` list. Positions in ångström."""
    top = Topology()
    chain = top.addChain("A")
    r1 = top.addResidue(r1_name, chain)
    r2 = top.addResidue(r2_name, chain)
    atoms_r1 = [top.addAtom(name, Element.getBySymbol(sym), r1)
                for name, sym in r1_atoms]
    atoms_r2 = [top.addAtom(name, Element.getBySymbol(sym), r2)
                for name, sym in r2_atoms]
    for i, j in (r1_bonds or []):
        top.addBond(atoms_r1[i], atoms_r1[j])
    for i, j in (r2_bonds or []):
        top.addBond(atoms_r2[i], atoms_r2[j])
    positions_nm = [tuple(x / 10.0 for x in p) for p in positions_a]
    return top, Quantity(positions_nm, nanometer)


def test_backbone_amide_hbond_not_reported_as_clash() -> None:
    """Classic backbone N-H..O=C hydrogen bond at 2.06 Å (α-helix i,i+4
    interaction distance) must be recognised as an H-bond and skipped,
    not reported as an ERROR clash. This is the primary 0.5.2
    regression case."""
    # Residue 1: carbonyl (C=O). Residue 2: amide (N-H) with the H
    # positioned 2.06 Å from the residue-1 O.
    top, pos = _make_two_residue_topology(
        r1_atoms=[("C", "C"), ("O", "O")],
        r2_atoms=[("N", "N"), ("H", "H")],
        positions_a=[
            (0.00, 0.00, 0.00),   # C
            (1.23, 0.00, 0.00),   # O (bonded to C)
            (3.29, 0.00, 0.00),   # N (2.06 Å from O)
            (2.20, 0.00, 0.00),   # H (bonded to N, 2.06 Å from O — H-bond!)
        ],
        r1_bonds=[(0, 1)],  # C=O
        r2_bonds=[(0, 1)],  # N-H
    )
    findings = clashes_python(top, pos)
    assert findings == [], (
        f"Expected zero findings (backbone amide H-bond), got: "
        f"{[str(f) for f in findings]}"
    )


def test_salt_bridge_arg_h_to_asp_o_not_reported() -> None:
    """ARG-guanidinium HH..carboxylate O at 1.79 Å (a strong salt
    bridge) must not be reported as a clash."""
    top, pos = _make_two_residue_topology(
        r1_atoms=[("NH1", "N"), ("HH11", "H")],   # ARG guanidinium
        r2_atoms=[("CG", "C"), ("OD1", "O")],      # ASP carboxylate
        positions_a=[
            (0.00, 0.00, 0.00),   # NH1
            (1.00, 0.00, 0.00),   # HH11 (bonded to NH1)
            (3.53, 0.00, 0.00),   # CG
            (2.79, 0.00, 0.00),   # OD1 (bonded to CG; 1.79 Å from HH11)
        ],
        r1_bonds=[(0, 1)],
        r2_bonds=[(0, 1)],
        r1_name="ARG",
        r2_name="ASP",
    )
    findings = clashes_python(top, pos)
    assert findings == []


def test_non_polar_h_close_to_o_still_flagged() -> None:
    """A non-polar H (bonded to C) close to an O must NOT be treated
    as an H-bond — polar-donor test protects genuine C-H..O clashes.
    Places the C-bonded H at 1.85 Å from a separate residue's O."""
    top, pos = _make_two_residue_topology(
        r1_atoms=[("C", "C"), ("H", "H")],       # C-H (aliphatic)
        r2_atoms=[("C", "C"), ("O", "O")],       # C=O
        positions_a=[
            (0.00, 0.00, 0.00),   # C
            (1.09, 0.00, 0.00),   # H (bonded to C, non-polar)
            (4.00, 0.00, 0.00),   # C
            (2.94, 0.00, 0.00),   # O (bonded to C; 1.85 Å from non-polar H)
        ],
        r1_bonds=[(0, 1)],
        r2_bonds=[(0, 1)],
    )
    findings = clashes_python(top, pos)
    # H..O overlap = 2.72 - 1.85 = 0.87 Å → still ERROR.
    assert len(findings) == 1
    assert findings[0].severity == Severity.ERROR


def test_mild_rotamer_hh_below_new_warning_floor() -> None:
    """H-H at 2.05 Å between two separate residues — overlap =
    2.40 - 2.05 = 0.35 Å, below the new 0.4 Å MolProbity floor.
    Under the old 0.2 Å floor this would have been WARNING; under
    0.5.1 it must be silent."""
    top = Topology()
    chain = top.addChain("X")
    r1 = top.addResidue("LG1", chain)
    r2 = top.addResidue("LG2", chain)
    top.addAtom("H1", Element.getBySymbol("H"), r1)
    top.addAtom("H1", Element.getBySymbol("H"), r2)
    positions_nm = [(0.0, 0.0, 0.0), (0.205, 0.0, 0.0)]
    pos = Quantity(positions_nm, nanometer)
    findings = clashes_python(top, pos)
    # 0.35 Å overlap — below WARNING floor, no findings expected.
    assert findings == []
