"""Regression: prepare must fix broken input geometry (coincident atoms).

The reported case: a C-terminal SER whose input file placed HG at the
same coordinates as OXT (0.001 Å apart — clearly a coordinate-corruption
bug in some upstream tool). If prepare passes this through unchanged,
OpenMM's `Modeller.addHydrogens` in a subsequent step re-places HG at
OXT's position (an OpenMM bug in the CSER template path). The user
observed HG at 1.7 Å from its own OG after `zbs`'s prepare step.

The workaround inside prepare: detect coincident (H, heavy-atom) pairs
in the same protein residue before stripping H's, and strip BOTH.
PDBFixer's `addMissingAtoms` re-adds the terminal atom (OXT) at the
correct position based on the C-terminal backbone; `addHydrogens` then
places HG correctly since OXT is no longer at its old broken position
during template matching.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

# All heavy imports live in the prepare pipeline — gate on OpenMM presence.
pytest.importorskip("openmm", reason="prepare needs OpenMM")
pytest.importorskip("pdbfixer", reason="prepare needs PDBFixer")

from dvbfixer.prepare import main as prepare_main  # noqa: E402


def _atoms(pdb: Path, chain: str, resid: str) -> dict[str, tuple[float, float, float]]:
    """Return ``{atom_name: (x, y, z)}`` for the given (chain, resid) in ``pdb``."""
    out: dict[str, tuple[float, float, float]] = {}
    for ln in pdb.read_text().splitlines():
        if not ln.startswith(("ATOM  ", "HETATM")):
            continue
        if ln[21] != chain or ln[22:26].strip() != resid:
            continue
        name = ln[12:16].strip()
        x = float(ln[30:38])
        y = float(ln[38:46])
        z = float(ln[46:54])
        out[name] = (x, y, z)
    return out


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def test_prepare_fixes_coincident_hg_and_oxt(fixtures_root: Path, tmp_workdir: Path) -> None:
    """prepare on `broken_SER/SER.pdb` must place HG within 1.2 Å of OG (a
    reasonable upper bound; correct O-H bond is ~0.97 Å).
    """
    broken = fixtures_root / "broken_SER" / "SER.pdb"
    if not broken.exists():
        pytest.skip(f"fixture missing: {broken}")
    out = tmp_workdir / "fixed.pdb"
    prepare_main([str(broken), "-o", str(out)])

    atoms = _atoms(out, "B", "126")
    assert "HG" in atoms, f"HG missing from output: {list(atoms)}"
    assert "OG" in atoms and "OXT" in atoms and "C" in atoms

    d_og_hg = _dist(atoms["OG"], atoms["HG"])
    d_c_oxt = _dist(atoms["C"], atoms["OXT"])
    d_hg_oxt = _dist(atoms["HG"], atoms["OXT"])

    assert d_og_hg < 1.2, (
        f"HG-OG distance {d_og_hg:.3f} Å is too large (broken geometry). "
        f"Expected ~0.97 Å."
    )
    assert d_c_oxt < 1.5, (
        f"C-OXT distance {d_c_oxt:.3f} Å is too large. Expected ~1.25 Å."
    )
    assert d_hg_oxt > 1.0, (
        f"HG and OXT coincident ({d_hg_oxt:.3f} Å apart) — the bug is back."
    )
