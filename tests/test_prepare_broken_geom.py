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


# Synthetic reproducer for the ACTUAL failure mode reported at 0.4.1:
# input has OXT PRESENT but HG MISSING on a C-terminal SER. OpenMM's
# addHydrogens then places the new HG on top of OXT. The 0.4.1 pre-strip
# couldn't catch this because the input has no coincident pair. Fixed
# in 0.4.2 by the post-addHydrogens repair pass in ffutils.geometry.
#
# Constructed in-test so the test doesn't depend on the `test/`
# fixtures directory (which is git-ignored).
_SER_INPUT_HG_MISSING_OXT_PRESENT = """\
ATOM   4731  N   SER B 125     -17.531  20.605 -24.080  1.00  0.00           N
ATOM   4732  H   SER B 125     -16.866  20.377 -24.816  1.00  0.00           H
ATOM   4733  CA  SER B 125     -18.451  21.721 -24.269  1.00  0.00           C
ATOM   4734  HA  SER B 125     -19.301  21.599 -23.599  1.00  0.00           H
ATOM   4735  C   SER B 125     -18.969  21.755 -25.702  1.00  0.00           C
ATOM   4736  O   SER B 125     -18.179  21.679 -26.655  1.00  0.00           O
ATOM   4737  CB  SER B 125     -17.748  23.028 -23.912  1.00  0.00           C
ATOM   4738  HB2 SER B 125     -17.280  22.939 -22.932  1.00  0.00           H
ATOM   4739  HB3 SER B 125     -16.994  23.237 -24.670  1.00  0.00           H
ATOM   4740  OG  SER B 125     -18.664  24.103 -23.849  1.00  0.00           O
ATOM   4741  HG  SER B 125     -18.389  24.725 -23.181  1.00  0.00           H
ATOM   4742  N   SER B 126     -20.246  21.928 -25.860  1.00  0.00           N
ATOM   4743  H   SER B 126     -20.762  22.238 -25.051  1.00  0.00           H
ATOM   4744  CA  SER B 126     -20.916  22.100 -27.152  1.00  0.00           C
ATOM   4745  HA  SER B 126     -20.275  21.751 -27.937  1.00  0.00           H
ATOM   4746  C   SER B 126     -21.114  23.596 -27.370  1.00  0.00           C
ATOM   4747  O   SER B 126     -20.149  24.343 -27.634  1.00  0.00           O
ATOM   4748  CB  SER B 126     -22.268  21.363 -27.174  1.00  0.00           C
ATOM   4749  HB2 SER B 126     -22.692  21.265 -28.187  1.00  0.00           H
ATOM   4750  HB3 SER B 126     -22.069  20.334 -26.888  1.00  0.00           H
ATOM   4751  OG  SER B 126     -23.160  21.914 -26.218  1.00  0.00           O
ATOM   4753  OXT SER B 126     -22.355  23.104 -27.085  1.00  0.00           O
TER
END
"""


def test_prepare_repairs_missing_hg_when_oxt_present(tmp_workdir: Path) -> None:
    """Actual reported failure mode: HG missing, OXT present.

    OpenMM's Modeller.addHydrogens places the new HG right on top of the
    existing OXT (a CSER template misbehaviour). The 0.4.2 post-addHydrogens
    repair pass in dvbfixer.ffutils.geometry detects the coincident /
    over-long H-parent distance and re-places HG at ~0.97 Å from OG in
    the linear-anti direction from CB.
    """
    in_pdb = tmp_workdir / "hg_missing.pdb"
    in_pdb.write_text(_SER_INPUT_HG_MISSING_OXT_PRESENT)
    out = tmp_workdir / "repaired.pdb"
    prepare_main([str(in_pdb), "-o", str(out)])

    atoms = _atoms(out, "B", "126")
    assert "HG" in atoms, "HG missing from repaired output"
    assert "OG" in atoms and "OXT" in atoms

    d_og_hg = _dist(atoms["OG"], atoms["HG"])
    d_hg_oxt = _dist(atoms["HG"], atoms["OXT"])

    assert d_og_hg < 1.2, (
        f"HG-OG distance {d_og_hg:.3f} Å is too large; post-check "
        f"didn't repair the misplaced hydrogen. Expected ~0.97 Å."
    )
    assert d_hg_oxt > 1.0, (
        f"HG still on top of OXT ({d_hg_oxt:.3f} Å) — post-check "
        f"didn't fire. The OpenMM CSER-template bug is back."
    )
