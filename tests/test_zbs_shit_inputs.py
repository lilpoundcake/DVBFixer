"""Slow end-to-end regression: the raw PDBs in ``test/shit/`` broke the
zbs pipeline with ``openmm.OpenMMException: Particle coordinate is NaN``.

Root cause (July 2026):
    Modeller.addHydrogens occasionally places every H on a methyl/
    methylene/NH3+ parent at exactly the same coordinate. The old
    ``repair_misplaced_hydrogens`` then relocated all coincident H's
    on the same parent to the same "linear-anti" target — still
    coincident, just at a different point. LJ 1/r^12 → NaN inside
    ``simulation.minimizeEnergy``.

Fix: ``repair_misplaced_hydrogens`` now uses proper sp3 tetrahedral
placement (see ``_tetrahedral_h_positions`` in ``ffutils/geometry.py``)
so distinct H's on the same parent land at distinct positions.

Runs the full zbs pipeline on each historically-failing input and
asserts (a) it doesn't crash, (b) the output has no coincident-atom
pairs under 0.5 Å. Marked ``slow`` because each pipeline invocation
takes ~1-2 minutes (Modeller + OpenMM minimize).

Run explicitly:
    pytest tests/test_zbs_shit_inputs.py -v -m slow
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("openmm", reason="zbs needs OpenMM")
pytest.importorskip("pdbfixer", reason="zbs needs PDBFixer")
try:
    import modeller  # noqa: F401
except Exception as exc:
    pytest.skip(f"zbs needs a licensed Sali-lab MODELLER: {exc}", allow_module_level=True)


SHIT_INPUTS = ["1EMV_original.pdb", "1FR2_original.pdb",
               "2VLN_original.pdb", "2VLQ_original.pdb"]


@pytest.mark.slow
@pytest.mark.parametrize("input_name", SHIT_INPUTS)
def test_zbs_completes_on_shit_input(
    input_name: str, shit_dir: Path, tmp_workdir: Path,
) -> None:
    """Full zbs pipeline must complete and emit a coincident-atom-free PDB.

    Regression for the tetrahedral-placement fix in
    ``ffutils/geometry.repair_misplaced_hydrogens``. Before that fix,
    Modeller's occasional "3 H's on same coord" output was amplified
    by our repair pass into 400-500 coincident-atom pairs, and minimize
    then crashed with a NaN force at startup.

    Uses ``--no-solvent --strip-heterogens`` per the session preference
    — the failure mode is independent of solvation and heterogens.
    """
    src = shit_dir / input_name
    if not src.exists():
        pytest.skip(f"fixture missing: {src}")

    # zbs writes side-by-side to the input, so copy into tmp_workdir first.
    import shutil
    local = tmp_workdir / input_name
    shutil.copy2(src, local)

    from dvbfixer.zbs import main as zbs_main
    out = tmp_workdir / f"{local.stem}_zbs.pdb"
    zbs_main([str(local), "-o", str(out),
              "--no-solvent", "--strip-heterogens", "-v"])

    assert out.exists(), f"zbs produced no output for {input_name}"

    # Coincident-atom check: any pair under 0.5 Å blows up LJ in a
    # downstream minimize and is a bug we want to fail loudly on.
    import numpy as np
    from openmm.app import PDBFile
    from openmm.unit import angstrom
    from scipy.spatial import cKDTree

    pdb = PDBFile(str(out))
    pos = np.array([list(pdb.positions[i].value_in_unit(angstrom))
                    for i in range(len(pdb.positions))])
    tree = cKDTree(pos)
    pairs = tree.query_pairs(0.5)
    if pairs:
        atoms = list(pdb.topology.atoms())
        offenders = []
        for i, j in list(pairs)[:5]:
            a1, a2 = atoms[i], atoms[j]
            d = float(np.linalg.norm(pos[i] - pos[j]))
            offenders.append(
                f"{a1.residue.chain.id}/{a1.residue.name}{a1.residue.id}/"
                f"{a1.name} <-> {a2.residue.chain.id}/{a2.residue.name}"
                f"{a2.residue.id}/{a2.name}: d={d:.4f} Å"
            )
        pytest.fail(
            f"{len(pairs)} coincident atom pair(s) < 0.5 Å in {out.name}. "
            f"First: {offenders}"
        )

    # Chirality check: the initial batch of broken inputs regressed
    # into 40+ D-Cα residues per chain because the Modeller improper
    # restraint had the wrong atom order. Enforce ≤ 5% of protein
    # residues in D so a resurgence gets caught (2VLQ occasionally
    # ships 1-2 stochastic Modeller edge cases which we accept as
    # noise).
    from dvbfixer.ffutils.geometry import find_d_residues
    offenders = find_d_residues(pdb.topology, pdb.positions)
    n_protein = sum(1 for r in pdb.topology.residues()
                    if any(a.element and a.element.symbol == "C"
                           and a.name == "CA" for a in r.atoms()))
    max_allowed = max(2, int(0.05 * n_protein))
    if len(offenders) > max_allowed:
        pytest.fail(
            f"{len(offenders)} D-Cα residues in {out.name} "
            f"(threshold {max_allowed} = 5% of {n_protein} protein "
            f"residues). Offenders (first 10): "
            + ", ".join(f"{c}/{n}{r}" for c, r, n, _ in offenders[:10])
        )
