"""Regression test for a real HETATM/protein resSeq-collision bug found
on `test/lipid/7x35_r_u.pdb`.

Chain A of this 3-chain viral capsid has no SEQRES records and only 267
ATOM residues, but its true (FASTA) sequence is 278 residues — an
11-residue gap that `model` fills in. The bound ligand `PLM` (palmitic
acid) got naively renumbered to chain A resSeq 268 by the FASTA-blind
standalone `renumber.py` step (one past the ATOM-only count) — which
collided exactly with the first gap-filled protein residue once `model`
expanded the chain to its full 278-residue length. `PLM`'s atom-naming
convention (`C1..C9, CA, CB, CC...`) then let `minimize`'s legacy
strip-and-splice position-restore merge (keyed by `(chain, resid,
atomname)`, no resname check) silently overwrite `PLM`'s `CA`/`CB`
atoms with the colliding protein residue's minimized coordinates —
"strange links between lipid atoms" (a bond to an atom now tens of
angstroms away) in any viewer.

Fixed in two places: the resSeq collision itself
(`model.renumber.build_resnum_mapping`, unit-tested directly in
`test_model_renumber.py`), and defense-in-depth in `minimize.pipeline`'s
restore-merge (now keyed by `(chain, resid, parent_resname, atomname)`).
This test exercises the real, full pipeline end to end.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("openmm", reason="zbs needs OpenMM")
pytest.importorskip("pdbfixer", reason="zbs needs PDBFixer")
pytest.importorskip("modeller", reason="zbs needs Sali-lab MODELLER")


@pytest.mark.slow
def test_zbs_no_lipid_protein_resseq_collision(lipid_dir: Path, tmp_workdir: Path) -> None:
    src = lipid_dir / "7x35_r_u.pdb"
    fasta = lipid_dir / "7x35_renamed.fasta"
    if not src.exists() or not fasta.exists():
        pytest.skip(f"fixture missing: {src} / {fasta}")

    import shutil
    local = tmp_workdir / src.name
    shutil.copy2(src, local)

    from dvbfixer.zbs import main as zbs_main
    out = tmp_workdir / f"{local.stem}_zbs.pdb"
    # NOTE: no --strip-heterogens here. The real bug reproduces on the
    # default keep-heterogens path (PLM goes through H-addition and
    # minimize's restraint/restore path, same as the user's original
    # report) — --strip-heterogens takes a different code path in
    # `prepare` that removes the ligand outright before minimize ever
    # runs, which doesn't exercise the collision bug at all.
    zbs_main([
        str(local), "-o", str(out), "--fasta", str(fasta),
        "--atom-naming", "standard", "--no-solvent", "-v",
    ])

    assert out.exists(), "zbs produced no output for the lipid-complex fixture"

    protein_resseqs_by_chain: dict[str, set[int]] = {}
    plm_atoms: dict[str, tuple[float, float, float]] = {}
    plm_chain = None
    for line in out.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        resname = line[17:20].strip()
        chain = line[21]
        resseq = int(line[22:26].strip())
        if resname == "PLM":
            plm_chain = chain
            plm_atoms[line[12:16].strip()] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54]),
            )
            continue
        protein_resseqs_by_chain.setdefault(chain, set()).add(resseq)

    assert plm_atoms, "PLM ligand missing from zbs output entirely"
    assert plm_chain is not None

    # The ligand's resSeq must not collide with ANY protein residue in its
    # own chain — the exact bug this test guards against.
    plm_resseq_lines = [
        int(ln[22:26].strip()) for ln in out.read_text().splitlines()
        if ln.startswith(("ATOM  ", "HETATM")) and ln[17:20].strip() == "PLM"
    ]
    assert len(set(plm_resseq_lines)) == 1
    plm_resseq = plm_resseq_lines[0]
    assert plm_resseq not in protein_resseqs_by_chain.get(plm_chain, set()), (
        f"PLM (chain {plm_chain}) resSeq {plm_resseq} collides with a "
        f"protein residue in the same chain"
    )

    # The ligand's own alkyl chain must be internally self-consistent —
    # consecutive-numbered carbons (C1..C9) should be within a normal
    # covalent-bond-ish distance of each other, not teleported across the
    # structure by a corrupted coordinate-restore merge.
    def _dist(a, b):
        return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5

    chain_atom_order = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
    present = [n for n in chain_atom_order if n in plm_atoms]
    assert len(present) >= 5, f"too few PLM chain atoms found: {list(plm_atoms)}"
    for a, b in zip(present, present[1:]):
        d = _dist(plm_atoms[a], plm_atoms[b])
        assert d < 3.0, (
            f"PLM {a}-{b} distance is {d:.1f} Å — the ligand's own "
            f"chain is no longer physically contiguous (coordinate "
            f"cross-contamination)"
        )
