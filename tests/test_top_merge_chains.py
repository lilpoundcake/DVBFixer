"""Regression test for the 2026-07-31 audit's `top --merge` fix.

`top.pipeline._merge_chains` built each merged `AtomEntry` without
`x`/`y`/`z`/`chain_id`/`orig_resseq`/`orig_resname` — those dataclass
fields silently fell back to their defaults (0.0, 0.0, 0.0, ' ', 0, ''),
so `--merge --pdb out.pdb` produced a PDB with every atom at the
origin. `resnr` was also copied unmodified per chain (not offset),
colliding residue numbers across originally-distinct chains.

Tests directly against `_merge_chains` (no FF/RTP setup needed) rather
than a full `dvbfixer top` CLI run, since the bug — and its fix — lives
entirely in this one function's `AtomEntry` construction.
"""
from __future__ import annotations

from dvbfixer.top.pipeline import AtomEntry, ChainTopology, _merge_chains


def _make_chain(name: str, n_atoms: int, chain_id: str) -> ChainTopology:
    ct = ChainTopology(name=name, nrexcl=3)
    for i in range(1, n_atoms + 1):
        ct.atoms.append(AtomEntry(
            index=i, atom_type="CT", resnr=i, resname="ALA", atomname="CA",
            cgnr=i, charge=0.0, mass=12.0,
            x=float(i), y=float(i) * 2, z=float(i) * 3,
            chain_id=chain_id, orig_resseq=i, orig_resname="ALA",
        ))
    ct.bonds = [(i, i + 1) for i in range(1, n_atoms)]
    return ct


def test_merge_chains_preserves_coordinates_and_identity() -> None:
    chain_a = _make_chain("A", 3, "A")
    chain_b = _make_chain("B", 2, "B")

    merged = _merge_chains([chain_a, chain_b])

    assert len(merged.atoms) == 5
    # No atom should have collapsed to the all-zero default.
    for atom in merged.atoms:
        assert (atom.x, atom.y, atom.z) != (0.0, 0.0, 0.0)
        assert atom.chain_id != ' '
        assert atom.orig_resname == "ALA"

    # Chain A's atoms keep their original coordinates/chain_id/orig_resseq.
    assert [a.x for a in merged.atoms[:3]] == [1.0, 2.0, 3.0]
    assert [a.chain_id for a in merged.atoms[:3]] == ["A", "A", "A"]
    assert [a.orig_resseq for a in merged.atoms[:3]] == [1, 2, 3]

    # Chain B's atoms keep THEIR original coordinates/chain_id too.
    assert [a.x for a in merged.atoms[3:]] == [1.0, 2.0]
    assert [a.chain_id for a in merged.atoms[3:]] == ["B", "B"]

    # Bond indices offset correctly (chain A: 1-2, 2-3; chain B's bond
    # 1-2 becomes 4-5 after chain A's 3 atoms).
    assert merged.bonds == [(1, 2), (2, 3), (4, 5)]


def test_merge_chains_offsets_resnr_across_chains() -> None:
    """resnr must not collide across originally-distinct chains — chain
    B's residue 1 must not still read as residue 1 once merged behind
    chain A's 3 residues."""
    chain_a = _make_chain("A", 3, "A")
    chain_b = _make_chain("B", 2, "B")

    merged = _merge_chains([chain_a, chain_b])

    resnrs = [a.resnr for a in merged.atoms]
    assert resnrs == sorted(resnrs), f"resnr must be monotonic, got {resnrs}"
    assert len(set(resnrs)) == len(resnrs), (
        f"resnr collided across chains: {resnrs}"
    )
