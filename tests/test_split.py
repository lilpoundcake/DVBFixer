"""Regression tests for `dvbfixer split`.

Uses the tracked `tests/fixtures/multistate.pdb` fixture —
it's a multi-MODEL PDB with 3 chains and no chain IDs, so it exercises
the multi-MODEL detection path.
"""

from __future__ import annotations

from pathlib import Path

from dvbfixer.molecule_chains import assign_unique_molecule_chains
from dvbfixer.split_chains import CHAIN_IDS, main


def _chain_ids(pdb: Path) -> set[str]:
    ids = set()
    for ln in pdb.read_text().splitlines():
        if ln.startswith(("ATOM  ", "HETATM")):
            cid = ln[21]
            if cid != " ":
                ids.add(cid)
    return ids


def _count_models(pdb: Path) -> int:
    return sum(1 for ln in pdb.read_text().splitlines() if ln.startswith("MODEL "))


def _count_ter(pdb: Path) -> int:
    return sum(1 for ln in pdb.read_text().splitlines() if ln.startswith("TER"))


def test_split_multistate_assigns_chain_ids(tmp_workdir: Path, multistate_pdb: Path) -> None:
    assert multistate_pdb.is_file(), f"tracked fixture missing: {multistate_pdb}"
    out = tmp_workdir / "split.pdb"
    main([str(multistate_pdb), "-o", str(out)])
    assert out.exists()
    ids = _chain_ids(out)
    # multistate PDB has 3 chains; after split each should have a letter ID.
    assert len(ids) >= 2, f"expected multiple chains, got {ids}"


def test_split_multistate_preserves_model_count(tmp_workdir: Path, multistate_pdb: Path) -> None:
    assert multistate_pdb.is_file(), f"tracked fixture missing: {multistate_pdb}"
    out = tmp_workdir / "split.pdb"
    main([str(multistate_pdb), "-o", str(out), "--unique-molecule-chains"])
    # multistate fixture is 11 MODELs.
    in_models = _count_models(multistate_pdb)
    out_models = _count_models(out)
    assert out_models == in_models, f"MODEL count changed: {in_models} -> {out_models}"


def test_split_inserts_ter_records(tmp_workdir: Path, multistate_pdb: Path) -> None:
    assert multistate_pdb.is_file(), f"tracked fixture missing: {multistate_pdb}"
    out = tmp_workdir / "split.pdb"
    main([str(multistate_pdb), "-o", str(out)])
    # Every chain boundary should get a TER; multi-MODEL multiplies.
    assert _count_ter(out) > 0


def test_unique_molecule_chains_reserve_assembly_ids_and_preserve_conect() -> None:
    lines = [
        "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, H, L\n",
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
        "HETATM   10  C1  LIG D 401       1.000   0.000   0.000  1.00  0.00           C  \n",
        "HETATM   11  C1  DRG D 402       2.000   0.000   0.000  1.00  0.00           C  \n",
        "CONECT   10   10\n",
        "END\n",
    ]
    output, count = assign_unique_molecule_chains(lines, CHAIN_IDS)
    assert count == 2
    heterogen_chains = [line[21] for line in output if line.startswith("HETATM")]
    assert heterogen_chains == ["B", "C"]
    assert not ({"H", "L"} & set(heterogen_chains))
    assert [line for line in output if line.startswith("CONECT")] == ["CONECT   10   10\n"]
