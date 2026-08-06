"""Regression tests for `dvbfixer split`.

Uses the tracked `tests/fixtures/multistate.pdb` fixture —
it's a multi-MODEL PDB with 3 chains and no chain IDs, so it exercises
the multi-MODEL detection path.
"""

from __future__ import annotations

from pathlib import Path

from dvbfixer.split_chains import main


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
    main([str(multistate_pdb), "-o", str(out)])
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
