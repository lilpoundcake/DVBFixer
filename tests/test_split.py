"""Regression tests for `dvbfixer split`.

Uses the tracked `tests/fixtures/multistate.pdb` fixture —
it's a multi-MODEL PDB with 3 chains and no chain IDs, so it exercises
the multi-MODEL detection path.
"""

from __future__ import annotations

from pathlib import Path

import pytest

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


@pytest.mark.parametrize("models", [1, 2])
@pytest.mark.parametrize("unique", [False, True])
def test_keep_heterogens_preserves_all_atoms_in_each_model(tmp_path: Path, models: int, unique: bool) -> None:
    atoms = [
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
        "HETATM   10  C1  LIG A 401       1.000   0.000   0.000  1.00  0.00           C  \n",
        "HETATM   11  O   HOH A 402       2.000   0.000   0.000  1.00  0.00           O  \n",
        "HETATM   12 ZN    ZN A 403       3.000   0.000   0.000  1.00  0.00          ZN  \n",
        "HETATM   13  C1  BUF A 404       4.000   0.000   0.000  1.00  0.00           C  \n",
    ]
    source = tmp_path / "input.pdb"
    source.write_text("".join(
        f"MODEL     {model:4d}\n" + "".join(atoms) + "ENDMDL\n"
        for model in range(1, models + 1)
    ) + "CONECT    1   10\nEND\n")
    output = tmp_path / "output.pdb"
    main([str(source), "-o", str(output), "--keep-heterogens"] +
         (["--unique-molecule-chains"] if unique else []))
    blocks = output.read_text().split("ENDMDL")[:models]
    for block in blocks:
        retained = [line for line in block.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
        assert len(retained) == len(atoms)
        # Only chain/residue identifiers may change; serials, names and
        # coordinates must survive, including solvent in every MODEL.
        assert [(line[:21], line[27:]) for line in retained] == [
            (line[:21], line.rstrip("\n")[27:]) for line in atoms
        ]
    assert "CONECT    1   10\n" in output.read_text()


@pytest.mark.parametrize("keep", [False, True])
def test_keep_heterogens_preserves_real_8ucd_molecules(tmp_path: Path, keep: bool) -> None:
    source = Path(__file__).parent / "fixtures" / "rename_mols" / "structure.pdb"
    output = tmp_path / "split.pdb"
    main([str(source), "-o", str(output), "--unique-molecule-chains"] +
         (["--keep-heterogens"] if keep else []))
    before = source.read_text().splitlines()
    after = output.read_text().splitlines()
    original = [line for line in before if line.startswith("HETATM")]
    retained = [line for line in after if line.startswith("HETATM")]
    assert len(original) == len(retained) == 807
    assert [(line[:21], line[27:]) for line in original] == [
        (line[:21], line[27:]) for line in retained
    ]
    assert len({line[21] for line in retained}) == 9
    assert [line for line in before if line.startswith("CONECT")] == [
        line for line in after if line.startswith("CONECT")
    ]
    # Compare atom identities, not just the unchanged CONECT text. The latter
    # can conceal corruption if coordinate serials were changed independently.
    def bond_identities(lines):
        atoms = {int(line[6:11]): (line[12:20], line[30:54])
                 for line in lines if line.startswith(("ATOM  ", "HETATM"))}
        bonds = []
        for line in lines:
            if not line.startswith("CONECT"):
                continue
            serials = [int(line[i:i + 5]) for i in range(6, len(line.rstrip()), 5)
                       if line[i:i + 5].strip()]
            bonds.extend((atoms[serials[0]], atoms[partner]) for partner in serials[1:])
        return bonds
    assert bond_identities(before) == bond_identities(after)

    # Exercise Model's actual bond restoration after its output serialization.
    from dvbfixer.model.pipeline import _renumber_atom_serials, restore_conect_records
    model_atoms = _renumber_atom_serials([
        line + "\n" for line in after if not line.startswith("CONECT")
    ])
    restored = restore_conect_records(model_atoms, [line + "\n" for line in after])
    assert bond_identities(restored) == bond_identities(before)
    helices_before = [line for line in before if line.startswith("HELIX ")]
    helices_after = [line for line in after if line.startswith("HELIX ")]
    assert len(helices_before) == len(helices_after) > 0
    atoms_before = [line for line in before if line.startswith("ATOM  ")]
    atoms_after = [line for line in after if line.startswith("ATOM  ")]
    residue_map = {
        (left[21], int(left[22:26]), left[26]): (right[21], int(right[22:26]), right[26])
        for left, right in zip(atoms_before, atoms_after)
    }
    for left, right in zip(helices_before, helices_after):
        for chain, start, end in [(19, 21, 25), (31, 33, 37)]:
            assert (right[chain], int(right[start:end]), right[end]) == residue_map[
                (left[chain], int(left[start:end]), left[end])
            ]


@pytest.mark.parametrize("models", [0, 2])
def test_keep_heterogens_accepts_solvent_only_input(tmp_path: Path, models: int) -> None:
    water = "HETATM   42  O   HOH A 402       2.000   0.000   0.000  1.00  0.00           O  \n"
    body = "".join(f"MODEL     {n:4d}\n{water}ENDMDL\n" for n in range(1, models + 1)) if models else water
    source = tmp_path / "water.pdb"
    source.write_text(body + "END\n")
    output = tmp_path / "retained.pdb"
    main([str(source), "-o", str(output), "--keep-heterogens", "--unique-molecule-chains"])
    assert output.read_text() == source.read_text()


def test_secondary_structure_drops_annotation_across_new_chain_boundary(capsys) -> None:
    from dvbfixer.split_chains import _remap_secondary_structure

    atoms = [
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
        "ATOM      2  CA  ALA A   2      30.000   0.000   0.000  1.00  0.00           C  \n",
    ]
    helix = list(" " * 80 + "\n")
    helix[:6] = "HELIX "
    for chain, start, end, number in [(19, 21, 25, 1), (31, 33, 37, 2)]:
        helix[chain] = "A"
        helix[start:end] = f"{number:4d}"
    header = "".join(helix)
    split_atoms = [atoms[0], atoms[1][:21] + "B" + atoms[1][22:]]
    result = _remap_secondary_structure([header, *atoms], [header, *split_atoms])
    assert result == split_atoms
    assert "omitted 1 HELIX/SHEET" in capsys.readouterr().out
