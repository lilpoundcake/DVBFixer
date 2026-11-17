"""Reducing-end oxygen transport through GLYCAM conversion and PDB writing."""

from io import StringIO
from pathlib import Path

import pytest

from dvbfixer.domain.force_field_naming import ForceFieldTarget
from dvbfixer.force_field_naming import NamingConversionRequest, convert_force_field_naming
from dvbfixer.glycam import _format_atom_line, convert_to_glycam


def _edges(text: str) -> set[tuple[int, int]]:
    edges = set()
    for line in text.splitlines():
        if line.startswith("CONECT"):
            serials = [int(line[index:index + 5]) for index in range(6, len(line), 5)]
            edges.update((serials[0], neighbor) for neighbor in serials[1:])
    return edges


def _source(path: Path, *, chain: str = "A", number: int = 1) -> str:
    atoms = [
        ("C1", "BGC", chain, number, 0.0),
        ("O1", "BGC", chain, number, 1.4),
        ("C2", "BGC", chain, number, -1.4),
        ("CA", "ALA", "B", 1, 20.0),
        ("CA", "ALA", "C", 219, 40.0),
    ]
    text = "REMARK synthetic reducing-end oxygen regression\n"
    for serial, (name, residue, atom_chain, resid, x) in enumerate(atoms, 1):
        text += _format_atom_line({
            "name": name, "chain": atom_chain, "resseq": resid,
            "x": x, "y": 0.0, "z": 0.0, "element": name[0],
        }, residue, serial)
    text += "TER\nCONECT    1    2    3\nCONECT    2    1\nEND\n"
    path.write_text(text)
    return text


def test_roh_openmm_roundtrip_naming(tmp_path: Path) -> None:
    pdbfile = pytest.importorskip("openmm.app").PDBFile
    source, converted = tmp_path / "input.pdb", tmp_path / "glycam.pdb"
    _source(source)
    convert_to_glycam(source, converted)
    # The converter has only one final TER. OpenMM introduces TER boundaries
    # between chains when writing the topology, as prepare/minimize do.
    assert sum(line.startswith("TER") for line in converted.read_text().splitlines()) == 1
    parsed = pdbfile(str(converted))
    output = StringIO()
    pdbfile.writeFile(parsed.topology, parsed.positions, output, keepIds=True)
    convert_force_field_naming(NamingConversionRequest(output.getvalue(), ForceFieldTarget.AMBER))
    assert [chain.id for chain in parsed.topology.chains()] == ["A", "B", "C"]
    residues = list(next(parsed.topology.chains()).residues())
    assert [(residue.name, residue.id) for residue in residues] == [("0GB", "1"), ("ROH", "2")]
    assert any(
        {left.residue.name, right.residue.name} == {"0GB", "ROH"}
        and {left.name, right.name} == {"C1", "O1"}
        for left, right in parsed.topology.bonds()
    )


@pytest.mark.parametrize("chain", ["A", "a"])
def test_roh_is_separate_local_residue_and_preserves_atoms_and_bonds(tmp_path: Path, chain: str) -> None:
    source, converted, repeated = (tmp_path / name for name in ("input.pdb", "glycam.pdb", "again.pdb"))
    original = _source(source, chain=chain)
    convert_to_glycam(source, converted)
    text = converted.read_text()
    records = [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    assert len(records) == 5  # O1 is moved, not added or deleted.
    assert [line[21] for line in records] == [chain, chain, chain, "B", "C"]
    cap = next(line for line in records if line[17:20] == "ROH")
    sugar = next(line for line in records if line[12:16].strip() == "C1")
    assert cap[21:27] != sugar[21:27]
    old_oxygen = original.splitlines()[2]
    assert cap[30:54] == old_oxygen[30:54]
    assert _edges(text) == {(1, 2), (1, 3), (2, 1), (3, 1)}
    assert (int(sugar[6:11]), int(cap[6:11])) in _edges(text)
    assert source.read_text() == original
    convert_to_glycam(converted, repeated)
    # CONECT record/neighbor order can change after serial remapping, but
    # neither coordinates, identities nor connectivity may change.
    assert _edges(repeated.read_text()) == _edges(text)
    assert [line for line in repeated.read_text().splitlines() if not line.startswith("CONECT")] == [
        line for line in text.splitlines() if not line.startswith("CONECT")
    ]


def test_no_roh_keeps_oxygen_in_sugar(tmp_path: Path) -> None:
    source, converted = tmp_path / "input.pdb", tmp_path / "glycam.pdb"
    _source(source)
    convert_to_glycam(source, converted, add_roh=False)
    text = converted.read_text()
    assert "ROH" not in text
    oxygen = next(line for line in text.splitlines() if line[12:16].strip() == "O1")
    assert oxygen[17:20] == "0GB"


def test_multiple_roots_reserve_numbers_against_all_chain_residues(tmp_path: Path) -> None:
    source, converted = tmp_path / "input.pdb", tmp_path / "glycam.pdb"
    text = _source(source)
    extra = ""
    for serial, name, residue, number, x in [
        (6, "C1", "BGC", 2, 100.0), (7, "O1", "BGC", 2, 101.4),
        (8, "O", "HOH", 3, 200.0),
    ]:
        extra += _format_atom_line({
            "name": name, "chain": "A", "resseq": number,
            "x": x, "y": 0.0, "z": 0.0, "element": name[0],
        }, residue, serial)
    # Keep input chains contiguous, with an existing heterogen after both sugars.
    lines = text.splitlines(keepends=True)
    source.write_text("".join(lines[:4]) + extra + "".join(lines[4:]))
    convert_to_glycam(source, converted)
    caps = [line for line in converted.read_text().splitlines() if line[17:20] == "ROH"]
    assert [int(line[22:26]) for line in caps] == [4, 5]


def test_roh_number_stays_within_pdb_columns(tmp_path: Path) -> None:
    source, converted = tmp_path / "input.pdb", tmp_path / "glycam.pdb"
    _source(source, number=9999)
    convert_to_glycam(source, converted)
    cap = next(line for line in converted.read_text().splitlines() if line[17:20] == "ROH")
    assert int(cap[22:26]) == 1
    assert cap[26] == " "
    assert float(cap[30:38]) == 1.4


def test_roh_does_not_inherit_sugar_insertion_code(tmp_path: Path) -> None:
    source, converted = tmp_path / "input.pdb", tmp_path / "glycam.pdb"
    text = _source(source)
    source.write_text("".join(
        line[:26] + "Z" + line[27:] if line.startswith("HETATM") else line
        for line in text.splitlines(keepends=True)
    ))
    convert_to_glycam(source, converted)
    records = [line for line in converted.read_text().splitlines() if line.startswith("HETATM")]
    assert [(line[17:20], line[22:27]) for line in records] == [
        ("0GB", "   1Z"), ("0GB", "   1Z"), ("ROH", "   2 "),
    ]
