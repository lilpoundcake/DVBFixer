import json

import pytest
from scripts.prepare_protein_only_reconstruction import prepare, protein_only


def atom(serial, name="CA", residue="ALA", chain="A", number=1, record="ATOM  ", icode=""):
    return (f"{record}{serial:5d} {name:^4s} {residue:>3s} {chain}{number:4d}{icode:1s}   "
            f"{1.:8.3f}{2.:8.3f}{3.:8.3f}{1.:6.2f}{20.:6.2f}          {name[0]:>2s}  \n")


def test_remove_nonprotein_before_reconstruction_but_preserve_protein_hetatm_and_case():
    protein = atom(1, chain="D")+atom(2, residue="MSE", chain="d", record="HETATM")
    water = atom(3, name="O", residue="HOH", number=201, record="HETATM")
    text = "HEADER    RETAIN SOURCE IDENTITY\n"+protein+"TER\n"+water+"TER\nEND\n"
    result, report = protein_only(text)
    assert result == "HEADER    RETAIN SOURCE IDENTITY\n"+protein+"TER\nEND\n"
    assert report["retained_chains"] == ["D", "d"]
    assert report["retained_protein_hetatm_atoms"] == 1
    assert len(report["removed_atoms"]) == 1
    assert report["reattach_removed_context"] is False


def test_keep_matching_anisou_and_drop_removed_atom_tensor():
    protein, water = atom(1), atom(2, name="O", residue="HOH", number=201, record="HETATM")
    def tensor(line):
        return "ANISOU"+line[6:27]+" "+"  10000"*6+"      "+line[76:]
    text = protein+tensor(protein)+water+tensor(water)+"END\n"
    result, report = protein_only(text)
    assert result == protein+tensor(protein)+"END\n"
    assert report["discarded_record_counts"]["ANISOU"] == 1


def test_do_not_detach_explicit_protein_ligand_conect():
    text = atom(1)+atom(2, name="C1", residue="LIG", number=201, record="HETATM")+"CONECT    1    2\nEND\n"
    with pytest.raises(ValueError, match="chemistry decision"):
        protein_only(text)


def test_preserve_protein_conect_but_drop_complete_dangling_records():
    text = atom(1)+atom(2, name="N", number=2)+"CONECT    1    2\nCONECT    1  999\nEND\n"
    result, report = protein_only(text)
    assert "CONECT    1    2\n" in result
    assert "999" not in result
    assert report["discarded_record_counts"]["CONECT"] == 1


def test_modified_unknown_atom_is_not_silently_deleted():
    with pytest.raises(ValueError, match="unclassified ATOM"):
        protein_only(atom(1, residue="UNK"))


def test_multimodel_and_serial_aliases_require_explicit_selection():
    with pytest.raises(ValueError, match="select one model"):
        protein_only("MODEL        1\n"+atom(1)+"ENDMDL\n")
    with pytest.raises(ValueError, match="duplicate atom serial"):
        protein_only(atom(1)+atom(1, name="N"))


def test_exclusive_preparation_preserves_source_and_requires_new_request(tmp_path):
    source = tmp_path / "source.pdb"
    original = atom(1)+atom(2, name="O", residue="HOH", number=201, record="HETATM")+"END\n"
    source.write_text(original)
    output = tmp_path / "prepared"
    prepare(source, output)
    report = json.loads((output / "preparation.json").read_text())
    assert source.read_text() == (output / "original.pdb").read_text() == original
    assert (output / "protein-only.pdb").read_text() == atom(1)+"END\n"
    assert report["fresh_request_required"] and not report["historical_requests_modified"]
    with pytest.raises(FileExistsError):
        prepare(source, output)
