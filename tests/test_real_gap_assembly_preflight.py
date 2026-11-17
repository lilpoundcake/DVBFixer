from types import SimpleNamespace

import pytest
from scripts import preflight_real_gap_assemblies as preflight
from scripts.mac_reconstruction_comparison import digest


def test_source_mutation_is_retained_refusal_not_an_adopted_input(tmp_path):
    source = tmp_path / "source.cif.gz"
    source.write_bytes(b"original public bytes")
    row = {"entry_id": "TEST", "source_path": str(source), "source_sha256": digest(source)}
    source.write_bytes(b"changed public bytes")
    result = preflight.inspect(row, tmp_path / "work")
    assert result["status"] == "assembly_inventory_refused"
    assert result["stage"] == "normalization"
    assert "changed" in result["error"]
    assert (tmp_path / "work/failure.txt").is_file()
    with pytest.raises(FileExistsError):
        preflight.inspect(row, tmp_path / "work")


def test_explicit_metal_link_is_not_automatically_called_covalent(tmp_path, monkeypatch):
    source = tmp_path / "source.cif.gz"
    source.write_bytes(b"retained original")
    lines = []
    for serial, name, xyz in ((1, "N", (0, 0, 0)), (2, "CA", (1, 0, 0)),
                              (3, "C", (2, 0, 0)), (4, "O", (3, 0, 0)), (5, "CB", (1, 1, 0))):
        element = name[0]
        lines.append(f"ATOM  {serial:5d} {name:^4} ALA a   1    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00          {element:>2}\n")
    lines.append("HETATM    6 MG    MG a 201       0.000   0.000   2.000  1.00 20.00          MG\n")
    lines.append("CONECT    1    6\n")

    def normalize(original, destination):
        destination.write_text("".join(lines))
        return SimpleNamespace(pdb_path=destination, chain_map={"a": "a"})

    monkeypatch.setattr(preflight, "normalize_structure", normalize)
    group = SimpleNamespace(chains=["a"], operators=lambda: [None])
    monkeypatch.setattr(preflight, "parse_biological_assemblies", lambda _: {"1": SimpleNamespace(groups=[group])})
    monkeypatch.setattr(preflight, "render_biological_assembly", lambda *args, **kwargs: lines)
    row = {"entry_id": "TEST", "source_path": str(source), "source_sha256": digest(source),
           "inventory": {"chains": [{"observed_author_axis": [{"author_chain": "a", "author_residue_number": "1", "insertion_code": ""}]}]}}
    result = preflight.inspect(row, tmp_path / "work")
    assert result["status"] == "assembly_inventory_complete", result
    assert result["protein_nonprotein_explicit_links"] == 1
    assert result["protein_nonprotein_link_details"][0]["chemical_bond_classification"].startswith("not established")
    assert any("coordination review" in reason for reason in result["strict_protein_route_obstacles"])
    assert result["all_original_context_still_retained"]
    assert not result["scientific_admission_complete"]
