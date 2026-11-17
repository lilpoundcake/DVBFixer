import pytest
from scripts.group_real_gap_candidates import components
from scripts.real_gap_source_inventory import inventory_mmcif, missing_intervals


def source():
    header = """data_TEST
_entry.id TEST
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code_can
1 'polypeptide(L)' ACDEFG
loop_
_struct_asym.id
_struct_asym.entity_id
a 1
A 1
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.occupancy
_atom_site.pdbx_PDB_model_num
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
"""
    atoms = [f"ATOM C a 1 {i} d {81+i} {'A' if i == 2 else '.'} {0.5 if i == 5 else 1} 1 {i} 0 0"
             for i in (1, 2, 5, 6)]
    atoms += ["ATOM C A 1 1 D 1 . 0 1 0 0 0", "ATOM C a 1 3 d 84 . 1 2 0 0 0"]
    return header + "\n".join(atoms) + "\n"


def metadata():
    return {"entry_id": "TEST", "chains": [
        {"label_chain": "a", "sequence": "ACDEFG", "gaps": missing_intervals(6, {1, 2, 5, 6})},
        {"label_chain": "A", "sequence": "ACDEFG", "gaps": missing_intervals(6, {1})}]}


def test_label_axis_gaps_not_greedy_alignment_and_case_icode_are_preserved():
    result = inventory_mmcif(source(), metadata())
    chains = {c["label_chain"]: c for c in result["chains"]}
    assert set(chains) == {"a", "A"}
    assert chains["a"]["natural_absent_coordinate_spans"] == [{"start": 2, "stop": 4, "length": 2, "kind": "internal"}]
    assert chains["a"]["observed_author_axis"][1] == {"label_seq_id": 2, "author_chain": "d",
        "author_residue_number": "83", "insertion_code": "A"}
    assert chains["a"]["residues_with_partial_occupancy"] == [5]
    assert chains["A"]["observed_residues_any_heavy_record"] == 1
    assert chains["A"]["observed_residues_positive_occupancy"] == 0
    assert result["multi_model_source_needs_explicit_policy"]
    assert not result["source_coordinate_changes"]
    assert not result["scientific_admission_complete"]


def test_source_change_is_not_silently_reconciled_to_metadata():
    changed = metadata()
    changed["chains"][0]["sequence"] = "AAAAAA"
    with pytest.raises(ValueError, match="sequence changed"):
        inventory_mmcif(source(), changed)
    with pytest.raises(ValueError, match="outside"):
        missing_intervals(6, {0, 7})


def test_shared_cluster_components_are_transitive_and_not_split_by_entry():
    rows = [{"entry_id": "one", "leakage_keys": ["cluster:a"]},
            {"entry_id": "two", "leakage_keys": ["cluster:a", "cluster:b"]},
            {"entry_id": "three", "leakage_keys": ["cluster:b"]},
            {"entry_id": "four", "leakage_keys": ["cluster:c"]}]
    assert sorted(components(rows)) == [["four"], ["one", "three", "two"]]


def test_whole_unobserved_chain_and_long_tail_are_not_truncated_or_conflated():
    assert missing_intervals(500, set()) == [{"start": 0, "stop": 500, "length": 500, "kind": "whole-chain"}]
    assert missing_intervals(500, set(range(1, 301))) == [{"start": 300, "stop": 500, "length": 200, "kind": "C-terminal"}]
