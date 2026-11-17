from copy import deepcopy

import gemmi
import pytest
from scripts.review_retained_construct_targets import target_facts


def fixture():
    block = gemmi.cif.Block("TEST")
    block.set_mmcif_category("_entity_poly.", {"entity_id": ["1"], "type": ["polypeptide(L)"],
        "pdbx_seq_one_letter_code_can": ["AAAA"]})
    block.set_mmcif_category("_atom_site.", {"label_asym_id": ["L", "L"], "label_entity_id": ["1", "1"],
        "label_seq_id": ["1", "3"], "type_symbol": ["C", "C"], "pdbx_PDB_model_num": ["1", "1"]})
    request = {"target_sequences": [{"chain": "a", "sequence": "AAAA"}], "gaps": [
        {"chain": "a", "target_interval": {"start": 1, "stop": 2}},
        {"chain": "a", "target_interval": {"start": 3, "stop": 4}}]}
    candidate = {"entity_id": "1", "source_label_chain": "L", "entity_annotations": [{"pdbx_fragment": "domain"}],
                 "database_alignments": [{"db_align_beg": "300", "db_align_end": "303"}]}
    construct = {"chains": [{"output_chain": "a", "source_label_candidates": [candidate],
        "gap_intervals_zero_based_half_open": [{"start": 1, "stop": 2}, {"start": 3, "stop": 4}]}]}
    return request, construct, block


def test_target_is_deposited_construct_not_expanded_to_database_coordinates():
    request, construct, block = fixture()
    original = deepcopy((request, construct))
    row, = target_facts(request, construct, block)
    assert (request, construct) == original
    assert row["output_chain"] == "a"
    assert row["target_length"] == 4
    assert row["all_original_gap_intervals"] == [{"start": 1, "stop": 2}, {"start": 3, "stop": 4}]
    source, = row["source_candidates"]
    assert source["deposited_observed_residues"] == 2
    assert source["raw_database_alignment_annotations"][0]["db_align_beg"] == "300"
    assert source["database_reference_outside_target_not_requested"]
    assert "scientific_admission" not in row  # Facts do not silently approve a task.


@pytest.mark.parametrize("change", ["sequence", "axis", "chain-case", "missing-candidate", "duplicate-candidate", "model", "entity", "unobserved-axis", "hydrogen-only", "unknown-chain"])
def test_discordant_construct_or_observed_axis_fails_closed(change):
    request, construct, block = fixture()
    atoms = block.get_mmcif_category("_atom_site.")
    if change == "sequence":
        request["target_sequences"][0]["sequence"] = "AAAG"
    elif change == "axis":
        request["gaps"][0]["target_interval"]["start"] = 0
    elif change == "chain-case":
        construct["chains"][0]["output_chain"] = "A"
    elif change == "missing-candidate":
        construct["chains"][0]["source_label_candidates"] = []
    elif change == "duplicate-candidate":
        construct["chains"][0]["source_label_candidates"] *= 2
    elif change == "model":
        atoms["pdbx_PDB_model_num"][1] = "2"
    elif change == "entity":
        atoms["label_entity_id"][1] = "2"
    elif change == "unobserved-axis":
        atoms["label_seq_id"][1] = "4"
    elif change == "hydrogen-only":
        atoms["type_symbol"][1] = "H"
    else:
        request["gaps"].append({"chain": "B", "target_interval": {"start": 0, "stop": 1}})
    block.set_mmcif_category("_atom_site.", atoms)
    with pytest.raises(ValueError):
        target_facts(request, construct, block)


def test_model_scope_is_not_silently_changed_to_a_different_model():
    request, construct, block = fixture()
    atoms = block.get_mmcif_category("_atom_site.")
    atoms["pdbx_PDB_model_num"] = ["2", "2"]
    block.set_mmcif_category("_atom_site.", atoms)
    with pytest.raises(ValueError, match="model1"):
        target_facts(request, construct, block)
