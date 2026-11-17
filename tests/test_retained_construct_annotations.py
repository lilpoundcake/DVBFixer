from types import SimpleNamespace

import gemmi
import pytest
from scripts.audit_retained_construct_annotations import chain_annotations, construct_facts

from dvbfixer.model.diffusion.contract import ResidueIdentity


def test_annotation_mapping_preserves_author_chain_case_and_raw_bounds():
    chain = {"entity_id": "1", "label_chain": "L", "author_chain_ids": ["a"]}
    references = [{"id": "ref", "entity_id": "1"}]
    alignment = {"align_id": "lower", "ref_id": "ref", "pdbx_strand_id": "a", "seq_align_beg": "-2"}
    other = {"align_id": "upper", "ref_id": "ref", "pdbx_strand_id": "A"}
    facts = chain_annotations(chain, [], references, [other, alignment],
                              [{"align_id": "lower", "details": "expression tag"}])
    assert facts["database_alignments"] == [alignment]
    assert facts["database_alignments"][0]["seq_align_beg"] == "-2"
    assert len(facts["deposited_sequence_differences"]) == 1


def fixture():
    observed = ResidueIdentity("a", "2", "A")
    request = SimpleNamespace(target_sequences=(SimpleNamespace(chain="a", sequence="AA"),),
        sequence_placements=(SimpleNamespace(chain="a", observed_residues=(observed,), observed_target_indices=(1,)),),
        gaps=(SimpleNamespace(chain="a", target_interval=SimpleNamespace(start=0, stop=1),
                             left_anchor=None, right_anchor=observed),))
    chain = {"sequence": "AA", "label_chain": "X", "entity_id": "1", "author_chain_ids": ["a"],
        "observed_author_axis": [{"author_residue_number": "2", "insertion_code": "A", "label_seq_id": 2}],
        "natural_absent_coordinate_spans": [{"start": 0, "stop": 1}]}
    return request, chain


def test_missing_annotations_remain_unknown_not_no_tags_or_admission():
    request, chain = fixture()
    facts = construct_facts(request, {"chains": [chain]}, gemmi.cif.read_string("data_test\n").sole_block())
    assert facts["gaps_checked"] == 1 and facts["terminal_gap_count"] == 1
    assert not facts["scientific_admission"] and not facts["new_tag_or_domain_classification"]
    assert facts["chains"][0]["source_label_candidates"][0]["database_alignments"] == []


def test_all_candidate_gap_axes_must_match_not_first_convenient_label():
    request, chain = fixture()
    incompatible = {**chain, "label_chain": "Y", "natural_absent_coordinate_spans": []}
    with pytest.raises(ValueError, match="absent-coordinate spans"):
        construct_facts(request, {"chains": [chain, incompatible]}, gemmi.cif.read_string("data_test\n").sole_block())
