from types import SimpleNamespace

from scripts.audit_real_gap_input_placement import placement_evidence

from dvbfixer.model.diffusion.contract import ResidueIdentity


def example(indices=(0, 3)):
    request = SimpleNamespace(target_sequences=[SimpleNamespace(chain="a", sequence="AAAA")],
        sequence_placements=[SimpleNamespace(chain="a", observed_target_indices=indices,
            observed_residues=[ResidueIdentity("a", "10"), ResidueIdentity("a", "13", "A")])])
    chain = {"label_chain": "ORIGINAL", "sequence": "AAAA", "observed_author_axis": [
        {"author_residue_number": "10", "insertion_code": "", "label_seq_id": 1},
        {"author_residue_number": "13", "insertion_code": "A", "label_seq_id": 4}]}
    return request, chain


def test_assembly_renamed_chain_can_match_exact_deposited_positions():
    request, chain = example()
    row, = placement_evidence(request, {"chains": [chain]})
    assert row["request_agrees_with_deposited_axis"]
    assert row["source_label_candidates"] == ["ORIGINAL"]


def test_equal_residue_sequence_does_not_authorize_wrong_repeat_placement():
    request, chain = example((0, 2))
    row, = placement_evidence(request, {"chains": [chain]})
    assert not row["request_agrees_with_deposited_axis"]


def test_ambiguous_original_chain_axes_are_not_guessed():
    request, chain = example()
    other = {**chain, "label_chain": "OTHER", "observed_author_axis": [
        {"author_residue_number": "10", "insertion_code": "", "label_seq_id": 2},
        {"author_residue_number": "13", "insertion_code": "A", "label_seq_id": 4}]}
    row, = placement_evidence(request, {"chains": [chain, other]})
    assert not row["unambiguous_deposited_position_evidence"]
    assert not row["request_agrees_with_deposited_axis"]
