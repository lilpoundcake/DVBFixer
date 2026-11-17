import io
import json

import pytest
from scripts import collect_real_gap_candidates as census


def entry():
    return {"rcsb_id": "TEST", "struct": {"title": "natural gaps"}, "exptl": [{"method": "X-RAY DIFFRACTION"}],
        "polymer_entities": [{"rcsb_id": "TEST_1", "entity_poly": {
            "type": "polypeptide(L)", "pdbx_seq_one_letter_code_can": "A"*180},
            "rcsb_cluster_membership": None, "polymer_entity_instances": [{
                "rcsb_id": "TEST.a", "rcsb_polymer_entity_instance_container_identifiers": {
                    "asym_id": "a", "auth_asym_id": "d"},
                "rcsb_polymer_instance_feature": [{"type": "UNOBSERVED_RESIDUE_XYZ", "feature_positions": [
                    {"beg_seq_id": 1, "end_seq_id": 120}, {"beg_seq_id": 150, "end_seq_id": 155}]}]}]}]}


def test_long_natural_tail_is_not_capped_and_chain_case_is_preserved():
    row = census.summarize(entry())
    assert row["longest_gap_metadata"] == 120
    assert row["natural_gap_count_metadata"] == 2
    chain = row["chains"][0]
    assert (chain["label_chain"], chain["author_chain"]) == ("a", "d")
    assert chain["gaps"][0] == {"start": 0, "stop": 120, "length": 120, "kind": "N-terminal"}
    assert not chain["cluster_assignment_available"]
    assert row["needs_family_grouping_before_final_split"]
    assert row["split_assignment"] is None
    assert not row["scientific_admission_complete"]


def test_feature_intervals_merge_without_double_counting_but_never_clip():
    features = [{"type": "UNOBSERVED_RESIDUE_XYZ", "feature_positions": [
        {"beg_seq_id": 2, "end_seq_id": 4}, {"beg_seq_id": 4, "end_seq_id": 6},
        {"beg_seq_id": 7, "end_seq_id": 8}]}]
    assert census.spans(features, 10) == [{"start": 1, "stop": 8, "length": 7, "kind": "internal"}]
    features[0]["feature_positions"].append({"beg_seq_id": 11, "end_seq_id": 12})
    with pytest.raises(ValueError, match="outside"):
        census.spans(features, 10)


def test_search_long_stratum_has_no_length_ceiling_and_complete_hits_are_required():
    query = census.search_query(100, None)
    missing_nodes = [n["parameters"] for n in query["query"]["nodes"]
                     if "unmodeled" in n["parameters"]["attribute"]]
    assert missing_nodes == [{"attribute": "rcsb_entry_info.deposited_unmodeled_polymer_monomer_count",
                              "operator": "greater_or_equal", "value": 100}]
    assert census.identifiers({"total_count": 2, "result_set": ["1HE8", "1I2M"]}) == ["1HE8", "1I2M"]
    with pytest.raises(ValueError, match="incomplete"):
        census.identifiers({"total_count": 3, "result_set": ["1HE8", "1I2M"]})


def test_receipt_reuse_verifies_exact_request_and_response(tmp_path, monkeypatch):
    (tmp_path / "http").mkdir()
    calls = []

    def fetch(request, timeout):
        calls.append(request.full_url)
        response = io.BytesIO(b'{"total_count": 1, "result_set": ["1HE8"]}')
        response.status = 200
        return response

    monkeypatch.setattr(census, "urlopen", fetch)
    payload = {"query": "public metadata only"}
    first = census.post(tmp_path, "one", census.SEARCH, payload)
    assert census.post(tmp_path, "one", census.SEARCH, payload) == first
    assert calls == [census.SEARCH]
    with pytest.raises(ValueError, match="changed"):
        census.post(tmp_path, "one", census.SEARCH, {"query": "different"})
    (tmp_path / "http/one.response.json").write_text(json.dumps({"different": True}))
    with pytest.raises(ValueError, match="changed"):
        census.post(tmp_path, "one", census.SEARCH, payload)


def test_partial_request_is_not_adopted_or_overwritten(tmp_path):
    (tmp_path / "http").mkdir()
    path = tmp_path / "http/partial.request.json"
    path.write_text("original partial request")
    with pytest.raises(ValueError, match="uncommitted"):
        census.post(tmp_path, "partial", census.SEARCH, {})
    assert path.read_text() == "original partial request"
