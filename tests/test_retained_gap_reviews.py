import json
from copy import deepcopy

import pytest
from scripts.build_retained_gap_reviews import build, gap_rows, scaffold_facts, unique
from scripts.mac_reconstruction_comparison import digest, save


def fixture():
    gaps = [{"chain": "a", "target_interval": {"start": 2, "stop": 4},
             "left_anchor": {}, "right_anchor": {}}]
    request = {"target_sequences": [{"chain": "a", "sequence": "AAAAAA"}],
               "sequence_placements": [{"chain": "a", "target_length": 6,
                                        "observed_target_indices": [0, 1, 4, 5]}], "gaps": gaps}
    case = {"entry_id": "TEST", "scope": "excluded_by_user", "component": 1,
            "gaps": 1, "questions": ["retain-explicit-exclusion-no-active-task"]}
    detail = {"entry_id": "TEST", "current_scope_status": "excluded_by_user",
              "gaps": [{"chain": "a", "start": 2, "stop": 4, "kind": "internal"}]}
    construct = {"entry_id": "TEST", "scope": "excluded_by_user", "chains": [
        {"output_chain": "a", "gap_intervals_zero_based_half_open": [{"start": 2, "stop": 4}]}]}
    differences = {"entry_id": "TEST", "scope": "excluded_by_user", "chains": [
        {"output_chain": "a", "gaps": [{"start": 2, "stop": 4, "terminal": False,
         "consensus_explicit_artifact_positions": [2], "possible_explicit_artifact_positions": [2]}]}]}
    return case, detail, construct, differences, request


def test_scaffold_reports_exact_runs_without_cutoff():
    assert scaffold_facts(20, [0, 1, 4, 5, 6, 19], 2, 4) == {
        "kind": "internal", "chain_target_residues": 20, "chain_observed_residues": 6,
        "chain_absent_residues": 14, "left_contiguous_observed_residues": 2,
        "right_contiguous_observed_residues": 3}
    assert scaffold_facts(2000, [0, 1999], 1, 1999)["kind"] == "internal"
    assert scaffold_facts(4, [2, 3], 0, 2)["kind"] == "N-terminal"
    assert scaffold_facts(4, [0, 1], 2, 4)["kind"] == "C-terminal"
    assert scaffold_facts(4, [], 0, 4)["kind"] == "whole-chain"


@pytest.mark.parametrize("observed,start,stop", [([0, 0, 4], 1, 4), ([0, 2, 4], 1, 4),
                                               ([0, 5], 2, 4), ([-1, 4], 0, 4)])
def test_invalid_scaffold_refused(observed, start, stop):
    with pytest.raises(ValueError):
        scaffold_facts(6, observed, start, stop)


def test_annotation_and_scope_are_not_automatic_admission_or_mask_trimming():
    inputs = fixture()
    original = deepcopy(inputs)
    row, = gap_rows(*inputs)
    assert inputs == original
    assert (row["start"], row["stop"]) == (2, 4)
    assert row["scope"] == "excluded_by_user"
    assert row["scientific_admission"] == "not-admitted"
    assert row["artifact_evidence"] == "partial-explicit-artifact-coverage"
    assert row["entry_review_questions"] == inputs[0]["questions"]


@pytest.mark.parametrize("consensus,possible,expected", [([], [], "no-mapped-artifact-annotation"),
    ([], [2], "candidate-dependent-artifact-coverage"), ([2, 3], [2, 3], "complete-explicit-artifact-coverage")])
def test_annotation_consensus_is_distinct_from_possibility(consensus, possible, expected):
    inputs = fixture()
    annotation = inputs[3]["chains"][0]["gaps"][0]
    annotation["consensus_explicit_artifact_positions"] = consensus
    annotation["possible_explicit_artifact_positions"] = possible
    assert gap_rows(*inputs)[0]["artifact_evidence"] == expected


@pytest.mark.parametrize("change", ["case", "chain", "axis", "partition", "duplicate", "anchor", "annotation", "denominator"])
def test_join_corruption_fails_closed(change):
    case, detail, construct, differences, request = fixture()
    if change == "case":
        differences["entry_id"] = "OTHER"
    elif change == "chain":
        differences["chains"][0]["output_chain"] = "A"
    elif change == "axis":
        construct["chains"][0]["gap_intervals_zero_based_half_open"][0]["start"] = 1
    elif change == "partition":
        request["sequence_placements"][0]["observed_target_indices"] = [0, 4, 5]
    elif change == "duplicate":
        request["gaps"] *= 2
    elif change == "anchor":
        request["gaps"][0]["left_anchor"] = None
    elif change == "annotation":
        differences["chains"][0]["gaps"][0]["possible_explicit_artifact_positions"] = [4]
    else:
        case["gaps"] = 2
    with pytest.raises(ValueError):
        gap_rows(case, detail, construct, differences, request)


def test_duplicate_entries_never_silently_overwrite():
    with pytest.raises(ValueError, match="Duplicate"):
        unique([{"entry_id": "X"}, {"entry_id": "X"}], "entry_id")


def retained_files(root):
    case, detail, construct, differences, request = fixture()
    request_path = root / "request.json"
    save(request_path, request)
    detail.update(request_path=str(request_path), request_sha256=digest(request_path))
    construct.update(request_path=str(request_path), request_sha256=digest(request_path))
    paths = [root / name for name in ("context.json", "construct.json", "difference.json")]
    for path, data in zip(paths, ({}, construct, differences), strict=True):
        save(path, data)
    case["evidence_paths"] = list(map(str, paths))
    detail_path = root / "detail.json"
    save(detail_path, detail)
    inventory_path, reviews_path, receipt_path = [root / name for name in ("inventory.json", "reviews.json", "receipt.json")]
    save(inventory_path, {"rows": [{"entry_id": "TEST", "entry_report": "detail.json"}], "bindings": {}})
    (root / "case-review.tsv").write_text("retained table\n")
    save(reviews_path, {"rows": [case], "summary": {"gaps": 1},
                       "bindings": {str(p): digest(p) for p in paths},
                       "TSV_sha256": digest(root / "case-review.tsv")})
    save(receipt_path, {"bindings": {str(detail_path): digest(detail_path)}})
    return reviews_path, inventory_path, receipt_path


def test_bound_build_retains_no_overwrite_and_table_digest(tmp_path):
    inputs = retained_files(tmp_path)
    output = tmp_path / "output"
    build(*inputs, output)
    data = json.loads((output / "gap-reviews.json").read_text())
    assert data["summary"]["gaps"] == 1
    assert data["summary"]["scientific_admissions_added"] == 0
    assert data["TSV_sha256"] == digest(output / "gap-review.tsv")
    assert data["bindings"][str(tmp_path / "request.json")] == digest(tmp_path / "request.json")
    with pytest.raises(FileExistsError):
        build(*inputs, output)


@pytest.mark.parametrize("name", ["detail.json", "request.json", "difference.json", "case-review.tsv"])
def test_bound_build_refuses_changed_input_before_publication(tmp_path, name):
    inputs = retained_files(tmp_path)
    path = tmp_path / name
    path.write_text(path.read_text() + "\n")
    output = tmp_path / "output"
    with pytest.raises(ValueError):
        build(*inputs, output)
    assert not output.exists()
