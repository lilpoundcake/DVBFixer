import hashlib
import json
from copy import deepcopy

import pytest
from scripts.mac_reconstruction_comparison import digest, save
from scripts.summarize_retained_gap_lengths import build, distribution, summarize


def fixture():
    source = {"entry_id": "TEST", "gap_index": 0, "chain": "a", "component": 1,
              "scope": "loop_scope_and_context_review_pending", "kind": "internal", "start": 2,
              "stop": 5, "length": 3, "scientific_admission": "not-admitted", "split_assignment": "unassigned"}
    second = {**source, "gap_index": 1, "chain": "A"}
    rows = [source, second]
    detail = {"current_scope_status": source["scope"], "target_residues": 20, "observed_residues": 14,
              "sequence_hashes_by_case_sensitive_chain": {"a": "same-sequence", "A": "same-sequence"},
              "gaps": [{key: row[key] for key in ("chain", "start", "stop", "length", "kind")} for row in rows]}
    return rows, {"TEST": detail}


def test_counts_are_gaps_not_residues_and_copies_not_independent_entries():
    rows, details = fixture()
    original = deepcopy((rows, details))
    data = summarize(rows, details)
    assert (rows, details) == original
    assert data["summary"]["entries"] == 1
    assert data["summary"]["gap_instances"] == 2
    assert data["summary"]["missing_residue_instances"] == 6
    assert data["summary"]["within_entry_sequence_interval_classes"] == 1
    assert set(data["entries"][0]["gap_instances_by_chain"]) == {"a", "A"}
    assert data["entries"][0]["internal_lengths"] == [3, 3]


def test_equivalent_sequence_axis_never_collapses_different_entries():
    rows, details = fixture()
    details["OTHER"] = deepcopy(details["TEST"])
    rows += [{**row, "entry_id": "OTHER"} for row in list(rows)]
    result = summarize(rows, details)
    assert result["summary"]["entries"] == 2
    assert result["summary"]["within_entry_sequence_interval_classes"] == 2


def test_excluded_oversized_gap_is_not_in_pending_distribution():
    rows, details = fixture()
    excluded = {**rows[0], "entry_id": "EXCLUDED", "scope": "excluded_by_user", "length": 812, "stop": 814}
    rows.append(excluded)
    details["EXCLUDED"] = {"current_scope_status": "excluded_by_user", "target_residues": 1000,
        "observed_residues": 188, "sequence_hashes_by_case_sensitive_chain": {"a": "other"}, "gaps": [excluded]}
    data = summarize(rows, details)
    pending = next(d for d in data["distributions"] if d["scope"] == rows[0]["scope"] and d["kind"] == "internal")
    assert pending["maximum"] == 3
    assert pending["gap_instances"] == 2
    assert data["summary"]["scope_entry_counts"]["excluded_by_user"] == 1


def test_empty_and_exact_histogram_boundaries():
    assert distribution([])["minimum"] is None
    lengths = [1, 8, 9, 24, 25, 49, 50, 99, 100, 1069]
    result = distribution(lengths)
    assert result["length_bins"] == {"1-8": 2, "9-24": 2, "25-49": 2, "50-99": 2, "100+": 2}
    assert result["median"] == 37
    assert result["missing_residue_instances"] == sum(lengths)


@pytest.mark.parametrize("value", [0, -1, 1.5, True])
def test_nonpositive_or_noninteger_length_refused(value):
    with pytest.raises(ValueError):
        distribution([value])


@pytest.mark.parametrize("mutation", ["duplicate", "length", "scope", "chain", "missing-entry", "gap-count", "residue-count", "admission"])
def test_inconsistent_bound_data_refused(mutation):
    rows, details = fixture()
    if mutation == "duplicate":
        rows.append(dict(rows[0]))
    elif mutation == "length":
        rows[0]["length"] = 4
    elif mutation == "scope":
        rows[0]["scope"] = "admitted"
    elif mutation == "chain":
        rows[0]["chain"] = "B"
    elif mutation == "missing-entry":
        details["OTHER"] = deepcopy(details["TEST"])
    elif mutation == "gap-count":
        rows.pop()
    elif mutation == "residue-count":
        details["TEST"]["observed_residues"] = 15
    else:
        rows[0]["scientific_admission"] = "admitted"
    with pytest.raises(ValueError):
        summarize(rows, details)


def retained_files(root):
    rows, details = fixture()
    (root / "entries").mkdir()
    request_path = root / "request.json"
    save(request_path, {"target_sequences": [{"chain": chain, "sequence": "A" * 10} for chain in ("a", "A")]})
    detail = details["TEST"]
    detail.update(entry_id="TEST", request_path=str(request_path), request_sha256=digest(request_path),
                  sequence_hashes_by_case_sensitive_chain={chain: hashlib.sha256(b"A" * 10).hexdigest() for chain in ("a", "A")})
    detail_path = root / "entries/TEST.json"
    save(detail_path, detail)
    (root / "gap-review.tsv").write_text("fixture\n")
    gap_path = root / "gap-reviews.json"
    save(gap_path, {"rows": rows, "summary": {"entries": 1, "gaps": 2},
                   "bindings": {str(p): digest(p) for p in (request_path, detail_path)},
                   "TSV_sha256": digest(root / "gap-review.tsv")})
    receipt_path = root / "readback.json"
    save(receipt_path, {"report_sha256": digest(gap_path)})
    return gap_path, receipt_path


def test_bound_build_and_no_overwrite(tmp_path):
    paths = retained_files(tmp_path)
    output = tmp_path / "output"
    build(*paths, output)
    report = json.loads((output / "length-report.json").read_text())
    assert report["summary"]["entries"] == 1
    assert report["summary"]["within_entry_sequence_interval_classes"] == 1
    for name, expected in report["tables_sha256"].items():
        assert digest(output / name) == expected
    with pytest.raises(FileExistsError):
        build(*paths, output)


@pytest.mark.parametrize("name", ["gap-reviews.json", "gap-review.tsv", "request.json", "entries/TEST.json"])
def test_changed_input_refused_before_publication(tmp_path, name):
    paths = retained_files(tmp_path)
    path = tmp_path / name
    path.write_text(path.read_text() + "\n")
    with pytest.raises(ValueError):
        build(*paths, tmp_path / "output")
    assert not (tmp_path / "output").exists()
