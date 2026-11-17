#!/usr/bin/env python3
"""Describe each retained gap's scaffold and annotations, without selecting tasks."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from scripts.build_retained_admission_reviews import checked
from scripts.mac_reconstruction_comparison import digest, save


def unique(rows, key):
    result = {}
    for row in rows:
        value = row[key]
        if value in result:
            raise ValueError(f"Duplicate {key}: {value}")
        result[value] = row
    return result


def scaffold_facts(length, observed, start, stop):
    """Exact sequence-axis occupancy, not a loop/domain or numeric admission rule."""
    if type(length) is not int or length < 1:
        raise ValueError("Invalid target length")
    if any(type(i) is not int for i in observed) or len(set(observed)) != len(observed):
        raise ValueError("Observed axis must contain distinct integer indices")
    observed = set(observed)
    if not observed <= set(range(length)) or not 0 <= start < stop <= length:
        raise ValueError("Invalid target interval or observation")
    if observed & set(range(start, stop)):
        raise ValueError("Declared gap overlaps observations")
    if (start > 0 and start - 1 not in observed) or (stop < length and stop not in observed):
        raise ValueError("Gap is not bounded by the declared observations")
    left = start
    while left > 0 and left - 1 in observed:
        left -= 1
    right = stop
    while right < length and right in observed:
        right += 1
    kind = ("whole-chain" if start == 0 and stop == length else
            "N-terminal" if start == 0 else "C-terminal" if stop == length else "internal")
    return {"kind": kind, "chain_target_residues": length,
            "chain_observed_residues": len(observed), "chain_absent_residues": length - len(observed),
            "left_contiguous_observed_residues": start - left,
            "right_contiguous_observed_residues": right - stop}


def gap_rows(case, detail, construct, differences, request):
    entry, scope = case["entry_id"], case["scope"]
    if any(data["entry_id"] != entry for data in (detail, construct, differences)):
        raise ValueError("Joined entry identity changed")
    if detail["current_scope_status"] != scope or any(data["scope"] != scope for data in (construct, differences)):
        raise ValueError("Joined scope changed")
    targets = unique(request["target_sequences"], "chain")
    placements = unique(request["sequence_placements"], "chain")
    if set(targets) != set(placements):
        raise ValueError("Target/placement chains differ")
    annotations = {}
    for chain in differences["chains"]:
        for gap in chain["gaps"]:
            key = (chain["output_chain"], gap["start"], gap["stop"])
            if key in annotations:
                raise ValueError("Duplicate annotated gap")
            annotations[key] = gap
    axes = [(c["output_chain"], g["start"], g["stop"]) for c in construct["chains"]
            for g in c["gap_intervals_zero_based_half_open"]]
    expected = [(g["chain"], g["target_interval"]["start"], g["target_interval"]["stop"]) for g in request["gaps"]]
    detailed = [(g["chain"], g["start"], g["stop"]) for g in detail["gaps"]]
    if (len(expected) != len(set(expected)) or len(axes) != len(set(axes)) or
            set(expected) != set(axes) or set(expected) != set(annotations) or detailed != expected or
            len(expected) != case["gaps"]):
        raise ValueError("Gap axes or denominators disagree")
    for chain, target in targets.items():
        placement = placements[chain]
        length = len(target["sequence"])
        if placement["target_length"] != length:
            raise ValueError("Target length changed")
        absent = [i for c, a, b in expected if c == chain for i in range(a, b)]
        observed = placement["observed_target_indices"]
        if (len(absent) != len(set(absent)) or len(observed) != len(set(observed)) or
                set(absent) & set(observed) or set(absent) | set(observed) != set(range(length))):
            raise ValueError("Gap and observed masks do not partition target")
    rows = []
    for index, (chain, start, stop) in enumerate(expected):
        facts = scaffold_facts(len(targets[chain]["sequence"]), placements[chain]["observed_target_indices"], start, stop)
        gap, annotation = request["gaps"][index], annotations[(chain, start, stop)]
        if (gap["left_anchor"] is not None) != (start > 0) or (gap["right_anchor"] is not None) != (stop < facts["chain_target_residues"]):
            raise ValueError("Anchor availability disagrees with target boundary")
        if annotation["terminal"] != (facts["kind"] != "internal") or detail["gaps"][index]["kind"] != facts["kind"]:
            raise ValueError("Gap boundary kind changed")
        consensus = annotation["consensus_explicit_artifact_positions"]
        possible = annotation["possible_explicit_artifact_positions"]
        if (len(consensus) != len(set(consensus)) or len(possible) != len(set(possible)) or
                not set(consensus) <= set(possible) <= set(range(start, stop))):
            raise ValueError("Invalid artifact annotation positions")
        coverage = ("complete-explicit-artifact-coverage" if len(consensus) == stop-start else
                    "partial-explicit-artifact-coverage" if consensus else
                    "candidate-dependent-artifact-coverage" if possible else "no-mapped-artifact-annotation")
        rows.append({"entry_id": entry, "component": case["component"], "scope": scope,
                     "gap_index": index, "chain": chain, "start": start, "stop": stop,
                     "length": stop-start, **facts, "artifact_evidence": coverage,
                     "consensus_artifact_residues": len(consensus), "possible_artifact_residues": len(possible),
                     "entry_review_questions": case["questions"],
                     "scientific_admission": "not-admitted", "split_assignment": "unassigned"})
    return rows


def build(reviews_path, inventory_path, receipt_path, output):
    reviews, inventory, receipt = map(checked, (reviews_path, inventory_path, receipt_path))
    if digest(reviews_path.parent / "case-review.tsv") != reviews["TSV_sha256"]:
        raise ValueError("Case review TSV changed")
    cases, entries = unique(reviews["rows"], "entry_id"), unique(inventory["rows"], "entry_id")
    if set(cases) != set(entries):
        raise ValueError("Retained entry denominators disagree")
    bindings = {str(p): digest(p) for p in (reviews_path, inventory_path, receipt_path)}
    rows = []
    for entry, case in cases.items():
        detail_path = inventory_path.parent / entries[entry]["entry_report"]
        if not detail_path.resolve().is_relative_to(inventory_path.parent.resolve()):
            raise ValueError("Entry report escapes inventory")
        if str(detail_path) not in receipt["bindings"] or digest(detail_path) != receipt["bindings"][str(detail_path)]:
            raise ValueError("Unbound entry report")
        detail = json.loads(detail_path.read_text())
        context_path, construct_path, difference_path = map(Path, case["evidence_paths"])
        for path in (context_path, construct_path, difference_path):
            if str(path) not in reviews["bindings"] or digest(path) != reviews["bindings"][str(path)]:
                raise ValueError("Unbound case evidence")
        construct, differences = [json.loads(p.read_text()) for p in (construct_path, difference_path)]
        request_path = Path(detail["request_path"])
        if (digest(request_path) != detail["request_sha256"] or
                str(request_path) != construct["request_path"] or digest(request_path) != construct["request_sha256"]):
            raise ValueError("Request changed")
        rows.extend(gap_rows(case, detail, construct, differences, json.loads(request_path.read_text())))
        for path in (detail_path, context_path, construct_path, difference_path, request_path):
            bindings[str(path)] = digest(path)
    if len(rows) != reviews["summary"]["gaps"]:
        raise ValueError("Global gap denominator changed")
    summary = {"entries": len(cases), "gaps": len(rows),
               "kind_counts": dict(Counter(r["kind"] for r in rows)),
               "boundary_artifact_counts": dict(Counter(f'{r["kind"]}|{r["artifact_evidence"]}' for r in rows)),
               "entries_with_internal_gaps": len({r["entry_id"] for r in rows if r["kind"] == "internal"}),
               "scientific_admissions_added": 0, "request_masks_changed": False,
               "new_sources_or_restoration_calculations": 0, "split_frozen": False,
               "interpretation": "Axis occupancy is not loop/domain admission; absent annotations are not proof of native sequence. Entry-wide questions still apply; no subset task is authorized."}
    output.mkdir(parents=True, exist_ok=False)
    with (output / "gap-review.tsv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows({**r, "entry_review_questions": ";".join(r["entry_review_questions"])} for r in rows)
    save(output / "gap-reviews.json", {"summary": summary, "rows": rows, "bindings": bindings,
         "builder_sha256": digest(Path(__file__)), "TSV_sha256": digest(output / "gap-review.tsv")})
    save(output / "summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reviews", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.reviews.resolve(), args.inventory.resolve(), args.receipt.resolve(), args.output.resolve())
