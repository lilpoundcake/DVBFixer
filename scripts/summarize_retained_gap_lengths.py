#!/usr/bin/env python3
"""Report gap lengths with explicit entry/copy denominators; never admit tasks."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

from scripts.build_retained_admission_reviews import checked
from scripts.mac_reconstruction_comparison import digest, save

SCOPES = ("loop_scope_and_context_review_pending", "small_control_only", "excluded_by_user")
KINDS = ("internal", "N-terminal", "C-terminal", "whole-chain")
BINS = (("1-8", 1, 8), ("9-24", 9, 24), ("25-49", 25, 49), ("50-99", 50, 99), ("100+", 100, None))


def distribution(lengths):
    if any(type(n) is not int or n < 1 for n in lengths):
        raise ValueError("Gap lengths must be positive integers")
    return {"gap_instances": len(lengths), "missing_residue_instances": sum(lengths),
            "minimum": min(lengths) if lengths else None, "median": median(lengths) if lengths else None,
            "maximum": max(lengths) if lengths else None,
            "length_bins": {name: sum(lo <= n and (hi is None or n <= hi) for n in lengths)
                            for name, lo, hi in BINS}}


def summarize(rows, details):
    grouped, seen = defaultdict(list), set()
    for row in rows:
        entry, index = row["entry_id"], row["gap_index"]
        if entry not in details or (entry, index) in seen or type(index) is not int or index < 0:
            raise ValueError("Duplicate/unmapped gap instance")
        seen.add((entry, index))
        detail = details[entry]
        if row["scope"] not in SCOPES or row["scope"] != detail["current_scope_status"] or row["kind"] not in KINDS:
            raise ValueError("Unknown or changed scope/kind")
        if index >= len(detail["gaps"]):
            raise ValueError("Gap index outside bound entry")
        source = detail["gaps"][index]
        if any(row[key] != source[key] for key in ("chain", "start", "stop", "length", "kind")):
            raise ValueError("Bound gap axis changed")
        if type(row["length"]) is not int or row["length"] != row["stop"] - row["start"] or row["length"] < 1:
            raise ValueError("Invalid gap length")
        if row["chain"] not in detail["sequence_hashes_by_case_sensitive_chain"]:
            raise ValueError("Unknown case-sensitive chain")
        if row["scientific_admission"] != "not-admitted" or row["split_assignment"] != "unassigned":
            raise ValueError("Unexpected admission/split mutation")
        grouped[entry].append(row)
    if grouped.keys() != details.keys():
        raise ValueError("Entry denominator changed")

    def equivalence(row):
        return (row["entry_id"], details[row["entry_id"]]["sequence_hashes_by_case_sensitive_chain"][row["chain"]],
                row["start"], row["stop"], row["kind"])

    entries = []
    for entry in sorted(grouped):
        gaps, detail = grouped[entry], details[entry]
        if {g["gap_index"] for g in gaps} != set(range(len(detail["gaps"]))):
            raise ValueError("Per-entry gap denominator changed")
        if sum(g["length"] for g in gaps) != detail["target_residues"] - detail["observed_residues"]:
            raise ValueError("Missing residues do not partition the bound target")
        lengths = {kind: [g["length"] for g in gaps if g["kind"] == kind] for kind in KINDS}
        entries.append({"entry_id": entry, "scope": detail["current_scope_status"],
            "component": gaps[0]["component"], "gap_instances": len(gaps),
            "within_entry_sequence_interval_classes": len({equivalence(g) for g in gaps}),
            "target_residue_instances": detail["target_residues"],
            "observed_residue_instances": detail["observed_residues"],
            "missing_residue_instances": sum(g["length"] for g in gaps),
            "longest_gap": max(g["length"] for g in gaps),
            "internal_lengths": lengths["internal"], "N_terminal_lengths": lengths["N-terminal"],
            "C_terminal_lengths": lengths["C-terminal"], "whole_chain_lengths": lengths["whole-chain"],
            "gap_instances_by_chain": {chain: [{"gap_index": g["gap_index"], "start": g["start"],
                "stop": g["stop"], "kind": g["kind"], "length": g["length"]} for g in gaps if g["chain"] == chain]
                for chain in sorted({g["chain"] for g in gaps})},
            "scientific_admission": "not-admitted", "split_assignment": "unassigned"})
    distributions = []
    for scope in ("all_retained", *SCOPES):
        for kind in ("all_kinds", *KINDS):
            subset = [g for g in rows if (scope == "all_retained" or g["scope"] == scope)
                      and (kind == "all_kinds" or g["kind"] == kind)]
            distributions.append({"scope": scope, "kind": kind, "entries": len({g["entry_id"] for g in subset}),
                "within_entry_sequence_interval_classes": len({equivalence(g) for g in subset}),
                **distribution([g["length"] for g in subset])})
    return {"entries": entries, "distributions": distributions,
            "summary": {"entries": len(entries), "gap_instances": len(rows),
                "scope_entry_counts": dict(Counter(e["scope"] for e in entries)),
                "within_entry_sequence_interval_classes": len({equivalence(g) for g in rows}),
                "missing_residue_instances": sum(g["length"] for g in rows),
                "bins_are_reporting_only_not_admission_thresholds": True,
                "equivalence_is_not_independence_or_identical_conditioning": True,
                "equivalence_definition": "Same entry, exact target sequence hash, interval and kind; chain copies remain separate tasks/records. Never collapsed across entries.",
                "scientific_admissions_or_mask_changes_or_new_sources_or_restoration_runs": 0}}


def build(gaps_path, receipt_path, output):
    receipt = json.loads(receipt_path.read_text())
    if digest(gaps_path) != receipt["report_sha256"]:
        raise ValueError("Gap report differs from its independent readback")
    gaps = checked(gaps_path)
    if digest(gaps_path.parent / "gap-review.tsv") != gaps["TSV_sha256"]:
        raise ValueError("Gap TSV changed")
    details, bindings = {}, {str(p): digest(p) for p in (gaps_path, receipt_path)}
    for name in gaps["bindings"]:
        path = Path(name)
        if path.parent.name != "entries" or path.suffix != ".json":
            continue
        detail = json.loads(path.read_text())
        if "sequence_hashes_by_case_sensitive_chain" not in detail:
            continue
        entry = detail["entry_id"]
        if entry in details:
            raise ValueError("Duplicate bound entry")
        request_path = Path(detail["request_path"])
        if digest(request_path) != detail["request_sha256"]:
            raise ValueError("Request changed")
        request = json.loads(request_path.read_text())
        targets = {t["chain"]: hashlib.sha256(t["sequence"].encode()).hexdigest() for t in request["target_sequences"]}
        if len(targets) != len(request["target_sequences"]) or targets != detail["sequence_hashes_by_case_sensitive_chain"]:
            raise ValueError("Case-sensitive sequence hashes changed")
        details[entry] = detail
        for p in (path, request_path):
            bindings[str(p)] = digest(p)
    result = summarize(gaps["rows"], details)
    if result["summary"]["entries"] != gaps["summary"]["entries"] or result["summary"]["gap_instances"] != gaps["summary"]["gaps"]:
        raise ValueError("Global denominator changed")
    output.mkdir(parents=True, exist_ok=False)
    with (output / "entry-gap-lengths.tsv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(result["entries"][0]), delimiter="\t")
        writer.writeheader()
        writer.writerows({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict)) else value
                          for key, value in row.items()} for row in result["entries"])
    with (output / "length-distributions.tsv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(result["distributions"][0]), delimiter="\t")
        writer.writeheader()
        writer.writerows({**row, "length_bins": json.dumps(row["length_bins"])} for row in result["distributions"])
    save(output / "length-report.json", {**result, "bindings": bindings, "reporter_sha256": digest(Path(__file__)),
         "tables_sha256": {name: digest(output / name) for name in ("entry-gap-lengths.tsv", "length-distributions.tsv")}})
    print(json.dumps(result["summary"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaps", type=Path, required=True)
    parser.add_argument("--readback", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.gaps.resolve(), args.readback.resolve(), args.output.resolve())
