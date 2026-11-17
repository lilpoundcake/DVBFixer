#!/usr/bin/env python3
"""Bind retained construct annotations to requested axes; never infer tags or admit loops."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import gemmi
from scripts.audit_real_gap_input_placement import placement_evidence
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import DiffusionRequest

SOURCE_LEDGERS = (
    "original-sources-attempt01",
    "source-expansion-attempt01/acquisition",
    "hard-original-sources-attempt01/acquisition",
    "medium-original-sources-attempt01/acquisition",
)


def category_rows(block, prefix):
    category = block.get_mmcif_category(prefix)
    if not category:
        return []
    lengths = {len(values) for values in category.values()}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent deposited category {prefix}")
    return [dict(zip(category, values, strict=True)) for values in zip(*category.values(), strict=True)]


def chain_annotations(chain, entities, references, alignments, differences):
    entity = chain["entity_id"]
    author_chains = set(chain["author_chain_ids"])
    refs = [row for row in references if row.get("entity_id") == entity]
    ref_ids = {row.get("id") for row in refs}
    spans = []
    for row in alignments:
        strands = {name.strip() for name in (row.get("pdbx_strand_id") or "").split(",") if name.strip()}
        if row.get("ref_id") in ref_ids and author_chains & strands:
            spans.append(row)
    alignment_ids = {row.get("align_id") for row in spans}
    return {"source_label_chain": chain["label_chain"], "source_author_chains": sorted(author_chains),
        "entity_id": entity, "entity_annotations": [row for row in entities if row.get("id") == entity],
        "database_references": refs, "database_alignments": spans,
        "deposited_sequence_differences": [row for row in differences if row.get("align_id") in alignment_ids],
        "reference_bounds_are_raw_deposited_fields_not_reinterpreted_author_or_label_axes": True,
        "absence_of_annotations_does_not_prove_absence_of_tags_or_construct_changes": True}


def construct_facts(request, inventory, block):
    placements = placement_evidence(request, inventory)
    if not all(row["request_agrees_with_deposited_axis"] for row in placements):
        raise ValueError("Request has ambiguous or discordant deposited placement")
    categories = {name: category_rows(block, name) for name in
                  ("_entity.", "_struct_ref.", "_struct_ref_seq.", "_struct_ref_seq_dif.")}
    source_chains = {chain["label_chain"]: chain for chain in inventory["chains"]}
    reports = []
    gaps_checked = 0
    for placement in placements:
        output_chain = placement["output_chain"]
        requested = [gap for gap in request.gaps if gap.chain == output_chain]
        candidates = []
        for label in placement["source_label_candidates"]:
            chain = source_chains[label]
            declared_spans = {(row["start"], row["stop"]) for row in chain["natural_absent_coordinate_spans"]}
            requested_spans = {(gap.target_interval.start, gap.target_interval.stop) for gap in requested}
            if requested_spans != declared_spans:
                raise ValueError("Requested gaps differ from all deposited absent-coordinate spans")
            candidates.append(chain_annotations(chain, categories["_entity."], categories["_struct_ref."],
                categories["_struct_ref_seq."], categories["_struct_ref_seq_dif."]))
        if not candidates:
            raise ValueError("No source-label candidates for requested output chain")
        reports.append({"output_chain": output_chain, "source_label_candidates": candidates,
            "source_label_assignment_is_unique": len(candidates) == 1,
            "all_candidate_deposited_gap_axes_match": True,
            "gap_intervals_zero_based_half_open": [{"start": gap.target_interval.start,
                "stop": gap.target_interval.stop, "terminal": gap.left_anchor is None or gap.right_anchor is None}
                for gap in requested],
            "construct_adjudication": "pending; annotations are evidence, not a loop-scale or tag decision"})
        gaps_checked += len(requested)
    return {"chains": reports, "gaps_checked": gaps_checked,
        "terminal_gap_count": sum(row["terminal"] for chain in reports for row in chain["gap_intervals_zero_based_half_open"]),
        "all_declared_gaps_match_deposited_absent_coordinate_axes": True,
        "new_tag_or_domain_classification": False, "scientific_admission": False}


def audit(inventory_path, sources_root, output):
    inventory = json.loads(inventory_path.read_text())
    for name, expected in inventory["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Retained input changed: {name}")
    bindings = {str(inventory_path): digest(inventory_path)}
    sources = {}
    for folder in SOURCE_LEDGERS:
        source_path = sources_root / folder / "source-results.json"
        complete_path = source_path.parent / "completed.json"
        complete = json.loads(complete_path.read_text())
        if digest(source_path) != complete["source_results_sha256"]:
            raise ValueError("Original-source inventory completion binding changed")
        for row in json.loads(source_path.read_text())["rows"]:
            if row["entry_id"] in sources:
                raise ValueError("Duplicate original source across retained ledgers")
            sources[row["entry_id"]] = row
        for path in (source_path, complete_path):
            bindings[str(path)] = digest(path)
    output.mkdir(parents=True, exist_ok=False)
    (output / "entries").mkdir()
    rows = []
    for row in inventory["rows"]:
        entry = row["entry_id"]
        detail_path = inventory_path.parent / row["entry_report"]
        if not detail_path.resolve().is_relative_to(inventory_path.parent.resolve()):
            raise ValueError("Entry report escapes retained review")
        detail = json.loads(detail_path.read_text())
        request_path = Path(detail["request_path"])
        if digest(request_path) != detail["request_sha256"]:
            raise ValueError("Retained request changed")
        request = DiffusionRequest.from_json(request_path.read_text())
        original = sources[entry]
        source_path = Path(original["source_path"])
        if digest(source_path) != original["source_sha256"]:
            raise ValueError("Original experimental mmCIF changed")
        block = gemmi.cif.read(str(source_path)).sole_block()
        facts = construct_facts(request, original["inventory"], block)
        save(output / "entries" / f"{entry}.json", {"entry_id": entry, "scope": row["current_scope_status"],
            "original_source_path": str(source_path), "original_source_sha256": digest(source_path),
            "request_path": str(request_path), "request_sha256": digest(request_path), **facts})
        chains = facts["chains"]
        candidates = [candidate for chain in chains for candidate in chain["source_label_candidates"]]
        rows.append({"entry_id": entry, "scope": row["current_scope_status"], "gaps_checked": facts["gaps_checked"],
            "terminal_gap_count": facts["terminal_gap_count"], "output_chains": len(chains),
            "chains_with_multiple_source_label_candidates": sum(not chain["source_label_assignment_is_unique"] for chain in chains),
            "has_entity_fragment_annotation": any(entity.get("pdbx_fragment") for c in candidates for entity in c["entity_annotations"]),
            "has_database_alignment": any(c["database_alignments"] for c in candidates),
            "has_deposited_sequence_difference": any(c["deposited_sequence_differences"] for c in candidates),
            "scientific_admission": False, "report": f"entries/{entry}.json"})
        for path in (detail_path, request_path, source_path):
            bindings[str(path)] = digest(path)
    summary = {"entries": len(rows), "gaps_checked": sum(row["gaps_checked"] for row in rows),
        "terminal_gaps": sum(row["terminal_gap_count"] for row in rows),
        "entries_with_terminal_gaps": sum(row["terminal_gap_count"] > 0 for row in rows),
        "entries_with_fragment_annotations": sum(row["has_entity_fragment_annotation"] for row in rows),
        "entries_with_database_alignments": sum(row["has_database_alignment"] for row in rows),
        "entries_with_deposited_sequence_differences": sum(row["has_deposited_sequence_difference"] for row in rows),
        "scope_counts": dict(Counter(row["scope"] for row in rows)),
        "requested_gap_axes_match_all_source_label_candidates": True,
        "scientific_admissions_or_new_sources_or_restoration_runs": 0,
        "terminal_tag_or_domain_semantics_established": False}
    save(output / "construct-audit.json", {"summary": summary, "rows": rows, "bindings": bindings,
        "auditor_sha256": digest(Path(__file__)),
        "reports_sha256": {row["report"]: digest(output / row["report"]) for row in rows}})
    save(output / "summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--sources-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.inventory.resolve(), args.sources_root.resolve(), args.output.resolve())
