#!/usr/bin/env python3
"""Map deposited difference annotations via their dictionary-defined entity label axis."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import gemmi
from scripts.audit_retained_construct_annotations import category_rows
from scripts.mac_reconstruction_comparison import digest, save

# Exact depositor descriptions, not motif inference or an automatic exclusion rule.
ARTIFACT_DESCRIPTIONS = frozenset({"expression tag", "cloning artifact"})
SEQ_NUM_DEFINITION = "https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_struct_ref_seq_dif.seq_num.html"


def text(value):
    return "" if value is None or value is False else str(value)


def map_difference(row, candidate, scheme):
    number = text(row.get("seq_num"))
    strand = text(row.get("pdbx_pdb_strand_id"))
    monomer = text(row.get("mon_id"))
    database_only_deletion = (text(row.get("details")).strip().casefold() == "deletion"
        and not number and not monomer and bool(text(row.get("db_mon_id")))
        and bool(text(row.get("pdbx_seq_db_seq_num"))))
    matches = []
    if number and monomer and strand in candidate["source_author_chains"]:
        for item in scheme:
            if (item.get("asym_id") == candidate["source_label_chain"]
                    and item.get("entity_id") == candidate["entity_id"]
                    and text(item.get("pdb_strand_id")) == strand
                    and text(item.get("seq_id")) == number
                    and text(item.get("mon_id")) == monomer):
                try:
                    position = int(text(item.get("seq_id")))
                except ValueError:
                    continue
                if position > 0:
                    matches.append(position-1)
    positions = sorted(set(matches))
    return {"annotation": row, "source_label_chain": candidate["source_label_chain"],
        "label_index_zero_based": positions[0] if len(matches) == 1 else None,
        "mapping_status": "not-applicable-database-deletion" if database_only_deletion else
                          "mapped-deposited-residue" if len(matches) == 1 else "unresolved",
        "scheme_matches": len(matches), "candidate_label_positions": positions,
        "mapping_basis": "struct_ref_seq_dif.seq_num points to entity_poly_seq.num; verify exact scheme seq_id/entity/label chain/author chain/monomer, never treat author number as label position",
        "dictionary_definition": SEQ_NUM_DEFINITION,
        "is_explicit_artifact_description": text(row.get("details")).strip().casefold() in ARTIFACT_DESCRIPTIONS}


def gap_evidence(gap, mapped_candidates):
    positions = set(range(gap["start"], gap["stop"]))
    candidates = []
    for rows in mapped_candidates:
        candidates.append({row["label_index_zero_based"] for row in rows
            if row["is_explicit_artifact_description"] and row["label_index_zero_based"] in positions})
    consensus = set.intersection(*candidates) if candidates else set()
    possible = set.union(*candidates) if candidates else set()
    return {**gap, "consensus_explicit_artifact_positions": sorted(consensus),
        "possible_explicit_artifact_positions": sorted(possible),
        "coverage": "complete-explicit-artifact-annotation" if consensus == positions and positions else
                    "partial-explicit-artifact-annotation" if consensus else
                    "candidate-dependent-artifact-annotation" if possible else "no-mapped-artifact-annotation",
        "no_annotation_is_not_proof_of_native_sequence": True,
        "request_mask_or_admission_changed": False}


def audit(construct_path, output):
    audit_data = json.loads(construct_path.read_text())
    for name, expected in audit_data["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Bound construct evidence changed: {name}")
    output.mkdir(parents=True, exist_ok=False)
    (output / "entries").mkdir()
    rows, bindings = [], {str(construct_path): digest(construct_path)}
    categories = Counter()
    mapped_count, unresolved_count, database_deletions = 0, 0, 0
    for entry in audit_data["rows"]:
        path = construct_path.parent / entry["report"]
        if digest(path) != audit_data["reports_sha256"][entry["report"]]:
            raise ValueError("Per-entry construct report changed")
        data = json.loads(path.read_text())
        source = Path(data["original_source_path"])
        if digest(source) != data["original_source_sha256"]:
            raise ValueError("Original deposited mmCIF changed")
        block = gemmi.cif.read(str(source)).sole_block()
        scheme = category_rows(block, "_pdbx_poly_seq_scheme.")
        chains, gap_rows = [], []
        for chain in data["chains"]:
            mapped_candidates = [[map_difference(row, candidate, scheme)
                for row in candidate["deposited_sequence_differences"]] for candidate in chain["source_label_candidates"]]
            mapped_count += sum(row["label_index_zero_based"] is not None for records in mapped_candidates for row in records)
            unresolved_count += sum(row["mapping_status"] == "unresolved" for records in mapped_candidates for row in records)
            database_deletions += sum(row["mapping_status"] == "not-applicable-database-deletion"
                                      for records in mapped_candidates for row in records)
            gaps = [gap_evidence(gap, mapped_candidates) for gap in chain["gap_intervals_zero_based_half_open"]]
            chains.append({"output_chain": chain["output_chain"], "mapped_candidate_annotations": mapped_candidates, "gaps": gaps})
            gap_rows.extend(gaps)
        categories.update(row["coverage"] for row in gap_rows)
        save(output / "entries" / f"{entry['entry_id']}.json", {"entry_id": entry["entry_id"], "scope": entry["scope"],
            "chains": chains, "scientific_admission": False, "original_request_unchanged": True})
        rows.append({"entry_id": entry["entry_id"], "scope": entry["scope"],
            "gap_annotation_counts": dict(Counter(row["coverage"] for row in gap_rows)),
            "gap_count": len(gap_rows), "has_consensus_artifact_overlap": any(row["consensus_explicit_artifact_positions"] for row in gap_rows),
            "has_completely_artifact_annotated_gap": any(row["coverage"] == "complete-explicit-artifact-annotation" for row in gap_rows),
            "report": f"entries/{entry['entry_id']}.json"})
        bindings[str(path)] = digest(path)
        bindings[str(source)] = digest(source)
    summary = {"entries": len(rows), "gaps": sum(row["gap_count"] for row in rows), "gap_annotation_counts": dict(categories),
        "entries_with_consensus_explicit_artifact_overlap": sum(row["has_consensus_artifact_overlap"] for row in rows),
        "entries_with_completely_artifact_annotated_gaps": sum(row["has_completely_artifact_annotated_gap"] for row in rows),
        "candidate_annotation_records_mapped_including_assembly_reuse": mapped_count,
        "candidate_annotation_records_unresolved_including_assembly_reuse": unresolved_count,
        "candidate_database_deletions_without_construct_position_including_assembly_reuse": database_deletions,
        "database_deletions_are_not_missing_coordinate_targets": True,
        "counts_are_nested_under_entries_not_independent_structures": True,
        "absence_of_mapping_not_no_tag_evidence": True, "automatic_exclusions_or_admissions": 0,
        "new_sources_or_restoration_runs": 0}
    save(output / "difference-audit.json", {"summary": summary, "rows": rows, "bindings": bindings,
        "reports_sha256": {row["report"]: digest(output / row["report"]) for row in rows},
        "auditor_sha256": digest(Path(__file__))})
    save(output / "summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--construct-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.construct_audit.resolve(), args.output.resolve())
