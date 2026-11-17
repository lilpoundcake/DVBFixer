#!/usr/bin/env python3
"""Source-only construct-target reviews; no backend admission or restoration."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import gemmi
from scripts.audit_retained_construct_annotations import category_rows
from scripts.build_retained_admission_reviews import checked
from scripts.build_retained_gap_reviews import unique
from scripts.mac_reconstruction_comparison import digest, save

DECISIONS = {
    "1AHW": "The deposited target already specifies the extracellular tissue-factor/Fab construct. Both tissue-factor copies have local internal absence plus short terminal absence within that deposited target, not a request to rebuild the full-length native membrane protein.",
    "1JPS": "The deposited target is the Fab/extracellular tissue-factor construct itself. All intended absent coordinates remain within those declared chains; database annotations are not expanded into native full-length antibody or membrane-protein targets.",
    "1KTZ": "The deposited target comprises TGF-beta3 and the receptor extracellular construct. The target is not expanded to the full precursor or full-length receptor; both internal gaps and terminal coordinate absences remain explicit parts of the original construct request.",
}


def target_facts(request, construct, block):
    targets = unique(request["target_sequences"], "chain")
    chains = unique(construct["chains"], "output_chain")
    if targets.keys() != chains.keys():
        raise ValueError("Request/construct chain denominator differs")
    if any(gap["chain"] not in targets for gap in request["gaps"]):
        raise ValueError("Gap references an unknown case-sensitive chain")
    entities = unique(category_rows(block, "_entity_poly."), "entity_id")
    atoms = category_rows(block, "_atom_site.")
    models = {row["pdbx_PDB_model_num"] for row in atoms}
    if models != {"1"}:
        raise ValueError("This bounded review requires a single deposited model1")
    result = []
    for chain, target in targets.items():
        cons = chains[chain]
        spans = sorted((g["target_interval"]["start"], g["target_interval"]["stop"])
                       for g in request["gaps"] if g["chain"] == chain)
        retained_spans = sorted((g["start"], g["stop"]) for g in cons["gap_intervals_zero_based_half_open"])
        if spans != retained_spans or len(spans) != len(set(spans)):
            raise ValueError("Construct gap axes changed")
        if any(not 0 <= start < stop <= len(target["sequence"]) for start, stop in spans):
            raise ValueError("Invalid target interval")
        requested_missing = {i for start, stop in spans for i in range(start, stop)}
        if len(requested_missing) != sum(stop - start for start, stop in spans):
            raise ValueError("Overlapping target masks")
        candidates = []
        seen_labels = set()
        for candidate in cons["source_label_candidates"]:
            entity_id, label = candidate["entity_id"], candidate["source_label_chain"]
            if label in seen_labels:
                raise ValueError("Repeated source-label candidate")
            seen_labels.add(label)
            entity = entities[entity_id]
            deposited = "".join(entity["pdbx_seq_one_letter_code_can"].split())
            if entity["type"] != "polypeptide(L)" or deposited != target["sequence"]:
                raise ValueError("Target differs from deposited canonical construct sequence")
            observed = set()
            for atom in atoms:
                if atom["label_asym_id"] != label or atom["type_symbol"] in {"H", "D"}:
                    continue
                if atom["label_entity_id"] != entity_id:
                    raise ValueError("Atom/entity identity mismatch")
                index = int(atom["label_seq_id"]) - 1
                if not 0 <= index < len(deposited):
                    raise ValueError("Observed label index outside deposited target")
                observed.add(index)
            if set(range(len(deposited))) - observed != requested_missing:
                raise ValueError("Request masks differ from deposited heavy-coordinate absence")
            candidates.append({"source_label_chain": label, "entity_id": entity_id,
                "deposited_target_length": len(deposited), "deposited_observed_residues": len(observed),
                "deposited_sequence_sha256": hashlib.sha256(deposited.encode()).hexdigest(),
                "entity_annotations": candidate["entity_annotations"],
                "raw_database_alignment_annotations": candidate["database_alignments"],
                "database_reference_outside_target_not_requested": True})
        if not candidates:
            raise ValueError("No compatible deposited chain candidate")
        result.append({"output_chain": chain, "target_length": len(target["sequence"]),
            "source_label_assignment_unique": len(candidates) == 1, "source_candidates": candidates,
            "all_source_candidates_match_deposited_sequence_and_absence": True,
            "all_original_gap_intervals": [{"start": start, "stop": stop} for start, stop in spans]})
    return result


def review(gaps_path, case_path, connections_path, output):
    gaps, cases, connections = map(checked, (gaps_path, case_path, connections_path))
    case_rows, connection_rows = [unique(data["rows"], "entry_id") for data in (cases, connections)]
    bindings = {str(path): digest(path) for path in (gaps_path, case_path, connections_path)}
    detail_paths = {}
    for name in gaps["bindings"]:
        path = Path(name)
        if path.parent.name != "entries" or path.suffix != ".json":
            continue
        data = json.loads(path.read_text())
        if "sequence_hashes_by_case_sensitive_chain" in data and data["entry_id"] in DECISIONS:
            if data["entry_id"] in detail_paths:
                raise ValueError("Duplicate bound entry detail")
            detail_paths[data["entry_id"]] = path
    if detail_paths.keys() != DECISIONS.keys():
        raise ValueError("Bounded review cases are missing")
    rows = []
    for entry, reason in DECISIONS.items():
        detail_path = detail_paths[entry]
        detail = json.loads(detail_path.read_text())
        case = case_rows[entry]
        if case["scope"] != "loop_scope_and_context_review_pending" or not detail["known_inspected_entries_in_component"]:
            raise ValueError("Reviewed case is not an existing exposed development candidate")
        context_path, construct_path, difference_path = map(Path, case["evidence_paths"])
        for path in (context_path, construct_path, difference_path):
            if digest(path) != cases["bindings"][str(path)]:
                raise ValueError("Joined evidence changed")
        context, construct = [json.loads(path.read_text()) for path in (context_path, construct_path)]
        if any(data["entry_id"] != entry for data in (context, construct, detail)):
            raise ValueError("Reviewed evidence entry identity changed")
        source_path, request_path = Path(construct["original_source_path"]), Path(detail["request_path"])
        if digest(source_path) != construct["original_source_sha256"] or digest(request_path) != detail["request_sha256"]:
            raise ValueError("Source/request custody changed")
        block = gemmi.cif.read(str(source_path)).sole_block()
        facts = target_facts(json.loads(request_path.read_text()), construct, block)
        connection_path = connections_path.parent / connection_rows[entry]["report"]
        if not connection_path.resolve().is_relative_to(connections_path.parent.resolve()):
            raise ValueError("Connection report escapes its root")
        connection_data = json.loads(connection_path.read_text())
        if connection_data["entry_id"] != entry:
            raise ValueError("Connection report belongs to another entry")
        connection = connection_data["actual_projection"]
        if not connection["resolved_retained_pair_custody_verified"] or connection["resolved_projection_boundary_pairs"]:
            raise ValueError("Unexpected explicit-link custody finding")
        names = context["omitted_atom_counts_by_residue_name"]
        if set(names) - {"HOH", "WAT"}:
            raise ValueError("This bounded review does not resolve nonwater chemistry")
        rows.append({"entry_id": entry, "reviewer": "OpenCode source-only construct review",
            "task_semantic_review": "supported-within-deposited-construct-not-full-native-sequence",
            "case_specific_basis": reason, "chains": facts,
            "all_original_gaps": detail["gaps"], "target_residue_instances": detail["target_residues"],
            "observed_residue_instances": detail["observed_residues"],
            "prepared_source_derivative": detail["prepared_source_derivative"],
            "full_context_omission_counts": names,
            "water_treatment": "No omitted atoms" if not names else "Original water context retained; endpoint water treatment remains to be declared, not assumed irrelevant",
            "no_annotation_does_not_prove_native_sequence_or_absence_of_tags": True,
            "database_only_deletions_not_coordinate_targets": case["database_deletion_annotations_not_coordinate_targets"],
            "known_development_component_exposure": detail["known_inspected_entries_in_component"],
            "recorded_preparation_scope_only_not_new_backend_test": detail["retained_preparation_profile_scope_only"],
            "remaining_gates": ["terminal completion is distinct from two-sided internal-loop reconstruction",
                "backend-specific full-request admission; do not bypass recorded terminal refusals",
                "complete-context endpoint assessment and declared water treatment where applicable",
                "prospective method/context/seed/resource manifest before any new structural-model call"],
            "original_effective_scope_unchanged": "loop_scope_and_context_review_pending",
            "scientific_admission": "not-admitted", "final_split_assignment": "unassigned",
            "new_source_or_restoration_or_optimizer_calls": 0,
            "request_mask_or_sequence_changed": False, "new_narrower_task_created": False})
        for path in (detail_path, context_path, construct_path, difference_path, source_path, request_path, connection_path):
            bindings[str(path)] = digest(path)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "construct-target-reviews.json", {"rows": rows, "bindings": bindings,
        "reviewer_sha256": digest(Path(__file__)), "case_specific_task_semantic_reviews": len(rows),
        "scientific_admissions_added": 0, "scope_overlay_changed": False,
        "new_structural_model_calls": 0, "limits": "Construct-target semantics only, not native-sequence identity, backend readiness or final scientific qualification"})
    print(json.dumps({"case_specific_task_semantic_reviews": len(rows), "scientific_admissions_added": 0}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaps", type=Path, required=True)
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--connections", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    review(args.gaps.resolve(), args.cases.resolve(), args.connections.resolve(), args.output.resolve())
