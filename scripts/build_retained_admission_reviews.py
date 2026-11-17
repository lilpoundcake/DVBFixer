#!/usr/bin/env python3
"""Join retained source facts into review questions; never approve, trim or resample."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save


def review_questions(scope, *, terminal_gaps, other_context, water_only, artifact_overlap,
                     complete_artifact_gap, unresolved_annotations, prepared, group_complete, exposed, database_deletions=0):
    if scope == "excluded_by_user":
        return ["retain-explicit-exclusion-no-active-task"]
    questions = []
    if scope == "small_control_only":
        questions.append("control-only-not-primary-cohort")
    if artifact_overlap:
        questions.append("decide-task-semantics-for-explicit-deposited-construct-artifacts")
    if complete_artifact_gap:
        questions.append("whole-gap-artifact-annotation-no-silent-mask-trimming")
    if terminal_gaps:
        questions.append("adjudicate-terminal-construct-not-just-gap-length")
    if unresolved_annotations:
        questions.append("resolve-unmapped-difference-annotations")
    if database_deletions:
        questions.append("retain-database-deletions-as-outside-deposited-coordinate-task")
    if other_context:
        questions.append("retain-and-assess-non-HOH-WAT-context-in-complete-endpoint")
    elif water_only:
        questions.append("declare-water-context-treatment")
    if prepared:
        questions.append("retain-prepared-source-reference-and-ownership")
    if not group_complete:
        questions.append("complete-family-grouping-before-split")
    if exposed:
        questions.append("known-development-component-not-untouched-final-evaluation")
    questions += ["adjudicate-loop-versus-domain-in-observed-scaffold", "backend-and-whole-context-admission-pending"]
    return questions


def checked(path):
    data = json.loads(path.read_text())
    for name, expected in data.get("bindings", {}).items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Bound retained evidence changed: {name}")
    for name, expected in data.get("reports_sha256", {}).items():
        report = path.parent / name
        if not report.resolve().is_relative_to(path.parent.resolve()) or digest(report) != expected:
            raise ValueError("Per-entry report escaped its root or changed")
    return data


def build(inventory_path, context_path, construct_path, differences_path, output):
    inventory, context, construct, differences = map(checked, (inventory_path, context_path, construct_path, differences_path))
    datasets = [{row["entry_id"]: row for row in data["rows"]} for data in (context, construct, differences)]
    entries = {row["entry_id"] for row in inventory["rows"]}
    if any(set(dataset) != entries for dataset in datasets) or len(entries) != len(inventory["rows"]):
        raise ValueError("Joined reports must cover exactly the same distinct retained entries")
    output.mkdir(parents=True, exist_ok=False)
    rows, bindings = [], {str(path): digest(path) for path in (inventory_path, context_path, construct_path, differences_path)}
    for original in inventory["rows"]:
        entry = original["entry_id"]
        ctx, cons, dif = [data[entry] for data in datasets]
        if len({original["current_scope_status"], ctx["scope"], cons["scope"], dif["scope"]}) != 1:
            raise ValueError("A joined report silently changed the entry scope")
        if len({original["gap_count"], ctx["gap_count"], cons["gaps_checked"], dif["gap_count"]}) != 1:
            raise ValueError("A joined report changed the gap denominator")
        dif_path = differences_path.parent / dif["report"]
        details = json.loads(dif_path.read_text())
        unresolved = sum(row.get("mapping_status", "unresolved" if row["label_index_zero_based"] is None else "mapped-deposited-residue") == "unresolved"
                         for chain in details["chains"]
                         for candidate in chain["mapped_candidate_annotations"] for row in candidate)
        database_deletions = sum(row.get("mapping_status") == "not-applicable-database-deletion" for chain in details["chains"]
                                 for candidate in chain["mapped_candidate_annotations"] for row in candidate)
        questions = review_questions(original["current_scope_status"], terminal_gaps=cons["terminal_gap_count"],
            other_context=ctx["outside_HOH_WAT_instances"] > 0, water_only=ctx["omitted_atoms"] > 0,
            artifact_overlap=dif["has_consensus_artifact_overlap"], complete_artifact_gap=dif["has_completely_artifact_annotated_gap"],
            unresolved_annotations=unresolved, database_deletions=database_deletions, prepared=original["prepared_source_derivative"],
            group_complete=original["group_membership_complete"], exposed=bool(original["known_inspected_entries_in_component"]))
        rows.append({"entry_id": entry, "scope": original["current_scope_status"], "component": original["union_component"],
            "gaps": original["gap_count"], "terminal_gaps": cons["terminal_gap_count"],
            "omitted_context_atoms": ctx["omitted_atoms"], "outside_HOH_WAT_instances": ctx["outside_HOH_WAT_instances"],
            "artifact_gap_annotation_counts": dif["gap_annotation_counts"], "unresolved_candidate_annotations": unresolved,
            "database_deletion_annotations_not_coordinate_targets": database_deletions,
            "questions": questions, "deposited_gap_axis_fact": "verified-for-all-source-label-candidates",
            "original_projection_observation_custody_fact": "verified",
            "actual_context_omission_declaration_fact": "verified",
            "scientific_admission": "not-admitted", "split_assignment": "unassigned",
            "reviewer": "", "decision_evidence": "",
            "evidence_paths": [str(context_path.parent / ctx["report"]), str(construct_path.parent / cons["report"]), str(dif_path)]})
        for path in rows[-1]["evidence_paths"]:
            bindings[path] = digest(Path(path))
    rows.sort(key=lambda row: (row["component"], row["entry_id"]))
    columns = tuple(rows[0])
    with (output / "case-review.tsv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows({**row, "questions": ";".join(row["questions"]),
                         "evidence_paths": ";".join(row["evidence_paths"])} for row in rows)
    summary = {"entries": len(rows), "gaps": sum(row["gaps"] for row in rows),
        "scope_counts": dict(Counter(row["scope"] for row in rows)),
        "question_counts": dict(Counter(q for row in rows for q in row["questions"])),
        "facts_verified_per_entry": ["deposited-gap-axis", "original-projection-custody", "actual-context-omissions"],
        "admissions_or_exclusions_added": 0, "split_frozen": False, "request_masks_changed": False,
        "new_sources_or_restoration_calculations": 0,
        "ordering": "existing conservative component then entry; not a performance/length ranking"}
    save(output / "case-reviews.json", {"summary": summary, "rows": rows, "bindings": bindings,
        "builder_sha256": digest(Path(__file__)), "TSV_sha256": digest(output / "case-review.tsv")})
    save(output / "summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--context-audit", type=Path, required=True)
    parser.add_argument("--construct-audit", type=Path, required=True)
    parser.add_argument("--difference-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.inventory.resolve(), args.context_audit.resolve(), args.construct_audit.resolve(),
          args.difference_audit.resolve(), args.output.resolve())
