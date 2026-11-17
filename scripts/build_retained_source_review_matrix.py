#!/usr/bin/env python3
"""Read-only all-source review facts; never a loop cutoff or scientific admission."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_real_gap_cohort_capacity import task_inventory
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"
SCOPE = ROOT / "current-loop-scope-attempt02/loop-scope-inventory.json"
CAPACITY = ROOT / "cohort-capacity-readback-attempt03/capacity-inventory.json"
GROUPS = ROOT / "union-candidate-groups-attempt01/grouped-census.json"
RECENT_EXPOSURE = frozenset({"1NON", "2BE1"})
OLDER_CONTROLS_OUTSIDE_CENSUS = ("1UBQ", "1FR2", "1A3N", "1PGA", "2CI2")


def review_topics(row, omitted, known_component_exposure):
    topics = []
    if row["current_scope_status"] == "excluded_by_user":
        return ["explicitly-excluded-retain-research-only"]
    if row["small_assembly_control"]:
        topics.append("small-control-not-primary-cohort")
    topics.append("loop-versus-domain-and-construct-adjudication")
    if row["terminal_construct_review_required"]:
        topics.append("terminal-construct-identity")
    if set(omitted) - {"HOH", "WAT"}:
        topics.append("omitted-context-other-than-HOH-WAT")
    if row["prepared_source_derivative"]:
        topics.append("prepared-source-provenance")
    if not row["group_membership_complete"]:
        topics.append("incomplete-family-membership")
    if known_component_exposure:
        topics.append("known-development-component-no-final-assignment")
    topics.append("whole-context-and-backend-admission")
    return topics


def worker(output):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("source review requires fresh guarded output")
    scope, capacity = json.loads(SCOPE.read_text()), json.loads(CAPACITY.read_text())
    if scope["source_capacity_sha256"] != digest(CAPACITY):
        raise ValueError("scope/capacity evidence changed")
    groups = {row["entry_id"]: row for row in json.loads(GROUPS.read_text())["rows"]}
    exposed = {}
    for group in groups.values():
        if group["previously_inspected"] or group["entry_id"] in RECENT_EXPOSURE:
            exposed.setdefault(group["conservative_entry_component"], []).append(group["entry_id"])
    capacity_rows = {row["entry_id"]: row for row in capacity["rows"]}
    report = output / "report"
    report.mkdir(exist_ok=False)
    (report / "entries").mkdir()
    bindings = {str(p): digest(p) for p in (SCOPE, CAPACITY, GROUPS)}
    rows, seen, total_gaps = [], set(), 0
    for row in scope["rows"]:
        entry = row["entry_id"]
        if entry in seen or entry not in groups or entry not in capacity_rows:
            raise ValueError("review sources have duplicate/unmapped entries")
        seen.add(entry)
        request_path = Path(row["request_path"])
        if digest(request_path) != row["request_sha256"] or capacity_rows[entry]["request_sha256"] != row["request_sha256"]:
            raise ValueError("source-bound request changed")
        result_path, projection_path = request_path.parent / "result.json", request_path.parent / "projection.json"
        result, projection = json.loads(result_path.read_text()), json.loads(projection_path.read_text())
        source = Path(result["source_input"])
        request = DiffusionRequest.from_json(request_path.read_text())
        if digest(source) != result["source_input_sha256"] or digest(source) != request.normalized_pdb.sha256:
            raise ValueError("prepared source differs from its request")
        actual = task_inventory(request)
        if any(actual[key] != row[key] for key in actual):
            raise ValueError("retained request inventory disagrees with source scope")
        group = groups[entry]
        if row["union_component"] != group["conservative_entry_component"]:
            raise ValueError("union component changed")
        points, _ = _coordinates(source.read_text())
        if not set(request.fixed_atoms) <= points.keys():
            raise ValueError("source omitted a fixed observation")
        gaps = []
        for index, gap in enumerate(request.gaps):
            anchors = []
            for residue in (gap.left_anchor, gap.right_anchor):
                if residue is not None:
                    atom = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, "CA")
                    if atom not in points:
                        raise ValueError("declared anchor lacks observed CA")
                    anchors.append(points[atom])
            gaps.append({"gap_index": index, **actual["gaps"][index],
                "left_anchor": asdict(gap.left_anchor) if gap.left_anchor else None,
                "right_anchor": asdict(gap.right_anchor) if gap.right_anchor else None,
                "anchor_CA_distance_angstrom": float(np.linalg.norm(anchors[0]-anchors[1])) if len(anchors) == 2 else None,
                "distance_or_length_does_not_establish_loop_admission": True})
        component_exposure = sorted(exposed.get(row["union_component"], []))
        omitted = projection["omitted_heterogen_atom_counts"]
        detailed = {"entry_id": entry, "title": group["title"], "current_scope_status": row["current_scope_status"],
            **actual, "gaps": gaps, "omitted_heterogen_atom_counts": omitted,
            "union_component": row["union_component"], "group_membership_complete": row["group_membership_complete"],
            "known_inspected_entries_in_component": component_exposure,
            "known_component_exposure_is_a_lower_bound_not_untouched_claim": True,
            "prepared_source_derivative": row["prepared_source_derivative"],
            "anchor_reference": "exact request-bound input; derivative if preparation changed source",
            "sequence_hashes_by_case_sensitive_chain": {t.chain: hashlib.sha256(t.sequence.encode()).hexdigest() for t in request.target_sequences},
            "retained_preparation_profile_scope_only": result.get("existing_public_scope"),
            "review_topics": review_topics(row, omitted, bool(component_exposure)),
            "scientific_test_admitted": False, "final_split_assignment": None, "new_restoration_calculations": False,
            "source_path": str(source), "source_sha256": digest(source),
            "request_path": str(request_path), "request_sha256": row["request_sha256"]}
        save(report / "entries" / f"{entry}.json", detailed)
        rows.append({key: detailed[key] for key in ("entry_id", "title", "current_scope_status", "target_residues",
            "observed_residues", "chains", "longest_gap", "union_component", "group_membership_complete",
            "prepared_source_derivative", "known_inspected_entries_in_component", "omitted_heterogen_atom_counts", "review_topics")}
            | {"gap_count": len(gaps), "entry_report": f"entries/{entry}.json"})
        total_gaps += len(gaps)
        for p in (request_path, result_path, projection_path, source):
            bindings[str(p)] = digest(p)
    if len(rows) != len(capacity_rows):
        raise ValueError("scope does not cover every retained ready source")
    from collections import Counter

    summary = {"ready_entries_reviewed": len(rows), "declared_gaps": total_gaps,
        "scope_counts": dict(Counter(row["current_scope_status"] for row in rows)),
        "entries_with_context_omission_review": sum(bool(set(row["omitted_heterogen_atom_counts"])-{"HOH", "WAT"}) for row in rows),
        "omission_rule": "components outside retained HOH/WAT whitelist need review; not a chemical classification",
        "entries_in_known_exposed_components": sum(bool(row["known_inspected_entries_in_component"]) for row in rows),
        "review_topic_counts": dict(Counter(topic for row in rows for topic in row["review_topics"])),
        "older_inspected_controls_missing_from_census": [entry for entry in OLDER_CONTROLS_OUTSIDE_CENSUS if entry not in groups],
        "exposure_overlay_not_complete_historical_leakage_qualification": True,
        "scientific_test_admitted": 0, "split_frozen": False, "recruitment_or_new_restoration_calculations": 0,
        "new_universal_numeric_gap_cutoff": None, "historical_quotas_reinstated": False}
    columns = ("entry_id", "current_scope_status", "target_residues", "observed_residues", "chains", "gap_count",
               "longest_gap", "union_component", "group_membership_complete", "prepared_source_derivative", "review_topics", "title")
    with (report / "review-matrix.tsv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({**row, "review_topics": ";".join(row["review_topics"])} for row in rows)
    save(report / "inventory.json", {"rows": rows, "summary": summary, "bindings": bindings,
        "reader_sha256": digest(Path(__file__)), "source_readiness_not_scientific_admission": True})
    save(report / "summary.json", summary)
    save(output / "result.json", {"status": "retained-source-review-written", **summary})
    print(json.dumps(summary), flush=True)


def declare(output):
    output.mkdir(parents=True, exist_ok=False)
    job = {"case_id": "retained-source-review", "case": str(SCOPE), "seed": 0,
        "command": [sys.executable, "-m", "scripts.build_retained_source_review_matrix", "worker", "--output",
                    str(output / "execution/retained-source-review-seed-0")]}
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(p): digest(p) for p in (Path(__file__), SCOPE, CAPACITY, GROUPS,
            PROJECT / "scripts/audit_real_gap_cohort_capacity.py", PROJECT / "scripts/run_mac_modeller_comparison.py")},
        "jobs": [job], "count": 1, "purpose": "read-only review facts for all194 retained ready inputs, not admission or modeling",
        "new_model_calls": 0, "new_sources": 0, "sequences_modified": False,
        "resource_guard": "existing serial CPU1800s/RSS8GiB/host pressure/evidence cap"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        worker(args.output.resolve())
