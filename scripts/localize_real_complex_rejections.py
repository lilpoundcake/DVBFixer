#!/usr/bin/env python3
"""Attribute retained final failures to internal, N-terminal and C-terminal gaps.

Counts gap x seed observations, NOT independently accepted per-gap candidates.
Uses detailed last-stage findings plus COMPLETE independent drift warnings;
the GUI's truncated finding previews and runner-authored warning copies are not
used for drift attribution. No coordinates or acceptance outcomes are changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path


def identity(value):
    return value["chain"], value["residue_number"], value["insertion_code"]


def main():
    from dvbfixer.model.diffusion.anchor_acceptance import _neighborhood
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.protein_context import independent_gap_requests

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    cohort = json.loads((args.root / "cohort-final-outcomes.json").read_text())
    stages = json.loads((args.root / "failure-analysis-2026-11-09/summary.json").read_text())
    stage_cases = {case["pdb"]: case for case in stages["cases"]}
    strata = defaultdict(Counter)
    combinations = Counter()
    candidates, gap_rows = [], []
    for case in cohort["rows"]:
        pdb = case["pdb"]
        raw = (args.root / case["job_record"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != case["job_record_sha256"]:
            raise ValueError("cohort job record changed")
        job = json.loads(raw)
        execution = args.root / ("1N8Z-managed-followup" if pdb == "1N8Z" else "managed-execution")
        workspace = execution / "data/projects" / job["workspaceId"]
        report_file = workspace / job["diffusion"]["report"]
        expected = next(r["sha256"] for r in job["provenance"]["outputs"] if r["file"] == job["diffusion"]["report"])
        if hashlib.sha256(report_file.read_bytes()).hexdigest() != expected:
            raise ValueError("published full seed outcome report changed")
        outcomes = {r["seed"]: r for r in json.loads(report_file.read_text())["outcomes"]}
        request = DiffusionRequest.from_json((args.root / "cases" / pdb / "request-before-preparation.json").read_text())
        detail = stage_cases[pdb]
        for gap, indexed in zip(request.gaps, detail["gaps"], strict=True):
            if (gap.chain != indexed["chain"] or len(gap.generated_residues) != indexed["length"]
                    or asdict(gap.generated_residues[0]) != indexed["first"]
                    or asdict(gap.generated_residues[-1]) != indexed["last"]):
                raise ValueError("stage report gap inventory differs from source-bound plan")
        expected_reports = {entry["archive_report"]: entry for entry in detail["seeds"]}
        reports = {}
        with tarfile.open(workspace / job["outputDir"] / "diagnostics.tar.gz", "r:gz") as archive:
            for member in archive:
                if member.name not in expected_reports:
                    continue
                entry = expected_reports[member.name]
                if not member.isfile() or member.size > 16_000_000:
                    raise ValueError("unbounded or nonregular stage report")
                data = archive.extractfile(member).read()
                if hashlib.sha256(data).hexdigest() != entry["report_sha256"]:
                    raise ValueError("archived stage report digest mismatch")
                if entry["seed"] in reports:
                    raise ValueError("duplicate archived stage report")
                reports[entry["seed"]] = json.loads(data)
        if set(reports) != {7, 101, 102, 103, 104}:
            raise ValueError("stage reports missing")
        generated = {identity(asdict(residue)): index for index, gap in enumerate(request.gaps) for residue in gap.generated_residues}
        anchors = defaultdict(set)
        labels = defaultdict(set)
        neighborhoods = []
        for index, single in enumerate(independent_gap_requests(request)):
            neighborhoods.append({identity(asdict(r)) for r in _neighborhood(single)})
            for anchor in (single.gaps[0].left_anchor, single.gaps[0].right_anchor):
                if anchor is not None:
                    anchors[identity(asdict(anchor))].add(index)
                    labels[(anchor.chain, anchor.residue_number + anchor.insertion_code)].add(index)
        per_gap = [{"pdb": pdb, **g, "geometry_problem_seeds": [], "anchor_drift_seeds": [],
                    "any_local_blocker_seeds": [], "geometry_categories": Counter(), "locations": Counter(),
                    "drift_findings": []} for g in detail["gaps"]]
        for entry in detail["seeds"]:
            seed = entry["seed"]
            last = entry["stages"][-1]
            # Re-read full stage details, verified by the previous archive audit.
            # Its per-gap summary is insufficient to distinguish outer-neighbor
            # issues from defects inside newly generated residues.
            report = reports[seed]
            assessment = report["attempts"][-1]["assessment"]
            if assessment is None:
                raise ValueError("cannot localize a missing final stage assessment")
            problems = defaultdict(list)
            unlocated = []
            for issue in assessment["issues"]:
                if not issue.get("actionable", True):
                    continue
                residues = {identity(r) for r in issue["residues"]}
                if not residues:
                    # Aggregate validator gates duplicate detailed findings and
                    # have no anatomical ownership. Keep them, never guess a gap.
                    unlocated.append(issue["category"])
                    continue
                generated_gaps = {generated[r] for r in residues if r in generated}
                affected = generated_gaps or {i for i, neighborhood in enumerate(neighborhoods) if residues & neighborhood}
                if not affected:
                    raise ValueError("actionable detailed finding outside all canonical gap neighborhoods")
                for index in affected:
                    gap_anchors = {r for r, indexes in anchors.items() if index in indexes}
                    if issue["category"] == "clash":
                        location = "observed_only_clash" if not generated_gaps else (
                            "inter_gap_clash" if len(generated_gaps) > 1 else (
                                "generated_observed_clash" if residues - generated.keys() else "within_generated_clash"))
                    elif generated_gaps:
                        location = "junction" if residues & gap_anchors else "inside_generated"
                    else:
                        location = "observed_anchor" if residues & gap_anchors else "observed_outer_neighbor"
                    problems[index].append({"category": issue["category"], "location": location, "message": issue["message"]})
            drift = defaultdict(list)
            pattern = re.compile(r"^Final validation: ([^:]):([^:]+): (?:backbone RMSD|Cα displacement|maximum heavy-atom displacement) = .+maximum allowed:")
            for warning in outcomes[seed]["validation"]["warnings"]:
                match = pattern.match(warning)
                if match:
                    affected = labels.get((match[1], match[2]))
                    if not affected:
                        raise ValueError("independent drift warning is not a declared gap anchor")
                    for index in affected:
                        drift[index].append(warning)
            if ("anchor-repair-per-anchor-movement-limit" in entry["gates"]) != bool(drift):
                raise ValueError("independent movement gate and anatomical attribution disagree")
            affected = set(problems) | set(drift)
            if entry["status"] == "accepted" and affected:
                raise ValueError("accepted candidate has a localized final blocker")
            if entry["status"] == "rejected" and not affected:
                raise ValueError("rejected candidate has no anatomically attributable blocker")
            kinds = sorted({detail["gaps"][index]["kind"] for index in affected})
            if entry["status"] == "rejected":
                group = "mixed_internal_and_terminal" if "internal" in kinds and len(kinds) > 1 else (
                    "internal_only" if kinds == ["internal"] else "terminal_only")
                combinations[group] += 1
            for index, gap in enumerate(per_gap):
                counter = strata[gap["kind"]]
                counter["gap_seed_observations"] += 1
                counter["geometry_problem_observations"] += index in problems
                counter["anchor_drift_observations"] += index in drift
                counter["any_local_blocker_observations"] += index in affected
                if index in problems:
                    gap["geometry_problem_seeds"].append(seed)
                    gap["geometry_categories"].update(p["category"] for p in problems[index])
                    gap["locations"].update(p["location"] for p in problems[index])
                if index in drift:
                    gap["anchor_drift_seeds"].append(seed)
                    gap["drift_findings"].extend(drift[index])
                if index in affected:
                    gap["any_local_blocker_seeds"].append(seed)
            candidates.append({"pdb": pdb, "seed": seed, "status": entry["status"],
                "gap_kinds_with_blockers": kinds, "geometry_gap_indices": sorted(problems),
                "drift_gap_indices": sorted(drift), "aggregate_unlocalized_duplicate_gates": unlocated,
                "last_stage": last["stage"]})
        for gap in per_gap:
            strata[gap["kind"]]["distinct_chain_specific_gaps"] += 1
            strata[gap["kind"]]["distinct_gaps_with_any_blocker"] += bool(gap["any_local_blocker_seeds"])
            gap_rows.append(gap)
    result = {"scope": "Final local blockers, including independently reported observed-anchor motion. No new sampling or minimization. Lack of a local blocker is NOT independent per-gap publication acceptance.",
        "strata": dict(strata), "rejected_candidate_localization": dict(combinations),
        "candidates": candidates, "gaps": gap_rows}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False)
        stream.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k in ("strata", "rejected_candidate_localization")}, indent=2))


if __name__ == "__main__":
    main()
