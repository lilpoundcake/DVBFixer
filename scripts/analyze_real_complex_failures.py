#!/usr/bin/env python3
"""Read-only stage and per-gap attribution for the retained ten-complex cohort.

Verifies archived report digests against run-summary.json. No extraction into
the run, coordinate mutation, refinement, sampling or acceptance-policy change.
Counts are overlapping: one candidate can fail several independent gates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
from collections import Counter
from pathlib import Path


def identity(value):
    return value["chain"], value["residue_number"], value["insertion_code"]


def stage_summary(attempt, generated, anchors):
    assessment = attempt.get("assessment")
    if assessment is None:
        return {"stage": attempt["stage"], "error": attempt.get("error"),
                "failure_kind": attempt.get("failure_kind"), "geometry_passed": None}
    issues = [i for i in assessment["issues"] if i.get("actionable", True)]
    counts = Counter(i["category"] for i in issues)
    by_gap = {}
    coupling = Counter()
    for issue in issues:
        residues = {identity(r) for r in issue["residues"]}
        affected = {generated[r] for r in residues if r in generated}
        if len(affected) > 1:
            coupling["different_generated_gaps"] += 1
        elif affected and residues - generated.keys():
            coupling["generated_and_observed"] += 1
        elif affected:
            coupling["same_generated_gap"] += 1
        else:
            coupling["observed_only"] += 1
        affected |= {g for r in residues for g in anchors.get(r, ())}
        for gap in affected:
            by_gap.setdefault(gap, Counter())[issue["category"]] += 1
    return {"stage": attempt["stage"], "geometry_passed": not issues,
            "admission_failures": assessment["admission_failures"],
            "issues": dict(counts), "issue_count": len(issues),
            "per_gap_issues": {str(k): dict(v) for k, v in by_gap.items()},
            "residue_coupling": dict(coupling), "error": attempt.get("error"),
            "restraint_residues": sorted({identity(r["atom"]) for r in attempt.get("restraints", [])})}


def analyze_case(root, case):
    job_path = root / case["job_record"]
    job = json.loads(job_path.read_text())
    execution = root / ("1N8Z-managed-followup" if case["pdb"] == "1N8Z" else "managed-execution")
    run = execution / "data/projects" / job["workspaceId"] / job["outputDir"]
    summary = json.loads((run / "run-summary.json").read_text())
    expected = {row["file"]: row for row in summary["diagnostics"]["files"]}
    reports, request = {}, None
    request_digest = None
    with tarfile.open(run / "diagnostics.tar.gz", "r:gz") as archive:
        for member in archive:
            match = re.search(r"/seed-(\d+)/refined/refinement-report\.json$", member.name)
            is_request = bool(re.fullmatch(r"\.dvbfixer-diffusion-run\.[^/]+/request\.json", member.name))
            if not match and not is_request:
                continue
            if not member.isfile() or member.size > 16_000_000:
                raise ValueError("report/request archive member is not a bounded regular file")
            data = archive.extractfile(member).read()
            digest = hashlib.sha256(data).hexdigest()
            if expected[member.name]["sha256"] != digest or expected[member.name]["bytes"] != len(data):
                raise ValueError("archived report/request differs from retained diagnostic manifest")
            value = json.loads(data)
            if is_request:
                if request is not None:
                    raise ValueError("ambiguous parent request")
                request, request_digest = value, digest
            else:
                seed = int(match[1])
                if seed in reports:
                    raise ValueError("duplicate seed refinement report")
                reports[seed] = (value, digest, member.name)
    if request is None or set(reports) != {7, 101, 102, 103, 104}:
        raise ValueError("parent request or five stage reports missing")
    gaps, generated, anchors = [], {}, {}
    for index, gap in enumerate(request["gaps"]):
        interval = gap["target_interval"]
        targets = {t["chain"]: t["sequence"] for t in request["target_sequences"]}
        gaps.append({"index": index, "chain": gap["chain"], "kind": gap["gap_kind"],
                     "length": interval["stop"] - interval["start"],
                     "sequence": targets[gap["chain"]][interval["start"]:interval["stop"]],
                     "first": gap["generated_residues"][0], "last": gap["generated_residues"][-1],
                     "left_anchor": gap["left_anchor"], "right_anchor": gap["right_anchor"]})
        for residue in gap["generated_residues"]:
            key = identity(residue)
            if key in generated:
                raise ValueError("overlapping gap residue ownership")
            generated[key] = index
        for anchor in (gap["left_anchor"], gap["right_anchor"]):
            if anchor is not None:
                anchors.setdefault(identity(anchor), set()).add(index)
    outcomes = {r["seed"]: r for r in job["diffusion"]["outcomes"]}
    rows = []
    for seed, (report, digest, archive_path) in sorted(reports.items()):
        stages = [stage_summary(a, generated, anchors) for a in report["attempts"]]
        rows.append({"seed": seed, "status": outcomes[seed]["status"],
                     "gates": outcomes[seed]["hardGateFailures"],
                     "final_findings": outcomes[seed]["findings"],
                     "omitted_final_findings": outcomes[seed]["omittedFindings"],
                     "stages": stages, "notices": report["notices"],
                     "archive_report": archive_path, "report_sha256": digest})
    return {"pdb": case["pdb"], "request_sha256": request_digest, "gaps": gaps, "seeds": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    cohort = json.loads((args.root / "cohort-final-outcomes.json").read_text())
    results = []
    for case in cohort["rows"]:
        row = analyze_case(args.root, case)
        (args.output / f"{case['pdb']}.json").write_text(json.dumps(row, indent=2, ensure_ascii=False) + "\n")
        results.append(row)
        print(row["pdb"], [(s["seed"], s["status"],
            [(a["stage"], a["geometry_passed"], a.get("issue_count"), a["error"]) for a in s["stages"]])
              for s in row["seeds"]], flush=True)
    rejected = [s for case in results for s in case["seeds"] if s["status"] == "rejected"]
    gates = Counter(g for seed in rejected for g in seed["gates"])
    result = {"scope": "digest-verified stage reports and actual gap ownership; no new sampling/refinement",
              "candidates": sum(len(c["seeds"]) for c in results), "rejected": len(rejected),
              "overlapping_rejection_gate_counts": dict(gates), "cases": results}
    (args.output / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "cases"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
