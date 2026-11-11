#!/usr/bin/env python3
"""Bounded post-refinement donor derivatives, retaining every rejected assembly.

Residual findings select which donor alternatives to try, NOT which generated
atoms may move: every trial still jointly minimizes ALL generated atoms. Retained-
donor cases honor user allow-lists; the RFD3 group is a separately labeled,
unreviewed control. No observed anchor is released or reference supplied.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import traceback
from dataclasses import asdict
from pathlib import Path

from scripts.refine_mac_reconstruction_smoke import run as refine

from dvbfixer.model.diffusion.contract import DiffusionRequest, ResidueIdentity
from dvbfixer.model.diffusion.gap_portfolio import (
    GapChoice,
    GapParent,
    _coordinates,
    audit_portfolio_materialization,
    plan_gap_repairs,
)
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment, RepairIssue
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(root, selection, review_path, output):
    base = root / "1he8-visual-compatible-selection"
    request = DiffusionRequest.from_json((base / "evaluation-mask.json").read_text())
    source = (base / "source.pdb").read_bytes().decode()
    baseline = (base / "selected.pdb").read_bytes().decode()
    protocol = json.loads((base / "protocol.json").read_text())
    if digest(review_path) != protocol["visual_review_sha256"]:
        raise ValueError("manual review binding changed")
    review = json.loads(review_path.read_text())
    metadata = json.loads((selection / "parents.json").read_text())
    parents, axes = [], []
    for item in metadata:
        folder = selection / f"parent-{item['seed']}"
        if digest(folder / "trace.json") != item["sampler_trace_sha256"] or digest(folder / "evidence.json") != item["evidence_sha256"]:
            raise ValueError("retained donor evidence changed")
        parent = GapParent(**item, pdb_text=(folder / "stage1.pdb").read_text())
        parents.append(parent)
        axes.append(_coordinates(parent.pdb_text)[0])
    initial = json.loads((base / "result.json").read_text())["selected"]["alternatives"]
    selected_seeds = [int(key.rsplit("-", 1)[1]) for key in initial]
    previous = root / "1he8-visual-all-generated-conditioned-cpu/status.json"
    findings = json.loads(previous.read_text())["after"]
    assessment = RepairAssessment(tuple(RepairIssue(item["category"], tuple(ResidueIdentity(**r) for r in item["residues"]),
        item["message"], item["actionable"]) for item in findings["issues"]), tuple(findings["admission_failures"]))
    if assessment.admission_failures:
        raise ValueError("admission failure cannot authorize donor/motion repair")
    scope = plan_gap_repairs(request, assessment)
    allowed = {g["gap_index"]: tuple(sorted(g["visually_approved_seeds"])) for g in review["gaps"]}
    by_seed = {p.seed:i for i,p in enumerate(parents)}
    owned = [tuple(a for a in request.generated_atoms if (a.chain, a.residue_number, a.insertion_code) in
                   {(r.chain, r.residue_number, r.insertion_code) for r in gap.generated_residues}) for gap in request.gaps]
    output.mkdir(exist_ok=False)
    cases = []
    for attempt in range(3):
        seeds = list(selected_seeds)
        for index in scope.gap_indices:
            alternatives = [seed for seed in allowed[index] if seed != selected_seeds[index]]
            if alternatives:
                seeds[index] = alternatives[attempt % len(alternatives)]
        if any(c["donor_seeds"] == seeds for c in cases) or seeds == selected_seeds:
            continue
        if any(seed not in allowed[index] for index, seed in enumerate(seeds)):
            raise ValueError("derivative violates manual allowed alternatives")
        xyz = {a:axes[by_seed[seed]][a] for index,seed in enumerate(seeds) for a in owned[index]}
        text = rewrite_generated_coordinates(parents[0].pdb_text, xyz)
        choices = tuple(GapChoice(i, by_seed[seed]) for i,seed in enumerate(seeds))
        audited = audit_portfolio_materialization(request, source, tuple(parents), choices, text)
        work = output / f"donor-derivative-{attempt}"
        work.mkdir()
        (work / "source-restored.pdb").write_text(text)
        (work / request.normalized_pdb.path).write_text(source)
        (work / "source.pdb").write_text(source)
        (work / "evaluation-mask.json").write_text(request.to_json())
        cases.append({"name":work.name, "donor_seeds":seeds, "work":str(work), "raw_sha256":digest(work / "source-restored.pdb"),
                      "raw_independent_audit":asdict(audited), "joint_group_origin":"independent retained motifs only"})
    context = root / "1he8-rfd3-context-assembly-recovered/joint-context"
    cases.append({"name":"rfd3-context", "work":str(context), "raw_sha256":digest(context / "source-restored.pdb"),
                  "donor_seeds":None, "joint_group_origin":"gaps 2/3 kept together from one RFD3 run; other eight retained donors"})
    declaration = {"input_source_sha256":request.normalized_pdb.sha256, "reference_baseline_raw_sha256":digest(base / "selected.pdb"),
        "prior_conditioned_result_sha256":digest(previous), "manual_review_sha256":digest(review_path),
        "donor_alternative_gap_scope":asdict(scope), "all_generated_atoms_mobile_in_every_case":len(request.generated_atoms),
        "all_observations_fixed":True, "platform":"CPU", "actual_threads":1, "refinement_seed":7,
        "iteration_caps":[100,75,75,1000], "restarts":0, "max_cases":4, "cases":cases,
        "ranking":"whole-assembly admission failures, actionable geometry findings, case name; no reference/confidence optimization",
        "artificial_surface_alarms_not_qualified_vetoes":True, "source_anchor_calibration_not_claimed":True,
        "old_cohort_outcomes_unchanged":True, "qualification":False}
    (output / "declaration.json").write_text(json.dumps(declaration, indent=2)+"\n")
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    results = []
    for case in cases:
        result = {"name":case["name"]}
        try:
            destination = output / f"{case['name']}-refinement"
            refine(Path(case["work"]), destination, platform="CPU", geometry_preconditioning=True, force_shared_refinement=True)
            result["status"] = json.loads((destination / "status.json").read_text())
            result["refined_pdb"] = str(destination / "refined.pdb")
        except Exception:
            result["failure"] = traceback.format_exc()
        results.append(result)
        (output / "results.json").write_text(json.dumps(results, indent=2)+"\n")
    (output / "completed.json").write_text(json.dumps({"cases_executed":len(results), "qualified":False})+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "selection", "review", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    run(args.root.absolute(), args.selection, args.review, args.output.absolute())
