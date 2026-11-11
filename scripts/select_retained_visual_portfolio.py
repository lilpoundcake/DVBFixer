#!/usr/bin/env python3
"""New constrained 1HE8 research assembly, never overwrite historical results.

Uses the user's hash-bound visual allow-lists plus independent scientific whole
assembly checks. Surface-probe alarms are retained diagnostics, not calibrated
vetoes. No resampling, observation movement or invented backend trace occurs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from dvbfixer.model.diffusion.compatible_portfolio import DonorAlternative, select_compatible
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import (
    GapChoice,
    GapParent,
    _coordinates,
    audit_portfolio_materialization,
)
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def digest(data):
    return hashlib.sha256(data).hexdigest()


def run(selection, review_path, output):
    manifest_path = review_path.parent / "manifest.json"
    review = json.loads(review_path.read_text())
    if digest(manifest_path.read_bytes()) != review["manifest_sha256"]:
        raise ValueError("visual-review manifest binding differs")
    for gap in review["gaps"]:
        if digest((review_path.parent / gap["file"]).read_bytes()) != gap["file_sha256"]:
            raise ValueError("visually reviewed coordinates changed")
    request = DiffusionRequest.from_json((selection / "request.json").read_text())
    source = (selection / "source.pdb").read_text()
    metadata = json.loads((selection / "parents.json").read_text())
    if json.loads(manifest_path.read_text())["donors"] != metadata:
        raise ValueError("visual-review donors differ from retained source custody")
    assess = assessment_adapter(request, source)
    parents, axes, counts = [], [], []
    for item in metadata:
        folder = selection / f"parent-{item['seed']}"
        for name, key in (("trace.json", "sampler_trace_sha256"), ("evidence.json", "evidence_sha256")):
            if digest((folder / name).read_bytes()) != item[key]:
                raise ValueError("retained parent custody digest differs")
        parent = GapParent(**item, pdb_text=(folder / "stage1.pdb").read_text())
        if parent.source_sha256 != digest(source.encode()) or parent.request_sha256 != digest(request.to_json().encode()):
            raise ValueError("retained parent request/source differs")
        assessment = assess(parent.pdb_text)
        if assessment.admission_failures:
            raise ValueError("retained donor changed observations/identity")
        counts.append([sum(issue.actionable and bool(set(issue.residues) &
            (set(gap.generated_residues) | {a for a in (gap.left_anchor, gap.right_anchor) if a is not None}))
            for issue in assessment.issues) for gap in request.gaps])
        parents.append(parent)
        axes.append(_coordinates(parent.pdb_text)[0])
    owned = [tuple(a for a in request.generated_atoms if (a.chain, a.residue_number, a.insertion_code) in
                   {(r.chain, r.residue_number, r.insertion_code) for r in gap.generated_residues}) for gap in request.gaps]
    alternatives, allowed, lookup = [], {}, {}
    for i, gap in enumerate(review["gaps"]):
        if gap["gap_index"] != i:
            raise ValueError("visual-review gap order differs")
        allowed[f"gap-{i:02d}"] = frozenset(f"gap-{i:02d}-seed-{seed}" for seed in gap["visually_approved_seeds"])
        for p, parent in enumerate(parents):
            key = f"gap-{i:02d}-seed-{parent.seed}"
            lineage = digest(json.dumps({"parent": metadata[p], "gap_index": i}, sort_keys=True).encode())
            alternatives.append(DonorAlternative(key, f"gap-{i:02d}", (i,), lineage, f"seed-{parent.seed}", counts[p][i]))
            lookup[key] = i, p

    def materialize(keys):
        xyz = {a: axes[p][a] for key in keys for i, p in (lookup[key],) for a in owned[i]}
        return rewrite_generated_coordinates(parents[0].pdb_text, xyz)

    output.mkdir(parents=True, exist_ok=False)
    (output / "protocol.json").write_text(json.dumps({"max_evaluations": 64, "beam_width": 4,
        "visual_review_sha256": digest(review_path.read_bytes()), "source_sha256": request.normalized_pdb.sha256,
        "pair_screen": "gap2/gap3 severe clash findings only; other pairs not pre-pruned; all full assemblies assessed",
        "surface_hits_are_not_vetoes": True, "new_sampling": False, "new_minimization": False,
        "joint_groups_not_mixed": "retained donors are independent gap motifs, not joint proposals",
        "script_sha256": digest(Path(__file__).read_bytes())}, indent=2) + "\n")
    conflicts, pair_reports = set(), []
    default = tuple(f"gap-{i:02d}-seed-{parents[0].seed}" for i in range(len(request.gaps)))
    for first in sorted(allowed["gap-02"]):
        for second in sorted(allowed["gap-03"]):
            keys = default[:2] + (first, second) + default[4:]
            assessment = assess(materialize(keys))
            left = set(request.gaps[2].generated_residues)
            right = set(request.gaps[3].generated_residues)
            bridges = [issue for issue in assessment.issues if issue.category == "clash" and issue.actionable
                       and set(issue.residues) & left and set(issue.residues) & right]
            if bridges:
                conflicts.add(frozenset((first, second)))
            pair_reports.append({"alternatives": (first, second), "direct_severe_clashes": [asdict(issue) for issue in bridges]})
    (output / "pair-screen.json").write_text(json.dumps(pair_reports, indent=2) + "\n")
    trial_number = 0
    def complete(keys):
        nonlocal trial_number
        text = materialize(keys)
        (output / f"trial-{trial_number:03d}.pdb").write_text(text)
        trial_number += 1
        return assess(text)
    result = select_compatible(tuple(alternatives), assess_complete=complete, allowed=allowed,
                               conflicts=frozenset(conflicts), max_evaluations=64, beam_width=4)
    if result.selected is not None:
        text = materialize(result.selected.alternatives)
        choices = tuple(GapChoice(*lookup[key]) for key in result.selected.alternatives)
        independent = audit_portfolio_materialization(request, source, tuple(parents), choices, text)
        (output / "selected.pdb").write_text(text)
        (output / "selected-independent-audit.json").write_text(json.dumps(asdict(independent), indent=2) + "\n")
        (output / "evaluation-mask.json").write_text(request.to_json())
        (output / "source.pdb").write_text(source)
        if request.normalized_pdb.path != "source.pdb":
            (output / request.normalized_pdb.path).write_text(source)
        (output / "source-restored.pdb").write_text(text)
    (output / "result.json").write_text(json.dumps(asdict(result), indent=2) + "\n")
    print(json.dumps({"evaluated": len(result.trials), "space": result.constrained_space_size,
                      "optimality": result.optimality_established,
                      "selected_score": result.selected.score if result.selected else None}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("selection", "review", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    run(args.selection, args.review, args.output)
