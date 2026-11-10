#!/usr/bin/env python3
"""Read-only, separate-process audit of private portfolio lineage and local repair."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path


def main():
    from dvbfixer.model.diffusion.contract import DiffusionRequest, SamplerTrace
    from dvbfixer.model.diffusion.gap_portfolio import (
        PORTFOLIO_KIND,
        PORTFOLIO_REVISION,
        GapChoice,
        GapParent,
        audit_portfolio_materialization,
        digest,
    )
    from dvbfixer.model.diffusion.refinement_evidence import (
        StagedRefinementEvidence,
        evidence_reference,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--portfolio", required=True, type=Path)
    parser.add_argument("--repair", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    request = DiffusionRequest.from_json((args.portfolio / "request.json").read_text())
    source = (args.portfolio / "source.pdb").read_text()
    selected = (args.portfolio / "result.pdb").read_text()
    selection = json.loads((args.portfolio / "result.json").read_text())
    if selection["kind"] != PORTFOLIO_KIND or selection["revision"] != PORTFOLIO_REVISION:
        raise ValueError("unsupported private portfolio contract")
    if digest(selected) != selection["selected"]["coordinate_sha256"]:
        raise ValueError("selected portfolio coordinate digest differs")
    parents = []
    for row in json.loads((args.portfolio / "parents.json").read_text()):
        directory = args.portfolio / f"parent-{row['seed']}"
        parent = GapParent(**row, pdb_text=(directory / "stage1.pdb").read_text())
        trace_text = (directory / "trace.json").read_text()
        evidence_text = (directory / "evidence.json").read_text()
        if digest(trace_text) != parent.sampler_trace_sha256 or digest(evidence_text) != parent.evidence_sha256:
            raise ValueError("portfolio parent custody differs")
        bound = evidence_reference(SamplerTrace.from_json(trace_text))
        if bound is None or bound.sha256 != parent.evidence_sha256:
            raise ValueError("parent trace does not bind the retained refinement evidence")
        evidence = StagedRefinementEvidence.from_json(evidence_text)
        evidence.validate_request(request, seed=parent.seed)
        if evidence.stage1 is None or evidence.stage1.sha256 != parent.coordinate_sha256:
            raise ValueError("trace/evidence does not bind this generated-only donor")
        parents.append(parent)
    choices = tuple(GapChoice(**c) for c in selection["selected"]["choices"])
    repaired = None
    if args.repair is not None:
        repaired = (args.repair / "result.pdb").read_text()
        operation = json.loads((args.repair / "repair-report.json").read_text())
        if (operation["operation"] != "generated-only-issue-local-refinement-v1"
                or operation["input_sha256"] != digest(selected) or operation["output_sha256"] != digest(repaired)):
            raise ValueError("local repair operation/input/output custody differs")
    assessment = audit_portfolio_materialization(request, source, tuple(parents), choices, selected,
        repaired_text=repaired)
    report = {"scope": "independent private portfolio geometry and custody, not public sampling qualification",
        "kind": PORTFOLIO_KIND, "request_sha256": digest(request.to_json()), "source_sha256": digest(source),
        "selected_sha256": digest(selected), "final_sha256": digest(repaired or selected),
        "passed": assessment.geometry_passed and not assessment.admission_failures,
        "assessment": asdict(assessment), "lineage": [
            {"gap_index": c.gap_index, "parent_candidate_id": parents[c.parent_index].candidate_id,
             "parent_seed": parents[c.parent_index].seed,
             "parent_coordinate_sha256": parents[c.parent_index].coordinate_sha256,
             "sampler_trace_sha256": parents[c.parent_index].sampler_trace_sha256,
             "evidence_sha256": parents[c.parent_index].evidence_sha256} for c in choices],
        "publication_qualified": False}
    with args.report.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: v for k, v in report.items() if k not in {"lineage", "assessment"}}), flush=True)


if __name__ == "__main__":
    main()
