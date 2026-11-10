#!/usr/bin/env python3
"""One declared issue-local repair of a private multi-parent portfolio."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path


def main():
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.gap_portfolio import (
        GapChoice,
        GapParent,
        audit_portfolio_materialization,
        digest,
        repair_portfolio_gaps,
    )
    from dvbfixer.model.diffusion.refinement_evidence import StagedRefinementEvidence

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--portfolio", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--refinement-seed", type=int, default=7)
    parser.add_argument("--max-iterations", type=int, default=250)
    parser.add_argument("--restart-count", type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    request = DiffusionRequest.from_json((args.portfolio / "request.json").read_text())
    source = (args.portfolio / "source.pdb").read_text()
    selected = (args.portfolio / "result.pdb").read_text()
    selection = json.loads((args.portfolio / "result.json").read_text())
    choices = tuple(GapChoice(**c) for c in selection["selected"]["choices"])
    parents = []
    for row in json.loads((args.portfolio / "parents.json").read_text()):
        directory = args.portfolio / f"parent-{row['seed']}"
        text = (directory / "stage1.pdb").read_text()
        trace = (directory / "trace.json").read_text()
        evidence_text = (directory / "evidence.json").read_text()
        parent = GapParent(**row, pdb_text=text)
        if digest(trace) != parent.sampler_trace_sha256 or digest(evidence_text) != parent.evidence_sha256:
            raise ValueError("portfolio parent trace/evidence custody changed")
        evidence = StagedRefinementEvidence.from_json(evidence_text)
        evidence.validate_request(request, seed=parent.seed)
        if evidence.stage1 is None or evidence.stage1.sha256 != parent.coordinate_sha256:
            raise ValueError("portfolio parent stage is not bound to its evidence")
        parents.append(parent)
    before = audit_portfolio_materialization(request, source, tuple(parents), choices, selected)
    protocol = {"kind": "experimental-multiparent-gap-portfolio", "scope": "one generated-only issue-local CPU repair; no new neural sample",
        "request_sha256": digest(request.to_json()), "source_sha256": digest(source),
        "input_sha256": digest(selected), "choices": [asdict(c) for c in choices],
        "refinement_seed": args.refinement_seed, "max_iterations": args.max_iterations,
        "restart_count": args.restart_count, "platform": "CPU", "anchor_mobility": False}
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    (args.output / "input.pdb").write_text(selected)
    import dvbfixer.model.diffusion.gap_portfolio as module

    for name in ("gap_portfolio.py", "boundary_refinement.py", "geometry_conditioning.py"):
        (args.output / name).write_bytes(Path(module.__file__).with_name(name).read_bytes())
    (args.output / "before-audit.json").write_text(json.dumps(asdict(before), indent=2) + "\n")
    try:
        result = repair_portfolio_gaps(request, source, selected, refinement_seed=args.refinement_seed,
            max_iterations=args.max_iterations, restart_count=args.restart_count)
        (args.output / "result.pdb").write_text(result.coordinate_text)
        report = {k: v for k, v in asdict(result).items() if k != "coordinate_text"}
        (args.output / "repair-report.json").write_text(json.dumps(report, indent=2) + "\n")
        audit = audit_portfolio_materialization(request, source, tuple(parents), choices, selected,
            repaired_text=result.coordinate_text)
        passed = audit.geometry_passed and not audit.admission_failures
        summary = {"passed_geometry_and_custody": passed, "publication_qualified": False,
            "scope": asdict(result.scope), "independent_assessment": asdict(audit),
            "output_sha256": digest(result.coordinate_text)}
        (args.output / "independent-audit.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary), flush=True)
    except Exception as exc:
        (args.output / "failure.json").write_text(json.dumps({"type": type(exc).__name__, "error": str(exc)}, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
