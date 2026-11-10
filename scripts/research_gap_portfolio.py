#!/usr/bin/env python3
"""Declare and test a private portfolio from verified old generated-only stages.

No neural sampling, no original-outcome rewrite, no schema-4 publication.
All donors are read from one exact common request and normalized source frame.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from dataclasses import asdict
from pathlib import Path, PurePosixPath


def main():
    from dvbfixer.model.diffusion.contract import DiffusionRequest, RunnerResult, SamplerTrace
    from dvbfixer.model.diffusion.gap_portfolio import GapParent, digest, search_gap_portfolio
    from dvbfixer.model.diffusion.refinement_evidence import (
        StagedRefinementEvidence,
        evidence_reference,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-evaluations", type=int, default=64)
    parser.add_argument("--beam-width", type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    protocol = {"scope": "post-hoc private portfolio, no sampling or original-outcome replacement",
        "case": args.case, "seeds": [7, 101, 102, 103, 104], "stage": "generated-only",
        "max_evaluations": args.max_evaluations, "beam_width": args.beam_width,
        "anchor_mobility": False, "kind": "experimental-multiparent-gap-portfolio"}
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    import dvbfixer.model.diffusion.gap_portfolio as module

    (args.output / "gap_portfolio.py").write_bytes(Path(module.__file__).read_bytes())
    (args.output / "run-script.py").write_bytes(Path(__file__).read_bytes())
    execution = args.cohort / "managed-execution"
    job = json.loads((execution / f"{args.case}-finished-job.json").read_text())
    run = execution / "data/projects" / job["workspaceId"] / job["outputDir"]
    manifest = {r["file"]: r for r in json.loads((run / "run-summary.json").read_text())["diagnostics"]["files"]}
    parents = []
    with tarfile.open(run / "diagnostics.tar.gz", "r:gz") as archive:
        names = archive.getnames()
        paths = [n for n in names if n.endswith("/request.json") and len(PurePosixPath(n).parts) == 2]
        if len(paths) != 1:
            raise ValueError("ambiguous parent request")
        base = PurePosixPath(paths[0]).parent

        def read(name):
            member = archive.getmember(name)
            if not member.isfile() or member.size > 16_000_000:
                raise ValueError("unbounded portfolio custody input")
            data = archive.extractfile(member).read()
            if len(data) != manifest[name]["bytes"] or hashlib.sha256(data).hexdigest() != manifest[name]["sha256"]:
                raise ValueError("portfolio custody differs from retained manifest")
            return data.decode("utf-8")

        request = DiffusionRequest.from_json(read(paths[0]))
        source = read((base / request.normalized_pdb.path).as_posix())
        (args.output / "request.json").write_text(request.to_json())
        (args.output / "source.pdb").write_text(source)
        for seed in protocol["seeds"]:
            prefix = base / f"candidates/protpardelle-1c-cc94-mps-v1/seed-{seed}"
            raw_result_text = read((prefix / "raw-result.json").as_posix())
            raw_result = RunnerResult.from_json(raw_result_text)
            if len(raw_result.candidates) != 1 or raw_result.candidates[0].seed != seed:
                raise ValueError("original runner seed differs from declaration")
            trace_text = read((prefix / "refined/sampler-trace.json").as_posix())
            reference = evidence_reference(SamplerTrace.from_json(trace_text))
            if reference is None:
                raise ValueError("portfolio parent lacks bound refinement evidence")
            evidence_text = read((base / reference.path).as_posix())
            if digest(evidence_text) != reference.sha256:
                raise ValueError("trace-bound evidence mismatch")
            evidence = StagedRefinementEvidence.from_json(evidence_text)
            evidence.validate_request(request, seed=seed)
            if evidence.stage1 is None:
                raise ValueError("portfolio donor lacks a completed generated-only stage")
            text = read((base / evidence.stage1.path).as_posix())
            if digest(text) != evidence.stage1.sha256:
                raise ValueError("portfolio donor coordinate mismatch")
            parent = GapParent(raw_result.candidates[0].candidate_id, seed, digest(request.to_json()),
                digest(source), digest(text), digest(trace_text), digest(evidence_text), text)
            parents.append(parent)
            directory = args.output / f"parent-{seed}"
            directory.mkdir()
            for name, data in (("stage1.pdb", text), ("trace.json", trace_text),
                               ("evidence.json", evidence_text), ("raw-result.json", raw_result_text)):
                (directory / name).write_text(data)
        (args.output / "parents.json").write_text(json.dumps([
            {k: v for k, v in asdict(p).items() if k != "pdb_text"} for p in parents], indent=2) + "\n")
    trials = args.output / "trials"
    trials.mkdir()
    counter = 0

    def retain(trial, text):
        nonlocal counter
        path = trials / f"trial-{counter:04d}"
        path.with_suffix(".pdb").write_text(text)
        path.with_suffix(".json").write_text(json.dumps(asdict(trial), indent=2) + "\n")
        counter += 1

    result = search_gap_portfolio(request, source, tuple(parents), max_evaluations=args.max_evaluations,
        beam_width=args.beam_width, retain_trial=retain)
    (args.output / "result.pdb").write_text(result.coordinate_text)
    report = {k: v for k, v in asdict(result).items() if k != "coordinate_text"}
    (args.output / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"case": args.case, "score": result.selected.score,
        "choices": [asdict(c) for c in result.selected.choices], "trials": len(result.trials),
        "repair_scope": asdict(result.repair_scope), "publication_qualified": False}), flush=True)


if __name__ == "__main__":
    main()
