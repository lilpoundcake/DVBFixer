#!/usr/bin/env python3
"""Post-hoc paired CPU iteration control on an unchanged saved stage-1 proposal.

Not a sampling repeat and never a replacement for any original cohort outcome.
Both arms reconstruct the same temporary FF input with the same hydrogen seed;
only max_iterations differs. No stochastic restarts or observed mobility.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path, PurePosixPath


def main():
    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--case", required=True)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    job = json.loads((args.root / "managed-execution" / f"{args.case}-finished-job.json").read_text())
    run = args.root / "managed-execution/data/projects" / job["workspaceId"] / job["outputDir"]
    summary = json.loads((run / "run-summary.json").read_text())
    manifest = {r["file"]: r for r in summary["diagnostics"]["files"]}
    with tarfile.open(run / "diagnostics.tar.gz", "r:gz") as archive:
        names = archive.getnames()
        requests = [n for n in names if len(PurePosixPath(n).parts) == 2 and n.endswith("/request.json")]
        if len(requests) != 1:
            raise ValueError("ambiguous archived request")
        base = PurePosixPath(requests[0]).parent

        def read(name):
            member = archive.getmember(name)
            if not member.isfile() or member.size > 16_000_000:
                raise ValueError("unbounded or nonregular control input")
            data = archive.extractfile(member).read()
            if hashlib.sha256(data).hexdigest() != manifest[name]["sha256"]:
                raise ValueError("control input digest differs from diagnostic manifest")
            return data

        request = DiffusionRequest.from_json(read(requests[0]).decode())
        source = read((base / request.normalized_pdb.path).as_posix()).decode()
        reports = [n for n in names if n.endswith(f"/seed-{args.seed}/refined/refinement-report.json")]
        if len(reports) != 1:
            raise ValueError("ambiguous stage report")
        report = json.loads(read(reports[0]))
        first = next(a for a in report["attempts"] if a["stage"] == "generated-only")
        proposal = read((base / first["coordinate_artifact"]["path"]).as_posix()).decode()
    (args.output / "stage1-input.pdb").write_text(proposal)
    (args.output / "source.pdb").write_text(source)
    (args.output / "request.json").write_text(request.to_json())
    assess = assessment_adapter(request, source)
    before = assess(proposal)
    protocol = {"scope": "post-hoc CPU diagnostic, no sampling or original-outcome replacement",
        "case": args.case, "seed": args.seed, "max_iterations": [250, 2500],
        "restart_count": 0, "platform": "CPU", "observed_mobility": False,
        "proposal_sha256": hashlib.sha256(proposal.encode()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "before_assessment": asdict(before)}
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    generated_residues = tuple(dict.fromkeys(r for gap in request.gaps for r in gap.generated_residues))
    results = []
    for iterations in (250, 2500):
        start = time.monotonic()
        try:
            result = refine_generated_region(proposal, generated_residues=generated_residues,
                generated_atoms=request.generated_atoms, max_iterations=iterations,
                restart_count=0, platform_name="CPU", random_seed=args.seed)
            text = rewrite_generated_coordinates(proposal, result.coordinates_angstrom)
            (args.output / f"iterations-{iterations}.pdb").write_text(text)
            assessment = assess(text)
            row = {"max_iterations": iterations, "geometry_passed": assessment.geometry_passed,
                "issues": dict(Counter(i.category for i in assessment.issues if i.actionable)),
                "admission_failures": assessment.admission_failures,
                "initial_energy_kj_mol": result.initial_energy_kj_mol,
                "final_energy_kj_mol": result.final_energy_kj_mol,
                "full_assessment": asdict(assessment)}
        except Exception as exc:
            row = {"max_iterations": iterations, "error_type": type(exc).__name__, "error": str(exc)}
        row["elapsed_seconds"] = round(time.monotonic() - start, 3)
        results.append(row)
        (args.output / f"iterations-{iterations}.json").write_text(json.dumps(row, indent=2) + "\n")
        print(json.dumps({k: v for k, v in row.items() if k != "full_assessment"}), flush=True)
    (args.output / "results.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
