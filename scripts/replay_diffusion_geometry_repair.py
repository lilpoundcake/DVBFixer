#!/usr/bin/env python3
"""Post-hoc v2 repair of digest-verified retained proposals, never new sampling.

Original archives/outcomes are read-only. Each invocation exclusively reserves
its output, records the request version change, and retains independent evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import time
from dataclasses import replace
from pathlib import Path, PurePosixPath


def main():
    from dvbfixer.model.diffusion.contract import DiffusionRequest, RunnerResult
    from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy, read_refinement_policy
    from dvbfixer.model.diffusion.runner_refinement import refine_runner_result
    from dvbfixer.model.diffusion.validate import build_validated_result

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--case", required=True)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    execution = args.root / "managed-execution"
    job = json.loads((execution / f"{args.case}-finished-job.json").read_text())
    run = execution / "data/projects" / job["workspaceId"] / job["outputDir"]
    manifest = {r["file"]: r for r in json.loads((run / "run-summary.json").read_text())["diagnostics"]["files"]}
    with tarfile.open(run / "diagnostics.tar.gz", "r:gz") as archive:
        names = archive.getnames()
        requests = [n for n in names if len(PurePosixPath(n).parts) == 2 and n.endswith("/request.json")]
        if len(requests) != 1:
            raise ValueError("ambiguous archived request")
        parent = PurePosixPath(requests[0]).parent

        def read(name):
            member = archive.getmember(name)
            if not member.isfile() or member.size > 16_000_000:
                raise ValueError("unbounded or nonregular replay input")
            data = archive.extractfile(member).read()
            expected = manifest[name]
            if len(data) != expected["bytes"] or hashlib.sha256(data).hexdigest() != expected["sha256"]:
                raise ValueError("replay artifact differs from retained manifest")
            return data

        original_bytes = read(requests[0])
        original = DiffusionRequest.from_json(original_bytes.decode())
        policy = read_refinement_policy(original.backend_options)
        if policy is None:
            raise ValueError("replay expects an explicit original staged policy")
        options = tuple(o for o in original.backend_options if not o.name.startswith("geometry_refinement."))
        updated = replace(original, backend_options=(*options,
            *RefinementPolicy(policy.relax_anchors, policy.bfactor_weighting).backend_options()))
        matches = [n for n in names if n.endswith(f"/seed-{args.seed}/raw-result.json")]
        if len(matches) != 1:
            raise ValueError("ambiguous raw result")
        raw_bytes = read(matches[0])
        raw_result = RunnerResult.from_json(raw_bytes.decode())
        candidate = raw_result.candidates[0]
        references = (original.normalized_pdb, candidate.coordinate_artifact, candidate.sampler_trace_artifact)
        for reference in references:
            relative = PurePosixPath(reference.path)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("unsafe input artifact path")
            data = read((parent / relative).as_posix())
            if hashlib.sha256(data).hexdigest() != reference.sha256:
                raise ValueError("reference digest mismatch")
            path = output / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("xb") as stream:
                stream.write(data)
    (output / "original-request.json").write_bytes(original_bytes)
    (output / "raw-result.json").write_bytes(raw_bytes)
    request_path = output / "request.json"
    request_path.write_text(updated.to_json())
    protocol = {"scope": "post-hoc retained-coordinate repair; no new sampling or original-outcome replacement",
        "case": args.case, "seed": args.seed, "original_policy": policy.revision,
        "new_policy": read_refinement_policy(updated.backend_options).revision,
        "original_request_sha256": hashlib.sha256(original_bytes).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "platform": "CPU"}
    import dvbfixer.model.diffusion.runner_refinement as refinement_module

    source_directory = Path(refinement_module.__file__).parent
    source_snapshot = output / "software-snapshot"
    source_snapshot.mkdir()
    protocol["software_sha256"] = {}
    for name in ("runner_refinement.py", "boundary_refinement.py", "geometry_conditioning.py", "refinement_policy.py"):
        data = (source_directory / name).read_bytes()
        (source_snapshot / name).write_bytes(data)
        protocol["software_sha256"][name] = hashlib.sha256(data).hexdigest()
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    start = time.monotonic()
    summary = {"error": "interrupted before result"}
    try:
        result = refine_runner_result(request_path, raw_result, output / "refined", platform_name="CPU")
        (output / "replay-result.json").write_text(result.to_json())
        validated = build_validated_result(updated, result, workspace=output)
        (output / "validated-result.json").write_text(validated.to_json())
        summary = {"postprocessing_failures": result.candidates[0].postprocessing_failures,
            "warnings": result.candidates[0].warnings,
            "independent_validation": [
                {"passed": v.passed, "hard_gate_failures": v.hard_gate_failures}
                for v in validated.validation_summaries]}
    except Exception as exc:
        summary = {"error_type": type(exc).__name__, "error": str(exc)}
        raise
    finally:
        summary["elapsed_seconds"] = time.monotonic() - start
        (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
