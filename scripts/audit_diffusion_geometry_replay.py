#!/usr/bin/env python3
"""Independent revalidation and bundle readback of completed post-hoc repairs.

Never substitutes these repairs for the original sampling cohort. All cases
must have finished before auditing; a missing result is not silently skipped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, replace
from pathlib import Path


def main():
    from dvbfixer.model.diffusion.contract import ArtifactReference, DiffusionRequest, RunnerResult
    from dvbfixer.model.diffusion.input_preparation import (
        PUBLICATION_ROLES,
        audit_preparation_payloads,
        preparation_options,
    )
    from dvbfixer.model.diffusion.pipeline import publish_diffusion_bundle
    from dvbfixer.model.diffusion.private_artifact import read_private_artifact
    from dvbfixer.model.diffusion.publication import read_published_bundle
    from dvbfixer.model.diffusion.refinement_policy import read_refinement_policy
    from dvbfixer.model.diffusion.validate import build_validated_result

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls", required=True, type=Path)
    parser.add_argument("--cohort", required=True, type=Path)
    args = parser.parse_args()
    controls = args.controls.resolve()
    declaration = json.loads((controls / "declaration.json").read_text())
    rows = []
    for case, seed in declaration["cases"]:
        workspace = controls / f"{case}-{seed}"
        original = DiffusionRequest.from_json((workspace / "original-request.json").read_text())
        request = DiffusionRequest.from_json((workspace / "request.json").read_text())
        original_options = tuple(o for o in original.backend_options if not o.name.startswith("geometry_refinement."))
        updated_options = tuple(o for o in request.backend_options if not o.name.startswith("geometry_refinement."))
        if original_options != updated_options or replace(request, backend_options=original.backend_options) != original:
            raise ValueError("repair replay modified more than the declared refinement policy")
        old_policy = read_refinement_policy(original.backend_options)
        new_policy = read_refinement_policy(request.backend_options)
        if old_policy.movement_limits != new_policy.movement_limits:
            raise ValueError("repair replay changed movement acceptance thresholds")
        result = RunnerResult.from_json((workspace / "replay-result.json").read_text())
        if result.candidates[0].seed != seed:
            raise ValueError("replay candidate seed differs from declaration")
        validated = build_validated_result(request, result, workspace=workspace)
        (workspace / "independent-revalidation.json").write_text(validated.to_json())
        summary = validated.validation_summaries[0]
        row = {"case": case, "seed": seed, "passed": summary.passed,
            "hard_gate_failures": summary.hard_gate_failures,
            "warnings": summary.warnings, "metrics": [asdict(m) for m in summary.metrics],
            "original_request_sha256": hashlib.sha256((workspace / "original-request.json").read_bytes()).hexdigest(),
            "request_sha256": hashlib.sha256(request.to_json().encode()).hexdigest(),
            "result_sha256": hashlib.sha256((workspace / "replay-result.json").read_bytes()).hexdigest()}
        job = json.loads((args.cohort / "managed-execution" / f"{case}-finished-job.json").read_text())
        baseline = next(o for o in job["diffusion"]["outcomes"] if o["seed"] == seed)
        row["unchanged_original_outcome"] = baseline["status"]
        if summary.passed:
            bundle = workspace / "bundle"
            payloads = None
            if preparation_options(request) is not None:
                original_bundle = (args.cohort / "managed-execution/data/projects" /
                    job["workspaceId"] / job["outputDir"] / "result").resolve()
                index = json.loads((original_bundle / "bundle.json").read_text())
                entry = next(c for c in index["candidates"] if c["seed"] == seed)

                def read(reference, root=original_bundle):
                    return read_private_artifact(root, ArtifactReference(**reference),
                        label="original preparation custody", max_bytes=16_000_000,
                        encoding="utf-8").encode("utf-8")

                provenance = json.loads(read(entry["provenance"]))
                payloads = {a["role"]: read({"path": a["path"], "sha256": a["sha256"]})
                    for a in provenance["artifacts"] if a["role"] in PUBLICATION_ROLES}
                # Audit the frozen original snapshot before explicitly binding
                # unchanged preparation bytes to the new post-hoc repair policy.
                audit_preparation_payloads(original, payloads)
                custody = workspace / "retained-preparation"
                custody.mkdir(exist_ok=True)
                for role, data in payloads.items():
                    path = custody / role
                    if path.exists():
                        if path.read_bytes() != data:
                            raise ValueError("retained preparation custody changed")
                    else:
                        path.write_bytes(data)
                row["preparation_snapshot_rebinding"] = {
                    "scope": "new post-hoc refinement policy only; original snapshot retained unchanged",
                    "original_sha256": hashlib.sha256(payloads["input-repair-request"]).hexdigest(),
                    "new_sha256": hashlib.sha256(request.to_json().encode()).hexdigest()}
                payloads["input-repair-request"] = request.to_json().encode()
                audit_preparation_payloads(request, payloads)
            if not bundle.exists():
                publish_diffusion_bundle(request, validated, workspace=workspace, destination_bundle=bundle,
                    input_preparation_artifacts=payloads)
            # Existing complete bundles are independently read, never adopted or
            # overwritten. A crash-left/incomplete reservation fails readback.
            published = read_published_bundle(bundle, request=request, source_root=workspace, work_parent=workspace)
            if published.selected().pdb_bytes != (workspace / result.candidates[0].coordinate_artifact.path).read_bytes():
                raise ValueError("published readback differs from independently validated coordinates")
            row["bundle_readback_passed"] = True
        else:
            row["bundle_readback_passed"] = None
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if k not in {"metrics", "warnings"}}), flush=True)
    report = {"scope": "independently revalidated post-hoc repairs and public bundle readback; original cohort unchanged",
        "cases": rows, "passed": sum(r["passed"] for r in rows), "total": len(rows)}
    with (controls / "independent-audit.json").open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
