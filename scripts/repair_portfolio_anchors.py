#!/usr/bin/env python3
"""Explicit private bounded-anchor escalation after a generated-only portfolio attempt.

Uses retained publication-evidence.json (including failed geometry) and its
workspace; never samples, mutates a frozen request, or adopts an old bundle.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path


def main():
    from dvbfixer.model.diffusion.portfolio_anchors import repair_portfolio_anchors
    from dvbfixer.model.diffusion.portfolio_contract import (
        PortfolioAnchorArtifacts,
        PortfolioAnchoredEvidence,
        PortfolioEvidence,
    )
    from dvbfixer.model.diffusion.portfolio_evidence import (
        inspect_portfolio_anchor_inputs,
        json_bytes,
        strict_json,
    )
    from dvbfixer.model.diffusion.portfolio_publication import (
        publish_portfolio_bundle,
        read_portfolio_bundle,
        resume_portfolio_bundle,
    )
    from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--refinement-seed", type=int, default=7)
    parser.add_argument("--max-iterations", type=int, default=250)
    parser.add_argument("--restart-count", type=int, default=4)
    parser.add_argument("--platform", choices=("CPU", "Reference"), default="CPU")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    out = args.output.resolve(strict=True)

    def write(name, data):
        from dvbfixer.model.diffusion.contract import ArtifactReference

        with (out / name).open("xb") as stream:
            stream.write(data)
        return ArtifactReference(name, hashlib.sha256(data).hexdigest())

    write("run-script.py", Path(__file__).read_bytes())
    write("protocol.json", json_bytes({"operation": "post-portfolio-bounded-anchors-v1",
        "refinement_seed": args.refinement_seed, "max_iterations": args.max_iterations,
        "restart_count": args.restart_count, "platform": args.platform, "new_sampling": False,
        "source_workspace": str(args.workspace.resolve()), "public_profile_enabled": False}))
    try:
        proof = PortfolioEvidence.from_dict(strict_json(read_private_file_bytes(
            args.workspace.resolve(strict=True), "publication-evidence.json", max_bytes=128_000)))
        inputs = inspect_portfolio_anchor_inputs(proof, args.workspace)
        write("before-audit.json", json_bytes(inputs.validation_record()))
        # Copy verified bytes only; original failed trials/bundles remain untouched.
        values = proof.to_dict()
        roles = {}
        for index, (role, data) in enumerate(inputs.payloads.items()):
            roles[role] = write(f"input-{index:03d}.payload", data)
        for role, reference in roles.items():
            if role.startswith("parent-"):
                _, index, name = role.split("-", 2)
                values["parents"][int(index)][name] = asdict(reference)
            elif role.startswith("input-"):
                name = {"input-original": "original", "input-prepared": "prepared",
                        "input-repair-report": "report", "input-repair-request": "request_snapshot"}[role]
                values["preparation"][name] = asdict(reference)
            else:
                values[role] = asdict(reference)
        source, selected, stage1 = (inputs.payloads[r].decode() for r in ("source", "selected", "candidate"))

        import dvbfixer.model.diffusion.portfolio_anchors as module

        snapshots = {}
        for role, name in (("adapter_source", "portfolio_anchors.py"), ("kernel_source", "boundary_refinement.py"),
                           ("conditioning_source", "geometry_conditioning.py"),
                           ("acceptance_source", "anchor_acceptance.py"), ("policy_source", "anchor_policy.py")):
            snapshots[role] = write("anchor-" + name, Path(module.__file__).with_name(name).read_bytes())

        def trial(index, text, report):
            write(f"trial-{index:03d}.pdb", text.encode())
            write(f"trial-{index:03d}.json", json_bytes(report))

        result = repair_portfolio_anchors(inputs.request, source, selected, stage1,
            refinement_seed=args.refinement_seed, max_iterations=args.max_iterations,
            restart_count=args.restart_count, platform_name=args.platform, retain_trial=trial)
        candidate = write("result.pdb", result.coordinate_text.encode())
        report = write("anchor-report.json", json_bytes(result.report))
        values.update(schema_version=2, candidate=asdict(candidate), anchor=asdict(PortfolioAnchorArtifacts(
            stage1=roles["candidate"], report=report, **snapshots)))
        anchored = PortfolioAnchoredEvidence.from_dict(values)
        write("publication-evidence.json", anchored.to_json().encode())
        passed = result.report["audit"]["summary"]["passed"]
        summary = {"passed_geometry_and_custody": passed, "publication_qualified": False,
            "new_sampling": False, "anchor_audit": result.report["audit"], "bundle_published": False}
        if passed and args.bundle is not None:
            publish_portfolio_bundle(anchored, out, args.bundle)
            published = read_portfolio_bundle(args.bundle)
            resume_portfolio_bundle(args.bundle, anchored)
            summary.update(bundle_published=True, exact_resume_passed=True,
                           dat_added_atoms=len(published.dat.added_atoms))
        write("summary.json", json_bytes(summary))
        print(json.dumps(summary), flush=True)
    except Exception as exc:
        write("failure.json", json_bytes({"type": type(exc).__name__, "error": str(exc)}))
        raise


if __name__ == "__main__":
    main()
