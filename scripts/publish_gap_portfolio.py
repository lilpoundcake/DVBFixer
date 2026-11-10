#!/usr/bin/env python3
"""Publish or exactly resume a private portfolio, without sampling/refinement.

Creates a new evidence workspace. Original run/cohort artifacts are read-only;
raw donor and source-preparation custody are restored from verified retention.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from dataclasses import asdict
from pathlib import Path, PurePosixPath


def main():
    from dvbfixer.model.diffusion.contract import ArtifactReference, DiffusionRequest, RunnerResult
    from dvbfixer.model.diffusion.input_preparation import PUBLICATION_ROLES, preparation_options
    from dvbfixer.model.diffusion.portfolio_contract import (
        PortfolioEvidence,
        PortfolioParentArtifacts,
        PortfolioPreparationArtifacts,
    )
    from dvbfixer.model.diffusion.portfolio_evidence import json_bytes, strict_json
    from dvbfixer.model.diffusion.portfolio_publication import (
        publish_portfolio_bundle,
        read_portfolio_bundle,
        resume_portfolio_bundle,
    )
    from dvbfixer.model.diffusion.private_artifact import (
        read_private_artifact_bytes,
        read_private_file_bytes,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--portfolio", type=Path, required=True)
    parser.add_argument("--repair", type=Path)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--resume", action="store_true", help="Read-only exact reuse of a complete matching bundle")
    args = parser.parse_args()
    if args.portfolio.is_symlink() or (args.repair is not None and args.repair.is_symlink()):
        raise ValueError("portfolio/repair inputs must be regular directories")
    portfolio = args.portfolio.resolve(strict=True)
    repair = args.repair.resolve(strict=True) if args.repair is not None else None
    args.work.mkdir(parents=True, exist_ok=False)
    workspace = args.work.resolve(strict=True)

    def read(root, name, max_bytes=16_000_000):
        return read_private_file_bytes(root, name, max_bytes=max_bytes)

    def write(name, data):
        with (workspace / name).open("xb") as stream:
            stream.write(data)
        return ArtifactReference(name, hashlib.sha256(data).hexdigest())

    request_bytes = read(portfolio, "request.json")
    request = DiffusionRequest.from_dict(strict_json(request_bytes))
    report_bytes = read(portfolio, "result.json")
    report = strict_json(report_bytes)
    protocol = strict_json(read(portfolio, "protocol.json"))
    if protocol["case"] != args.case:
        raise ValueError("requested case differs from retained selection protocol")
    values = {"request": write("request.json", request_bytes), "source": write("source.pdb", read(portfolio, "source.pdb")),
        "selected": write("selected.pdb", read(portfolio, "result.pdb")),
        "candidate": write("candidate.pdb", read(repair or portfolio, "result.pdb")),
        "search_protocol": write("protocol.json", read(portfolio, "protocol.json")),
        "search_report": write("search.json", report_bytes),
        "selection_source": write("selection.py", read(portfolio, "gap_portfolio.py", 256_000)),
        "choices": list(report["selected"]["choices"])}
    cohort_root = args.cohort.resolve(strict=True)
    cohort_index = strict_json(read(cohort_root, "cohort-final-outcomes.json", 1_000_000))
    retained_case = next(row for row in cohort_index["rows"] if row["pdb"] == args.case)
    job_reference = ArtifactReference(retained_case["job_record"], retained_case["job_record_sha256"])
    job = strict_json(read_private_artifact_bytes(cohort_root, job_reference, max_bytes=16_000_000))
    execution = (cohort_root / job_reference.path).parent
    run = execution / "data/projects" / job["workspaceId"] / job["outputDir"]
    if not run.resolve(strict=True).is_relative_to(execution):
        raise ValueError("retained case run escapes its execution root")
    manifest = {row["file"]: row for row in json.loads((run / "run-summary.json").read_text())["diagnostics"]["files"]}
    parent_entries = []
    with tarfile.open(run / "diagnostics.tar.gz", "r:gz") as archive:
        paths = [name for name in archive.getnames() if name.endswith("/request.json") and len(PurePosixPath(name).parts) == 2]
        if len(paths) != 1:
            raise ValueError("ambiguous retained parent request")
        base = PurePosixPath(paths[0]).parent

        def archived(name):
            member = archive.getmember(name)
            if not member.isfile() or member.size > 16_000_000:
                raise ValueError("unbounded/nonregular retained portfolio artifact")
            data = archive.extractfile(member).read()
            if len(data) != manifest[name]["bytes"] or hashlib.sha256(data).hexdigest() != manifest[name]["sha256"]:
                raise ValueError("retained portfolio artifact manifest mismatch")
            return data

        if archived(paths[0]) != request_bytes:
            raise ValueError("portfolio common request differs from frozen cohort custody")
        for i, row in enumerate(strict_json(read(portfolio, "parents.json"))):
            seed = row["seed"]
            if type(seed) is not int or seed < 0:
                raise ValueError("portfolio parent seed must be an integer")
            parent_dir = f"parent-{seed}"
            raw_result_bytes = read(portfolio, f"{parent_dir}/raw-result.json")
            prefix = base / f"candidates/protpardelle-1c-cc94-mps-v1/seed-{seed}"
            if archived((prefix / "raw-result.json").as_posix()) != raw_result_bytes:
                raise ValueError("portfolio parent runner result differs from original custody")
            result = RunnerResult.from_dict(strict_json(raw_result_bytes))
            if len(result.candidates) != 1:
                raise ValueError("portfolio donor needs one original candidate")
            raw_reference = result.candidates[0].coordinate_artifact
            raw_bytes = archived((base / raw_reference.path).as_posix())
            if hashlib.sha256(raw_bytes).hexdigest() != raw_reference.sha256:
                raise ValueError("original raw donor digest mismatch")
            references = {"runner_result": write(f"parent-{i}-runner.json", raw_result_bytes),
                          "raw": write(f"parent-{i}-raw.pdb", raw_bytes)}
            for role, name, expected in (("coordinates", "stage1.pdb", "coordinate_sha256"),
                                        ("sampler_trace", "trace.json", "sampler_trace_sha256"),
                                        ("staged_evidence", "evidence.json", "evidence_sha256")):
                data = read(portfolio, f"{parent_dir}/{name}")
                if hashlib.sha256(data).hexdigest() != row[expected]:
                    raise ValueError("portfolio parent evidence differs from declared custody")
                references[role] = write(f"parent-{i}-{name}", data)
            parent_entries.append(PortfolioParentArtifacts(row["candidate_id"], seed, **references))
    values["parents"] = tuple(parent_entries)
    if repair is not None:
        values.update(repair_report=write("repair-report.json", read(repair, "repair-report.json")),
            repair_selector_source=write("repair-source.py", read(repair, "gap_portfolio.py", 256_000)),
            repair_kernel_source=write("kernel-source.py", read(repair, "boundary_refinement.py", 256_000)),
            repair_conditioning_source=write("conditioning-source.py", read(repair, "geometry_conditioning.py", 256_000)))
    if preparation_options(request) is not None:
        original_bundle = (run / "result").resolve(strict=True)
        index = strict_json(read(original_bundle, "bundle.json", 128_000))
        candidate = next(c for c in index["candidates"] if c["seed"] == parent_entries[0].seed)
        provenance = strict_json(read_private_artifact_bytes(original_bundle,
            ArtifactReference(**candidate["provenance"]), max_bytes=16_000_000))
        retained = {}
        for artifact in provenance["artifacts"]:
            role = artifact["role"]
            if role in PUBLICATION_ROLES:
                data = read_private_artifact_bytes(original_bundle,
                    ArtifactReference(artifact["path"], artifact["sha256"]), max_bytes=16_000_000)
                retained[role] = write(role + ".payload", data)
        values["preparation"] = PortfolioPreparationArtifacts(retained["input-original"], retained["input-prepared"],
            retained["input-repair-report"], retained["input-repair-request"])
    # Strict decoding verifies choices/booleans/version before any publication.
    proof = PortfolioEvidence.from_dict({"schema_version": 1,
        **{key: ([asdict(p) for p in value] if key == "parents" else
                   asdict(value) if hasattr(value, "__dataclass_fields__") else value)
           for key, value in values.items()}})
    write("publication-evidence.json", proof.to_json().encode())
    if args.resume:
        published = resume_portfolio_bundle(args.bundle, proof)
    else:
        publish_portfolio_bundle(proof, workspace, args.bundle)
        published = read_portfolio_bundle(args.bundle)
        assert resume_portfolio_bundle(args.bundle, proof).candidate_bytes == published.candidate_bytes
    summary = {"kind": published.index.kind, "resume_key": published.index.resume_key,
        "candidate_sha256": published.evidence.candidate.sha256, "dat_added_atoms": len(published.dat.added_atoms),
        "preparation_verified": published.audit.preparation_verified, "readback_passed": True,
        "exact_resume_passed": True, "new_sampling": False, "new_refinement": False,
        "scientific_qualification": False, "public_profile_enabled": False, "resumed": args.resume}
    write("summary.json", json_bytes(summary))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
