"""Replay an archived raw PDB through staged refinement, validation and resume.

No neural sampler is executed. Trace evidence is explicitly incomplete and the
profile/provenance identify archival coordinate replay, not sampler qualification.
Original source, masks and seeds remain unchanged. Sampler-only context plans
remain in archived-request.json, not the explicitly non-sampling replay request.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from dataclasses import replace
from importlib import import_module
from pathlib import Path

from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    BackendOption,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResult,
)
from dvbfixer.model.diffusion.pipeline import publish_diffusion_bundle
from dvbfixer.model.diffusion.publication import read_published_bundle
from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy
from dvbfixer.model.diffusion.runner_refinement import refine_runner_result
from dvbfixer.model.diffusion.trace import build_sampler_trace, write_sampler_trace
from dvbfixer.model.diffusion.validate import build_validated_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--raw", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--uniform", action="store_true")
    parser.add_argument("--no-relax-anchors", action="store_true")
    args = parser.parse_args()
    archived_bytes = args.request.read_bytes()
    request = DiffusionRequest.from_json(archived_bytes.decode())
    source_bytes = args.source.read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("source does not match the archived request")
    if request.candidate_count != 1:
        raise ValueError("replay requires a single-seed archived request")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = root / request.normalized_pdb.path
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_bytes(source_bytes)
    profile = "archived-staged-refinement-replay-v1"
    (root / "archived-request.json").write_bytes(archived_bytes)
    request = replace(request, protein_context=None, gap_contexts=None, backend_options=(
        *(option for option in request.backend_options
          if option.name != "profile" and not option.name.startswith("geometry_refinement.")),
        BackendOption("profile", profile),
        BackendOption("archived_request_sha256", hashlib.sha256(archived_bytes).hexdigest()),
        *RefinementPolicy(not args.no_relax_anchors, not args.uniform).backend_options(),
    ))
    request_path = root / "request.json"
    request_path.write_text(request.to_json())
    raw_bytes = args.raw.read_bytes()
    raw_path = root / "archived-raw.pdb"
    raw_path.write_bytes(raw_bytes)
    engine = "builtin://archived-coordinate-replay"
    fixed = request.fixed_atoms
    axis = (*fixed, *request.generated_atoms)
    trace = build_sampler_trace(profile=profile, engine_repository=engine, engine_revision="v1",
        patch_identity="no-neural-sampling", atom_order=axis, fixed_atoms=request.fixed_atoms,
        represented_fixed_atoms=fixed, device="cpu", fallback_disabled=True,
        denoising_update_count=None, final_fixed_coordinate_restoration=True,
        sampler_evidence_complete=False)
    trace_ref = write_sampler_trace(root, root / "archived-trace.json", trace)
    candidate = RunnerCandidate(candidate_id="archive-replay", seed=request.seeds[0],
        coordinate_artifact=ArtifactReference(raw_path.name, hashlib.sha256(raw_bytes).hexdigest()),
        sampler_trace_artifact=trace_ref, generated_atoms=request.generated_atoms,
        generated_residues=tuple(r for gap in request.gaps for r in gap.generated_residues),
        raw_backend_score=None, score_provenance="archived-raw:no-new-neural-sampling")
    result = RunnerResult(schema_version=request.schema_version, status=DiffusionStatus.SUCCESS, candidates=(candidate,),
        runner_diagnostics=RunnerDiagnostics(exit_code=0, timed_out=False),
        backend_provenance=BackendProvenance(backend="archived-coordinate-replay", runner_protocol_version=2,
            engine_repository=engine, engine_revision="v1", device="cpu"))
    # Capture the exact working-tree implementation before native execution;
    # a repository commit alone cannot identify uncommitted experimental code.
    import_module("dvbfixer.model.diffusion.staged_runtime")

    snapshot = root / "code"
    snapshot.mkdir()
    digests = {}
    paths = [Path(__file__), *(Path(module.__file__) for name, module in tuple(sys.modules.items())
        if name.startswith("dvbfixer.model.diffusion.") and getattr(module, "__file__", None))]
    for code_path in paths:
        shutil.copy2(code_path, snapshot / code_path.name)
        digests[code_path.name] = hashlib.sha256(code_path.read_bytes()).hexdigest()
    (root / "launch.json").write_text(json.dumps({"code_sha256": digests,
        "raw_sha256": hashlib.sha256(raw_bytes).hexdigest(), "source_sha256": request.normalized_pdb.sha256,
        "request_sha256": hashlib.sha256(request.to_json().encode()).hexdigest(),
        "platform": "CPU", "solvent": False, "sampler_qualification": False}, indent=2) + "\n")
    refined = refine_runner_result(request_path, result, root / "refined", platform_name="CPU")
    (root / "runner-result.json").write_text(refined.to_json())
    validated = build_validated_result(request, refined, workspace=root)
    publish_diffusion_bundle(request, validated, workspace=root, destination_bundle=root / "bundle")
    resumed = read_published_bundle(root / "bundle", request=request, source_root=root, work_parent=root)
    report = {"status": resumed.status, "sampler_qualification": False,
        "raw_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "summary": json.loads(validated.to_json())["validation_summaries"][0]}
    (root / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    shutil.copy2(__file__, root / "replay-script.py")
    print(json.dumps(report, ensure_ascii=False))
    raise SystemExit(0 if resumed.status == "success" else 1)


if __name__ == "__main__":
    main()
