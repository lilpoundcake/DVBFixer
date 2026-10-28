"""Shared generated-region refinement for external diffusion runners."""

from __future__ import annotations

import hashlib
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from openmm import OpenMMException

from dvbfixer.ffutils.geometry import ChiralityError
from dvbfixer.model.diffusion.boundary_refinement import (
    BOUNDARY_REFINEMENT_MAX_ITERATIONS,
    BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM,
    BOUNDARY_REFINEMENT_RESTART_COUNT,
    BOUNDARY_REFINEMENT_REVISION,
    BoundaryRefinementError,
    refine_generated_region,
)
from dvbfixer.model.diffusion.contract import (
    MAX_TRACE_BYTES,
    ArtifactReference,
    AtomIdentity,
    BackendOption,
    DiffusionRequest,
    RunnerCandidate,
    RunnerResourceMetrics,
    RunnerResult,
    SamplerTrace,
)
from dvbfixer.model.diffusion.private_artifact import read_private_artifact


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_private_artifact(
    workspace: Path,
    artifact: ArtifactReference,
    *,
    label: str,
    max_bytes: int,
    encoding: str,
) -> str:
    return read_private_artifact(workspace, artifact, label=label, max_bytes=max_bytes, encoding=encoding)


def _atom_identity(line: str) -> AtomIdentity:
    return AtomIdentity(
        line[21:22],
        line[22:26].strip(),
        line[26:27].strip(),
        line[12:16].strip(),
    )


def rewrite_generated_coordinates(
    pdb_text: str,
    coordinates: dict[AtomIdentity, np.ndarray],
) -> str:
    """Replace only generated coordinate columns in a candidate PDB."""
    expected = set(coordinates)
    seen: set[AtomIdentity] = set()
    output: list[str] = []
    for line in pdb_text.splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            output.append(line)
            continue
        identity = _atom_identity(line)
        if identity not in expected:
            output.append(line)
            continue
        if identity in seen:
            raise ValueError(f"candidate contains duplicate generated atom: {identity}")
        xyz = np.asarray(coordinates[identity], dtype=np.float64)
        fields = tuple(f"{float(value):8.3f}" for value in xyz)
        if xyz.shape != (3,) or not np.isfinite(xyz).all() or any(
            len(field) != 8 for field in fields
        ):
            raise ValueError(f"invalid refined coordinate for {identity}")
        output.append(line[:30] + "".join(fields) + line[54:])
        seen.add(identity)
    if seen != expected:
        raise ValueError(f"candidate omits {len(expected - seen)} generated atoms")
    return "".join(output)


def refine_runner_result(
    request_path: Path,
    raw_result: RunnerResult,
    output_dir: Path,
    *,
    platform_name: str,
) -> RunnerResult:
    """Refine one runner candidate while preserving its sampling provenance."""
    request_path = request_path.resolve()
    workspace = request_path.parent
    output_dir = output_dir.resolve()
    if not output_dir.is_relative_to(workspace):
        raise ValueError("refinement output directory must remain inside the request workspace")
    request = DiffusionRequest.from_json(request_path.read_text(encoding="utf-8"))
    if len(raw_result.candidates) != 1:
        raise ValueError("production refinement requires one candidate")
    raw_candidate = raw_result.candidates[0]
    raw_text = _read_private_artifact(
        workspace,
        raw_candidate.coordinate_artifact,
        label="coordinate artifact",
        max_bytes=1_000_000_000,
        encoding="ascii",
    )
    raw_trace = SamplerTrace.from_json(
        _read_private_artifact(
            workspace,
            raw_candidate.sampler_trace_artifact,
            label="sampler trace",
            max_bytes=MAX_TRACE_BYTES,
            encoding="utf-8",
        )
    )
    from dvbfixer.model.diffusion.refinement_policy import read_refinement_policy

    policy = read_refinement_policy(request.backend_options)
    if policy is not None:
        return _refine_staged(request, workspace, raw_result, raw_text, raw_trace, output_dir, platform_name, policy)
    movable_residues = tuple(dict.fromkeys(
        residue for gap in request.gaps for residue in gap.movable_junction_residues
    ))
    short_gaps = all(len(gap.generated_residues) <= 5 for gap in request.gaps)

    start = time.perf_counter()
    refined_text = raw_text
    refinement_platform = platform_name
    try:
        # One immutable-source / mutable-generated mask across all gaps avoids
        # misclassifying another generated region's D centre as a fixed defect.
        generated_residues = tuple(dict.fromkeys(
            residue for gap in request.gaps for residue in gap.generated_residues
        ))
        refinement_input = (_protein_only_refinement_text(refined_text)
            if request.heterogen_contexts else refined_text)
        refinement = refine_generated_region(
            refinement_input, generated_residues=generated_residues,
            generated_atoms=request.generated_atoms,
            max_iterations=500 if short_gaps else BOUNDARY_REFINEMENT_MAX_ITERATIONS,
            restart_count=BOUNDARY_REFINEMENT_RESTART_COUNT,
            perturbation_angstrom=0.25 if short_gaps else BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM,
            platform_name=platform_name, random_seed=request.seeds[0],
        )
        refinement_platform = refinement.platform
        refined_text = rewrite_generated_coordinates(refined_text, refinement.coordinates_angstrom)
    except (BoundaryRefinementError, ChiralityError, OpenMMException, ValueError) as exc:
        elapsed = time.perf_counter() - start
        rejected = replace(
            raw_candidate,
            warnings=(
                *raw_candidate.warnings,
                "Localized OpenMM boundary refinement failed; the raw sampler "
                f"candidate was retained for inspection. {type(exc).__name__}: {exc}",
            ),
            postprocessing_failures=(
                *raw_candidate.postprocessing_failures,
                "localized-openmm-boundary-refinement-failed",
            ),
        )
        raw_metrics = raw_result.resource_metrics
        return replace(
            raw_result,
            candidates=(rejected,),
            resource_metrics=RunnerResourceMetrics(
                wall_time_seconds=(raw_metrics.wall_time_seconds or 0.0) + elapsed,
                model_load_seconds=raw_metrics.model_load_seconds,
                peak_ram_bytes=raw_metrics.peak_ram_bytes,
                peak_vram_bytes=raw_metrics.peak_vram_bytes,
            ),
            message=(
                "localized OpenMM boundary refinement failed; raw sampler "
                f"candidate retained: {type(exc).__name__}: {exc}"
            ),
        )
    output_dir.mkdir(parents=True, exist_ok=False)
    refined_path = output_dir / "candidate.pdb"
    refined_path.write_text(refined_text, encoding="ascii")
    trace_path = output_dir / "sampler-trace.json"
    trace = replace(
        raw_trace,
        refinement_mode="localized-openmm-boundary-refinement",
        refinement_parameters=(
            BackendOption("platform", refinement_platform),
            BackendOption("gap_count", str(len(request.gaps))),
            BackendOption(
                "max_iterations",
                str(500 if short_gaps else BOUNDARY_REFINEMENT_MAX_ITERATIONS),
            ),
            BackendOption("restart_count", str(BOUNDARY_REFINEMENT_RESTART_COUNT)),
            BackendOption(
                "perturbation_angstrom",
                str(0.25 if short_gaps else BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM),
            ),
            *(BackendOption(name, str(value)) for name in ("initial_energy_kj_mol", "final_energy_kj_mol")
                if (value := getattr(refinement, name, None)) is not None),
        ),
        localized_refinement_residues=movable_residues,
        final_heavy_coordinate_operations=(
            *raw_trace.final_heavy_coordinate_operations,
            "localized-openmm-boundary-refinement",
        ),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=(raw_trace.resource_metrics.wall_time_seconds or 0.0)
            + (time.perf_counter() - start),
            model_load_seconds=raw_trace.resource_metrics.model_load_seconds,
            peak_ram_bytes=raw_trace.resource_metrics.peak_ram_bytes,
            peak_vram_bytes=raw_trace.resource_metrics.peak_vram_bytes,
        ),
    )
    trace_path.write_text(trace.to_json(), encoding="utf-8")
    relative_path = refined_path.resolve().relative_to(workspace).as_posix()
    candidate = RunnerCandidate(
        candidate_id=f"{raw_candidate.candidate_id}-refined",
        seed=raw_candidate.seed,
        coordinate_artifact=ArtifactReference(relative_path, _sha256(refined_path)),
        sampler_trace_artifact=ArtifactReference(
            trace_path.resolve().relative_to(workspace).as_posix(),
            _sha256(trace_path),
        ),
        generated_atoms=raw_candidate.generated_atoms,
        generated_residues=raw_candidate.generated_residues,
        raw_backend_score=raw_candidate.raw_backend_score,
        score_provenance=f"{BOUNDARY_REFINEMENT_REVISION}:unranked",
        warnings=raw_candidate.warnings,
        postprocessing_failures=raw_candidate.postprocessing_failures,
    )
    raw_metrics = raw_result.resource_metrics
    elapsed = time.perf_counter() - start
    platform_label = refinement.platform
    return replace(
        raw_result,
        candidates=(candidate,),
        backend_provenance=replace(
            raw_result.backend_provenance,
            backend=(
                f"{raw_result.backend_provenance.backend}"
                f"+{platform_label.lower()}-boundary-refinement"
            ),
            deterministic_flags=(
                *raw_result.backend_provenance.deterministic_flags,
                f"refinement-platform={platform_label}",
                f"refinement={BOUNDARY_REFINEMENT_REVISION}",
                f"refinement-restart-count={BOUNDARY_REFINEMENT_RESTART_COUNT}",
            ),
            known_nondeterministic_operations=(
                raw_result.backend_provenance.known_nondeterministic_operations
                if platform_label == "Reference"
                else (
                    *raw_result.backend_provenance.known_nondeterministic_operations,
                    f"OpenMM {platform_label} minimization",
                )
            ),
        ),
        resource_metrics=RunnerResourceMetrics(
            wall_time_seconds=(raw_metrics.wall_time_seconds or 0.0) + elapsed,
            model_load_seconds=raw_metrics.model_load_seconds,
            peak_ram_bytes=raw_metrics.peak_ram_bytes,
            peak_vram_bytes=raw_metrics.peak_vram_bytes,
        ),
    )


def _protein_only_refinement_text(pdb_text: str) -> str:
    """Remove retained chemistry from the temporary protein FF refinement input."""
    return "".join(
        line
        for line in pdb_text.splitlines(keepends=True)
        if not line.startswith(("HETATM", "CONECT", "LINK  "))
    )


def _refine_staged(request, workspace, raw_result, raw_text, raw_trace, output_dir, platform_name, policy):
    """Versioned path; historical requests continue through the unchanged refiner."""
    import json
    from dataclasses import asdict

    from dvbfixer.model.diffusion.refinement_evidence import StagedRefinementEvidence, bind_evidence
    from dvbfixer.model.diffusion.refinement_policy import POLICY_REVISION
    from dvbfixer.model.diffusion.refinement_stages import run_refinement_stages
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter, validate_staged_evidence

    start = time.perf_counter()
    if not output_dir.resolve().is_relative_to(workspace):
        raise ValueError("staged output directory escapes request workspace")
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "raw.pdb").write_text(raw_text, encoding="ascii")
    source_text = _read_private_artifact(workspace, request.normalized_pdb,
        label="source structure", max_bytes=1_000_000_000, encoding="ascii")
    source_assess = assessment_adapter(request, source_text)
    allowed = frozenset()
    notices = []
    platforms = []
    generated_residues = tuple(dict.fromkeys(r for gap in request.gaps for r in gap.generated_residues))
    short = all(len(gap.generated_residues) <= 5 for gap in request.gaps)

    def minimize_native(text, restraints):
        nonlocal allowed
        allowed = frozenset(item.atom for item in restraints)
        for reason in sorted({item.fallback_reason for item in restraints if item.fallback_reason}):
            if policy.bfactor_weighting:
                notices.append(f"Observed-anchor B weighting uses uniform fallback: {reason}.")
        result = refine_generated_region(
            _protein_only_refinement_text(text) if request.heterogen_contexts else text,
            generated_residues=generated_residues, generated_atoms=request.generated_atoms,
            observed_anchor_restraints=restraints,
            max_iterations=500 if short else BOUNDARY_REFINEMENT_MAX_ITERATIONS,
            restart_count=BOUNDARY_REFINEMENT_RESTART_COUNT,
            perturbation_angstrom=0.25 if short else BOUNDARY_REFINEMENT_PERTURBATION_ANGSTROM,
            platform_name=platform_name, random_seed=raw_result.candidates[0].seed,
        )
        platforms.append(result.platform)
        rewritten = rewrite_generated_coordinates(text, result.coordinates_angstrom)
        (output_dir / ("anchor-attempt.pdb" if restraints else "generated-attempt.pdb")).write_text(
            rewritten, encoding="ascii")
        return rewritten

    def minimize(text, restraints):
        try:
            return minimize_native(text, restraints)
        except (OpenMMException, ChiralityError) as exc:
            # Neither a configuration error nor an arbitrary engine exception
            # is proof of the specifically typed numerical recovery condition.
            raise RuntimeError(f"{type(exc).__name__}: {exc}") from exc

    staged = run_refinement_stages(raw_text, source_text, request,
        assess=lambda text: source_assess(text, allowed), minimize=minimize,
        relax_anchors=policy.relax_anchors,
        bfactor_source="auto" if policy.bfactor_weighting else "disabled",
        notify=lambda line: notices.append(line.removeprefix("WARNING: ")))

    def write_text(name, text):
        path = output_dir / name
        path.write_text(text, encoding="ascii")
        return ArtifactReference(path.relative_to(workspace).as_posix(), _sha256(path))

    raw_ref = write_text("raw.pdb", raw_text)
    first_ref = None
    reports = []
    for index, attempt in enumerate(staged.attempts):
        report = asdict(attempt)
        text = report.pop("output_text")
        if text is not None:
            reference = write_text(f"stage-{index}.pdb", text)
            report["coordinate_artifact"] = asdict(reference)
            if index == 1:
                first_ref = reference
        reports.append(report)
    candidate = raw_result.candidates[0]
    selected_stage = staged.selected_stage
    coordinate_text = staged.coordinate_text
    refinement_failed = selected_stage is None
    if refinement_failed:
        # Diagnostic publication is not scientific acceptance. Keep the last
        # completed processing result even when its geometry is rejected; raw
        # and stage-1 coordinates remain digest-bound evidence in the bundle.
        completed = next((attempt for attempt in reversed(staged.attempts)
            if attempt.stage != "raw" and attempt.output_text is not None), None)
        notices.extend(f"Refinement {attempt.stage} failed: {attempt.error}"
            for attempt in staged.attempts if attempt.error is not None)
        if completed is not None:
            selected_stage = completed.stage
            coordinate_text = completed.output_text
            notices.append(f"Geometry-gated refinement failed; last completed {selected_stage} "
                "candidate retained for inspection, not accepted. Raw coordinates are retained separately.")
        else:
            notices.append("Geometry-gated refinement failed; no processed coordinates were produced. "
                "Raw candidate retained for inspection, not accepted.")
    # Failure/cancellation diagnostics remain in the run even when no candidate
    # qualifies. Nonfinite diagnostics are not serialized as invalid JSON.
    (output_dir / "refinement-report.json").write_text(json.dumps(
        {"policy": POLICY_REVISION, "attempts": reports, "notices": notices},
        ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if selected_stage is None:
        for notice in dict.fromkeys(notices):
            print(f"WARNING: {notice}")
        return replace(raw_result, candidates=(replace(candidate,
            warnings=(*candidate.warnings, *notices),
            postprocessing_failures=(*candidate.postprocessing_failures, "geometry-gated-refinement-failed")),),
            resource_metrics=replace(raw_result.resource_metrics,
                wall_time_seconds=(raw_result.resource_metrics.wall_time_seconds or 0.0) + time.perf_counter() - start))
    assert coordinate_text is not None
    selected_ref = write_text("candidate.pdb", coordinate_text)
    first = staged.attempts[1] if len(staged.attempts) > 1 else None
    evidence = StagedRefinementEvidence(
        request_sha256=hashlib.sha256(request.to_json().encode()).hexdigest(), raw=raw_ref,
        selected=selected_ref, selected_stage=selected_stage, stage1=first_ref,
        stage1_input_sha256=first.input_sha256 if first is not None else None,
        stage1_failure=first.numerical_failure if first is not None else None,
        stage2_input_sha256=(staged.attempts[2].input_sha256
            if selected_stage == "generated-and-anchors" else None),
        report=ArtifactReference((output_dir / "refinement-report.json").relative_to(workspace).as_posix(),
            _sha256(output_dir / "refinement-report.json")),
    )
    audit = validate_staged_evidence(request, evidence, source_text=source_text, selected=selected_ref,
        load_text=lambda reference: _read_private_artifact(workspace, reference,
            label="staged coordinate", max_bytes=1_000_000_000, encoding="ascii"))
    notices.extend(audit.warnings)
    for notice in dict.fromkeys(notices):
        print(f"WARNING: {notice}")
    evidence_path = output_dir / "refinement-evidence.json"
    evidence_path.write_text(evidence.to_json(), encoding="utf-8")
    evidence_ref = ArtifactReference(evidence_path.relative_to(workspace).as_posix(), _sha256(evidence_path))
    elapsed = time.perf_counter() - start
    trace = bind_evidence(replace(raw_trace,
        refinement_mode="none" if selected_stage == "raw" else POLICY_REVISION,
        refinement_parameters=(BackendOption("selected_stage", selected_stage),
            BackendOption("selection_reason", "last-completed-rejected-stage" if refinement_failed
                else "geometry-accepted-stage")),
        localized_refinement_residues=() if selected_stage == "raw" else (
            *generated_residues,
            *(staged.selected_anchors if selected_stage == "generated-and-anchors" else ())),
        final_heavy_coordinate_operations=raw_trace.final_heavy_coordinate_operations + (
            () if selected_stage == "raw" else (POLICY_REVISION,)),
        resource_metrics=replace(raw_trace.resource_metrics,
            wall_time_seconds=(raw_trace.resource_metrics.wall_time_seconds or 0.0) + elapsed),
    ), evidence_ref)
    trace_path = output_dir / "sampler-trace.json"
    trace_path.write_text(trace.to_json(), encoding="utf-8")
    candidate = replace(candidate, candidate_id=f"{candidate.candidate_id}-staged",
        coordinate_artifact=selected_ref,
        sampler_trace_artifact=ArtifactReference(trace_path.relative_to(workspace).as_posix(), _sha256(trace_path)),
        warnings=tuple(dict.fromkeys((*candidate.warnings, *notices))),
        postprocessing_failures=(*candidate.postprocessing_failures,
            *(("geometry-gated-refinement-failed",) if refinement_failed else ()),
            *audit.hard_gate_failures),
        score_provenance=f"{POLICY_REVISION}:unranked")
    provenance = raw_result.backend_provenance
    if platforms:
        provenance = replace(provenance, deterministic_flags=(*provenance.deterministic_flags,
            f"refinement={POLICY_REVISION}", f"refinement-platform={platforms[-1]}"),
            known_nondeterministic_operations=tuple(dict.fromkeys((
                *provenance.known_nondeterministic_operations,
                *((f"OpenMM {platforms[-1]} minimization",) if platforms[-1] != "Reference" else ()),
            ))))
    return replace(raw_result, candidates=(candidate,), backend_provenance=provenance,
        resource_metrics=replace(raw_result.resource_metrics,
            wall_time_seconds=(raw_result.resource_metrics.wall_time_seconds or 0.0) + elapsed))
