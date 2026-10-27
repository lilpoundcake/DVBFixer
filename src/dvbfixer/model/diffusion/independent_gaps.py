"""Corrected local sampling per gap, with source-bound aggregate evidence."""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    RunnerResourceMetrics,
    RunnerResult,
    SamplerTrace,
    atom_identity_digest,
    fixed_mask_digest,
)
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.protein_context import independent_gap_requests
from dvbfixer.model.diffusion.trace import validate_sampler_trace_context, write_sampler_trace


def sample_independent_gaps(
    request_path: Path, output: Path,
    sample: Callable[[Path, Path], RunnerResult],
) -> RunnerResult:
    """Sample immutable source-only views, then merge only generated coordinates.

    Never condition a later gap on a previously generated prediction. Every
    sub-run proves the existing corrected-profile numerical/position contract.
    Full-complex refinement and independent validation follow the merged result.
    """
    request_path = request_path.resolve()
    workspace = request_path.parent
    request = DiffusionRequest.from_json(request_path.read_text())
    if request.gap_contexts is None or request.candidate_count != 1:
        raise ValueError("independent gap sampling requires prepared contexts and one seed")
    output.mkdir(parents=True, exist_ok=False)
    source = _read(workspace, request.normalized_pdb)
    fixed = set(request.fixed_atoms)
    source_records = {_identity(line): line for line in source.splitlines()
        if line.startswith(("ATOM  ", "HETATM"))}
    if not fixed <= set(source_records):
        raise ValueError("source omits requested fixed atoms")
    results: list[RunnerResult] = []
    traces = []
    coordinates = {}
    names = {}
    flags = []
    for index, single in enumerate(independent_gap_requests(request)):
        single_path = workspace / f"gap-{index}-seed-{request.seeds[0]}.json"
        single_path.write_text(single.to_json())
        directory = output / f"gap-{index}"
        directory.mkdir()
        result = sample(single_path, directory)
        if results:
            for field in ("engine_repository", "engine_revision", "checkpoint_sha256", "device", "precision", "framework", "framework_version"):
                if getattr(result.backend_provenance, field) != getattr(results[0].backend_provenance, field):
                    raise ValueError("independent gap runs changed model/runtime identity")
        if len(result.candidates) != 1:
            raise ValueError("independent gap sampler must return one candidate")
        candidate = result.candidates[0]
        if candidate.postprocessing_failures or candidate.seed != request.seeds[0]:
            raise ValueError("independent gap sampler returned incompatible provenance")
        trace = SamplerTrace.from_json(_read(workspace, candidate.sampler_trace_artifact))
        validate_sampler_trace_context(trace, single, result.backend_provenance)
        text = _read(workspace, candidate.coordinate_artifact)
        expected = set(single.generated_atoms)
        found = set()
        fixed_seen = set()
        for line in text.splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = _identity(line)
            if atom in source_records:
                if atom in fixed_seen or line != source_records[atom]:
                    raise ValueError("independent gap sampler changed an immutable source atom")
                fixed_seen.add(atom)
            if atom not in expected:
                if atom not in source_records:
                    raise ValueError("independent gap sampler added unrequested atoms")
                continue
            if atom in found or atom in coordinates:
                raise ValueError("independent gap coordinates overlap or repeat atom identities")
            found.add(atom)
            coordinates[atom] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            names[atom] = line[17:20].strip()
        if found != expected:
            raise ValueError("independent gap sampler omitted generated atoms")
        if fixed_seen != set(source_records):
            raise ValueError("independent gap sampler omitted immutable source atoms")
        for flag in result.backend_provenance.deterministic_flags:
            if flag.startswith(("cc89-position-features-sha256=", "cc89-first-layer-rope-calls=", "cc94-")):
                flags.append(f"gap-{index}:{flag}")
        results.append(result)
        traces.append(trace)
    if set(coordinates) != set(request.generated_atoms):
        raise ValueError("independent gap merge changed the generated mask")
    path = output / "candidate.pdb"
    path.write_text(materialize_candidate_pdb(source, request.generated_atoms,
        [coordinates[atom] for atom in request.generated_atoms],
        tuple(names[atom] for atom in request.generated_atoms), request), encoding="ascii")
    axis = tuple(dict.fromkeys(atom for trace in traces for atom in trace.atom_order))
    selected = {atom for trace in traces for atom in trace.represented_fixed_atoms}
    represented = tuple(atom for atom in request.fixed_atoms if atom in selected)
    metrics = RunnerResourceMetrics(
        wall_time_seconds=sum(result.resource_metrics.wall_time_seconds or 0 for result in results),
        model_load_seconds=sum(result.resource_metrics.model_load_seconds or 0 for result in results),
        peak_ram_bytes=max((result.resource_metrics.peak_ram_bytes or 0 for result in results), default=0),
        peak_vram_bytes=max((result.resource_metrics.peak_vram_bytes or 0 for result in results), default=0),
    )
    evidence = []
    for item in traces:
        assert item.protein_context_evidence is not None
        evidence.append(item.protein_context_evidence)
    trace = replace(traces[0], atom_order=axis, atom_order_sha256=atom_identity_digest(axis),
        fixed_mask_sha256=fixed_mask_digest(axis, represented), represented_fixed_atoms=represented,
        protein_context_evidence=None, gap_context_evidence=tuple(evidence),
        gap_atom_orders=tuple(trace.atom_order for trace in traces), conditioning_contexts=(), resource_metrics=metrics)
    reference = write_sampler_trace(workspace, output / "sampler-trace.json", trace)
    candidate = replace(results[0].candidates[0], coordinate_artifact=ArtifactReference(
        path.relative_to(workspace).as_posix(), hashlib.sha256(path.read_bytes()).hexdigest()),
        generated_atoms=request.generated_atoms,
        generated_residues=tuple(residue for gap in request.gaps for residue in gap.generated_residues),
        sampler_trace_artifact=reference)
    first = results[0]
    provenance = replace(first.backend_provenance,
        backend=first.backend_provenance.backend + "+independent-local-gaps-v1",
        deterministic_flags=tuple(flag for flag in first.backend_provenance.deterministic_flags
            if not flag.startswith(("cc89-position-features-sha256=", "cc89-first-layer-rope-calls=", "cc94-")))
            + tuple(flags) + (f"independent-local-gap-invocations={len(results)}",))
    validate_sampler_trace_context(trace, request, provenance)
    return replace(first, candidates=(candidate,), backend_provenance=provenance, resource_metrics=metrics)


def _read(workspace: Path, reference: ArtifactReference) -> str:
    # Reuse the core bounded/private artifact reader, including symlink and hash checks.
    from dvbfixer.model.diffusion.private_artifact import read_private_artifact

    return read_private_artifact(workspace, reference, label="independent gap artifact",
        max_bytes=16_000_000, encoding="utf-8")


def _identity(line: str) -> AtomIdentity:
    return AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
