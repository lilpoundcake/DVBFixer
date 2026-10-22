"""Source-bound per-gap execution and composite provenance regressions."""

from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    BackendOption,
    BackendProvenance,
    DiffusionContractError,
    DiffusionRequest,
    DiffusionStatus,
    GapKind,
    GapRegion,
    ResidueIdentity,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResult,
    SamplerTrace,
    TargetInterval,
)
from dvbfixer.model.diffusion.independent_gaps import sample_independent_gaps
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.protein_context import (
    MPS_CC89_PROFILE,
    context_evidence,
    independent_gap_requests,
    prepare_mps_cc89_request,
    sampling_entities,
    source_context_position_features,
    validate_context_capability,
    validate_context_source,
)
from dvbfixer.model.diffusion.scope import canonical_heavy_atom_identities
from dvbfixer.model.diffusion.trace import (
    build_sampler_trace,
    validate_sampler_trace_context,
    write_sampler_trace,
)
from tests.test_protein_context import context_request


def prepared():
    request, pdb = context_request()
    terminal = tuple(ResidueIdentity("D", str(i)) for i in (6, 7))
    request = replace(request,
        target_sequences=(replace(request.target_sequences[0], sequence="AAAAAAA"), request.target_sequences[1]),
        sequence_placements=(replace(request.sequence_placements[0], target_length=7), request.sequence_placements[1]),
        gaps=(*request.gaps, GapRegion("D", TargetInterval(5, 7), ResidueIdentity("D", "5"), None,
            terminal, (ResidueIdentity("D", "5"), *terminal), gap_kind=GapKind.C_TERMINAL)),
        generated_atoms=(*request.generated_atoms, *(atom for residue in terminal
            for atom in canonical_heavy_atom_identities(residue, "ALA"))),
        backend_options=(BackendOption("profile", MPS_CC89_PROFILE),))
    return prepare_mps_cc89_request(request, pdb, "auto"), pdb


def test_independent_contexts_cover_each_gap_without_fabricating_other_missing_regions():
    request, pdb = prepared()
    assert request.gap_contexts is not None and len(request.gap_contexts) == 2
    assert DiffusionRequest.from_json(request.to_json()) == request
    validate_context_source(request, pdb)
    validate_context_capability(request, "protpardelle")
    for single in independent_gap_requests(request):
        assert prepare_mps_cc89_request(single, pdb, "auto") == single
        assert single.fixed_atoms == request.fixed_atoms
        assert not set(single.generated_atoms) & {item.atom for item in single.protein_context.fixed_conditioning}
        represented = {residue for entity in sampling_entities(single) for residue in entity.residues}
        other = {residue for gap in request.gaps if gap != single.gaps[0] for residue in gap.generated_residues}
        assert not represented & other
    assert {span.chain for plan in request.gap_contexts for span in plan.spans} == {"D", "d"}


def test_independent_context_coordinates_are_recomputed_not_trusted():
    request, pdb = prepared()
    plan = request.gap_contexts[0]
    item = replace(plan.fixed_conditioning[0], coordinate_angstrom=(999.0, 0.0, 0.0))
    tampered = replace(request, gap_contexts=(replace(plan, fixed_conditioning=(item, *plan.fixed_conditioning[1:])), request.gap_contexts[1]))
    with pytest.raises(ValueError, match="coordinates do not match"):
        validate_context_source(tampered, pdb)


@pytest.mark.parametrize("mutation", [None, "duplicate-fixed", "missing-generated"])
def test_samples_merge_only_generated_atoms_and_prove_each_invocation(tmp_path: Path, mutation):
    request, pdb = prepared()
    (tmp_path / "input.pdb").write_bytes(pdb)
    path = tmp_path / "request.json"
    path.write_text(request.to_json())
    calls = []

    def sample(path, output):
        single = DiffusionRequest.from_json(path.read_text())
        calls.append(single)
        selected = tuple(item.atom for item in single.protein_context.fixed_conditioning)
        axis = (*single.generated_atoms, *reversed(selected))
        patch = "sha256:627891e28d5055cb0d903f542af695569d7133b0480cab8f25d154dc8ccc78c9;cc89-chain-phase=neutral-v1;cc89-position-features=source-physical-v1"
        trace = build_sampler_trace(profile=MPS_CC89_PROFILE, engine_repository="test://cc89", engine_revision="revision",
            patch_identity=patch, atom_order=axis, fixed_atoms=single.fixed_atoms, represented_fixed_atoms=selected,
            device="mps", fallback_disabled=True, denoising_update_count=500, final_fixed_coordinate_restoration=True,
            protein_context_evidence=context_evidence(single, axis))
        ref = write_sampler_trace(tmp_path, output / "trace.json", trace)
        coords = np.asarray([[int(atom.residue_number) * 3, 20, i] for i, atom in enumerate(single.generated_atoms)], dtype=float)
        candidate = output / "candidate.pdb"
        candidate.write_text(materialize_candidate_pdb(pdb.decode(), single.generated_atoms, coords,
            ("ALA",) * len(coords), single))
        if mutation == "duplicate-fixed":
            candidate.write_text(candidate.read_text() + pdb.decode().splitlines(keepends=True)[0])
        if mutation == "missing-generated":
            first = single.generated_atoms[0]
            candidate.write_text("".join(line for line in candidate.read_text().splitlines(keepends=True)
                if not (line.startswith("ATOM") and line[21] == first.chain
                    and line[22:26].strip() == first.residue_number and line[12:16].strip() == first.atom_name)))
        positions, chains = source_context_position_features(single)
        digest = hashlib.sha256(np.asarray((positions, chains), dtype="<i8").tobytes()).hexdigest()
        provenance = BackendProvenance("cc89", 4, "test://cc89", "revision", device="mps",
            deterministic_flags=(f"cc89-position-features-sha256={digest}", "cc89-first-layer-rope-calls=1000"))
        return RunnerResult(schema_version=4, status=DiffusionStatus.SUCCESS,
            candidates=(RunnerCandidate("sample", 7, ArtifactReference(candidate.relative_to(tmp_path).as_posix(),
                hashlib.sha256(candidate.read_bytes()).hexdigest()), ref, single.generated_atoms,
                single.gaps[0].generated_residues, None, "unranked"),),
            runner_diagnostics=RunnerDiagnostics(0, False), backend_provenance=provenance)

    if mutation is not None:
        with pytest.raises(ValueError, match="immutable source atom|omitted generated atoms"):
            sample_independent_gaps(path, tmp_path / "raw", sample)
        return
    result = sample_independent_gaps(path, tmp_path / "raw", sample)
    assert len(calls) == 2
    candidate = result.candidates[0]
    trace = SamplerTrace.from_json((tmp_path / candidate.sampler_trace_artifact.path).read_text())
    validate_sampler_trace_context(trace, request, result.backend_provenance)
    merged = (tmp_path / candidate.coordinate_artifact.path).read_text()
    for line in pdb.decode().splitlines():
        if line.startswith("ATOM"):
            assert line in merged.splitlines()
    assert set(candidate.generated_atoms) == set(request.generated_atoms)
    assert len(trace.gap_context_evidence) == 2
    flags = tuple(flag for flag in result.backend_provenance.deterministic_flags if not flag.startswith("gap-1:"))
    with pytest.raises(DiffusionContractError, match="provenance mismatch"):
        validate_sampler_trace_context(trace, request, replace(result.backend_provenance, deterministic_flags=flags))
    with pytest.raises(DiffusionContractError, match="actual sampled atom axis"):
        replace(trace, gap_atom_orders=None)


def test_chirality_is_a_minimum_admission_gate(monkeypatch):
    import dvbfixer.model.diffusion.boundary_refinement as boundary

    monkeypatch.setattr(boundary, "_state_is_finite", lambda _state: True)
    monkeypatch.setattr(boundary, "find_d_residues", lambda *_args: [("C", "13", "LEU", -0.1)])
    class State:
        def getPositions(self):
            return None
    assert not boundary._acceptable_minimum(None, State())
    monkeypatch.setattr(boundary, "find_d_residues", lambda *_args: [])
    assert boundary._acceptable_minimum(None, State())
    assert not boundary._acceptable_minimum(None, None)


def test_real_5j13_prepares_all_ten_gaps_without_sampling_or_dropping_tails(tmp_path):
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request

    fixture = Path(__file__).parent / "fixtures/protein_context/5j13_selected.pdb"
    source = tmp_path / "input.pdb"
    source.write_text("".join(line for line in fixture.read_text().splitlines(keepends=True)
        if not line.startswith(("HETATM", "CONECT", "LINK  "))))
    request = build_cli_diffusion_request(source, _target_sequences(source, None), seeds=(7,), profile=MPS_CC89_PROFILE)
    assert len(request.gaps) == 10
    assert sum(len(gap.generated_residues) for gap in request.gaps) == 124
    prepared_request = prepare_mps_cc89_request(request, source.read_bytes(), "auto")
    assert prepared_request.fixed_atoms == request.fixed_atoms
    assert prepared_request.generated_atoms == request.generated_atoms
    for single in independent_gap_requests(prepared_request):
        validate_context_source(single, source.read_bytes())
        validate_context_capability(single, "protpardelle")
        assert len(source_context_position_features(single)[0]) <= 512
