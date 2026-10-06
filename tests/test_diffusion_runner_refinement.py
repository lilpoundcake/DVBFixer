from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import dvbfixer.model.diffusion.runner_refinement as runner_refinement
from dvbfixer.model.diffusion.boundary_refinement import BoundaryRefinementError
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    GapRegion,
    ResidueIdentity,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResult,
    SamplerTrace,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.runner_refinement import (
    refine_runner_result,
    rewrite_generated_coordinates,
)
from dvbfixer.model.diffusion.trace import build_sampler_trace, write_sampler_trace


def _atom_line(serial: int, residue_number: int, atom_name: str, x: float) -> str:
    return (
        f"ATOM  {serial:5d}  {atom_name:<3} ALA A{residue_number:4d}    "
        f"{x:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00  0.00           C\n"
    )


def test_rewrite_generated_coordinates_preserves_non_coordinate_text() -> None:
    fixed = _atom_line(1, 1, "CA", 1.0)
    generated = _atom_line(2, 2, "CA", 4.0)
    text = "HEADER preserved\n" + fixed + generated + "END\n"

    rewritten = rewrite_generated_coordinates(
        text,
        {AtomIdentity("A", "2", "", "CA"): np.asarray([7.0, 8.0, 9.0])},
    )

    original_lines = text.splitlines(keepends=True)
    rewritten_lines = rewritten.splitlines(keepends=True)
    assert rewritten_lines[:2] == original_lines[:2]
    assert rewritten_lines[2][:30] == original_lines[2][:30]
    assert rewritten_lines[2][30:54] == "   7.000   8.000   9.000"
    assert rewritten_lines[2][54:] == original_lines[2][54:]
    assert rewritten_lines[3] == original_lines[3]


def test_rewrite_generated_coordinates_rejects_missing_atoms() -> None:
    coordinates = {AtomIdentity("A", "2", "", "CA"): np.asarray([7.0, 8.0, 9.0])}
    with pytest.raises(ValueError, match="omits 1 generated atoms"):
        rewrite_generated_coordinates("END\n", coordinates)


def test_private_coordinate_read_rejects_digest_links_and_oversize(tmp_path: Path) -> None:
    coordinate = tmp_path / "candidate.pdb"
    data = b"END\n"
    coordinate.write_bytes(data)
    valid = ArtifactReference("candidate.pdb", hashlib.sha256(data).hexdigest())

    with pytest.raises(ValueError, match="digest mismatch"):
        runner_refinement._read_private_artifact(
            tmp_path,
            ArtifactReference("candidate.pdb", "0" * 64),
            label="coordinate artifact",
            max_bytes=100,
            encoding="ascii",
        )

    linked = tmp_path / "linked.pdb"
    linked.hardlink_to(coordinate)
    with pytest.raises(ValueError, match="private regular file"):
        runner_refinement._read_private_artifact(
            tmp_path,
            valid,
            label="coordinate artifact",
            max_bytes=100,
            encoding="ascii",
        )
    linked.unlink()

    symlink = tmp_path / "symlink.pdb"
    symlink.symlink_to(coordinate.name)
    with pytest.raises(ValueError, match="contains a symlink"):
        runner_refinement._read_private_artifact(
            tmp_path,
            ArtifactReference("symlink.pdb", valid.sha256),
            label="coordinate artifact",
            max_bytes=100,
            encoding="ascii",
        )

    with pytest.raises(ValueError, match="size limit"):
        runner_refinement._read_private_artifact(
            tmp_path,
            valid,
            label="coordinate artifact",
            max_bytes=3,
            encoding="ascii",
        )


def test_reference_refinement_preserves_and_updates_sampler_trace(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixed_residue = ResidueIdentity("A", "1")
    generated_residue = ResidueIdentity("A", "2")
    right_residue = ResidueIdentity("A", "3")
    fixed_atoms = (
        AtomIdentity("A", "1", "", "CA"),
        AtomIdentity("A", "3", "", "CA"),
    )
    generated_atom = AtomIdentity("A", "2", "", "CA")
    raw_text = (
        _atom_line(1, 1, "CA", 1.0)
        + _atom_line(2, 2, "CA", 4.0)
        + _atom_line(3, 3, "CA", 7.0)
        + "END\n"
    )
    raw_path = tmp_path / "raw.pdb"
    raw_path.write_text(raw_text, encoding="ascii")
    request = DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("raw.pdb", hashlib.sha256(raw_text.encode()).hexdigest()),
        target_sequences=(TargetSequence("A", "AAA"),),
        sequence_placements=(
            SequencePlacement("A", 3, (0, 2), (fixed_residue, right_residue)),
        ),
        gaps=(
            GapRegion(
                "A",
                TargetInterval(1, 2),
                fixed_residue,
                right_residue,
                (generated_residue,),
                (fixed_residue, generated_residue, right_residue),
            ),
        ),
        fixed_atoms=fixed_atoms,
        generated_atoms=(generated_atom,),
        retained_explicit_links=(),
        candidate_count=1,
        seeds=(7,),
    )
    request_path = tmp_path / "request.json"
    request_path.write_text(request.to_json(), encoding="utf-8")
    raw_trace = build_sampler_trace(
        profile="test",
        engine_repository="builtin://test",
        engine_revision="revision",
        patch_identity="patch",
        atom_order=(*fixed_atoms, generated_atom),
        fixed_atoms=fixed_atoms,
        device="cpu",
        fallback_disabled=True,
        denoising_update_count=0,
        final_fixed_coordinate_restoration=True,
    )
    trace_artifact = write_sampler_trace(tmp_path, tmp_path / "raw.trace.json", raw_trace)
    raw_result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=(
            RunnerCandidate(
                candidate_id="raw",
                seed=7,
                coordinate_artifact=ArtifactReference(
                    "raw.pdb",
                    hashlib.sha256(raw_text.encode()).hexdigest(),
                ),
                sampler_trace_artifact=trace_artifact,
                generated_atoms=(generated_atom,),
                generated_residues=(generated_residue,),
                raw_backend_score=None,
                score_provenance="test",
            ),
        ),
        runner_diagnostics=RunnerDiagnostics(0, False),
        backend_provenance=BackendProvenance(
            backend="test",
            runner_protocol_version=4,
            engine_repository="builtin://test",
            engine_revision="revision",
        ),
    )
    monkeypatch.setattr(
        runner_refinement,
        "refine_generated_region",
        lambda *_args, **_kwargs: SimpleNamespace(
            coordinates_angstrom={generated_atom: np.asarray([5.0, 6.0, 7.0])},
            platform="Reference",
        ),
    )

    refined = refine_runner_result(
        request_path,
        raw_result,
        tmp_path / "refined",
        platform_name="Reference",
    )
    refined_reference = refined.candidates[0].sampler_trace_artifact
    refined_trace = SamplerTrace.from_json(
        (tmp_path / refined_reference.path).read_text(encoding="utf-8")
    )

    assert refined_reference.sha256 != trace_artifact.sha256
    assert refined_trace.engine_revision == raw_trace.engine_revision
    assert refined_trace.atom_order_sha256 == raw_trace.atom_order_sha256
    assert refined_trace.sampler_evidence_complete == raw_trace.sampler_evidence_complete
    assert refined_trace.refinement_mode == "localized-openmm-boundary-refinement"
    assert refined_trace.localized_refinement_residues == request.gaps[0].movable_junction_residues
    assert refined.backend_provenance.known_nondeterministic_operations == ()

    def fail_refinement(*_args: object, **_kwargs: object) -> None:
        raise BoundaryRefinementError("injected scientific refinement failure")

    monkeypatch.setattr(
        runner_refinement,
        "refine_generated_region",
        fail_refinement,
    )
    rejected = refine_runner_result(
        request_path,
        raw_result,
        tmp_path / "failed-refinement",
        platform_name="Reference",
    )

    assert rejected.candidates[0].coordinate_artifact == raw_result.candidates[0].coordinate_artifact
    assert rejected.candidates[0].postprocessing_failures == (
        "localized-openmm-boundary-refinement-failed",
    )
    assert "raw sampler candidate retained" in rejected.message
    assert not (tmp_path / "failed-refinement").exists()
