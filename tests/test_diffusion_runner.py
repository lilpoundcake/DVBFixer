"""Tests for the contained external diffusion runner boundary."""

from __future__ import annotations

import hashlib
import json
import sys
import textwrap
from pathlib import Path

import pytest

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.runner import (
    DIFFUSION_RUNNER_PROTOCOL_VERSION,
    DiffusionRunnerError,
    RunnerLimits,
    prepare_diffusion_workspace,
    run_diffusion_runner,
    run_prepared_diffusion_runner,
)


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _request(input_bytes: bytes) -> DiffusionRequest:
    residues = tuple(ResidueIdentity("D", str(number)) for number in range(10, 17))
    generated_residues = residues[1:6]
    generated_atoms = tuple(
        AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, "CA")
        for residue in generated_residues
    )
    return DiffusionRequest(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("input/normalized.pdb", _digest(input_bytes)),
        target_sequences=(TargetSequence("D", "FSGSKSG"),),
        sequence_placements=(
            SequencePlacement(
                chain="D",
                target_length=7,
                observed_target_indices=(0, 6),
                observed_residues=(residues[0], residues[-1]),
            ),
        ),
        gaps=(
            GapRegion(
                chain="D",
                target_interval=TargetInterval(1, 6),
                left_anchor=residues[0],
                right_anchor=residues[-1],
                generated_residues=generated_residues,
                movable_junction_residues=(residues[0], *generated_residues, residues[-1]),
            ),
        ),
        fixed_atoms=(
            AtomIdentity("D", "10", "", "CA"),
            AtomIdentity("D", "16", "", "CA"),
        ),
        generated_atoms=generated_atoms,
        retained_explicit_links=(),
        candidate_count=2,
        seeds=(7, 11),
    )


def _write_input(source_root: Path, data: bytes) -> None:
    path = source_root / "input" / "normalized.pdb"
    path.parent.mkdir(parents=True)
    path.write_bytes(data)


def _write_runner(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "runner.py"
    path.write_text(textwrap.dedent(body), encoding="utf-8")
    return path


def _command(script: Path, *arguments: str) -> tuple[str, ...]:
    return (sys.executable, str(script), *arguments)


def _valid_runner(tmp_path: Path) -> Path:
    return _write_runner(
        tmp_path,
        f"""
        import hashlib
        import json
        from pathlib import Path

        request = json.loads(Path("request.json").read_text())
        input_path = Path(request["normalized_pdb"]["path"])
        candidate_path = Path("candidates/candidate-0001.pdb")
        candidate_path.parent.mkdir(parents=True)
        candidate_path.write_bytes(input_path.read_bytes())
        digest = hashlib.sha256(candidate_path.read_bytes()).hexdigest()
        result = {{
            "schema_version": {DIFFUSION_SCHEMA_VERSION},
            "status": "success",
            "candidates": [{{
                "candidate_id": "candidate-0001",
                "seed": request["seeds"][0],
                "coordinate_artifact": {{
                    "path": str(candidate_path),
                    "sha256": digest,
                }},
                "generated_atoms": request["generated_atoms"],
                "generated_residues": request["gaps"][0]["generated_residues"],
                "raw_backend_score": 1.0,
                "score_provenance": "fake-runner:test-score-v1",
                "warnings": [],
            }}],
            "runner_diagnostics": {{
                "exit_code": None,
                "timed_out": False,
                "stdout": "runner supplied diagnostics are ignored",
                "stderr": "",
            }},
            "backend_provenance": {{
                "backend": "fake",
                "runner_protocol_version": {DIFFUSION_RUNNER_PROTOCOL_VERSION},
                "engine_repository": "https://example.invalid/fake",
                "engine_revision": "test-revision",
                "checkpoint_sha256": "",
                "environment_hash": "",
            }},
            "message": "",
        }}
        Path("result.json").write_text(json.dumps(result, sort_keys=True))
        print("fake runner complete")
        """,
    )


def test_runner_validates_manifests_digests_and_uses_actual_diagnostics(
    tmp_path: Path,
) -> None:
    source_root = tmp_path / "source"
    source_root.mkdir()
    input_bytes = b"HEADER    FAKE INPUT\nEND\n"
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)
    workspace = tmp_path / "workspace"

    result = run_diffusion_runner(
        request,
        _command(_valid_runner(tmp_path)),
        source_root=source_root,
        workspace=workspace,
    )

    assert result.candidates[0].coordinate_artifact.sha256 == _digest(input_bytes)
    assert (workspace / result.candidates[0].coordinate_artifact.path).read_bytes() == input_bytes
    assert result.runner_diagnostics.exit_code == 0
    assert result.runner_diagnostics.timed_out is False
    assert result.runner_diagnostics.stdout == "fake runner complete\n"
    assert result.runner_diagnostics.stderr == ""
    assert (workspace / "request.json").is_file()
    assert (workspace / "result.json").is_file()


def test_builtin_fake_runner_is_deterministic_and_external(tmp_path: Path) -> None:
    input_bytes = b"HEADER    FAKE INPUT\nEND\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    first = prepare_diffusion_workspace(
        request,
        source_root=source_root,
        workspace=tmp_path / "first-workspace",
    )
    first_result = run_prepared_diffusion_runner(
        first,
        (sys.executable, "-m", "dvbfixer.model.diffusion.fake_runner"),
    )
    second = prepare_diffusion_workspace(
        request,
        source_root=source_root,
        workspace=tmp_path / "second-workspace",
    )
    second_result = run_prepared_diffusion_runner(
        second,
        (sys.executable, "-m", "dvbfixer.model.diffusion.fake_runner"),
    )

    assert first_result == second_result
    assert len(first_result.candidates) == request.candidate_count
    assert [candidate.raw_backend_score for candidate in first_result.candidates] == [7.0, 11.0]
    assert all(
        candidate.coordinate_artifact.sha256 == request.normalized_pdb.sha256
        for candidate in first_result.candidates
    )
    assert first_result.backend_provenance.backend == "dvbfixer-fake"


def test_runner_rejects_modified_prepared_request(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)
    prepared = prepare_diffusion_workspace(
        request,
        source_root=source_root,
        workspace=tmp_path / "workspace",
    )
    raw = json.loads(prepared.request_manifest.read_text(encoding="utf-8"))
    raw["candidate_count"] = 1
    raw["seeds"] = [7]
    prepared.request_manifest.write_text(json.dumps(raw), encoding="utf-8")

    with pytest.raises(DiffusionRunnerError, match="does not match"):
        run_prepared_diffusion_runner(
            prepared,
            (sys.executable, "-m", "dvbfixer.model.diffusion.fake_runner"),
        )


def test_runner_rejects_missing_or_malformed_result_manifest(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    missing = _write_runner(tmp_path, "print('no result')")
    with pytest.raises(DiffusionRunnerError, match="did not create") as missing_error:
        run_diffusion_runner(
            request,
            _command(missing),
            source_root=source_root,
            workspace=tmp_path / "missing-workspace",
        )
    assert missing_error.value.diagnostics is not None
    assert missing_error.value.diagnostics.exit_code == 0

    malformed = _write_runner(
        tmp_path,
        """
        from pathlib import Path
        Path("result.json").write_text("{not-json")
        """,
    )
    with pytest.raises(DiffusionRunnerError, match="invalid result.json"):
        run_diffusion_runner(
            request,
            _command(malformed),
            source_root=source_root,
            workspace=tmp_path / "malformed-workspace",
        )


def test_runner_rejects_digest_mismatch_and_candidate_symlink(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    bad_digest = _valid_runner(tmp_path)
    script_text = bad_digest.read_text(encoding="utf-8").replace(
        '"sha256": digest,',
        '"sha256": "0" * 64,',
    )
    bad_digest.write_text(script_text, encoding="utf-8")
    with pytest.raises(DiffusionRunnerError, match="SHA-256 mismatch"):
        run_diffusion_runner(
            request,
            _command(bad_digest),
            source_root=source_root,
            workspace=tmp_path / "digest-workspace",
        )

    symlink_runner = _write_runner(
        tmp_path,
        f"""
        import hashlib
        import json
        from pathlib import Path

        request = json.loads(Path("request.json").read_text())
        candidate = Path("candidate.pdb")
        candidate.symlink_to(request["normalized_pdb"]["path"])
        result = {{
            "schema_version": {DIFFUSION_SCHEMA_VERSION},
            "status": "success",
            "candidates": [{{
                "candidate_id": "candidate-0001",
                "seed": request["seeds"][0],
                "coordinate_artifact": {{
                    "path": str(candidate),
                    "sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
                }},
                "generated_atoms": request["generated_atoms"],
                "generated_residues": request["gaps"][0]["generated_residues"],
                "raw_backend_score": 1.0,
                "score_provenance": "fake",
                "warnings": [],
            }}],
            "runner_diagnostics": {{
                "exit_code": 0,
                "timed_out": False,
                "stdout": "",
                "stderr": "",
            }},
            "backend_provenance": {{
                "backend": "fake",
                "runner_protocol_version": {DIFFUSION_RUNNER_PROTOCOL_VERSION},
                "engine_repository": "https://example.invalid/fake",
                "engine_revision": "test-revision",
                "checkpoint_sha256": "",
                "environment_hash": "",
            }},
            "message": "",
        }}
        Path("result.json").write_text(json.dumps(result))
        """,
    )
    with pytest.raises(DiffusionRunnerError, match="contains a symlink"):
        run_diffusion_runner(
            request,
            _command(symlink_runner),
            source_root=source_root,
            workspace=tmp_path / "symlink-workspace",
        )


def test_runner_rejects_oversized_manifest_artifact_and_hard_link(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    oversized_manifest = _write_runner(
        tmp_path,
        """
        from pathlib import Path
        Path("result.json").write_text("x" * 1000)
        """,
    )
    with pytest.raises(DiffusionRunnerError, match="manifest exceeds the size limit"):
        run_diffusion_runner(
            request,
            _command(oversized_manifest),
            source_root=source_root,
            workspace=tmp_path / "manifest-workspace",
            limits=RunnerLimits(max_manifest_bytes=100),
        )

    oversized_artifact = _valid_runner(tmp_path)
    larger_candidate = oversized_artifact.read_text(encoding="utf-8").replace(
        "candidate_path.write_bytes(input_path.read_bytes())",
        'candidate_path.write_bytes(b"candidate-bytes")',
    )
    oversized_artifact.write_text(larger_candidate, encoding="utf-8")
    with pytest.raises(DiffusionRunnerError, match="artifact exceeds the size limit"):
        run_diffusion_runner(
            request,
            _command(oversized_artifact),
            source_root=source_root,
            workspace=tmp_path / "artifact-workspace",
            limits=RunnerLimits(max_artifact_bytes=len(input_bytes)),
        )

    hard_link = _write_runner(
        tmp_path,
        f"""
        import hashlib
        import json
        import os
        from pathlib import Path

        request = json.loads(Path("request.json").read_text())
        candidate = Path("candidate.pdb")
        os.link(request["normalized_pdb"]["path"], candidate)
        result = {{
            "schema_version": {DIFFUSION_SCHEMA_VERSION},
            "status": "success",
            "candidates": [{{
                "candidate_id": "candidate-0001",
                "seed": request["seeds"][0],
                "coordinate_artifact": {{
                    "path": str(candidate),
                    "sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
                }},
                "generated_atoms": request["generated_atoms"],
                "generated_residues": request["gaps"][0]["generated_residues"],
                "raw_backend_score": 1.0,
                "score_provenance": "fake",
                "warnings": [],
            }}],
            "runner_diagnostics": {{
                "exit_code": 0,
                "timed_out": False,
                "stdout": "",
                "stderr": "",
            }},
            "backend_provenance": {{
                "backend": "fake",
                "runner_protocol_version": {DIFFUSION_RUNNER_PROTOCOL_VERSION},
                "engine_repository": "https://example.invalid/fake",
                "engine_revision": "test-revision",
                "checkpoint_sha256": "",
                "environment_hash": "",
            }},
            "message": "",
        }}
        Path("result.json").write_text(json.dumps(result))
        """,
    )
    with pytest.raises(DiffusionRunnerError, match="must not be hard-linked"):
        run_diffusion_runner(
            request,
            _command(hard_link),
            source_root=source_root,
            workspace=tmp_path / "hard-link-workspace",
        )


def test_runner_environment_is_minimal_and_overrides_are_allowlisted(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)
    environment_probe = _write_runner(
        tmp_path,
        f"""
        import hashlib
        import json
        import os
        from pathlib import Path

        request = json.loads(Path("request.json").read_text())
        if "AWS_ACCESS_KEY_ID" in os.environ or "SSH_AUTH_SOCK" in os.environ:
            raise SystemExit(7)
        if os.environ.get("LANG") != "test-locale":
            raise SystemExit(8)
        candidate = Path("candidate.pdb")
        candidate.write_bytes(Path(request["normalized_pdb"]["path"]).read_bytes())
        result = {{
            "schema_version": {DIFFUSION_SCHEMA_VERSION},
            "status": "success",
            "candidates": [{{
                "candidate_id": "candidate-0001",
                "seed": request["seeds"][0],
                "coordinate_artifact": {{
                    "path": str(candidate),
                    "sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
                }},
                "generated_atoms": request["generated_atoms"],
                "generated_residues": request["gaps"][0]["generated_residues"],
                "raw_backend_score": 1.0,
                "score_provenance": "fake",
                "warnings": [],
            }}],
            "runner_diagnostics": {{
                "exit_code": 0,
                "timed_out": False,
                "stdout": "",
                "stderr": "",
            }},
            "backend_provenance": {{
                "backend": "fake",
                "runner_protocol_version": {DIFFUSION_RUNNER_PROTOCOL_VERSION},
                "engine_repository": "https://example.invalid/fake",
                "engine_revision": "test-revision",
                "checkpoint_sha256": "",
                "environment_hash": "",
            }},
            "message": "",
        }}
        Path("result.json").write_text(json.dumps(result))
        """,
    )

    result = run_diffusion_runner(
        request,
        _command(environment_probe),
        source_root=source_root,
        workspace=tmp_path / "workspace",
        environment={"LANG": "test-locale"},
    )

    assert result.status.value == "success"


def test_runner_kills_pipe_owning_descendants_after_parent_exit(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)
    spawning_runner = _write_runner(
        tmp_path,
        """
        import subprocess
        import sys

        subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)"],
        )
        """,
    )

    with pytest.raises(DiffusionRunnerError, match="did not create"):
        run_diffusion_runner(
            request,
            _command(spawning_runner),
            source_root=source_root,
            workspace=tmp_path / "workspace",
            limits=RunnerLimits(timeout_seconds=1.0),
        )


def test_runner_times_out_and_bounds_stdout_and_stderr(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    noisy = _write_runner(
        tmp_path,
        """
        import sys
        import time

        sys.stdout.write("O" * 10000)
        sys.stdout.flush()
        sys.stderr.write("E" * 10000)
        sys.stderr.flush()
        time.sleep(10)
        """,
    )
    with pytest.raises(DiffusionRunnerError, match="timed out") as error:
        run_diffusion_runner(
            request,
            _command(noisy),
            source_root=source_root,
            workspace=tmp_path / "timeout-workspace",
            limits=RunnerLimits(timeout_seconds=0.1, max_output_bytes=128),
        )

    diagnostics = error.value.diagnostics
    assert diagnostics is not None
    assert diagnostics.timed_out is True
    assert len(diagnostics.stdout.encode()) <= 128
    assert len(diagnostics.stderr.encode()) <= 128
    assert "truncated" in diagnostics.stdout
    assert "truncated" in diagnostics.stderr


def test_runner_rejects_incompatible_protocol_and_input_mutation(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    incompatible = _valid_runner(tmp_path)
    script_text = incompatible.read_text(encoding="utf-8").replace(
        f'"runner_protocol_version": {DIFFUSION_RUNNER_PROTOCOL_VERSION},',
        f'"runner_protocol_version": {DIFFUSION_RUNNER_PROTOCOL_VERSION + 1},',
    )
    incompatible.write_text(script_text, encoding="utf-8")
    with pytest.raises(DiffusionRunnerError, match="incompatible protocol version"):
        run_diffusion_runner(
            request,
            _command(incompatible),
            source_root=source_root,
            workspace=tmp_path / "protocol-workspace",
        )

    modifying = _valid_runner(tmp_path)
    script_text = modifying.read_text(encoding="utf-8").replace(
        'candidate_path.parent.mkdir(parents=True)',
        'input_path.chmod(0o600)\n'
        'input_path.write_text("modified")\n'
        'candidate_path.parent.mkdir(parents=True)',
    )
    modifying.write_text(script_text, encoding="utf-8")
    with pytest.raises(DiffusionRunnerError, match="modified its staged input"):
        run_diffusion_runner(
            request,
            _command(modifying),
            source_root=source_root,
            workspace=tmp_path / "mutation-workspace",
        )


def test_runner_accepts_symlinked_system_interpreter(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)

    result = run_diffusion_runner(
        _request(input_bytes),
        (sys.executable, "-m", "dvbfixer.model.diffusion.fake_runner"),
        source_root=source_root,
        workspace=tmp_path / "workspace",
    )

    assert result.backend_provenance.backend == "dvbfixer-fake"


def test_runner_rejects_missing_executable_existing_workspace_and_credentials(
    tmp_path: Path,
) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    with pytest.raises(DiffusionRunnerError, match="not found"):
        run_diffusion_runner(
            request,
            (str(tmp_path / "missing-runner"),),
            source_root=source_root,
            workspace=tmp_path / "missing-executable-workspace",
        )

    workspace = tmp_path / "existing-workspace"
    workspace.mkdir()
    with pytest.raises(DiffusionRunnerError, match="must not already exist"):
        run_diffusion_runner(
            request,
            _command(_valid_runner(tmp_path)),
            source_root=source_root,
            workspace=workspace,
        )

    with pytest.raises(DiffusionRunnerError, match="not allowlisted"):
        run_diffusion_runner(
            request,
            _command(_valid_runner(tmp_path)),
            source_root=source_root,
            workspace=tmp_path / "credential-workspace",
            environment={"MODEL_API_TOKEN": "do-not-forward"},
        )


def test_runner_rejects_unrequested_candidate_seed(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    unrequested = _valid_runner(tmp_path)
    script_text = unrequested.read_text(encoding="utf-8").replace(
        '"seed": request["seeds"][0],',
        '"seed": 999,',
    )
    unrequested.write_text(script_text, encoding="utf-8")
    with pytest.raises(DiffusionRunnerError, match="unrequested candidate seed"):
        run_diffusion_runner(
            request,
            _command(unrequested),
            source_root=source_root,
            workspace=tmp_path / "unrequested-seed-workspace",
        )


def test_runner_rejects_duplicate_candidate_seeds(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)
    duplicate_seeds = _write_runner(
        tmp_path,
        f"""
        import hashlib
        import json
        from pathlib import Path

        request = json.loads(Path("request.json").read_text())
        candidates = []
        for index in range(2):
            candidate = Path(f"candidate-{{index}}.pdb")
            candidate.write_bytes(Path(request["normalized_pdb"]["path"]).read_bytes())
            candidates.append({{
                "candidate_id": f"candidate-{{index}}",
                "seed": request["seeds"][0],
                "coordinate_artifact": {{
                    "path": str(candidate),
                    "sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
                }},
                "generated_atoms": request["generated_atoms"],
                "generated_residues": request["gaps"][0]["generated_residues"],
                "raw_backend_score": float(index),
                "score_provenance": "fake",
                "warnings": [],
            }})
        result = {{
            "schema_version": {DIFFUSION_SCHEMA_VERSION},
            "status": "success",
            "candidates": candidates,
            "runner_diagnostics": {{
                "exit_code": 0,
                "timed_out": False,
                "stdout": "",
                "stderr": "",
            }},
            "backend_provenance": {{
                "backend": "fake",
                "runner_protocol_version": {DIFFUSION_RUNNER_PROTOCOL_VERSION},
                "engine_repository": "https://example.invalid/fake",
                "engine_revision": "test-revision",
            }},
            "message": "",
        }}
        Path("result.json").write_text(json.dumps(result))
        """,
    )

    with pytest.raises(DiffusionRunnerError, match="duplicate candidate seeds"):
        run_diffusion_runner(
            request,
            _command(duplicate_seeds),
            source_root=source_root,
            workspace=tmp_path / "duplicate-seed-workspace",
        )


def test_runner_rejects_result_symlink_and_out_of_contract_masks(tmp_path: Path) -> None:
    input_bytes = b"END\n"
    source_root = tmp_path / "source"
    source_root.mkdir()
    _write_input(source_root, input_bytes)
    request = _request(input_bytes)

    outside_result = tmp_path / "outside-result.json"
    outside_result.write_text("{}", encoding="utf-8")
    result_symlink = _write_runner(
        tmp_path,
        f"""
        from pathlib import Path
        Path("result.json").symlink_to({str(outside_result)!r})
        """,
    )
    with pytest.raises(DiffusionRunnerError, match="regular result.json"):
        run_diffusion_runner(
            request,
            _command(result_symlink),
            source_root=source_root,
            workspace=tmp_path / "result-symlink-workspace",
        )

    changed_mask = _valid_runner(tmp_path)
    script_text = changed_mask.read_text(encoding="utf-8").replace(
        '"generated_atoms": request["generated_atoms"],',
        '"generated_atoms": request["generated_atoms"][:-1],',
    )
    changed_mask.write_text(script_text, encoding="utf-8")
    with pytest.raises(DiffusionRunnerError, match="changed the generated atom mask"):
        run_diffusion_runner(
            request,
            _command(changed_mask),
            source_root=source_root,
            workspace=tmp_path / "mask-workspace",
        )
