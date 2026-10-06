"""Tests for fail-closed local diffusion adapter preflight."""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import pytest

from dvbfixer.model.diffusion.preflight import (
    DIFFUSION_PROFILES,
    AdapterAvailability,
    AdapterPreflightCode,
    AdapterPreflightSpec,
    LocalCheckpoint,
    RunnerPreflightFacts,
    RunnerPreflightReport,
    assess_adapter_preflight,
    invoke_runner_preflight,
)
from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION


def _codes(spec: AdapterPreflightSpec, availability: AdapterAvailability) -> set[str]:
    return {
        issue.code.value
        for issue in assess_adapter_preflight(
            spec,
            availability=availability,
        ).issues
    }


def test_preflight_accepts_local_pinned_artifacts_without_downloading(
    tmp_path: Path,
) -> None:
    checkpoint_bytes = b"pinned-checkpoint"
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(checkpoint_bytes)
    expected = hashlib.sha256(checkpoint_bytes).hexdigest()
    spec = AdapterPreflightSpec(
        runner=sys.executable,
        protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
        checkpoint=LocalCheckpoint(checkpoint, expected),
        image="engine@sha256:immutable",
        requires_cuda=True,
    )

    report = assess_adapter_preflight(
        spec,
        availability=AdapterAvailability(
            available_images=frozenset({"engine@sha256:immutable"}),
            cuda_devices=("cuda:0",),
        ),
    )

    assert report.passed
    assert report.runner_path
    assert report.checkpoint_sha256 == expected


def test_preflight_reports_each_adapter_specific_blocker(tmp_path: Path) -> None:
    missing_checkpoint = tmp_path / "missing.pt"
    spec = AdapterPreflightSpec(
        runner=str(tmp_path / "missing-runner"),
        protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION + 1,
        checkpoint=LocalCheckpoint(missing_checkpoint, "0" * 64),
        image="missing@sha256:image",
        requires_cuda=True,
    )

    assert _codes(spec, AdapterAvailability()) == {
        AdapterPreflightCode.MISSING_RUNNER.value,
        AdapterPreflightCode.INCOMPATIBLE_PROTOCOL.value,
        AdapterPreflightCode.MISSING_IMAGE.value,
        AdapterPreflightCode.MISSING_CHECKPOINT.value,
        AdapterPreflightCode.MISSING_CUDA.value,
    }


def test_preflight_rejects_checkpoint_digest_mismatch_and_symlink(
    tmp_path: Path,
) -> None:
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"unexpected")
    mismatch = AdapterPreflightSpec(
        runner=sys.executable,
        protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
        checkpoint=LocalCheckpoint(checkpoint, "0" * 64),
    )
    assert _codes(mismatch, AdapterAvailability()) == {
        AdapterPreflightCode.CHECKPOINT_DIGEST_MISMATCH.value,
    }

    link = tmp_path / "checkpoint-link.pt"
    link.symlink_to(checkpoint)
    symlinked = AdapterPreflightSpec(
        runner=sys.executable,
        protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
        checkpoint=LocalCheckpoint(link, hashlib.sha256(b"unexpected").hexdigest()),
    )
    assert _codes(symlinked, AdapterAvailability()) == {
        AdapterPreflightCode.MISSING_CHECKPOINT.value,
    }


def test_preflight_metadata_validation_is_fail_closed(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="64 hexadecimal"):
        LocalCheckpoint(tmp_path / "checkpoint.pt", "short")
    with pytest.raises(ValueError, match="hexadecimal"):
        LocalCheckpoint(tmp_path / "checkpoint.pt", "z" * 64)
    with pytest.raises(ValueError, match="runner"):
        AdapterPreflightSpec("", DIFFUSION_RUNNER_PROTOCOL_VERSION)


def _preflight_runner(tmp_path: Path, body: str) -> Path:
    runner = tmp_path / "runner.py"
    runner.write_text(f"#!{sys.executable}\n{body}\n", encoding="utf-8")
    runner.chmod(0o755)
    return runner


def _successful_response(profile_name: str = "protenix-v1-cuda") -> RunnerPreflightReport:
    profile = DIFFUSION_PROFILES[profile_name]
    is_cuda = profile_name == "protenix-v1-cuda"
    return RunnerPreflightReport(
        profile=profile_name,
        runner_protocol_version=profile.runner_protocol_version,
        contract_schema_version=profile.contract_schema_version,
        facts=RunnerPreflightFacts(
            platform="linux" if is_cuda else "darwin",
            architecture="x86_64" if is_cuda else "arm64",
            python_version="3.13" if is_cuda else "3.12",
            engine_revision=profile.source_revision,
            patch_identity=profile.patch_identity,
            environment_identity=(
                "linux-amd64;python=3.13;torch=2.13.0;cuda=12.9"
                if is_cuda
                else "macos-arm64;python=3.12;torch=2.6.0"
            ),
            checkpoint_sha256=profile.expected_checkpoint_sha256,
            framework_version="2.13.0" if is_cuda else "2.6.0",
            accelerator_available=True,
            effective_device="cuda:0" if is_cuda else "mps",
            fallback_disabled=True,
            sampling_platform=profile.sampling_platform,
            refinement_platform=profile.refinement_platform,
        ),
    )


def test_runner_preflight_strict_round_trip_and_success(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    response = _successful_response()
    assert RunnerPreflightReport.from_json(response.to_json()) == response
    runner = _preflight_runner(
        tmp_path,
        "import os\n"
        "assert 'PRIVATE_CREDENTIAL' not in os.environ\n"
        f"print({response.to_json()!r}, end='')",
    )
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"not read by mock")
    monkeypatch.setenv("PRIVATE_CREDENTIAL", "must-not-reach-runner")

    observed = invoke_runner_preflight(
        profile="protenix-v1-cuda",
        runner=str(runner),
        checkpoint=checkpoint,
    )

    assert observed.passed


@pytest.mark.parametrize(
    ("body", "timeout", "code"),
    [
        ("print('not-json')", 1.0, "malformed-response"),
        ("print('x' * 70000)", 1.0, "oversized-response"),
        ("import time; time.sleep(2)", 0.01, "runner-timeout"),
        ("raise SystemExit(7)", 1.0, "runner-failed"),
    ],
)
def test_runner_preflight_transport_failures_are_stable(
    tmp_path: Path,
    body: str,
    timeout: float,
    code: str,
) -> None:
    runner = _preflight_runner(tmp_path, body)
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"checkpoint")

    report = invoke_runner_preflight(
        profile="protenix-v1-cuda",
        runner=str(runner),
        checkpoint=checkpoint,
        timeout_seconds=timeout,
    )

    assert [issue.code.value for issue in report.issues] == [code]
    assert "not-json" not in report.to_json()


def test_runner_preflight_rejects_protocol_and_digest_mismatch(tmp_path: Path) -> None:
    response = _successful_response()
    bad = RunnerPreflightReport(
        profile=response.profile,
        runner_protocol_version=response.runner_protocol_version + 1,
        contract_schema_version=response.contract_schema_version,
        facts=RunnerPreflightFacts(
            **{
                **response.facts.to_dict(),
                "checkpoint_sha256": "0" * 64,
            }
        ),
    )
    runner = _preflight_runner(tmp_path, f"print({bad.to_json()!r}, end='')")
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"checkpoint")

    report = invoke_runner_preflight(
        profile="protenix-v1-cuda",
        runner=str(runner),
        checkpoint=checkpoint,
    )

    assert {issue.code.value for issue in report.issues} == {
        "incompatible-protocol",
        "checkpoint-digest-mismatch",
    }


def test_runner_preflight_parser_rejects_unknown_fields_and_unbounded_issues() -> None:
    raw = json.loads(_successful_response().to_json())
    raw["private_path"] = "/secret/checkpoint"
    with pytest.raises(ValueError, match="fields"):
        RunnerPreflightReport.from_json(json.dumps(raw))

    raw.pop("private_path")
    raw["issues"] = [{"code": "missing-cuda", "message": "x"}] * 33
    raw["passed"] = False
    with pytest.raises(ValueError, match="issues"):
        RunnerPreflightReport.from_json(json.dumps(raw))


@pytest.mark.parametrize("timeout", [0.0, -1.0, float("nan"), float("inf")])
def test_runner_preflight_rejects_invalid_timeout(tmp_path: Path, timeout: float) -> None:
    with pytest.raises(ValueError, match="positive finite"):
        invoke_runner_preflight(
            profile="protenix-v1-cuda",
            runner=str(tmp_path / "runner"),
            checkpoint=tmp_path / "checkpoint.pt",
            timeout_seconds=timeout,
        )


def test_runner_preflight_timeout_cleans_up_inherited_output_pipes(tmp_path: Path) -> None:
    response = _successful_response()
    runner = _preflight_runner(
        tmp_path,
        "import subprocess, sys\n"
        "subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(10)'])\n"
        f"print({response.to_json()!r}, end='')",
    )
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"checkpoint")
    started = time.monotonic()

    report = invoke_runner_preflight(
        profile="protenix-v1-cuda",
        runner=str(runner),
        checkpoint=checkpoint,
        timeout_seconds=0.2,
    )

    assert time.monotonic() - started < 2.0
    assert report.passed


def test_runner_preflight_allows_launcher_owned_checkpoint(tmp_path: Path) -> None:
    response = _successful_response()
    runner = _preflight_runner(
        tmp_path,
        "import sys\n"
        "assert '--checkpoint' not in sys.argv\n"
        f"print({response.to_json()!r}, end='')",
    )

    report = invoke_runner_preflight(
        profile="protenix-v1-cuda",
        runner=str(runner),
    )

    assert report.passed
