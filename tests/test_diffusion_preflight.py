"""Tests for fail-closed local diffusion adapter preflight."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

from dvbfixer.model.diffusion.preflight import (
    AdapterAvailability,
    AdapterPreflightCode,
    AdapterPreflightSpec,
    LocalCheckpoint,
    assess_adapter_preflight,
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
