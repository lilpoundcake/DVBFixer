"""Fail-closed local preflight primitives for external diffusion adapters."""

from __future__ import annotations

import hashlib
import os
import shutil
import stat
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION


class AdapterPreflightCode(StrEnum):
    """Stable machine-readable adapter preflight failure categories."""

    MISSING_RUNNER = "missing-runner"
    MISSING_IMAGE = "missing-image"
    MISSING_CHECKPOINT = "missing-checkpoint"
    CHECKPOINT_DIGEST_MISMATCH = "checkpoint-digest-mismatch"
    MISSING_CUDA = "missing-cuda"
    INCOMPATIBLE_PROTOCOL = "incompatible-protocol"


@dataclass(frozen=True, slots=True)
class LocalCheckpoint:
    """A local-only checkpoint requirement with an immutable expected digest."""

    path: Path
    sha256: str

    def __post_init__(self) -> None:
        if len(self.sha256) != 64:
            raise ValueError("checkpoint SHA-256 must contain 64 hexadecimal characters")
        try:
            int(self.sha256, 16)
        except ValueError as exc:
            raise ValueError("checkpoint SHA-256 must be hexadecimal") from exc


@dataclass(frozen=True, slots=True)
class AdapterAvailability:
    """Runtime facts supplied by an engine-specific launcher or container host."""

    available_images: frozenset[str] = frozenset()
    cuda_devices: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class AdapterPreflightSpec:
    """Requirements that must be satisfied before launching an adapter."""

    runner: str
    protocol_version: int
    checkpoint: LocalCheckpoint | None = None
    image: str = ""
    requires_cuda: bool = False

    def __post_init__(self) -> None:
        if not self.runner:
            raise ValueError("adapter runner must not be empty")
        if self.protocol_version <= 0:
            raise ValueError("adapter protocol version must be positive")


@dataclass(frozen=True, slots=True)
class AdapterPreflightIssue:
    code: AdapterPreflightCode
    message: str


@dataclass(frozen=True, slots=True)
class AdapterPreflightReport:
    issues: tuple[AdapterPreflightIssue, ...]
    runner_path: str = ""
    checkpoint_sha256: str = ""

    @property
    def passed(self) -> bool:
        return not self.issues


def assess_adapter_preflight(
    spec: AdapterPreflightSpec,
    *,
    availability: AdapterAvailability | None = None,
) -> AdapterPreflightReport:
    """Check local immutable adapter requirements without downloading artifacts."""
    runtime = availability or AdapterAvailability()
    issues: list[AdapterPreflightIssue] = []
    runner_path = _resolve_runner(spec.runner)
    if not runner_path:
        issues.append(
            AdapterPreflightIssue(
                AdapterPreflightCode.MISSING_RUNNER,
                "external diffusion runner executable is unavailable",
            )
        )
    if spec.protocol_version != DIFFUSION_RUNNER_PROTOCOL_VERSION:
        issues.append(
            AdapterPreflightIssue(
                AdapterPreflightCode.INCOMPATIBLE_PROTOCOL,
                "adapter protocol version does not match the DVBFixer runner protocol",
            )
        )
    if spec.image and spec.image not in runtime.available_images:
        issues.append(
            AdapterPreflightIssue(
                AdapterPreflightCode.MISSING_IMAGE,
                "required immutable diffusion image is unavailable",
            )
        )
    checkpoint_sha256 = ""
    if spec.checkpoint is not None:
        checkpoint_sha256, checkpoint_issue = _checkpoint_digest(spec.checkpoint)
        if checkpoint_issue is not None:
            issues.append(checkpoint_issue)
    if spec.requires_cuda and not runtime.cuda_devices:
        issues.append(
            AdapterPreflightIssue(
                AdapterPreflightCode.MISSING_CUDA,
                "adapter requires a CUDA device but none was reported",
            )
        )
    return AdapterPreflightReport(
        issues=tuple(issues),
        runner_path=runner_path,
        checkpoint_sha256=checkpoint_sha256,
    )


def _resolve_runner(runner: str) -> str:
    if os.sep in runner or (os.altsep is not None and os.altsep in runner):
        path = Path(runner).expanduser()
        try:
            resolved = path.resolve(strict=True)
        except FileNotFoundError:
            return ""
        if not resolved.is_file() or not os.access(resolved, os.X_OK):
            return ""
        return str(resolved)
    return shutil.which(runner) or ""


def _checkpoint_digest(
    checkpoint: LocalCheckpoint,
) -> tuple[str, AdapterPreflightIssue | None]:
    path = checkpoint.path.expanduser()
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except FileNotFoundError:
        return "", AdapterPreflightIssue(
            AdapterPreflightCode.MISSING_CHECKPOINT,
            "required local checkpoint is unavailable",
        )
    except OSError:
        return "", AdapterPreflightIssue(
            AdapterPreflightCode.MISSING_CHECKPOINT,
            "required local checkpoint could not be opened safely",
        )
    try:
        file_stat = os.fstat(descriptor)
        if not stat.S_ISREG(file_stat.st_mode):
            return "", AdapterPreflightIssue(
                AdapterPreflightCode.MISSING_CHECKPOINT,
                "required local checkpoint is not a regular file",
            )
        digest = hashlib.sha256()
        while chunk := os.read(descriptor, 1024 * 1024):
            digest.update(chunk)
    finally:
        os.close(descriptor)
    actual = digest.hexdigest()
    if actual != checkpoint.sha256.lower():
        return actual, AdapterPreflightIssue(
            AdapterPreflightCode.CHECKPOINT_DIGEST_MISMATCH,
            "local checkpoint SHA-256 does not match the pinned digest",
        )
    return actual, None
