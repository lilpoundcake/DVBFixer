"""Contained subprocess boundary for experimental diffusion runners."""

from __future__ import annotations

import hashlib
import math
import os
import shutil
import signal
import stat
import subprocess
import threading
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import BinaryIO

from dvbfixer.model.diffusion.contract import (
    MAX_TRACE_BYTES,
    ArtifactReference,
    DiffusionContractError,
    DiffusionRequest,
    DiffusionStatus,
    RunnerDiagnostics,
    RunnerResult,
    SamplerTrace,
)
from dvbfixer.model.diffusion.trace import validate_sampler_trace_context

DIFFUSION_RUNNER_PROTOCOL_VERSION = 4
REQUEST_MANIFEST = "request.json"
RESULT_MANIFEST = "result.json"
_TRUNCATION_MARKER = b"\n...[output truncated by DVBFixer]"
_PIPE_DRAIN_TIMEOUT_SECONDS = 1.0
_ENVIRONMENT_ALLOWLIST = (
    "CUDA_VISIBLE_DEVICES",
    "LANG",
    "LC_ALL",
    "LC_CTYPE",
    "LD_LIBRARY_PATH",
    "NVIDIA_DRIVER_CAPABILITIES",
    "NVIDIA_VISIBLE_DEVICES",
    "PATH",
    "PYTHONNOUSERSITE",
    "PYTHONPATH",
    "PYTORCH_ENABLE_MPS_FALLBACK",
    "PYTORCH_MPS_FAST_MATH",
    "PYTORCH_MPS_HIGH_WATERMARK_RATIO",
    "PYTORCH_MPS_LOW_WATERMARK_RATIO",
    "PYTORCH_MPS_PREFER_METAL",
    "SYSTEMROOT",
    "TZ",
)


class DiffusionRunnerError(RuntimeError):
    """Raised when runner preflight, execution, or artifact validation fails."""

    def __init__(
        self,
        message: str,
        *,
        diagnostics: RunnerDiagnostics | None = None,
    ) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


@dataclass(frozen=True, slots=True)
class RunnerLimits:
    """Resource bounds enforced around one external runner invocation."""

    timeout_seconds: float = 300.0
    max_output_bytes: int = 1_000_000
    max_manifest_bytes: int = 10_000_000
    max_artifact_bytes: int = 1_000_000_000
    max_workspace_bytes: int = 2_000_000_000

    def __post_init__(self) -> None:
        if not math.isfinite(self.timeout_seconds) or self.timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be a positive finite number")
        if self.max_output_bytes <= 0:
            raise ValueError("max_output_bytes must be positive")
        if self.max_manifest_bytes <= 0:
            raise ValueError("max_manifest_bytes must be positive")
        if self.max_artifact_bytes <= 0:
            raise ValueError("max_artifact_bytes must be positive")
        if self.max_workspace_bytes <= 0:
            raise ValueError("max_workspace_bytes must be positive")


@dataclass(frozen=True, slots=True)
class PreparedDiffusionRun:
    """Manifest paths created for an already isolated runner workspace."""

    request: DiffusionRequest
    workspace: Path
    request_manifest: Path


class _BoundedCapture:
    def __init__(self, limit: int) -> None:
        self._limit = limit
        self._data = bytearray()
        self._truncated = False

    def drain(self, stream: BinaryIO, stop: threading.Event) -> None:
        descriptor = stream.fileno()
        os.set_blocking(descriptor, False)
        while not stop.is_set():
            try:
                chunk = os.read(descriptor, 65536)
            except BlockingIOError:
                stop.wait(0.01)
                continue
            if not chunk:
                return
            remaining = self._limit - len(self._data)
            if remaining > 0:
                self._data.extend(chunk[:remaining])
            if len(chunk) > remaining:
                self._truncated = True

    def text(self) -> str:
        data = bytes(self._data)
        if self._truncated:
            if self._limit >= len(_TRUNCATION_MARKER):
                data = data[: self._limit - len(_TRUNCATION_MARKER)] + _TRUNCATION_MARKER
            else:
                data = _TRUNCATION_MARKER[: self._limit]
        return data.decode("utf-8", errors="replace")


def prepare_diffusion_workspace(
    request: DiffusionRequest,
    *,
    source_root: Path,
    workspace: Path,
) -> PreparedDiffusionRun:
    """Copy the normalized input and request manifest into a new workspace."""
    root = _require_directory(source_root, "source_root")
    run_root = _create_workspace(workspace)

    input_source = _validated_artifact_path(root, request.normalized_pdb)
    input_destination = run_root.joinpath(*Path(request.normalized_pdb.path).parts)
    input_destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(input_source, input_destination)
    input_destination.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)

    request_path = run_root / REQUEST_MANIFEST
    request_path.write_text(request.to_json(), encoding="utf-8")
    return PreparedDiffusionRun(
        request=request,
        workspace=run_root,
        request_manifest=request_path,
    )


def run_prepared_diffusion_runner(
    prepared: PreparedDiffusionRun,
    command: Sequence[str | os.PathLike[str]],
    *,
    limits: RunnerLimits | None = None,
    environment: Mapping[str, str] | None = None,
) -> RunnerResult:
    """Launch an external adapter for a workspace prepared by DVBFixer."""
    active_limits = limits or RunnerLimits()
    argv = _preflight_command(command)
    run_root = _require_directory(prepared.workspace, "runner workspace")
    if prepared.request_manifest != run_root / REQUEST_MANIFEST:
        raise DiffusionRunnerError("prepared request manifest is outside its runner workspace")
    if prepared.request_manifest.is_symlink() or not prepared.request_manifest.is_file():
        raise DiffusionRunnerError("prepared request manifest is not a regular request.json file")
    try:
        staged_request = DiffusionRequest.from_json(
            _read_bounded_text(prepared.request_manifest, active_limits.max_manifest_bytes)
        )
    except (DiffusionContractError, OSError, UnicodeError) as exc:
        raise DiffusionRunnerError(f"prepared request.json manifest is invalid: {exc}") from exc
    if staged_request != prepared.request:
        raise DiffusionRunnerError("prepared request.json manifest does not match the request")
    _validate_staged_input(
        run_root,
        prepared.request.normalized_pdb,
        max_bytes=active_limits.max_artifact_bytes,
    )

    diagnostics = _execute(
        argv,
        cwd=run_root,
        limits=active_limits,
        environment=_runner_environment(run_root, environment),
    )
    if diagnostics.timed_out:
        raise DiffusionRunnerError(
            f"external diffusion runner timed out after {active_limits.timeout_seconds:g} seconds",
            diagnostics=diagnostics,
        )
    if diagnostics.exit_code != 0:
        raise DiffusionRunnerError(
            f"external diffusion runner exited with code {diagnostics.exit_code}",
            diagnostics=diagnostics,
        )

    _validate_staged_input(
        run_root,
        prepared.request.normalized_pdb,
        max_bytes=active_limits.max_artifact_bytes,
    )
    result_path = run_root / RESULT_MANIFEST
    if result_path.is_symlink() or not result_path.is_file():
        raise DiffusionRunnerError(
            "external diffusion runner did not create a regular result.json manifest",
            diagnostics=diagnostics,
        )
    try:
        result = RunnerResult.from_json(
            _read_bounded_text(result_path, active_limits.max_manifest_bytes)
        )
    except (DiffusionContractError, OSError, UnicodeError) as exc:
        raise DiffusionRunnerError(
            f"external diffusion runner produced an invalid result.json manifest: {exc}",
            diagnostics=diagnostics,
        ) from exc

    if result.backend_provenance.runner_protocol_version != DIFFUSION_RUNNER_PROTOCOL_VERSION:
        raise DiffusionRunnerError(
            "external diffusion runner uses an incompatible protocol version: "
            f"{result.backend_provenance.runner_protocol_version}; "
            f"expected {DIFFUSION_RUNNER_PROTOCOL_VERSION}",
            diagnostics=diagnostics,
        )

    _validate_result_artifacts(
        run_root,
        prepared.request,
        result,
        max_artifact_bytes=active_limits.max_artifact_bytes,
    )
    return replace(result, runner_diagnostics=diagnostics)


def run_diffusion_runner(
    request: DiffusionRequest,
    command: Sequence[str | os.PathLike[str]],
    *,
    source_root: Path,
    workspace: Path,
    limits: RunnerLimits | None = None,
    environment: Mapping[str, str] | None = None,
) -> RunnerResult:
    """Prepare a dedicated workspace and run one trusted backend adapter.

    ``request.normalized_pdb.path`` is resolved below ``source_root``, checked
    against its declared digest, and copied into ``workspace``. The runner is
    launched without a shell and communicates only through fixed
    ``request.json`` and ``result.json`` manifests in its current directory.
    Candidate artifacts are returned only after path, file-type, identity-mask,
    and SHA-256 validation. This protocol boundary is not an OS/filesystem
    sandbox and cannot prevent the executable from writing elsewhere.
    """
    prepared = prepare_diffusion_workspace(
        request,
        source_root=source_root,
        workspace=workspace,
    )
    return run_prepared_diffusion_runner(
        prepared,
        command,
        limits=limits,
        environment=environment,
    )


def _preflight_command(command: Sequence[str | os.PathLike[str]]) -> tuple[str, ...]:
    if not command:
        raise DiffusionRunnerError("external diffusion runner command must not be empty")
    argv = tuple(os.fspath(part) for part in command)
    if any("\x00" in part for part in argv):
        raise DiffusionRunnerError("external diffusion runner command contains a NUL byte")

    executable = argv[0]
    if os.sep in executable or (os.altsep is not None and os.altsep in executable):
        executable_path = Path(executable).expanduser()
        try:
            executable_path = executable_path.resolve(strict=True)
        except FileNotFoundError as exc:
            raise DiffusionRunnerError(
                "external diffusion runner executable was not found"
            ) from exc
        if not executable_path.is_file():
            raise DiffusionRunnerError("external diffusion runner executable was not found")
        if not os.access(executable_path, os.X_OK):
            raise DiffusionRunnerError("external diffusion runner executable is not executable")
        executable = str(executable_path)
    else:
        resolved = shutil.which(executable)
        if resolved is None:
            raise DiffusionRunnerError("external diffusion runner executable was not found on PATH")
        executable = resolved
    return (executable, *argv[1:])


def _require_directory(path: Path, field_name: str) -> Path:
    if path.is_symlink() or not path.is_dir():
        raise DiffusionRunnerError(f"{field_name} must be a regular directory, not a symlink")
    return path.resolve(strict=True)


def _create_workspace(workspace: Path) -> Path:
    absolute = workspace.expanduser().absolute()
    if absolute.is_symlink() or absolute.exists():
        raise DiffusionRunnerError("runner workspace must not already exist")
    if absolute.parent.is_symlink() or not absolute.parent.is_dir():
        raise DiffusionRunnerError("runner workspace parent must be a regular directory")
    absolute.mkdir(mode=0o700)
    return absolute.resolve(strict=True)


def _reject_reserved_artifact_path(relative_path: str) -> None:
    parts = Path(relative_path).parts
    if not parts:
        raise DiffusionRunnerError("artifact path must name a contained file")
    first = parts[0]
    if first in {REQUEST_MANIFEST, RESULT_MANIFEST, ".home", ".tmp", ".cache"}:
        raise DiffusionRunnerError(f"artifact path uses reserved runner name {first!r}")


def _read_bounded_text(path: Path, max_bytes: int) -> str:
    descriptor = _open_regular_file(path, label="manifest", max_bytes=max_bytes)
    try:
        with os.fdopen(descriptor, "rb") as handle:
            return handle.read().decode("utf-8")
    except UnicodeDecodeError as exc:
        raise DiffusionRunnerError(f"manifest is not valid UTF-8: {path.name}") from exc


def _validated_artifact_path(
    root: Path,
    artifact: ArtifactReference,
    *,
    max_bytes: int | None = None,
) -> Path:
    _reject_reserved_artifact_path(artifact.path)
    candidate = root.joinpath(*Path(artifact.path).parts)
    current = root
    for part in Path(artifact.path).parts:
        current = current / part
        if current.is_symlink():
            raise DiffusionRunnerError(f"artifact path contains a symlink: {artifact.path}")
    try:
        resolved = candidate.resolve(strict=True)
    except FileNotFoundError as exc:
        raise DiffusionRunnerError(f"artifact does not exist: {artifact.path}") from exc
    if not resolved.is_relative_to(root):
        raise DiffusionRunnerError(f"artifact escapes its workspace: {artifact.path}")
    descriptor = _open_regular_file(
        candidate,
        label="artifact",
        display_name=artifact.path,
        max_bytes=max_bytes,
    )
    try:
        actual_digest = _sha256_descriptor(descriptor)
    finally:
        os.close(descriptor)
    if actual_digest != artifact.sha256.lower():
        raise DiffusionRunnerError(f"artifact SHA-256 mismatch: {artifact.path}")
    return candidate


def _validate_staged_input(
    root: Path,
    artifact: ArtifactReference,
    *,
    max_bytes: int,
) -> None:
    try:
        _validated_artifact_path(root, artifact, max_bytes=max_bytes)
    except DiffusionRunnerError as exc:
        raise DiffusionRunnerError(
            f"external diffusion runner modified its staged input: {exc}"
        ) from exc


def _validate_result_artifacts(
    root: Path,
    request: DiffusionRequest,
    result: RunnerResult,
    *,
    max_artifact_bytes: int,
) -> None:
    if len(result.candidates) > request.candidate_count:
        raise DiffusionRunnerError("external diffusion runner returned too many candidates")

    candidate_ids = [candidate.candidate_id for candidate in result.candidates]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise DiffusionRunnerError("external diffusion runner returned duplicate candidate IDs")

    seeds = [candidate.seed for candidate in result.candidates]
    if len(seeds) != len(set(seeds)):
        raise DiffusionRunnerError("external diffusion runner returned duplicate candidate seeds")
    if not set(seeds) <= set(request.seeds):
        raise DiffusionRunnerError("external diffusion runner returned an unrequested candidate seed")
    if result.status is DiffusionStatus.SUCCESS and len(result.candidates) != request.candidate_count:
        raise DiffusionRunnerError(
            "external diffusion runner did not return the requested candidate count"
        )

    paths = [
        path
        for candidate in result.candidates
        for path in (
            candidate.coordinate_artifact.path,
            candidate.sampler_trace_artifact.path,
        )
    ]
    if len(paths) != len(set(paths)):
        raise DiffusionRunnerError("external diffusion runner returned duplicate artifact paths")

    expected_residues = {
        residue
        for gap in request.gaps
        for residue in gap.generated_residues
    }
    expected_atoms = set(request.generated_atoms)
    for candidate in result.candidates:
        if set(candidate.generated_residues) != expected_residues:
            raise DiffusionRunnerError(
                f"candidate {candidate.candidate_id!r} changed the generated residue mask"
            )
        if set(candidate.generated_atoms) != expected_atoms:
            raise DiffusionRunnerError(
                f"candidate {candidate.candidate_id!r} changed the generated atom mask"
            )
        if candidate.coordinate_artifact.path == request.normalized_pdb.path:
            raise DiffusionRunnerError(
                f"candidate {candidate.candidate_id!r} aliases the staged input artifact"
            )
        _validated_artifact_path(
            root,
            candidate.coordinate_artifact,
            max_bytes=max_artifact_bytes,
        )
        trace_path = _validated_artifact_path(
            root,
            candidate.sampler_trace_artifact,
            max_bytes=min(max_artifact_bytes, MAX_TRACE_BYTES),
        )
        try:
            trace = SamplerTrace.from_json(
                _read_bounded_text(trace_path, MAX_TRACE_BYTES)
            )
        except (DiffusionContractError, OSError, UnicodeError) as exc:
            raise DiffusionRunnerError(
                f"candidate {candidate.candidate_id!r} has an invalid sampler trace: {exc}"
            ) from exc
        try:
            validate_sampler_trace_context(
                trace,
                request,
                result.backend_provenance,
            )
        except DiffusionContractError as exc:
            raise DiffusionRunnerError(
                f"candidate {candidate.candidate_id!r} has inconsistent sampler trace: {exc}"
            ) from exc


def _runner_environment(
    workspace: Path,
    additions: Mapping[str, str] | None,
) -> dict[str, str]:
    environment = {
        name: os.environ[name]
        for name in _ENVIRONMENT_ALLOWLIST
        if name in os.environ
    }
    if additions:
        rejected = sorted(name for name in additions if name not in _ENVIRONMENT_ALLOWLIST)
        if rejected:
            raise DiffusionRunnerError(
                "runner environment overrides are not allowlisted: "
                + ", ".join(rejected)
            )
        environment.update(additions)

    home = workspace / ".home"
    temporary = workspace / ".tmp"
    cache = workspace / ".cache"
    for directory in (home, temporary, cache):
        directory.mkdir(mode=0o700)
    environment.update(
        {
            "HOME": str(home),
            "TMPDIR": str(temporary),
            "XDG_CACHE_HOME": str(cache),
            "DVBFIXER_DIFFUSION_PROTOCOL_VERSION": str(
                DIFFUSION_RUNNER_PROTOCOL_VERSION
            ),
            "DVBFIXER_DIFFUSION_REQUEST": REQUEST_MANIFEST,
            "DVBFIXER_DIFFUSION_RESULT": RESULT_MANIFEST,
        }
    )
    return environment


def _open_regular_file(
    path: Path,
    *,
    label: str,
    max_bytes: int | None,
    display_name: str | None = None,
) -> int:
    name = display_name or path.name
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise DiffusionRunnerError(f"{label} could not be opened safely: {name}") from exc
    file_stat = os.fstat(descriptor)
    if not stat.S_ISREG(file_stat.st_mode):
        os.close(descriptor)
        raise DiffusionRunnerError(f"{label} is not a regular file: {name}")
    if file_stat.st_nlink != 1:
        os.close(descriptor)
        raise DiffusionRunnerError(f"{label} must not be hard-linked: {name}")
    if max_bytes is not None and file_stat.st_size > max_bytes:
        os.close(descriptor)
        raise DiffusionRunnerError(f"{label} exceeds the size limit: {name}")
    return descriptor


def _execute(
    argv: tuple[str, ...],
    *,
    cwd: Path,
    limits: RunnerLimits,
    environment: Mapping[str, str],
) -> RunnerDiagnostics:
    try:
        process = subprocess.Popen(
            argv,
            cwd=cwd,
            env=dict(environment),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
        )
    except OSError as exc:
        raise DiffusionRunnerError(
            f"external diffusion runner could not be started: {exc.strerror or type(exc).__name__}"
        ) from exc

    assert process.stdout is not None
    assert process.stderr is not None
    stdout_capture = _BoundedCapture(limits.max_output_bytes)
    stderr_capture = _BoundedCapture(limits.max_output_bytes)
    stop_readers = threading.Event()
    threads = (
        threading.Thread(
            target=stdout_capture.drain,
            args=(process.stdout, stop_readers),
            daemon=True,
        ),
        threading.Thread(
            target=stderr_capture.drain,
            args=(process.stderr, stop_readers),
            daemon=True,
        ),
    )
    for thread in threads:
        thread.start()

    timed_out = False
    workspace_exceeded = False
    deadline = time.monotonic() + limits.timeout_seconds
    try:
        while True:
            if _workspace_size_exceeds(cwd, limits.max_workspace_bytes):
                workspace_exceeded = True
                _terminate_process_group(process)
                exit_code = process.wait()
                break
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                timed_out = True
                _terminate_process_group(process)
                exit_code = process.wait()
                break
            try:
                exit_code = process.wait(timeout=min(remaining, 0.1))
                workspace_exceeded = _workspace_size_exceeds(
                    cwd,
                    limits.max_workspace_bytes,
                )
                break
            except subprocess.TimeoutExpired:
                continue
    finally:
        _terminate_process_group(process, include_descendants=True)
        for thread in threads:
            thread.join(timeout=_PIPE_DRAIN_TIMEOUT_SECONDS)
        readers_alive = any(thread.is_alive() for thread in threads)
        if readers_alive:
            timed_out = True
        stop_readers.set()
        for thread in threads:
            thread.join()
        process.stdout.close()
        process.stderr.close()

    diagnostics = RunnerDiagnostics(
        exit_code=exit_code,
        timed_out=timed_out,
        stdout=stdout_capture.text(),
        stderr=stderr_capture.text(),
    )
    if workspace_exceeded:
        raise DiffusionRunnerError(
            "external diffusion runner exceeded the workspace size limit",
            diagnostics=diagnostics,
        )
    return diagnostics


def _workspace_size_exceeds(root: Path, limit: int) -> bool:
    """Return early once regular files below ``root`` exceed ``limit`` bytes."""
    total = 0
    pending = [root]
    while pending:
        directory = pending.pop()
        try:
            entries = tuple(os.scandir(directory))
        except OSError:
            return True
        for entry in entries:
            try:
                if entry.is_dir(follow_symlinks=False):
                    pending.append(Path(entry.path))
                elif entry.is_file(follow_symlinks=False):
                    total += entry.stat(follow_symlinks=False).st_size
                    if total > limit:
                        return True
            except OSError:
                return True
    return False


def _terminate_process_group(
    process: subprocess.Popen[bytes],
    *,
    include_descendants: bool = False,
) -> None:
    if process.poll() is not None and not include_descendants:
        return
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except (AttributeError, ProcessLookupError, PermissionError):
        if process.poll() is None:
            process.kill()


def _sha256(path: Path) -> str:
    descriptor = _open_regular_file(path, label="artifact", max_bytes=None)
    try:
        return _sha256_descriptor(descriptor)
    finally:
        os.close(descriptor)


def _sha256_descriptor(descriptor: int) -> str:
    digest = hashlib.sha256()
    while chunk := os.read(descriptor, 1024 * 1024):
        digest.update(chunk)
    return digest.hexdigest()
