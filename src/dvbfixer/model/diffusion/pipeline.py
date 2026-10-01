"""Internal fail-closed orchestration and atomic diffusion bundle publication."""

from __future__ import annotations

import ctypes
import errno
import hashlib
import json
import os
import shutil
import stat
import sys
import tempfile
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from dvbfixer.ffutils.dat import AddedAtom, DatRecord, ResidueSummary
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    DiffusionCandidate,
    DiffusionRequest,
    DiffusionResult,
    DiffusionStatus,
    ResidueIdentity,
    ValidationSummary,
)
from dvbfixer.model.diffusion.provenance import (
    PublicationArtifact,
    build_provenance_manifest,
)
from dvbfixer.model.diffusion.runner import (
    RunnerLimits,
    run_diffusion_runner,
)
from dvbfixer.model.diffusion.scope import (
    DiffusionScopeError,
    ScopeAdmission,
    assess_diffusion_scope,
)
from dvbfixer.model.diffusion.validate import (
    ValidationThresholds,
    build_validated_result,
)

BUNDLE_INDEX_SCHEMA_VERSION = 1
_RENAME_NOREPLACE = 1
_RENAME_EXCL = 4
_AT_FDCWD = -100


class DiffusionPipelineError(RuntimeError):
    """Raised when internal orchestration or publication fails closed."""


@dataclass(frozen=True, slots=True)
class DiffusionPipelineOutcome:
    result: DiffusionResult | None
    admission: ScopeAdmission | None
    published_bundle: Path | None
    message: str = ""

    @property
    def status(self) -> DiffusionStatus:
        if self.admission is not None and not self.admission.supported:
            return DiffusionStatus.UNSUPPORTED
        if self.result is None:
            return DiffusionStatus.FAILED
        return self.result.status


def run_diffusion_pipeline(
    request: DiffusionRequest,
    command: Sequence[str | os.PathLike[str]],
    *,
    source_root: Path,
    work_parent: Path,
    destination_bundle: Path,
    limits: RunnerLimits | None = None,
    environment: Mapping[str, str] | None = None,
    thresholds: ValidationThresholds | None = None,
    repository_root: Path | None = None,
) -> DiffusionPipelineOutcome:
    """Run one internal diffusion request and publish only validated outputs."""
    source_bytes = _read_source_artifact(
        source_root,
        request.normalized_pdb,
        max_bytes=(limits or RunnerLimits()).max_artifact_bytes,
    )
    try:
        admission = assess_diffusion_scope(request, source_bytes)
    except DiffusionScopeError as exc:
        raise DiffusionPipelineError(str(exc)) from exc
    if not admission.supported:
        return DiffusionPipelineOutcome(
            result=None,
            admission=admission,
            published_bundle=None,
            message=", ".join(admission.reasons),
        )

    work_parent = _require_directory(work_parent, "work_parent")
    workspace = Path(
        tempfile.mkdtemp(
            prefix=".dvbfixer-diffusion-run.",
            dir=work_parent,
        )
    )
    shutil.rmtree(workspace)
    try:
        runner_result = run_diffusion_runner(
            request,
            command,
            source_root=source_root,
            workspace=workspace,
            limits=limits,
            environment=environment,
        )
        result = build_validated_result(
            request,
            runner_result,
            workspace=workspace,
            thresholds=thresholds,
        )
        if result.status is not DiffusionStatus.SUCCESS:
            return DiffusionPipelineOutcome(
                result=result,
                admission=admission,
                published_bundle=None,
                message=result.message,
            )
        publish_diffusion_bundle(
            request,
            result,
            workspace=workspace,
            destination_bundle=destination_bundle,
            repository_root=repository_root,
        )
        return DiffusionPipelineOutcome(
            result=result,
            admission=admission,
            published_bundle=destination_bundle,
            message=result.message,
        )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def publish_diffusion_bundle(
    request: DiffusionRequest,
    result: DiffusionResult,
    *,
    workspace: Path,
    destination_bundle: Path,
    repository_root: Path | None = None,
    before_commit: Callable[[Path], None] | None = None,
) -> Path:
    """Publish selected PDB/.dat/provenance sets with one directory rename."""
    if result.status is not DiffusionStatus.SUCCESS:
        raise DiffusionPipelineError("only a successful validated result can be published")
    if len(result.candidates) != len(result.validation_summaries):
        raise DiffusionPipelineError("candidate and validation counts do not match")

    root = _require_directory(workspace, "validation workspace")
    destination = destination_bundle.expanduser().absolute()
    parent = _require_directory(destination.parent, "publication parent")
    if destination.exists() or destination.is_symlink():
        raise DiffusionPipelineError("diffusion destination bundle already exists")

    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=parent,
        )
    )
    committed = False
    try:
        bundle_candidates: list[dict[str, object]] = []
        for candidate, summary in zip(
            result.candidates,
            result.validation_summaries,
        ):
            bundle_candidates.append(
                _stage_candidate(
                    request,
                    candidate,
                    summary,
                    result,
                    workspace=root,
                    staging=staging,
                    repository_root=repository_root,
                )
            )
        index_path = staging / "bundle.json"
        _write_file_fsync(
            index_path,
            (
                json.dumps(
                    {
                        "schema_version": BUNDLE_INDEX_SCHEMA_VERSION,
                        "status": "success",
                        "candidates": bundle_candidates,
                    },
                    ensure_ascii=False,
                    sort_keys=True,
                    indent=2,
                )
                + "\n"
            ).encode("utf-8"),
        )
        _fsync_directory(staging)
        if before_commit is not None:
            before_commit(staging)
        _rename_bundle(staging, destination)
        committed = True
        _fsync_directory(parent)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        if committed:
            shutil.rmtree(destination, ignore_errors=True)
            try:
                _fsync_directory(parent)
            except OSError:
                pass
        raise
    return destination


def _stage_candidate(
    request: DiffusionRequest,
    candidate: DiffusionCandidate,
    summary: ValidationSummary,
    result: DiffusionResult,
    *,
    workspace: Path,
    staging: Path,
    repository_root: Path | None,
) -> dict[str, object]:
    candidate_bytes = _read_workspace_artifact(
        workspace,
        candidate.coordinate_artifact,
    )
    stem = _safe_stem(candidate.candidate_id)
    pdb_name = f"{stem}.pdb"
    dat_name = f"{stem}.dat"
    provenance_name = f"{stem}.diffusion.json"

    pdb_path = staging / pdb_name
    _write_file_fsync(pdb_path, candidate_bytes)
    dat_record = _candidate_dat_record(
        candidate_bytes,
        set(candidate.generated_residues),
    )
    dat_bytes = (
        json.dumps(
            dat_record.to_dict(),
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
        )
        + "\n"
    ).encode("utf-8")
    dat_path = staging / dat_name
    _write_file_fsync(dat_path, dat_bytes)

    artifacts = (
        PublicationArtifact("pdb", pdb_name, _sha256_bytes(candidate_bytes)),
        PublicationArtifact("dat", dat_name, _sha256_bytes(dat_bytes)),
    )
    manifest = build_provenance_manifest(
        request,
        candidate,
        summary,
        result.backend_provenance,
        result.runner_diagnostics,
        artifacts,
        resource_metrics=result.resource_metrics,
        repository_root=repository_root,
    )
    manifest_bytes = manifest.to_json().encode("utf-8")
    manifest_path = staging / provenance_name
    _write_file_fsync(manifest_path, manifest_bytes)
    return {
        "candidate_id": candidate.candidate_id,
        "seed": candidate.seed,
        "pdb": {
            "path": pdb_name,
            "sha256": artifacts[0].sha256,
        },
        "dat": {
            "path": dat_name,
            "sha256": artifacts[1].sha256,
        },
        "provenance": {
            "path": provenance_name,
            "sha256": _sha256_bytes(manifest_bytes),
        },
    }


def _candidate_dat_record(
    candidate_bytes: bytes,
    generated_residues: set[ResidueIdentity],
) -> DatRecord:
    try:
        text = candidate_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise DiffusionPipelineError("candidate PDB is not valid UTF-8") from exc

    added_atoms: list[AddedAtom] = []
    residue_summary: dict[str, ResidueSummary] = {}
    seen: set[AtomIdentity] = set()
    for line_number, line in enumerate(text.splitlines(), 1):
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 54:
            raise DiffusionPipelineError(
                f"candidate PDB has a truncated coordinate on line {line_number}"
            )
        identity = AtomIdentity(
            line[21],
            line[22:26].strip(),
            line[26].strip(),
            line[12:16].strip(),
        )
        residue = ResidueIdentity(
            identity.chain,
            identity.residue_number,
            identity.insertion_code,
        )
        if residue not in generated_residues:
            continue
        if identity in seen:
            raise DiffusionPipelineError("candidate PDB has duplicate generated atom identity")
        seen.add(identity)
        resname = line[17:20].strip()
        element = line[76:78].strip() if len(line) >= 78 else ""
        if not element:
            element = next(
                (
                    character
                    for character in identity.atom_name
                    if character.isalpha()
                ),
                "",
            )
        element = element.upper()
        added_atoms.append(
            {
                "chain": identity.chain,
                "resid": identity.residue_number,
                "icode": identity.insertion_code,
                "resname": resname,
                "atom": identity.atom_name,
                "element": element,
            }
        )
        key = (
            f"{identity.chain}/{resname}"
            f"{identity.residue_number}{identity.insertion_code}"
        )
        bucket = residue_summary.setdefault(
            key,
            {"heavy": 0, "hydrogen": 0},
        )
        if element == "H":
            bucket["hydrogen"] += 1
        else:
            bucket["heavy"] += 1

    represented = {
        ResidueIdentity(atom["chain"], atom["resid"], atom["icode"])
        for atom in added_atoms
    }
    if represented != generated_residues:
        raise DiffusionPipelineError(
            "candidate PDB does not represent every generated residue in DatRecord"
        )
    return DatRecord(
        description=(
            "Diffusion gap-reconstruction data. Every atom in a generated "
            "residue is recorded as added."
        ),
        added_atoms=added_atoms,
        residue_summary=residue_summary,
    )


def _read_source_artifact(
    root: Path,
    artifact: ArtifactReference,
    *,
    max_bytes: int,
) -> bytes:
    source_root = _require_directory(root, "source_root")
    return _read_contained_artifact(
        source_root,
        artifact,
        max_bytes=max_bytes,
    )


def _read_workspace_artifact(
    workspace: Path,
    artifact: ArtifactReference,
) -> bytes:
    return _read_contained_artifact(
        workspace,
        artifact,
        max_bytes=None,
    )


def _read_contained_artifact(
    root: Path,
    artifact: ArtifactReference,
    *,
    max_bytes: int | None,
) -> bytes:
    path = root.joinpath(*Path(artifact.path).parts)
    current = root
    for part in Path(artifact.path).parts:
        current = current / part
        if current.is_symlink():
            raise DiffusionPipelineError(
                f"artifact path contains a symlink: {artifact.path}"
            )
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise DiffusionPipelineError(
            f"artifact does not exist: {artifact.path}"
        ) from exc
    if not resolved.is_relative_to(root):
        raise DiffusionPipelineError(
            f"artifact escapes its workspace: {artifact.path}"
        )
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise DiffusionPipelineError(
            f"artifact could not be opened safely: {artifact.path}"
        ) from exc
    try:
        file_stat = os.fstat(descriptor)
        if not stat.S_ISREG(file_stat.st_mode) or file_stat.st_nlink != 1:
            raise DiffusionPipelineError(
                f"artifact is not a private regular file: {artifact.path}"
            )
        if max_bytes is not None and file_stat.st_size > max_bytes:
            raise DiffusionPipelineError(
                f"artifact exceeds the size limit: {artifact.path}"
            )
        data = bytearray()
        digest = hashlib.sha256()
        while chunk := os.read(descriptor, 1024 * 1024):
            data.extend(chunk)
            digest.update(chunk)
    finally:
        os.close(descriptor)
    if digest.hexdigest() != artifact.sha256.lower():
        raise DiffusionPipelineError(
            f"artifact SHA-256 mismatch: {artifact.path}"
        )
    return bytes(data)


def _rename_bundle(staging: Path, destination: Path) -> None:
    if os.name != "posix":
        raise DiffusionPipelineError(
            "atomic no-replace bundle publication is unavailable on this platform"
        )
    if sys.platform == "darwin":
        _rename_bundle_macos(staging, destination)
        return
    _rename_bundle_linux(staging, destination)


def _rename_bundle_linux(staging: Path, destination: Path) -> None:
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        renameat2 = libc.renameat2
    except (AttributeError, OSError) as exc:
        raise DiffusionPipelineError(
            "atomic no-replace bundle publication requires renameat2"
        ) from exc
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    error_number = _call_rename(
        renameat2,
        (
            _AT_FDCWD,
            os.fsencode(staging),
            _AT_FDCWD,
            os.fsencode(destination),
            _RENAME_NOREPLACE,
        ),
    )
    _raise_rename_error(error_number, destination)


def _rename_bundle_macos(staging: Path, destination: Path) -> None:
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        renamex_np = libc.renamex_np
    except (AttributeError, OSError) as exc:
        raise DiffusionPipelineError(
            "atomic no-replace bundle publication requires renamex_np"
        ) from exc
    renamex_np.argtypes = (
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renamex_np.restype = ctypes.c_int
    error_number = _call_rename(
        renamex_np,
        (
            os.fsencode(staging),
            os.fsencode(destination),
            _RENAME_EXCL,
        ),
    )
    _raise_rename_error(error_number, destination)


def _call_rename(function: object, arguments: tuple[object, ...]) -> int:
    ctypes.set_errno(0)
    result = function(*arguments)  # type: ignore[operator]
    return 0 if result == 0 else ctypes.get_errno()


def _raise_rename_error(error_number: int, destination: Path) -> None:
    if error_number == 0:
        return
    if error_number == errno.EEXIST:
        raise DiffusionPipelineError("diffusion destination bundle already exists")
    if error_number in {errno.ENOSYS, errno.EINVAL, errno.ENOTSUP}:
        raise DiffusionPipelineError(
            "atomic no-replace bundle publication is unavailable on this filesystem"
        )
    raise OSError(
        error_number,
        os.strerror(error_number),
        os.fspath(destination),
    )


def _write_file_fsync(path: Path, data: bytes) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(path, flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        path.unlink(missing_ok=True)
        raise


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _safe_stem(candidate_id: str) -> str:
    if not candidate_id or candidate_id in {".", ".."}:
        raise DiffusionPipelineError("candidate ID is not a safe bundle stem")
    if any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-" for character in candidate_id):
        raise DiffusionPipelineError("candidate ID is not a safe bundle stem")
    return candidate_id


def _require_directory(path: Path, field_name: str) -> Path:
    if path.is_symlink() or not path.is_dir():
        raise DiffusionPipelineError(
            f"{field_name} must be a regular directory, not a symlink"
        )
    return path.resolve(strict=True)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()
