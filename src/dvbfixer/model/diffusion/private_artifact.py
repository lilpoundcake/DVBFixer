"""Bounded private-artifact reading without importing a refinement runtime."""
import hashlib
import os
import stat
from pathlib import Path

from dvbfixer.model.diffusion.contract import ArtifactReference


def read_private_artifact(workspace: Path, artifact: ArtifactReference, *, label: str,
                          max_bytes: int, encoding: str) -> str:
    path = workspace.joinpath(*Path(artifact.path).parts)
    current = workspace
    for part in Path(artifact.path).parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f"raw {label} path contains a symlink")
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise ValueError(f"raw {label} does not exist") from exc
    if not resolved.is_relative_to(workspace):
        raise ValueError(f"raw {label} escapes the request workspace")
    info = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError(f"raw {label} is not a private regular file")
    if info.st_size > max_bytes:
        raise ValueError(f"raw {label} exceeds the size limit")
    data = path.read_bytes()
    if len(data) > max_bytes:
        raise ValueError(f"raw {label} exceeds the size limit")
    if hashlib.sha256(data).hexdigest() != artifact.sha256:
        raise ValueError(f"raw {label} digest mismatch")
    try:
        return data.decode(encoding)
    except UnicodeDecodeError as exc:
        raise ValueError(f"raw {label} is not valid {encoding}") from exc


def read_private_artifact_bytes(workspace: Path, artifact: ArtifactReference, *, max_bytes: int) -> bytes:
    """Bounded digest-checked binary evidence; never follow artifact symlinks."""
    data = read_private_file_bytes(workspace, artifact.path, max_bytes=max_bytes)
    if hashlib.sha256(data).hexdigest() != artifact.sha256.lower():
        raise ValueError("private evidence digest mismatch")
    return data


def read_private_file_bytes(workspace: Path, relative_path: str, *, max_bytes: int) -> bytes:
    """Bounded private-file reading for an index whose digest is not yet known."""
    ArtifactReference(relative_path, "0" * 64)
    current = workspace
    for part in Path(relative_path).parts:
        current /= part
        if current.is_symlink():
            raise ValueError("private evidence path contains a symlink")
    if not current.resolve(strict=True).is_relative_to(workspace.resolve(strict=True)):
        raise ValueError("private evidence escapes workspace")
    descriptor = os.open(current, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1 or metadata.st_size > max_bytes:
            raise ValueError("private evidence must be a bounded single-linked regular file")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            data = stream.read(max_bytes + 1)
        if len(data) > max_bytes:
            raise ValueError("private evidence exceeds byte bound")
        return data
    finally:
        os.close(descriptor)
