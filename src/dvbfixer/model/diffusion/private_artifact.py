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
