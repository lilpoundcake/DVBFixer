"""Manifest-last directory commit for POSIX filesystems without rename flags.

The exclusive directory reservation owns the destination name. Payloads are
durable and single-linked before bundle.json is linked into place last. Readers
must require that index; directory existence alone never establishes readiness.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

INDEX = "bundle.json"


def _identity(metadata: os.stat_result) -> tuple[int, int]:
    return metadata.st_dev, metadata.st_ino


def commit_manifest_last(staging: Path, destination: Path) -> None:
    """Consume a flat private staging directory without replacing any entry.

    A crash before the index commit can leave a reserved, incomplete directory.
    It is deliberately not reusable automatically. A competing publisher must
    fail rather than adopt, remove or overwrite that directory.
    """
    entries = {entry.name: entry.lstat() for entry in staging.iterdir()}
    if INDEX not in entries or any(
        not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1
        for metadata in entries.values()
    ):
        raise OSError("bundle staging must contain a private index and flat regular files")
    # mkdir, unlike an exists()/rename() pair, reserves the name atomically even
    # on CephFS. Do not catch FileExistsError here: this directory is not ours.
    destination.mkdir(mode=0o700)
    reserved = destination.lstat()
    descriptor = None
    linked: list[str] = []
    try:
        descriptor = os.open(destination, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        if _identity(os.fstat(descriptor)) != _identity(reserved):
            raise OSError("bundle destination changed after reservation")
        for name in entries:
            if name != INDEX:
                _transfer(staging, name, descriptor, linked)
        # All payload links/removals precede the only readiness marker. The
        # staged bytes were fsynced by the publisher before this helper is used.
        os.fsync(descriptor)
        if _identity(destination.lstat()) != _identity(reserved):
            raise OSError("bundle destination changed before index commit")
        _transfer(staging, INDEX, descriptor, linked)
        os.fsync(descriptor)
        staging.rmdir()
    except BaseException:
        try:
            _rollback(destination, reserved, descriptor, entries, linked)
        except OSError:
            # Preserve the original error. If index removal fails, rollback
            # stops before removing payloads, so it cannot leave a ready shell.
            pass
        raise
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _transfer(staging: Path, name: str, descriptor: int, linked: list[str]) -> None:
    # link() is atomic/exclusive for files: neither payloads nor the final index
    # can overwrite an entry introduced by another process. Unlink the staging
    # name immediately; bundle readers reject multiply-linked artifacts.
    os.link(staging / name, name, dst_dir_fd=descriptor, follow_symlinks=False)
    linked.append(name)
    (staging / name).unlink()


def _rollback(
    destination: Path,
    reserved: os.stat_result,
    descriptor: int | None,
    entries: dict[str, os.stat_result],
    linked: list[str],
) -> None:
    if descriptor is not None:
        # Remove the readiness marker first, and only entries we actually own.
        # A replaced/injected index belongs to someone else: preserve everything.
        try:
            index = os.stat(INDEX, dir_fd=descriptor, follow_symlinks=False)
        except FileNotFoundError:
            index = None
        if index is not None and (INDEX not in linked or _identity(index) != _identity(entries[INDEX])):
            return
        for name in reversed(linked):
            try:
                metadata = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
            except FileNotFoundError:
                continue
            if _identity(metadata) == _identity(entries[name]):
                os.unlink(name, dir_fd=descriptor)
        os.fsync(descriptor)
    try:
        current = destination.lstat()
    except FileNotFoundError:
        return
    if _identity(current) == _identity(reserved):
        # rmdir cannot erase an unrelated entry injected into our directory.
        destination.rmdir()
