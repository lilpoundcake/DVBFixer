#!/usr/bin/env python3
"""Install a private launcher for one separately provisioned diffusion backend."""

from __future__ import annotations

import argparse
import hashlib
import shlex
import stat
import sys
from dataclasses import dataclass
from pathlib import Path

from dvbfixer.model.diffusion.runtime import diffusion_runner_directory


@dataclass(frozen=True, slots=True)
class Backend:
    runner: str
    launcher: str
    checkpoint_sha256: str
    engine_python_subdir: str
    mps: bool = False


BACKENDS = {
    "protpardelle": Backend(
        "deploy/protpardelle-1c/production_runner.py",
        "dvbfixer-diffusion-protpardelle-1c-mps",
        "dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483",
        "src",
        mps=True,
    ),
    "protenix": Backend(
        "deploy/protenix-v1/production_runner.py",
        "dvbfixer-diffusion-protenix-v1-cuda",
        "2b7d5a8b30494514fc47fd2271a16260528cdba170ba09cc112fdecd8f85ec04",
        ".",
    ),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Register an installed diffusion backend with DVBFixer",
    )
    parser.add_argument("model", choices=tuple(BACKENDS))
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument(
        "--engine-root",
        type=Path,
        help="Pinned, patched engine checkout (needed when it is not installed cleanly)",
    )
    args = parser.parse_args()

    repository = Path(__file__).resolve().parent.parent
    backend = BACKENDS[args.model]
    runner = repository / backend.runner
    checkpoint = args.checkpoint.expanduser().resolve()
    if not runner.is_file():
        parser.error(f"runner is missing from this DVBFixer checkout: {runner}")
    if checkpoint.is_symlink() or not checkpoint.is_file():
        parser.error(f"checkpoint is not a regular file: {checkpoint}")
    if _sha256(checkpoint) != backend.checkpoint_sha256:
        parser.error(f"checkpoint does not match the frozen {args.model} model")

    environment = ["export PYTHONNOUSERSITE=1"]
    if args.engine_root is not None:
        engine_root = args.engine_root.expanduser().resolve()
        engine_python_path = engine_root / backend.engine_python_subdir
        if not engine_python_path.is_dir():
            parser.error(f"engine Python source directory is missing: {engine_python_path}")
        environment.append(f"export PYTHONPATH={shlex.quote(str(engine_python_path))}")
    if backend.mps:
        environment.extend((
            "export PYTORCH_ENABLE_MPS_FALLBACK=0",
            "unset PYTORCH_MPS_FAST_MATH PYTORCH_MPS_PREFER_METAL",
            "unset PYTORCH_MPS_HIGH_WATERMARK_RATIO PYTORCH_MPS_LOW_WATERMARK_RATIO",
        ))
    command = " ".join(
        shlex.quote(value)
        for value in (
            sys.executable,
            str(runner),
            "--checkpoint",
            str(checkpoint),
        )
    )
    contents = "\n".join(("#!/bin/sh", *environment, f'exec {command} "$@"', ""))
    destination = diffusion_runner_directory() / backend.launcher
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".tmp")
    temporary.write_text(contents, encoding="utf-8")
    temporary.chmod(stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    temporary.replace(destination)
    print(f"Installed {args.model} diffusion backend: {destination}")


if __name__ == "__main__":
    main()
