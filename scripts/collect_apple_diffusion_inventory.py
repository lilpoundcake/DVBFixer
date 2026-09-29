#!/usr/bin/env python3
"""Collect fail-closed Apple Silicon diffusion runtime provenance."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

MPS_ENVIRONMENT_KEYS = (
    "PYTORCH_ENABLE_MPS_FALLBACK",
    "PYTORCH_MPS_FAST_MATH",
    "PYTORCH_MPS_PREFER_METAL",
    "PYTORCH_MPS_HIGH_WATERMARK_RATIO",
    "PYTORCH_MPS_LOW_WATERMARK_RATIO",
)
PACKAGE_NAMES = (
    "torch",
    "numpy",
    "MDAnalysis",
    "biopython",
    "biotite",
    "einops",
    "hydra-core",
    "omegaconf",
    "scipy",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _command(*argv: str) -> str:
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"inventory command failed ({' '.join(argv)}): {detail}")
    return completed.stdout.strip()


def _sysctl(name: str) -> str:
    return _command("/usr/sbin/sysctl", "-n", name)


def _package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for name in PACKAGE_NAMES:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def collect(environment_lock: Path, repository: Path) -> dict[str, Any]:
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise RuntimeError("Apple diffusion inventory requires native macOS arm64 Python")
    if os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK") != "0":
        raise RuntimeError("set PYTORCH_ENABLE_MPS_FALLBACK=0 before starting Python")
    if not environment_lock.is_file():
        raise FileNotFoundError(f"environment lock does not exist: {environment_lock}")

    import torch

    if not torch.backends.mps.is_built() or not torch.backends.mps.is_available():
        raise RuntimeError("PyTorch MPS is not built and available")
    revision = _command("git", "-C", str(repository), "rev-parse", "HEAD")
    status = _command("git", "-C", str(repository), "status", "--porcelain")
    if status:
        raise RuntimeError("repository must be clean before freezing Apple inventory")

    disk = shutil.disk_usage(repository)
    recommended = getattr(torch.mps, "recommended_max_memory", None)
    return {
        "schema_version": 1,
        "host": {
            "architecture": platform.machine(),
            "chip": _sysctl("machdep.cpu.brand_string"),
            "physical_memory_bytes": int(_sysctl("hw.memsize")),
            "macos": platform.mac_ver()[0],
            "platform": platform.platform(),
            "xcode_select_path": _command("/usr/bin/xcode-select", "-p"),
            "free_disk_bytes": disk.free,
        },
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
        },
        "torch": {
            "version": torch.__version__,
            "mps_built": bool(torch.backends.mps.is_built()),
            "mps_available": bool(torch.backends.mps.is_available()),
            "mps_device_count": int(torch.mps.device_count()),
            "mps_recommended_max_memory_bytes": (
                int(recommended()) if callable(recommended) else None
            ),
        },
        "mps_environment": {key: os.environ.get(key) for key in MPS_ENVIRONMENT_KEYS},
        "packages": _package_versions(),
        "environment_lock": {
            "path": str(environment_lock.resolve()),
            "sha256": _sha256(environment_lock),
        },
        "repository": {
            "path": str(repository.resolve()),
            "revision": revision,
            "clean": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--environment-lock", required=True, type=Path)
    parser.add_argument("--repository", default=Path.cwd(), type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"inventory output already exists: {args.output}")
    value = collect(args.environment_lock.resolve(), args.repository.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=args.output.parent,
        prefix=f".{args.output.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(args.output)
    print(f"{args.output}  sha256={_sha256(args.output)}")


if __name__ == "__main__":
    main()
