#!/usr/bin/env python3
"""Compare paired CPU and MPS Protpardelle operator-smoke coordinates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import tempfile
from pathlib import Path
from typing import Any

import numpy as np


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> tuple[dict[str, Any], np.ndarray]:
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("status") != "passed":
        raise ValueError(f"operator smoke did not pass: {path}")
    coordinate_path = path.parent / str(summary["coordinate_artifact"])
    if _sha256(coordinate_path) != summary["coordinate_artifact_sha256"]:
        raise ValueError(f"coordinate artifact digest mismatch: {coordinate_path}")
    coordinates = np.load(coordinate_path, allow_pickle=False)
    if list(coordinates.shape) != summary["shape"] or not np.all(np.isfinite(coordinates)):
        raise ValueError(f"coordinate artifact is invalid: {coordinate_path}")
    return summary, coordinates


def compare(cpu_path: Path, mps_path: Path) -> dict[str, Any]:
    cpu, cpu_coordinates = _load(cpu_path)
    mps, mps_coordinates = _load(mps_path)
    if cpu["device"].split(":", 1)[0] != "cpu" or mps["device"].split(":", 1)[0] != "mps":
        raise ValueError("comparison requires CPU first and MPS second")
    keys = ("config_sha256", "checkpoint_sha256", "length", "steps", "seed", "shape")
    mismatches = [key for key in keys if cpu[key] != mps[key]]
    if mismatches:
        raise ValueError(f"operator-smoke inputs differ: {mismatches}")
    displacement = np.linalg.norm(mps_coordinates.astype(np.float64) - cpu_coordinates, axis=-1)
    return {
        "schema_version": 1,
        "cpu_summary_sha256": _sha256(cpu_path),
        "mps_summary_sha256": _sha256(mps_path),
        "coordinate_rmsd": math.sqrt(float(np.mean(np.square(displacement)))),
        "coordinate_maximum_displacement": float(np.max(displacement)),
        "coordinate_mean_displacement": float(np.mean(displacement)),
        "coordinate_byte_identical": bool(np.array_equal(cpu_coordinates, mps_coordinates)),
        "cpu_sample_seconds": cpu["sample_seconds"],
        "mps_sample_seconds": mps["sample_seconds"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cpu_summary", type=Path)
    parser.add_argument("mps_summary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"comparison output already exists: {args.output}")
    result = compare(args.cpu_summary, args.mps_summary)
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
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(args.output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
