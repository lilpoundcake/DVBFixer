"""Bounded native single-point control for a context worker's exact full XYZ."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--identities", type=Path, required=True)
    parser.add_argument("--xtb", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    coordinates = np.load(args.reference, allow_pickle=False)
    identities = json.loads(args.identities.read_text())["input_atom_identities"]
    if coordinates.shape != (len(identities), 3) or not np.isfinite(coordinates).all():
        raise ValueError("Reference coordinates do not match identities")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "input.xyz").write_text(f"{len(identities)}\nExact full-context reference\n" + "".join(
        f"{identity[-1]} {x:.12f} {y:.12f} {z:.12f}\n" for identity, (x, y, z) in zip(identities, coordinates)))
    command = [str(args.xtb.resolve()), "input.xyz", "--gfnff", "--sp", "--norestart", "--verbose"]
    env = {**os.environ, **{k: str(args.threads) for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]}}
    report = {"command": command, "threads": args.threads,
              "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
              "binary_sha256": hashlib.sha256(args.xtb.read_bytes()).hexdigest()}
    started = time.monotonic()
    try:
        run = subprocess.run(command, cwd=args.output, env=env, capture_output=True, timeout=args.timeout)
        report["returncode"] = run.returncode
        (args.output / "stdout.raw").write_bytes(run.stdout)
        (args.output / "stderr.raw").write_bytes(run.stderr)
    except subprocess.TimeoutExpired as exc:
        report["timeout"] = True
        (args.output / "stdout.raw").write_bytes(exc.stdout or b"")
        (args.output / "stderr.raw").write_bytes(exc.stderr or b"")
    report["seconds"] = time.monotonic() - started
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
