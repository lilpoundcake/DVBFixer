"""Passive optimizer callbacks for private research, never a convergence gate."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from openmm import MinimizationReporter


class PassiveReporter(MinimizationReporter):
    def __init__(self, output: Path):
        super().__init__()
        output.mkdir(parents=True, exist_ok=False)
        self.output = output
        self.calls = 0
        self.last_iteration = None

    def report(self, iteration, x, grad, args):
        self.calls += 1
        self.last_iteration = int(iteration)
        coordinates = np.asarray(x).reshape(-1, 3)
        gradient = np.asarray(grad).reshape(-1, 3)
        if coordinates.shape != gradient.shape:
            raise ValueError("optimizer coordinates/gradient axes differ")
        stats = {str(key): float(value) if math.isfinite(float(value)) else None for key, value in args.items()}
        row = {"callback": self.calls, "optimizer_iteration_index": int(iteration), "stats": stats,
               "objective_gradient_component_rms": float(np.sqrt(np.mean(gradient**2))) if np.isfinite(gradient).all() else None,
               "coordinates_finite": bool(np.isfinite(coordinates).all()),
               "gradient_finite": bool(np.isfinite(gradient).all()), "requested_termination": False}
        with (self.output / "callbacks.jsonl").open("a") as stream:
            stream.write(json.dumps(row, allow_nan=False)+"\n")
        if self.calls <= 5 or self.calls % 50 == 0:
            with (self.output / f"callback-{self.calls:05d}.npz").open("xb") as stream:
                np.savez_compressed(stream, positions_nm=coordinates, objective_gradient=gradient)
        # In particular, never use a geometry finding to terminate the optimizer.
        return False
