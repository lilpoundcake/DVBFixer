"""MDAnalysis input-boundary checks, before expensive trajectory analysis."""

from __future__ import annotations

from contextlib import suppress
from pathlib import Path
from typing import Any


class TrajectoryInputError(ValueError):
    """A dependency or topology/trajectory pair cannot be used for analysis."""


def load_universe(topology: str | Path, trajectory: str | Path | None = None) -> Any:
    try:
        import MDAnalysis as mda
    except Exception as exc:
        raise TrajectoryInputError(
            f"MDAnalysis could not be imported ({type(exc).__name__}: {exc}). "
            "Install/recreate Standard with compatible NumPy/SciPy/MDAnalysis binaries."
        ) from exc
    topology = Path(topology)
    trajectory = Path(trajectory) if trajectory is not None else None
    for path in (topology, trajectory):
        if path is not None and not path.is_file():
            raise TrajectoryInputError(f"Analysis input file not found: {path}")
    universe = None
    try:
        paths = [str(topology)] + ([str(trajectory)] if trajectory is not None else [])
        universe = mda.Universe(*paths)
        if len(universe.atoms) == 0:
            raise ValueError("topology contains no atoms")
        reader = universe.trajectory
        if reader.n_atoms != len(universe.atoms):
            raise ValueError(f"topology has {len(universe.atoms)} atoms, trajectory has {reader.n_atoms}")
        reader[0]
        return universe
    except Exception as exc:
        if universe is not None:
            # A partially initialized reader must not mask the original input error.
            with suppress(Exception):
                universe.trajectory.close()
        advice = ("This MDAnalysis version may not support the TPR writer version. "
                  "Use a supported matching TPR or a matching GRO/PDB topology, preserving atom order. "
                  if topology.suffix.lower() == ".tpr" else "Use a topology and trajectory with matching atom counts/order. ")
        raise TrajectoryInputError(
            f"Cannot load analysis topology {topology}"
            + (f" with trajectory {trajectory}" if trajectory is not None else "")
            + f" ({type(exc).__name__}: {exc}). {advice}"
        ) from exc
