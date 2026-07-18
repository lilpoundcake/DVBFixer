"""``dvbfixer minimize`` — energy minimisation with selective restraints.

Phase 2.1 split the flat 1933-line ``minimize.py`` into a package:

- :mod:`dvbfixer.minimize.cli` — argparse only; safe to import without
  pulling in OpenMM. Defines the ``parse_args`` function plus the shared
  module-level constants (``DEFAULT_*``, ``BACKBONE_NAMES``).
- :mod:`dvbfixer.minimize.refine` — xtb + OpenBabel refinement passes,
  the Kabsch glycan-tracking helper, and the heterogen-subsystem
  extract/splice utilities.
- :mod:`dvbfixer.minimize.pipeline` — the ``minimize()`` and ``main()``
  entry points, restraint tiers, the .dat handoff, and the OpenMM
  Simulation wiring.

The public API — ``main``, ``parse_args``, ``refine_with_xtb``,
``refine_with_obminimize``, and the concrete pipeline helpers — is
re-exported here so ``from dvbfixer.minimize import main`` continues to
work from ``cli.py`` and ``zbs.py``.
"""

from __future__ import annotations

from dvbfixer.minimize.cli import (
    BACKBONE_NAMES,
    DEFAULT_FF,
    DEFAULT_MAX_ITER,
    DEFAULT_PADDING,
    DEFAULT_PH,
    DEFAULT_RESTRAINT_K,
    DEFAULT_WEAK_K,
    parse_args,
)
from dvbfixer.minimize.pipeline import (
    build_restraint_force,
    load_dat,
    main,
    minimize,
    resolve_new_atom_indices,
    strip_solvent,
)
from dvbfixer.minimize.refine import refine_with_obminimize, refine_with_xtb

__all__ = [
    "BACKBONE_NAMES",
    "DEFAULT_FF",
    "DEFAULT_MAX_ITER",
    "DEFAULT_PADDING",
    "DEFAULT_PH",
    "DEFAULT_RESTRAINT_K",
    "DEFAULT_WEAK_K",
    "build_restraint_force",
    "load_dat",
    "main",
    "minimize",
    "parse_args",
    "refine_with_obminimize",
    "refine_with_xtb",
    "resolve_new_atom_indices",
    "strip_solvent",
]
