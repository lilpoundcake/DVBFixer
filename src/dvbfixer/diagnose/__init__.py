"""``dvbfixer diagnose`` — standalone structure diagnostic tool.

Inspired by BioLuminate's Protein Report widget and MolProbity's
all-atom validation. Emits a plain-text per-residue findings report;
never mutates the input. Users pick the appropriate dvbfixer tool
(``prepare`` / ``minimize`` / ``protonate`` / ``pull``) to actually
fix any issue.

Layout:

- :mod:`dvbfixer.diagnose.cli` — argparse only, safe to import
  without OpenMM / PDBFixer.
- :mod:`dvbfixer.diagnose.report` — :class:`Finding` dataclass,
  :class:`Severity` enum, plain-text formatter.
- :mod:`dvbfixer.diagnose.structural` — missing atoms/residues,
  coincident atoms, misplaced H, altLoc conflicts, chain breaks.
- :mod:`dvbfixer.diagnose.chemistry` — valence, bond geometry,
  cis-peptides, chirality.
- :mod:`dvbfixer.diagnose.steric` — all-atom clash detection
  (Python NeighborList + optional MolProbity ``probe`` binary).
- :mod:`dvbfixer.diagnose.pipeline` — orchestrator + ``main()``.
"""

from __future__ import annotations

from dvbfixer.diagnose.cli import parse_args
from dvbfixer.diagnose.pipeline import main

__all__ = ["main", "parse_args"]
