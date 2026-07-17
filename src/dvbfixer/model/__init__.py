"""``dvbfixer model`` — Modeller-based loop / gap rebuilding.

Phase 2.2 of the revision plan is progressive: this initial split
carves out the argparse layer as a mechanical first step and keeps
the numbering + Modeller-run logic in :mod:`~dvbfixer.model.pipeline`
so the change stays reviewable. Follow-up work can extract
``renumber.py`` (``build_resnum_mapping`` and its private helpers) and
``modeller_run.py`` (``run_modeller``, ``_PinnedLoopModel``,
``_fix_terminal_alignment``, ``_explain_modeller_error``) — the
boundaries are already clear inside pipeline.py.

Current layout:

- :mod:`dvbfixer.model.cli` — argparse only + shared constants
  (``AA3TO1``, ``WATER_RESNAMES``). Safe to import without Modeller.
- :mod:`dvbfixer.model.pipeline` — everything else: SEQRES / FASTA
  parsing, numbering, Modeller invocation, PDB post-processing, ``main()``.

The public API is re-exported so ``from dvbfixer.model import main``
keeps working from ``cli.py``, ``zbs.py``, and downstream tests.
"""

from __future__ import annotations

from dvbfixer.model.cli import AA3TO1, WATER_RESNAMES, parse_args
from dvbfixer.model.modeller_run import parse_pir_sequence
from dvbfixer.model.pipeline import (
    build_model_dat,
    main,
    remove_water_lines,
)

__all__ = [
    "AA3TO1",
    "WATER_RESNAMES",
    "build_model_dat",
    "main",
    "parse_args",
    "parse_pir_sequence",
    "remove_water_lines",
]
