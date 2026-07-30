"""``dvbfixer top`` — GROMACS topology from PDB/GRO.

Current layout:

- :mod:`dvbfixer.top.cli` — argparse only + FF-directory constants
  (``FF_DIR``, ``FF_CHOICES``). Safe to import without rtp_parser.
- :mod:`dvbfixer.top.ff_data` — RTP/ARN/R2B/TDB loading,
  ``_dedup_atomtypes``, ``_GLYCAN_LINKAGE_PARAMS``, ``ION_PARAMS``.
- :mod:`dvbfixer.top.writers` — ``_write_moleculetype``, ``write_top``,
  ``write_pdb``, ``_read_ff_content``, ``_write_water_topology``.
- :mod:`dvbfixer.top.acpype` — the ``--acpype`` mode dispatcher.
- :mod:`dvbfixer.top.pipeline` — chain/glycan/glycolipid topology
  assembly (``build_chain``, ``build_glycan_chain``,
  ``build_glycolipid_chain``) and ``main``. Still the largest file
  (~2900 lines) — extracting the topology-builder functions above into
  a ``rtp_build.py`` is the one remaining follow-up from the original
  split.

Public API re-exported so ``from dvbfixer.top import main`` keeps
working from ``cli.py`` and every downstream test.
"""

from __future__ import annotations

from dvbfixer.top.cli import FF_CHOICES, FF_DIR, parse_args
from dvbfixer.top.pipeline import main

__all__ = ["FF_CHOICES", "FF_DIR", "main", "parse_args"]
