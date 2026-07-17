"""``dvbfixer top`` — GROMACS topology from PDB/GRO.

Phase 2.4 of the revision plan is progressive: this initial split
carves out the argparse layer + FF-directory constants as a mechanical
first step and keeps the RTP parser, topology builder, and output
writers in :mod:`~dvbfixer.top.pipeline` so the change stays
reviewable.

Current layout:

- :mod:`dvbfixer.top.cli` — argparse only + FF-directory constants
  (``FF_DIR``, ``FF_CHOICES``). Safe to import without rtp_parser.
- :mod:`dvbfixer.top.pipeline` — everything else: RTP/ARN/R2B/TDB
  loading, chain/glycan/glycolipid topology assembly, output writers,
  and ``main``.

Follow-up (queued):

- ``ff_data.py`` — RTP/ARN/R2B/TDB loading, ``_dedup_atomtypes``,
  ``_GLYCAN_LINKAGE_PARAMS``, ``ION_PARAMS``.
- ``rtp_build.py`` — ``build_chain``, ``build_glycan_chain``,
  ``build_glycolipid_chain``, chain topology assembly.
- ``writers.py`` — ``_write_moleculetype``, ``write_top``,
  ``write_pdb``, ``_read_ff_content``, ``_write_water_topology``.
- ``acpype.py`` — the ``--acpype`` mode dispatcher currently living
  as a 300+ line branch inside ``main``.

Public API re-exported so ``from dvbfixer.top import main`` keeps
working from ``cli.py`` and every downstream test.
"""

from __future__ import annotations

from dvbfixer.top.cli import FF_CHOICES, FF_DIR, parse_args
from dvbfixer.top.pipeline import main

__all__ = ["FF_CHOICES", "FF_DIR", "main", "parse_args"]
