"""``dvbfixer top`` — GROMACS topology from PDB/GRO.

Current layout:

- :mod:`dvbfixer.top.cli` — argparse only + FF-directory constants
  (``FF_DIR``, ``FF_CHOICES``). Safe to import without rtp_parser.
- :mod:`dvbfixer.top.ff_data` — RTP/ARN/R2B/TDB loading,
  ``_dedup_atomtypes``, ``_GLYCAN_LINKAGE_PARAMS``, ``ION_PARAMS``.
- :mod:`dvbfixer.top.writers` — ``_write_moleculetype``, ``write_top``,
  ``write_pdb``, ``_read_ff_content``, ``_write_water_topology``.
- :mod:`dvbfixer.top.acpype` — the ``--acpype`` mode dispatcher.
- :mod:`dvbfixer.top.types` — the shared dataclasses (``PDBResidue``,
  ``PDBChain``, ``AtomEntry``, ``ChainTopology``) used by both
  ``pipeline`` and ``topology_builder``; split out to avoid a circular
  import between the two.
- :mod:`dvbfixer.top.glycan` — glycan/glycolipid link detection
  (``detect_glycan_links``, ``build_glycan_trees``, ``_is_ceramide``),
  run once from ``main()`` before any topology building starts.
- :mod:`dvbfixer.top.topology_builder` — ``TopologyBuilder``, the
  per-chain topology assembler (``build_chain``, ``build_glycan_chain``,
  ``build_glycolipid_chain``), plus its private atom-name-matching
  helpers.
- :mod:`dvbfixer.top.pipeline` — CLI orchestration (``main``): PDB
  reading/writing, chain classification, and wiring the modules above
  together. No longer the largest file in the package now that
  ``topology_builder``/``glycan``/``types`` are split out.

Public API re-exported so ``from dvbfixer.top import main`` keeps
working from ``cli.py`` and every downstream test.
"""

from __future__ import annotations

from dvbfixer.top.cli import FF_CHOICES, FF_DIR, parse_args
from dvbfixer.top.pipeline import main

__all__ = ["FF_CHOICES", "FF_DIR", "main", "parse_args"]
