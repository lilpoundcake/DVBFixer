"""Shared PDB utilities for dvbfixer.

Package layout (Phase 1.4):

- :mod:`dvbfixer.pdbutils.inference` — CONECT-record inference (OpenBabel
  bond perception, distance fallback, standard-AA filter, and the three
  domain overrides for SS bonds / glycosidic linkages / glycosylation).
- :mod:`dvbfixer.pdbutils.io` — line-level PDB read/write helpers
  (``build_serial_map``, ``remap_conect_records``, ``append_before_end``).

The public API is re-exported here so every existing caller
(``from dvbfixer.pdbutils import X``) keeps working across the split.
"""

from __future__ import annotations

from dvbfixer.pdbutils.inference import (
    _has_any_conect,
    _materialise_inferred_pdb,
    _strip_existing_conect,
    infer_conect_records,
    write_conect_block,
)
from dvbfixer.pdbutils.io import (
    append_before_end,
    build_serial_map,
    remap_conect_records,
)

__all__ = [
    "infer_conect_records",
    "write_conect_block",
    "_materialise_inferred_pdb",
    "_has_any_conect",
    "_strip_existing_conect",
    "build_serial_map",
    "remap_conect_records",
    "append_before_end",
]
