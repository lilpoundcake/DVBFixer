"""``dvbfixer prepare`` — PDBFixer wrapper with heterogen-H and mutations.

Phase 2.3 of the revision plan is progressive: this initial split
carves out the argparse layer as a mechanical first step and keeps
the rest of the logic in :mod:`~dvbfixer.prepare.pipeline` so the
change stays reviewable.

Current layout:

- :mod:`dvbfixer.prepare.cli` — argparse only + shared constants
  (``DEFAULT_PH``, ``GLYCOSYLATED_RESIDUES``, ``SUGAR_RESNAMES``). Safe
  to import without OpenMM / PDBFixer.
- :mod:`dvbfixer.prepare.pipeline` — everything else: glycan detection,
  mutation parsing, deletion / substitution cleanup, glycoprotein input
  preprocessing, PDBFixer invocation, heterogen-H addition, and ``main``.

Follow-up (not this commit): pipeline.py still has natural boundaries
for extraction — ``mutations.py`` (``parse_mutations``,
``apply_deletions_to_pdb_text``, and the glycan-walk BFS + SSBOND repair
helpers) and ``glycan.py``
(``find_glycosylated_atoms_with_sugar``,
``rename_glycosylated_protein_residues``,
``remove_extra_glycan_hydrogens``, ``add_heterogen_h_via_rdkit`` /
``_via_openbabel``, ``_build_glycan_trees``).

Public API re-exported so ``from dvbfixer.prepare import main`` keeps
working from ``cli.py``, ``zbs.py``, and downstream tests.
"""

from __future__ import annotations

from dvbfixer.prepare.cli import (
    DEFAULT_PH,
    GLYCOSYLATED_RESIDUES,
    SUGAR_RESNAMES,
    parse_args,
)
from dvbfixer.prepare.pipeline import (
    build_dat,
    main,
    parse_mutations,
    run_pdbfixer,
    write_dat,
)

__all__ = [
    "DEFAULT_PH",
    "GLYCOSYLATED_RESIDUES",
    "SUGAR_RESNAMES",
    "build_dat",
    "main",
    "parse_args",
    "parse_mutations",
    "run_pdbfixer",
    "write_dat",
]
