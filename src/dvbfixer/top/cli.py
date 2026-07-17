"""CLI parsing + FF-directory constants for ``dvbfixer top``.

Split out of the flat ``top.py`` in Phase 2.4 of the revision plan.
Argparse only — safe to import without OpenMM / rtp_parser. Shared
constants (``FF_DIR``, ``FF_CHOICES``) live here so the pipeline and
its future extracted submodules (``rtp_build`` / ``ff_data`` /
``writers`` / ``acpype``) can pick them up without a circular dep.
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Bundled GROMACS force-field directories. ``FF_DIR`` points at the
# top-level ``FF/`` folder shipped with the package; ``FF_CHOICES``
# maps the ``--ff`` short name to its subdirectory.
FF_DIR = Path(__file__).parent.parent.parent.parent / "FF"
FF_CHOICES = {
    "amber": "amber99sb-ildn-lipid21.ff",
    "charmm": "charmm36_ljpme-jul2022.ff",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dvbfixer top",
        description="Generate GROMACS topology files from PDB",
    )
    parser.add_argument("input", help="Input PDB file")
    parser.add_argument("-o", "--output", help="Output .top file (default: topol.top)")
    parser.add_argument("--ff", choices=["amber", "charmm"], default="amber",
                        help="Force field (default: amber)")
    parser.add_argument("--ff-dir", help="Custom force field directory")
    parser.add_argument("--water", default="tip3p",
                        choices=["tip3p", "spc", "spce", "tip4p", "tip4pew", "opc"],
                        help="Water model (default: tip3p). With --ff charmm only "
                             "tip3p/spc/spce are accepted; OPC/TIP4P/TIP4P-Ew are "
                             "not parametrized for CHARMM36 ions.")
    parser.add_argument("--ion-set", default="auto", dest="ion_set",
                        choices=["auto", "jc-tip3p", "jc-spce", "jc-tip4pew",
                                 "lm-hfe-opc", "lm-iod-opc", "dang-legacy"],
                        help="Ion LJ parameter set (default: auto, picks the set "
                             "matched to the water model). Ignored with --ff charmm.")
    parser.add_argument("--ignh", action="store_true",
                        help="Ignore hydrogens in input PDB")
    parser.add_argument("--keep-all-hydrogens", dest="keep_all_hydrogens",
                        action="store_true",
                        help="Do not remove any hydrogen atoms from the "
                             "input (default OFF: HO1/HO2/HO3/HO4/HO6 at "
                             "glycosidic linkage sites are stripped and "
                             "their charge is redistributed onto the linked "
                             "O — matches the CHARMM RTP template for a "
                             "formed glycosidic bond). Use when the input "
                             "has a free reducing end that must keep its H, "
                             "or when charges must round-trip untouched. "
                             "WARNING: at a real glycosidic linkage this "
                             "produces an over-valent O (H + neighbour C "
                             "bonded to the same O); grompp may complain "
                             "and the resulting energy is chemically wrong. "
                             "Use only when you know why.")
    parser.add_argument("--no-infer-conect", dest="no_infer_conect",
                        action="store_true",
                        help="Skip automatic CONECT inference. By default "
                             "missing SS/glycosidic/glycosylation bonds are "
                             "perceived from coordinates before topology build.")
    parser.add_argument("--ss", action="append", default=[],
                        help="Disulfide bond: CHAIN1:NUM1:CHAIN2:NUM2 (repeatable)")
    parser.add_argument("--his", action="append", default=[],
                        help="HIS protonation: CHAIN:NUM:STATE (HIE/HID/HIP, repeatable)")
    parser.add_argument("--protonate", default=None,
                        help="Protonate residues. \"all\" protonates every ASP->ASPP, "
                             "GLU->GLUP, HIS->HSP. Comma-separated list protonates "
                             "specific residues: CHAIN:NUM[:STATE],... "
                             "(e.g. --protonate all, --protonate H:66,K:50:GLUP).")
    parser.add_argument("--merge", action="store_true",
                        help="Merge all chains into single moleculetype")
    parser.add_argument("--pdb", help="Output PDB file with topology-matched atom names")
    parser.add_argument("--acpype", action="store_true",
                        help="Use ACPYPE pipeline (AMBER14+GLYCAM -> ParmEd -> GROMACS). "
                             "Handles mixed 1-4 scaling via [ pairs_nb ].")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")
    return parser.parse_args(argv)
