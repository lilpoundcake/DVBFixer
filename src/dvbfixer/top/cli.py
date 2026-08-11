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
    io = parser.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB file")
    io.add_argument("-o", "--output", help="Output .top file (default: topol.top)")
    io.add_argument("--pdb", help="Output PDB file with topology-matched atom names")

    ff = parser.add_argument_group("Force field / solvation")
    ff.add_argument("--ff", choices=["amber", "charmm"], default="amber",
                    help="Force field (default: amber)")
    ff.add_argument("--ff-dir", help="Custom force field directory")
    ff.add_argument("--water", default="tip3p",
                    choices=["tip3p", "spc", "spce", "tip4p", "tip4pew", "opc"],
                    help="Water model (default: tip3p). With --ff charmm only "
                         "tip3p/spc/spce are accepted; OPC/TIP4P/TIP4P-Ew are "
                         "not parametrized for CHARMM36 ions.")
    ff.add_argument("--ion-set", default="auto", dest="ion_set",
                    choices=["auto", "jc-tip3p", "jc-spce", "jc-tip4pew",
                             "lm-hfe-opc", "lm-iod-opc", "dang-legacy"],
                    help="Ion LJ parameter set (default: auto, picks the set "
                         "matched to the water model). Ignored with --ff charmm.")

    protonation = parser.add_argument_group("Protonation / bonds")
    protonation.add_argument("--ss", action="append", default=[],
                             help="Disulfide bond: CHAIN1:NUM1:CHAIN2:NUM2 (repeatable)")
    protonation.add_argument("--his", action="append", default=[],
                             help="HIS protonation: CHAIN:NUM:STATE (HIE/HID/HIP, repeatable)")
    protonation.add_argument("--protonate", default=None,
                             help="Protonate residues. \"all\" protonates every ASP->ASPP, "
                                  "GLU->GLUP, HIS->HSP. Comma-separated list protonates "
                                  "specific residues: CHAIN:NUM[:STATE],... "
                                  "(e.g. --protonate all, --protonate H:66,K:50:GLUP).")

    content = parser.add_argument_group("Content / behaviour")
    content.add_argument("--ignh", action="store_true",
                         help="Ignore hydrogens in input PDB")
    content.add_argument("--keep-all-hydrogens", dest="keep_all_hydrogens",
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
    content.add_argument("--no-infer-conect", dest="no_infer_conect",
                         action="store_true",
                         help="Skip automatic CONECT inference. By default "
                              "missing SS/glycosidic/glycosylation bonds are "
                              "perceived from coordinates before topology build.")
    content.add_argument("--merge", action="store_true",
                         help="Merge all chains into single moleculetype")

    mode = parser.add_argument_group("Pipeline mode")
    mode.add_argument("--acpype", action="store_true",
                      help="Use ACPYPE pipeline (AMBER14+GLYCAM -> ParmEd -> GROMACS). "
                           "Handles mixed 1-4 scaling via [ pairs_nb ].")

    diag = parser.add_argument_group("Diagnostics")
    diag.add_argument("-v", "--verbose", action="store_true",
                      help="Verbose output")

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(parser)
    return parser.parse_args(argv)
