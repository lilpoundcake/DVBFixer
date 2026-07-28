"""CLI parsing + shared constants for ``dvbfixer prepare``.

Split out of the flat ``prepare.py`` in Phase 2.3 of the revision plan.
Argparse only — safe to import without OpenMM / PDBFixer. Shared
constants (``DEFAULT_PH``, ``GLYCOSYLATED_RESIDUES``, ``SUGAR_RESNAMES``)
live here so the pipeline and its future extracted submodules can pick
them up without a circular dep.
"""

from __future__ import annotations

import argparse

DEFAULT_PH = 7.0

# Residues that form glycosidic bonds through their sidechain donor
# (ASN via ND2, SER via OG, THR via OG1).
GLYCOSYLATED_RESIDUES = {"ASN", "SER", "THR"}

# Known sugar residue names across the three force fields:
#   - PDB-style 3-char codes from the RCSB Chemical Component Dictionary
#   - CHARMM-GUI 4-char codes (BGLC, AMAN, BGAL, BGLCNA, ...)
#   - GLYCAM 3-char codes are detected separately via is_glycam_sugar()
SUGAR_RESNAMES = {
    # PDB 3-char
    "NAG", "NDG", "BMA", "MAN", "FUC", "FUL", "GAL", "BGC", "GLC", "SIA",
    "NGA", "A2G", "AFU", "AMA", "BGA", "BGL", "XYS", "XYP", "RIB", "GCU",
    "IDS", "RAM", "NAN",
    # CHARMM-GUI 4-char
    "BGLC", "AGLC", "BMAN", "AMAN", "BGAL", "AGAL", "BGLCNA", "BGALNA",
    "AFUC", "BFUC", "ANE5AC", "BNE5AC", "BXYL", "AXYL", "BIDOA", "BGLCA",
    "AGLCA",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="dvbfixer prepare",
        description="Fix missing atoms and residues in a PDB structure using PDBFixer. "
        "Writes a .dat file recording added atoms for selective restraints "
        "during minimization.",
    )
    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB file")
    io.add_argument("-o", "--output", help="Output PDB file (default: <input>_prepared.pdb)")
    io.add_argument("--dat", help="Restraint data file path (default: <output>.dat)")

    ff = p.add_argument_group("Force field / pH")
    ff.add_argument("--backend", choices=["tleap-reduce", "legacy"],
                    default="legacy",
                    help="Prep backend. 'legacy' (default): PDBFixer + "
                         "Modeller.addHydrogens; handles glycans, ligands, "
                         "heterogens and covalent-HETATM links. The chirality "
                         "invariant is enforced by minimize's post-phase-2 "
                         "unconditional force-reflect (0.7.4+), so legacy "
                         "prep's D-Cα risk is neutralised downstream. "
                         "'tleap-reduce': opt-in deterministic AmberTools + "
                         "MolProbity pipeline (tleap for heavy atoms, reduce "
                         "for H). Pure-protein only — rejects non-canonical "
                         "residues. Use when you specifically want L-only "
                         "heavy atoms produced by tleap itself.")
    ff.add_argument("--ph", type=float, default=DEFAULT_PH,
                    help=f"pH for adding hydrogens (default: {DEFAULT_PH})")
    ff.add_argument("--ff", nargs="+", default=["auto"],
                    help="Force field selection for heterogen-H addition. "
                         "Accepts a short name (auto, amber, amber+glycam, "
                         "charmm, ...) or an explicit list of OpenMM XML "
                         "paths. Default: 'auto' — detect from residue names "
                         "in the input. See docs/force-fields.md. Only "
                         "consulted when heterogen-H addition runs; the "
                         "protein-only PDBFixer path is unaffected.")

    content = p.add_argument_group("Content selection")
    content.add_argument("--keep-water", action="store_true",
                         help="Keep crystallographic waters")
    content.add_argument("--strip-heterogens", dest="keep_heterogens",
                         action="store_false", default=True,
                         help="Remove heterogens (sugars, ligands, ions) before processing "
                              "(protein-only mode). Default: keep heterogens.")
    content.add_argument("--no-heterogen-h", dest="heterogen_h",
                         action="store_false", default=True,
                         help="Skip hydrogen addition for heterogens (sugars/ligands).")
    content.add_argument("--rename", action="store_true",
                         help="Rename non-canonical residues (AMBER/CHARMM) to standard names before processing")
    content.add_argument("--no-infer-conect", dest="no_infer_conect",
                         action="store_true",
                         help="Skip automatic CONECT inference. By default dvbfixer "
                              "perceives missing CONECT records (SS, glycosidic, "
                              "glycosylation) so glycoprotein flows work even on "
                              "inputs without CONECT.")

    mutations = p.add_argument_group("Mutations")
    mutations.add_argument("--mutate", action="append", default=[],
                           metavar="CHAIN:RESNUM:NEW_AA",
                           help="Mutate a residue (e.g. A:39:ALA, A:83:HIP). Use "
                                "CHAIN:RESNUM:del to DELETE a residue (e.g. H:446:del). "
                                "Insertion codes are supported in RESNUM (e.g. H:100A:del). "
                                "Can be used multiple times.")

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument("-v", "--verbose", action="store_true",
                      help="Print detailed progress")

    return p.parse_args(argv)
