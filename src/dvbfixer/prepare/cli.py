"""CLI parsing + shared constants for ``dvbfixer prepare``.

Split out of the flat ``prepare.py`` in Phase 2.3 of the revision plan.
Argparse only — safe to import without OpenMM / PDBFixer. Shared
constants (``DEFAULT_PH``, ``GLYCOSYLATED_RESIDUES``, ``SUGAR_RESNAMES``)
live here so the pipeline and its future extracted submodules can pick
them up without a circular dep.
"""

from __future__ import annotations

import argparse

from dvbfixer.residue_registry import CHARMM_SUGAR_RESNAMES, PDB_SUGAR_RESNAMES

DEFAULT_PH = 7.0

# Residues that form glycosidic bonds through their sidechain donor
# (ASN via ND2, SER via OG, THR via OG1).
GLYCOSYLATED_RESIDUES = {"ASN", "SER", "THR"}

# Known sugar residue names across the three force fields:
#   - PDB-style 3-char codes from the RCSB Chemical Component Dictionary
#   - CHARMM-GUI 4-char codes (BGLC, AMAN, BGAL, BGLCNA, ...)
#   - GLYCAM 3-char codes are detected separately via is_glycam_sugar()
SUGAR_RESNAMES = set(PDB_SUGAR_RESNAMES | CHARMM_SUGAR_RESNAMES)


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
    ff.add_argument("--atom-naming", choices=["gromacs", "standard"],
                    default="gromacs",
                    help="Atom-naming convention for the output PDB. "
                         "'gromacs' (default): GROMACS amber99sb-ildn "
                         "shifts (HB3→HB1 keeping HB2, HZ3→HZ1 on LYN, "
                         "O/OXT→OC2/OC1, H→HN for CHARMM). 'standard': "
                         "IUPAC/AMBER-native names (HB2/HB3, HZ1/HZ2/HZ3, "
                         "O/OXT, plain H).")
    ff.add_argument("--propka", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="Run PROPKA3 for pKa-driven AMBER variant "
                         "renames (ASH/GLH/HIP/CYM/LYN/CYX). Default ON. "
                         "Pass --no-propka to skip; variants then come "
                         "from --mutate + input HD1/HE2 atoms only "
                         "(0.7.5/0.7.6 behaviour).")
    ff.add_argument("--protassign", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="Run MolProbity Reduce for HIS tautomer "
                         "(HID vs HIE) + ASN/GLN flip detection. "
                         "Default ON.")
    ff.add_argument("--his-default", choices=["HIE", "HID"], default="HIE",
                    help="Default HIS tautomer when PROPKA says neutral "
                         "AND Reduce didn't place either HD1 or HE2 "
                         "(rare — deprotonated HIS). Default: HIE.")
    ff.add_argument("--cys-ss-pka", type=float, default=99.99,
                    help="PROPKA pKa threshold above which CYS is "
                         "assumed to be in a disulfide bond and renamed "
                         "to CYX (default: 99.99, matching PROPKA's sentinel). "
                         "Explicit CONECT-detected SS pairs override PROPKA regardless.")
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
    content.add_argument(
        "--smiles", action="append", default=[], metavar="RESNAME=SMILES",
        help="Use SMILES chemistry when adding H to an isolated small-molecule "
             "residue (for example --smiles 'LIG=CC(=O)[O-]'). Applies to "
             "every matching residue and may be repeated. Optional; unmapped "
             "heterogens keep the existing automatic RDKit/OpenBabel path.",
    )
    content.add_argument("--rename", action="store_true",
                         help="Rename non-canonical residues (AMBER/CHARMM) to standard names before processing")
    content.add_argument("--no-infer-conect", dest="no_infer_conect",
                         action="store_true",
                         help="Skip automatic CONECT inference. By default dvbfixer "
                              "perceives missing CONECT records (SS, glycosidic, "
                              "glycosylation) so glycoprotein flows work even on "
                              "inputs without CONECT.")
    content.add_argument(
        "--cap-termini", action="store_true",
        help="Add a neutral ACE N-cap and NME C-cap to protein chains. "
             "By default every protein chain is capped; use --cap-chain "
             "to restrict the selection.",
    )
    content.add_argument(
        "--cap-chain", action="append", default=[], metavar="CHAIN",
        help="Protein chain to cap (repeatable). Use '_' for a blank chain ID. "
             "Only meaningful with --cap-termini.",
    )

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

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(p, batch=True)
    return p.parse_args(argv)
