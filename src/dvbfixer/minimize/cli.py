"""CLI parsing for ``dvbfixer minimize``.

Split out of the flat ``minimize.py`` in Phase 2.1 of the revision plan.
Argparse only — no side effects on import. The concrete pipeline lives
in :mod:`dvbfixer.minimize.pipeline`.
"""

from __future__ import annotations

import argparse

DEFAULT_FF = "auto"
DEFAULT_PH = 7.0
DEFAULT_PADDING = 1.0  # nm
DEFAULT_RESTRAINT_K = 100.0  # kcal/mol/A^2 for original atoms
DEFAULT_WEAK_K = 5.0  # kcal/mol/A^2 for added atoms (backbone only)
DEFAULT_MAX_ITER = 1000

BACKBONE_NAMES = {"N", "CA", "C", "O", "CB"}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="dvbfixer minimize",
        description="Energy-minimize a PDB structure with OpenMM. Uses selective "
        "restraints: original atoms are restrained, newly added atoms (from "
        "PDBFixer .dat file) are free to relax.",
    )
    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB file")
    io.add_argument("-o", "--output", help="Output minimized PDB (default: <input>_minimized.pdb)")
    io.add_argument("--dat", help="Restraint data file from 'dvbfixer prepare' (default: <input>.dat)")

    ff = p.add_argument_group("Force field")
    ff.add_argument("--ph", type=float, default=DEFAULT_PH,
                    help=f"pH for hydrogen addition if needed (default: {DEFAULT_PH})")
    ff.add_argument("--ff", nargs="+", default=[DEFAULT_FF],
                    help="Force field selection. Accepts a short name "
                         "(auto, amber, amber+glycam, charmm, ...) or an "
                         "explicit list of OpenMM XML paths. Default: 'auto' — "
                         "detect from residue names in the input. "
                         "See docs/force-fields.md.")
    ff.add_argument("--parametrize-ligands", action="store_true",
                    help="For each heterogen residue that lacks a template in "
                         "the resolved --ff, run GAFF2 + AM1-BCC parametrisation "
                         "via antechamber/parmchk2 and register the resulting "
                         "GAFF template with OpenMM before createSystem. "
                         "Requires openmmforcefields + openff-toolkit + "
                         "AmberTools (antechamber, parmchk2). Cached to "
                         "~/.cache/dvbfixer/lig_params/ (override with "
                         "$DVBFIXER_LIG_CACHE). See docs/force-fields.md.")
    ff.add_argument("--atom-naming", choices=["gromacs", "standard"],
                    default="gromacs",
                    help="Atom-naming convention for the output PDB. "
                         "'gromacs' (default): rewrite atom names to "
                         "GROMACS amber99sb-ildn conventions (HB2/HB3 → "
                         "HB1/HB2, HZ3 → HZ1 on LYN, O → OC2, OXT → OC1, "
                         "H → HN on CHARMM). Matches `pdb2gmx -ff "
                         "amber99sb-ildn` expectations. 'standard': keep "
                         "IUPAC/AMBER-native names (HB2/HB3, HZ1/HZ2/HZ3, "
                         "O/OXT, plain H) — matches ff14SB / most PDB "
                         "downloaders / VMD.")

    physics = p.add_argument_group("Physics / restraints")
    physics.add_argument("--padding", type=float, default=DEFAULT_PADDING,
                         help=f"Solvent padding in nm (default: {DEFAULT_PADDING})")
    physics.add_argument("--no-solvent", action="store_true",
                         help="Minimize in vacuum (no solvent box)")
    physics.add_argument("--restraint-k", type=float, default=DEFAULT_RESTRAINT_K,
                         help=f"Restraint force constant for original atoms in kcal/mol/A^2 (default: {DEFAULT_RESTRAINT_K})")
    physics.add_argument("--weak-k", type=float, default=DEFAULT_WEAK_K,
                         help=f"Restraint force constant for added backbone atoms in kcal/mol/A^2 (default: {DEFAULT_WEAK_K})")
    physics.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER,
                         help=f"Max minimization iterations per phase (default: {DEFAULT_MAX_ITER})")

    content = p.add_argument_group("Content selection")
    content.add_argument("--strip-heterogens", dest="keep_heterogens",
                         action="store_false", default=True,
                         help="Strip heterogens before minimization, restore coords after "
                              "(protein-only mode). Default: minimize the whole system.")
    content.add_argument("--rebuild-h", action="store_true",
                         help="Strip and re-add hydrogens via OpenMM (default: keep existing)")
    content.add_argument("--rename", action="store_true",
                         help="Rename non-canonical residues (AMBER/CHARMM) to standard names before processing")
    content.add_argument("--no-infer-conect", dest="no_infer_conect",
                         action="store_true",
                         help="Skip automatic CONECT inference (default: infer missing "
                              "SS / glycosidic / glycosylation bonds before minimize).")

    refine = p.add_argument_group("Refinement (post-OpenMM)")
    refine.add_argument("--xtb-refine", action="store_true",
                        help="After OpenMM minimization, run xtb GFN-FF universal "
                             "force field as a refinement pass. Auto-parametrizes any "
                             "organic molecule (sugars, ligands) without templates. "
                             "Requires `xtb` binary in PATH. Slower but higher quality.")
    refine.add_argument("--xtb-cycles", type=int, default=200,
                        help="Max xtb optimization cycles (default: 200)")
    refine.add_argument("--obminimize-refine", action="store_true",
                        help="After OpenMM minimization, run OpenBabel obminimize "
                             "(MMFF94) as a refinement pass. Auto-typing for any "
                             "organic molecule. Faster than xtb, lower quality.")
    refine.add_argument("--obminimize-ff", default="UFF",
                        choices=["MMFF94", "MMFF94s", "UFF", "GAFF", "Ghemical"],
                        help="OpenBabel force field for --obminimize-refine (default: UFF — "
                             "handles N-glycosidic linkages correctly; MMFF94s mistypes the "
                             "anomeric C as sp2 giving 120° angles instead of 109°)")
    refine.add_argument("--obminimize-steps", type=int, default=500,
                        help="OpenBabel minimization steps (default: 500)")
    refine.add_argument("--refine-heterogens-only", action="store_true",
                        help="Restrict xtb/obminimize refinement to heterogen "
                             "residues (protein heavy atoms frozen). Refines only "
                             "the ligand's INTERNAL geometry — the protein-ligand "
                             "INTERFACE (contacts, H-bonds) is NOT relaxed and "
                             "any pre-existing clash there will persist. Use "
                             "without this flag for whole-system refinement when "
                             "the interface matters (whole-system xtb auto-switches "
                             "to heterogens-only above ~5000 atoms for performance).")

    runtime = p.add_argument_group("Runtime")
    runtime.add_argument("--platform", choices=["CPU", "CUDA", "OpenCL", "Reference"],
                         help="OpenMM platform (default: auto-select fastest)")
    runtime.add_argument("-v", "--verbose", action="store_true",
                         help="Print detailed progress")

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(p, batch=True)
    return p.parse_args(argv)
