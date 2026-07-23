"""CLI parsing + shared constants for ``dvbfixer model``.

Split out of the flat ``model.py`` in Phase 2.2 of the revision plan.
Argparse only — no OpenMM / Modeller imports at module scope. Shared
tables (``AA3TO1``, ``WATER_RESNAMES``) live here so both the numbering
helpers and the pipeline can import them without a circular dep.
"""

from __future__ import annotations

import argparse

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Nonstandard / protonation variants
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
    "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K",
}

WATER_RESNAMES = {"HOH", "WAT", "TIP3", "TIP", "SOL", "T3P", "T4P", "T5P"}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="dvbfixer model",
        description="Rebuild missing loops and gaps in a PDB structure using Modeller. "
        "Identifies gaps from SEQRES vs ATOM records (or a provided FASTA), "
        "then uses Modeller's loop modeling to fill them."
    )
    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB file (must contain SEQRES or use --fasta)")
    io.add_argument("-o", "--output", help="Output PDB file (default: <input>_model.pdb)")
    io.add_argument(
        "--fasta", help="FASTA file with complete sequence(s). Headers must encode "
        "chain IDs: '>chain_X', '>PDBID_X', or '>X'. Mapping is by chain ID, "
        "not file order. Use instead of SEQRES."
    )

    modelling = p.add_argument_group("Modelling parameters")
    modelling.add_argument(
        "-n", "--num-models", type=int, default=1,
        help="Number of initial models to generate (default: 1)"
    )
    modelling.add_argument(
        "--num-loops", type=int, default=2,
        help="Number of loop refinement models per initial model (default: 2)"
    )
    modelling.add_argument(
        "--num-output", type=int, default=1, dest="num_output",
        help="Number of top-ranked candidate models to save (default: 1; "
             "ceiling: num_models × num_loops). Output PDBs are sorted "
             "ascending by Modeller's molpdf score (best first). With "
             "--num-output > 1, output filenames get a _N suffix: "
             "<stem>_model_1.pdb, _model_2.pdb, ... (and matching .dat). "
             "With --num-output 1 the filename is unchanged from today."
    )
    modelling.add_argument(
        "--md-level", choices=["none", "fast", "slow", "very_slow", "slow_large"],
        default="fast",
        help="MD refinement level for loop modeling (default: fast)"
    )
    modelling.add_argument(
        "--pin-input", dest="pin_input",
        action=argparse.BooleanOptionalAction, default=True,
        help="During Modeller's loop refinement MD, allow only the "
             "newly-modeled gap residues to move — no ±flank margin. "
             "Default ON: prevents flanking-residue drift; the initial "
             "automodel CG still runs on all atoms so any input close "
             "contacts get relaxed normally. Pass --no-pin-input for the "
             "legacy LoopModel behaviour (gap ±~3 residue flank mobile)."
    )
    modelling.add_argument(
        "--no-terminal", action="store_true",
        help="Do not model missing N/C terminal residues (only rebuild internal gaps)"
    )

    content = p.add_argument_group("Content selection")
    content.add_argument(
        "--keep-water", action="store_true",
        help="Keep water molecules (HOH, WAT, TIP3, SOL) in output (default: remove)"
    )
    content.add_argument(
        "--strip-heterogens", dest="keep_heterogens",
        action="store_false", default=True,
        help="Remove all HETATM records (ligands, sugars, ions, cofactors) "
             "before Modeller runs. Waters are also removed unless "
             "--keep-water is passed. Off by default — Modeller usually "
             "benefits from heterogen context for loop refinement, but "
             "in some cases (bad ligand geometry, ambiguous CONECT) they "
             "cause artifacts.",
    )

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument(
        "--keep-workdir", action="store_true",
        help="Keep the Modeller working directory (for debugging)"
    )
    diag.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print Modeller progress"
    )

    return p.parse_args(argv)
