"""dvbfixer — unified CLI for PDB structure preparation tools."""

from __future__ import annotations

import sys
from importlib import import_module

COMMANDS: dict[str, str] = {
    "split": "Split chains empirically or extract PDB biological assemblies",
    "renumber": "Renumber residues using FASTA or SEQRES alignment",
    "model": "Rebuild missing loops/gaps with Modeller",
    "pull": "Pull atoms together to form a bond (geometry-only)",
    "prepare": "Fix missing atoms/residues with PDBFixer",
    "minimize": "Energy-minimize with OpenMM using selective restraints",
    "protonate": "Set protonation states using PROPKA3 pKa predictions",
    "rename": "Rename non-canonical residues (AMBER/CHARMM) to standard names",
    "top": "Generate GROMACS .itp/.top topology files from PDB",
    "transplant": "Transplant molecules from donor PDB to acceptor PDB",
    "puppet": "Strip PDB to backbone-only polyglycine model",
    "convert": "Convert between PDB/AMBER/GLYCAM and CHARMM naming (sugars + protonation variants)",
    "conect": "Add inferred CONECT records (SS bonds, glycosidic links, glycosylation)",
    "cluster": "Cluster glycan conformations from MD trajectory",
    "parametrize": "Parametrize small molecules with GAFF2 + AM1-BCC/RESP",
    "homology": "Multi-template homology modeling with Modeller",
    "diagnose": "Report structure-quality issues (missing atoms, clashes, valence, ...)",
    "doctor": "Report installed backends, executables, and OpenMM platforms",
    "zbs": "Full pipeline: renumber -> model -> prepare -> minimize",
}


def print_help() -> None:
    print("dvbfixer — PDB structure preparation tools\n")
    print("Usage: dvbfixer <command> [options]\n")
    print("Commands:")
    for cmd, desc in COMMANDS.items():
        print(f"  {cmd:<12s}  {desc}")
    print("\n  --version     Show version")
    print("\nBatch mode (runs a command independently on each structure;")
    print("continues after per-file failures by default):")
    print("  --input-dir DIR      Process every .pdb/.ent structure in DIR")
    print("  --output-dir DIR     Batch output directory")
    print("  --recursive          Include subdirectories")
    print("  --fail-fast          Stop after the first failed structure")
    print("\nRun 'dvbfixer <command> --help' for command-specific options.")


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print_help()
        sys.exit(0)

    if sys.argv[1] == "--version":
        from dvbfixer import __version__
        print(f"dvbfixer {__version__}")
        sys.exit(0)

    command = sys.argv[1]
    argv = sys.argv[2:]

    from dvbfixer.batch import extract_batch_options

    batch_options, argv = extract_batch_options(argv)

    if command == "glycam":
        # `glycam` is the legacy name; `convert` is preferred. The module
        # filename stays as `glycam.py` for now to keep imports stable.
        print("[deprecated] 'dvbfixer glycam' is now 'dvbfixer convert'. "
              "The old name still works but please update scripts.",
              file=sys.stderr)
    if command not in COMMANDS and command != "glycam":
        print(f"Unknown command: {command}\n", file=sys.stderr)
        print_help()
        sys.exit(1)

    module_name = {
        "split": "split_chains",
        "convert": "glycam",
        "glycam": "glycam",
    }.get(command, command)
    cmd_main = import_module(f"dvbfixer.{module_name}").main

    if batch_options.input_dir:
        from dvbfixer.batch import run_directory

        run_directory(command, cmd_main, batch_options, argv)
    else:
        if batch_options.output_dir or batch_options.recursive or batch_options.fail_fast:
            raise SystemExit("--output-dir, --recursive, and --fail-fast require --input-dir")
        cmd_main(argv)


if __name__ == "__main__":
    main()
