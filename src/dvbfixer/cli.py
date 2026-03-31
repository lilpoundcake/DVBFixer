"""dvbfixer — unified CLI for PDB structure preparation tools."""

import sys


COMMANDS = {
    "split": "Split chains empirically (distance + numbering)",
    "renumber": "Renumber residues using SEQRES alignment",
    "model": "Rebuild missing loops/gaps with Modeller",
    "pull": "Pull atoms together to form a bond (geometry-only)",
    "prepare": "Fix missing atoms/residues with PDBFixer",
    "minimize": "Energy-minimize with OpenMM using selective restraints",
    "protonate": "Set protonation states using PROPKA3 pKa predictions",
    "rename": "Rename non-canonical residues (AMBER/CHARMM) to standard names",
    "top": "Generate GROMACS .itp/.top topology files from PDB",
    "transplant": "Transplant molecules from donor PDB to acceptor PDB",
    "puppet": "Strip PDB to backbone-only polyglycine model",
    "glycam": "Convert glycan PDB to GLYCAM force field nomenclature",
    "cluster": "Cluster glycan conformations from MD trajectory",
    "parametrize": "Parametrize small molecules with GAFF2 + AM1-BCC/RESP",
    "homology": "Multi-template homology modeling with Modeller",
    "zbs": "Full pipeline: renumber -> model -> prepare -> minimize -> protonate",
}


def print_help():
    print("dvbfixer — PDB structure preparation tools\n")
    print("Usage: dvbfixer <command> [options]\n")
    print("Commands:")
    for cmd, desc in COMMANDS.items():
        print(f"  {cmd:<12s}  {desc}")
    print(f"\n  --version     Show version")
    print(f"\nRun 'dvbfixer <command> --help' for command-specific options.")


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print_help()
        sys.exit(0)

    if sys.argv[1] == "--version":
        from dvbfixer import __version__
        print(f"dvbfixer {__version__}")
        sys.exit(0)

    command = sys.argv[1]
    argv = sys.argv[2:]

    if command == "split":
        from dvbfixer.split_chains import main as cmd_main
    elif command == "renumber":
        from dvbfixer.renumber import main as cmd_main
    elif command == "model":
        from dvbfixer.model import main as cmd_main
    elif command == "pull":
        from dvbfixer.pull import main as cmd_main
    elif command == "rename":
        from dvbfixer.rename import main as cmd_main
    elif command == "top":
        from dvbfixer.top import main as cmd_main
    elif command == "prepare":
        from dvbfixer.prepare import main as cmd_main
    elif command == "minimize":
        from dvbfixer.minimize import main as cmd_main
    elif command == "protonate":
        from dvbfixer.protonate import main as cmd_main
    elif command == "transplant":
        from dvbfixer.transplant import main as cmd_main
    elif command == "puppet":
        from dvbfixer.puppet import main as cmd_main
    elif command == "glycam":
        from dvbfixer.glycam import main as cmd_main
    elif command == "cluster":
        from dvbfixer.cluster import main as cmd_main
    elif command == "parametrize":
        from dvbfixer.parametrize import main as cmd_main
    elif command == "homology":
        from dvbfixer.homology import main as cmd_main
    elif command == "zbs":
        from dvbfixer.zbs import main as cmd_main
    else:
        print(f"Unknown command: {command}\n", file=sys.stderr)
        print_help()
        sys.exit(1)

    cmd_main(argv)


if __name__ == "__main__":
    main()
