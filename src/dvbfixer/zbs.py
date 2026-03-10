"""Full pipeline: renumber -> model -> prepare -> minimize -> protonate -> minimize.

Runs the complete PDB preparation workflow in sequence, passing output
of each step as input to the next. Two minimize passes: first with standard
protonation to get good heavy-atom positions, then protonate renames residues
to AMBER names (HIE/HID/HIP, ASH, GLH, CYM, CYX, LYN), then second minimize
re-adds hydrogens matching the correct protonation states.
"""

import argparse
import sys
from pathlib import Path


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer zbs",
        description="Run the full preparation pipeline: "
        "renumber -> model -> prepare -> minimize -> protonate -> minimize.",
    )
    p.add_argument("input", help="Input PDB file (must contain SEQRES)")
    p.add_argument("-o", "--output", help="Final output PDB file (default: <input>_zbs.pdb)")
    p.add_argument("--ph", type=float, default=7.0,
                   help="pH for protonation and hydrogen addition (default: 7.0)")
    p.add_argument("--ff", nargs='+',
                   default=['amber19/protein.ff19SB.xml', 'amber19/tip3p.xml'],
                   help="Force field XML files for minimization")

    # renumber
    p.add_argument("--skip-renumber", action="store_true",
                   help="Skip the renumber step")

    # model
    p.add_argument("--skip-model", action="store_true",
                   help="Skip the model step")
    p.add_argument("--no-terminal", action="store_true",
                   help="Do not model missing N/C terminal residues")
    p.add_argument("--num-loops", type=int, default=2,
                   help="Number of loop models (default: 2)")
    p.add_argument("--md-level", choices=["none", "fast", "slow", "very_slow", "slow_large"],
                   default="fast", help="Modeller MD refinement level (default: fast)")
    p.add_argument("--fasta", help="FASTA file with complete sequence(s) for model step")

    # prepare
    p.add_argument("--skip-prepare", action="store_true",
                   help="Skip the prepare step")
    p.add_argument("--keep-heterogens", action="store_true",
                   help="Keep all heterogens during prepare step")
    p.add_argument("--mutate", action="append", default=[],
                   metavar="CHAIN:RESNUM:NEW_AA",
                   help="Mutate a residue during prepare step (can be used multiple times)")

    # minimize
    p.add_argument("--skip-minimize", action="store_true",
                   help="Skip the minimize step")
    p.add_argument("--no-solvent", action="store_true",
                   help="Minimize in vacuum (no solvent box)")
    p.add_argument("--rebuild-h", action="store_true",
                   help="Strip and re-add hydrogens via OpenMM during minimization (default: keep existing)")
    p.add_argument("--restraint-k", type=float, default=100.0,
                   help="Restraint force constant for original atoms (default: 100)")
    p.add_argument("--max-iter", type=int, default=1000,
                   help="Max minimization iterations per phase (default: 1000)")
    p.add_argument("--platform", choices=["CPU", "CUDA", "OpenCL", "Reference"],
                   help="OpenMM platform (default: auto)")

    # protonate
    p.add_argument("--skip-protonate", action="store_true",
                   help="Skip the protonate step")
    p.add_argument("--no-hydrogens", action="store_true",
                   help="Only rename residues in protonate step, do not add/fix H atoms")
    p.add_argument("--keep-water", action="store_true",
                   help="Keep water molecules in output (default: remove)")

    p.add_argument("--keep-interim", action="store_true",
                   help="Keep all intermediate files (default: only final output)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print detailed progress for all steps")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    if args.output:
        final_output = Path(args.output)
    else:
        final_output = input_path.with_stem(input_path.stem + "_zbs")

    current = str(input_path)
    step_num = 0
    interim_files = []  # files to clean up unless --keep-interim

    def step_output(name):
        path = str(input_path.with_stem(input_path.stem + f"_{name}"))
        interim_files.append(path)
        # Also track .dat files that steps may produce
        interim_files.append(str(Path(path).with_suffix('.dat')))
        return path

    # 1. Renumber
    if not args.skip_renumber:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: RENUMBER")
        print(f"{'='*60}")
        from dvbfixer.renumber import main as renumber_main
        out = step_output("renum")
        renumber_argv = [current, "-o", out]
        if args.keep_water:
            renumber_argv.append("--keep-water")
        if args.verbose:
            renumber_argv.append("-v")
        renumber_main(renumber_argv)
        current = out

    # 2. Model
    if not args.skip_model:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: MODEL")
        print(f"{'='*60}")
        from dvbfixer.model import main as model_main
        out = step_output("model")
        model_argv = [current, "-o", out,
                      "--num-loops", str(args.num_loops),
                      "--md-level", args.md_level]
        if args.no_terminal:
            model_argv.append("--no-terminal")
        if args.fasta:
            model_argv.extend(["--fasta", args.fasta])
        if args.verbose:
            model_argv.append("-v")
        model_main(model_argv)
        current = out

    # 3. Prepare
    if not args.skip_prepare:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: PREPARE")
        print(f"{'='*60}")
        from dvbfixer.prepare import main as prepare_main
        out = step_output("prepared")
        prepare_argv = [current, "-o", out, "--ph", str(args.ph)]
        if args.keep_heterogens:
            prepare_argv.append("--keep-heterogens")
        for mut in args.mutate:
            prepare_argv.extend(["--mutate", mut])
        if args.verbose:
            prepare_argv.append("-v")
        prepare_main(prepare_argv)
        current = out

    # 4. Minimize (first pass — standard protonation, good heavy-atom positions)
    if not args.skip_minimize:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: MINIMIZE (pass 1)")
        print(f"{'='*60}")
        from dvbfixer.minimize import main as minimize_main
        out = step_output("minimized")
        minimize_argv = [current, "-o", out,
                         "--ph", str(args.ph),
                         "--ff"] + args.ff + [
                         "--restraint-k", str(args.restraint_k),
                         "--max-iter", str(args.max_iter)]
        if args.no_solvent:
            minimize_argv.append("--no-solvent")
        if args.rebuild_h:
            minimize_argv.append("--rebuild-h")
        if args.platform:
            minimize_argv.extend(["--platform", args.platform])
        if args.verbose:
            minimize_argv.append("-v")
        minimize_main(minimize_argv)
        current = out

    # 5. Protonate (rename only — text-based, no hydrogen addition)
    if not args.skip_protonate:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: PROTONATE (rename residues)")
        print(f"{'='*60}")
        from dvbfixer.protonate import main as protonate_main
        out = step_output("prot")
        protonate_argv = [current, "-o", out,
                          "--ph", str(args.ph),
                          "--no-hydrogens"]  # always rename-only in pipeline
        if args.verbose:
            protonate_argv.extend(["-v", "--summary"])
        protonate_main(protonate_argv)
        current = out

    # 6. Minimize (pass 2 — re-adds H with correct AMBER protonation templates)
    if not args.skip_minimize and not args.skip_protonate:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: MINIMIZE (pass 2 — correct protonation)")
        print(f"{'='*60}")
        out = step_output("minimized2")
        minimize_argv = [current, "-o", out,
                         "--ph", str(args.ph),
                         "--ff"] + args.ff + [
                         "--restraint-k", str(args.restraint_k),
                         "--max-iter", str(args.max_iter),
                         "--rebuild-h"]  # must rebuild H with correct AMBER protonation
        if args.no_solvent:
            minimize_argv.append("--no-solvent")
        if args.platform:
            minimize_argv.extend(["--platform", args.platform])
        if args.verbose:
            minimize_argv.append("-v")
        minimize_main(minimize_argv)
        current = out

    # 7. Re-apply AMBER residue names (minimize reverts them to standard names)
    if not args.skip_protonate:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: PROTONATE (re-apply AMBER names)")
        print(f"{'='*60}")
        out = step_output("final")
        protonate_argv = [current, "-o", out,
                          "--ph", str(args.ph),
                          "--no-hydrogens"]
        if args.verbose:
            protonate_argv.append("-v")
        protonate_main(protonate_argv)
        current = out

    # Copy last output to final destination
    if str(current) != str(final_output):
        import shutil
        shutil.copy2(current, final_output)

    # Clean up interim files
    if not args.keep_interim:
        removed = 0
        for f in interim_files:
            p = Path(f)
            if p.exists() and p.resolve() != final_output.resolve():
                p.unlink()
                removed += 1
        if removed and args.verbose:
            print(f"Cleaned up {removed} interim file(s)")

    print(f"\n{'='*60}")
    print(f"Pipeline complete: {final_output}")
    print(f"{'='*60}")
