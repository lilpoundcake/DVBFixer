"""Full pipeline: renumber -> model -> prepare -> minimize -> protonate -> minimize.

Runs the complete PDB preparation workflow in sequence, passing output
of each step as input to the next. Two minimize passes:

  pass 1: relax heavy atoms using prepare's approximate hydrogens
  pass 2: refine hydrogens placed by protonate's variant-aware addHydrogens

Between the two minimize passes, `protonate` runs the FULL pipeline
(PROPKA + MolProbity Reduce for HIS tautomers + ASN/GLN flip detection +
variant-aware OpenMM addHydrogens). No `--no-hydrogens` flag is passed to
protonate — that flag leaves existing H atoms in wrong positions after the
variant renames.

Pass 2 minimize does NOT force `--rebuild-h`: hydrogens placed by
protonate are already correct via the variant-aware addHydrogens; the
minimize just relaxes them. minimize.py preserves AMBER variant names on
output via its `_input_variants` capture-restore path, so no final "re-apply
names" step is needed.
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
    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB file (must contain SEQRES)")
    io.add_argument("-o", "--output", help="Final output PDB file (default: <input>_zbs.pdb)")

    ff = p.add_argument_group("Force field")
    ff.add_argument("--ph", type=float, default=7.0,
                    help="pH for protonation and hydrogen addition (default: 7.0)")
    ff.add_argument("--ff", nargs='+', default=['auto'],
                    help="Force field selection forwarded to prepare / "
                         "minimize / protonate. Accepts a short name (auto, "
                         "amber, amber+glycam, charmm, ...) or an explicit "
                         "list of OpenMM XML paths. Default: 'auto'. "
                         "See docs/force-fields.md.")
    ff.add_argument("--parametrize-ligands", action="store_true",
                    help="Forward --parametrize-ligands to both minimize "
                         "passes (GAFF2 + AM1-BCC for unknown ligands via "
                         "antechamber). See docs/force-fields.md.")

    skip = p.add_argument_group("Pipeline skip flags")
    skip.add_argument("--skip-renumber", action="store_true",
                      help="Skip the renumber step")
    skip.add_argument("--skip-model", action="store_true",
                      help="Skip the model step")
    skip.add_argument("--skip-prepare", action="store_true",
                      help="Skip the prepare step")
    skip.add_argument("--skip-minimize", action="store_true",
                      help="Skip the minimize step")
    skip.add_argument("--skip-protonate", action="store_true",
                      help="Skip the protonate step")

    model_grp = p.add_argument_group("Model step (Modeller)")
    model_grp.add_argument("--fasta", help="FASTA file with complete sequence(s) for model step")
    model_grp.add_argument("--no-terminal", action="store_true",
                           help="Do not model missing N/C terminal residues")
    model_grp.add_argument("--num-loops", type=int, default=2,
                           help="Number of loop models (default: 2)")
    model_grp.add_argument("--md-level", choices=["none", "fast", "slow", "very_slow", "slow_large"],
                           default="fast", help="Modeller MD refinement level (default: fast)")
    model_grp.add_argument("--num-output", type=int, default=1, dest="num_output",
                           help="Save top-N candidate models from Modeller (default: 1). "
                                "With N>1, the model step writes <stem>_model_1.pdb, "
                                "..._2.pdb, ... zbs picks _1 (best) for downstream.")
    model_grp.add_argument("--pin-input", dest="pin_input",
                           action=argparse.BooleanOptionalAction, default=True,
                           help="During Modeller's loop refinement MD, allow only "
                                "the gap residues to move (no ±flank margin). "
                                "Default ON — pass --no-pin-input to restore the "
                                "legacy LoopModel behaviour (gap ±~3 residue flank "
                                "mobile).")

    prep = p.add_argument_group("Prepare step")
    prep.add_argument("--strip-heterogens", dest="keep_heterogens",
                      action="store_false", default=True,
                      help="Strip heterogens before processing (protein-only pipeline). "
                           "Default: keep heterogens through prepare and minimize the whole system.")
    prep.add_argument("--no-heterogen-h", dest="heterogen_h",
                      action="store_false", default=True,
                      help="Skip hydrogen addition for heterogens in prepare "
                           "(default: add H to heterogens BioLuminate-style).")
    prep.add_argument("--mutate", action="append", default=[],
                      metavar="CHAIN:RESNUM:NEW_AA",
                      help="Mutate a residue during prepare step (can be used multiple times)")
    prep.add_argument("--rename", action="store_true",
                      help="Canonicalise non-standard residue names before prepare/minimize.")

    minz = p.add_argument_group("Minimize step (OpenMM)")
    minz.add_argument("--no-solvent", action="store_true",
                      help="Minimize in vacuum (no solvent box)")
    minz.add_argument("--rebuild-h", action="store_true",
                      help="Force --rebuild-h on both minimize passes (default: off; "
                           "pass 2 already has correct H from protonate)")
    minz.add_argument("--restraint-k", type=float, default=100.0,
                      help="Restraint force constant for original atoms (default: 100)")
    minz.add_argument("--max-iter", type=int, default=1000,
                      help="Max minimization iterations per phase (default: 1000)")
    minz.add_argument("--refine", choices=["none", "xtb", "obminimize"], default="none",
                      help="Post-minimize refinement pass in minimize step 2 "
                           "(default: none). 'xtb' uses GFN-FF, 'obminimize' uses "
                           "OpenBabel UFF.")
    minz.add_argument("--refine-heterogens-only", action="store_true",
                      help="Restrict --refine pass to heterogen residues "
                           "(protein backbone frozen). Only meaningful with "
                           "--refine != none.")

    prot = p.add_argument_group("Protonate step")
    prot.add_argument("--no-propka", dest="propka",
                      action="store_false", default=True,
                      help="Skip PROPKA3 in the protonate step. Reduce "
                           "(--protassign) becomes the only source of HIS "
                           "tautomer picks and ASN/GLN flip detection. "
                           "Combining --no-propka with --no-protassign is an "
                           "error.")
    prot.add_argument("--no-protassign", dest="protassign",
                      action="store_false", default=True,
                      help="Skip MolProbity Reduce (HIS tautomer / ASN-GLN flip "
                           "detection) in protonate. Default: run Reduce, matches "
                           "standalone `dvbfixer protonate` default since Jun 2026.")

    general = p.add_argument_group("Pipeline behaviour")
    general.add_argument("--keep-water", action="store_true",
                         help="Keep water molecules in output (default: remove)")
    general.add_argument("--no-infer-conect", dest="no_infer_conect",
                         action="store_true",
                         help="Skip automatic CONECT inference in prepare/minimize/"
                              "protonate. Default: infer missing CONECT bonds "
                              "(SS/glycosidic/glycosylation) from coordinates.")
    general.add_argument("--keep-interim", action="store_true",
                         help="Keep all intermediate files (default: only final output)")
    general.add_argument("--dry-run", action="store_true",
                         help="Print the planned pipeline steps + output "
                              "filenames without running anything. Useful "
                              "when many skip flags are in play.")
    general.add_argument("--align-to-input", dest="align_to_input",
                         action=argparse.BooleanOptionalAction, default=True,
                         help="After every pipeline step, Kabsch-align the output "
                              "back to the ORIGINAL input on protein backbone "
                              "atoms. Prevents accumulated rigid-body drift so "
                              "residue-by-residue comparisons in a viewer line up. "
                              "Default ON — pass --no-align-to-input for the "
                              "legacy behaviour (each step's output in its own "
                              "frame).")

    runtime = p.add_argument_group("Runtime")
    runtime.add_argument("--platform", choices=["CPU", "CUDA", "OpenCL", "Reference"],
                         help="OpenMM platform (default: auto)")
    runtime.add_argument("-v", "--verbose", action="store_true",
                         help="Print detailed progress for all steps")

    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    try:
        _run_pipeline(args, input_path)
    except _lazy_import_chirality_error() as e:
        print(f"\nERROR: chirality guard tripped in the zbs pipeline.\n"
              f"{e}\n"
              f"Downstream MD would produce nonsense on a D-Cα residue; "
              f"fix the upstream model / prepare / protonate step or "
              f"rebuild the affected residue by hand.",
              file=sys.stderr)
        sys.exit(2)


def _lazy_import_chirality_error():
    """Deferred import — keeps ``dvbfixer zbs --help`` free of the
    ffutils import chain (which pulls OpenMM)."""
    from dvbfixer.ffutils.geometry import ChiralityError
    return ChiralityError


def _print_dry_run(args, input_path, final_output):
    """Print the planned pipeline steps + outputs without running anything."""
    print("Planned zbs pipeline:")
    step = 1

    def _line(name, out_suffix, notes=""):
        nonlocal step
        out = str(input_path.with_stem(input_path.stem + out_suffix))
        line = f"  {step}. {name:12s} → {Path(out).name}"
        if notes:
            line += f"   ({notes})"
        print(line)
        step += 1

    if not args.skip_renumber:
        _line("renumber", "_renum",
              "SEQRES-based residue renumbering"
              + (", keep water" if args.keep_water else ""))
    if not args.skip_model:
        notes = f"Modeller LoopModel, --md-level {args.md_level}"
        if args.fasta:
            notes += f", --fasta {args.fasta}"
        if not args.pin_input:
            notes += ", --no-pin-input"
        _line("model", "_model", notes)
    if not args.skip_prepare:
        notes = "backend: legacy (PDBFixer + Modeller.addHydrogens)"
        if not args.keep_heterogens:
            notes += ", --strip-heterogens"
        if args.mutate:
            notes += f", {len(args.mutate)} mutation(s)"
        _line("prepare", "_prepared", notes)
    if not args.skip_minimize:
        notes = "OpenMM ff14SB"
        if args.no_solvent:
            notes += ", --no-solvent"
        if args.rebuild_h:
            notes += ", --rebuild-h"
        if args.refine != "none":
            notes += f", --refine {args.refine}"
        _line("minimize", "_minimized", notes)

    print(f"  final       → {final_output.name}")
    if args.verbose:
        interim_kept = "keep" if args.keep_interim else "delete"
        print(f"  (interim files: {interim_kept})")


def _run_pipeline(args, input_path):
    if args.output:
        final_output = Path(args.output)
    else:
        final_output = input_path.with_stem(input_path.stem + "_zbs")

    if getattr(args, "dry_run", False):
        _print_dry_run(args, input_path, final_output)
        return

    current = str(input_path)
    step_num = 0
    interim_files = []  # files to clean up unless --keep-interim

    def step_output(name):
        path = str(input_path.with_stem(input_path.stem + f"_{name}"))
        interim_files.append(path)
        # Also track .dat files that steps may produce
        interim_files.append(str(Path(path).with_suffix('.dat')))
        return path

    def _maybe_align(path):
        """Kabsch-superpose `path` onto the original input in-place.
        No-op if --no-align-to-input, or the file didn't get written,
        or the align helper couldn't match enough atoms."""
        if not args.align_to_input:
            return
        if not Path(path).exists():
            return
        try:
            from dvbfixer.align import kabsch_align_pdb
            kabsch_align_pdb(path, str(input_path), path,
                             selection='backbone', verbose=args.verbose)
        except Exception as e:
            print(f"  WARNING: align-to-input failed on {path} ({e}); "
                  f"leaving unaligned")

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
        _maybe_align(out)
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
                      "--md-level", args.md_level,
                      "--num-output", str(args.num_output)]
        if args.no_terminal:
            model_argv.append("--no-terminal")
        if args.fasta:
            model_argv.extend(["--fasta", args.fasta])
        if not args.pin_input:
            model_argv.append("--no-pin-input")
        # Propagate the same heterogen policy zbs was invoked with, so
        # Modeller doesn't refine a loop against ligand/glycan context
        # the user has explicitly asked to strip downstream.
        if not args.keep_heterogens:
            model_argv.append("--strip-heterogens")
        if args.keep_water:
            model_argv.append("--keep-water")
        if args.verbose:
            model_argv.append("-v")
        model_main(model_argv)
        # With --num-output > 1, model writes _model_1.pdb, _model_2.pdb, ...
        # zbs picks the best (_1) for the downstream pipeline; other candidates
        # remain on disk for the user to inspect.
        if args.num_output > 1:
            multi_out = str(input_path.with_stem(input_path.stem + "_model_1"))
            if Path(multi_out).exists():
                interim_files.append(multi_out)
                interim_files.append(str(Path(multi_out).with_suffix('.dat')))
                out = multi_out
        _maybe_align(out)
        current = out

    # 3. Prepare
    if not args.skip_prepare:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: PREPARE")
        print(f"{'='*60}")
        from dvbfixer.prepare import main as prepare_main
        out = step_output("prepared")
        prepare_argv = [current, "-o", out, "--ph", str(args.ph),
                        "--ff"] + args.ff
        if not args.keep_heterogens:
            prepare_argv.append("--strip-heterogens")
        if not args.heterogen_h:
            prepare_argv.append("--no-heterogen-h")
        for mut in args.mutate:
            prepare_argv.extend(["--mutate", mut])
        if args.rename:
            prepare_argv.append("--rename")
        if args.no_infer_conect:
            prepare_argv.append("--no-infer-conect")
        if args.verbose:
            prepare_argv.append("-v")
        prepare_main(prepare_argv)
        _maybe_align(out)
        current = out

    # 4. Minimize — single pass. The new tleap-reduce prepare backend
    # already ran PROPKA + Reduce + variant H patching, so there's no
    # protonate step and no second minimize. Force ff14SB in the
    # minimize --ff so it matches tleap's ff14SB templates.
    _minimize_ff = list(args.ff)
    if _minimize_ff == ["auto"]:
        _minimize_ff = ["amber14-all.xml", "amber14/tip3pfb.xml"]
    if not args.skip_minimize:
        step_num += 1
        print(f"\n{'='*60}")
        print(f"Step {step_num}: MINIMIZE")
        print(f"{'='*60}")
        from dvbfixer.minimize import main as minimize_main
        out = step_output("minimized")
        minimize_argv = [current, "-o", out,
                         "--ph", str(args.ph),
                         "--ff"] + _minimize_ff + [
                         "--restraint-k", str(args.restraint_k),
                         "--max-iter", str(args.max_iter)]
        if args.no_solvent:
            minimize_argv.append("--no-solvent")
        if args.rebuild_h:
            minimize_argv.append("--rebuild-h")
        if not args.keep_heterogens:
            minimize_argv.append("--strip-heterogens")
        if args.rename:
            minimize_argv.append("--rename")
        if args.no_infer_conect:
            minimize_argv.append("--no-infer-conect")
        if args.parametrize_ligands:
            minimize_argv.append("--parametrize-ligands")
        if args.refine == "xtb":
            minimize_argv.append("--xtb-refine")
        elif args.refine == "obminimize":
            minimize_argv.append("--obminimize-refine")
        if args.refine != "none" and args.refine_heterogens_only:
            minimize_argv.append("--refine-heterogens-only")
        if args.platform:
            minimize_argv.extend(["--platform", args.platform])
        if args.verbose:
            minimize_argv.append("-v")
        minimize_main(minimize_argv)
        _maybe_align(out)
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
