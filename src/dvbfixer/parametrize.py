"""Parametrize small molecules with GAFF2 force field for GROMACS MD.

Pipeline: antechamber (atom types + charges) → parmchk2 (missing params)
→ tleap (AMBER topology) → ParmEd (AMBER → GROMACS conversion).

Charge methods:
  bcc  — AM1-BCC (default, fast, no QM needed)
  resp — RESP (requires Gaussian .log file with HF/6-31G* ESP)
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


# Input format detection from file extension
_FORMAT_MAP = {
    '.pdb': 'pdb',
    '.mol2': 'mol2',
    '.sdf': 'sdf',
    '.mol': 'mdl',
}


def _run_cmd(cmd, verbose=False, cwd=None):
    """Run a subprocess command. Returns (returncode, stdout, stderr)."""
    if verbose:
        print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=cwd)
    if verbose and result.stdout.strip():
        for line in result.stdout.strip().split('\n'):
            print(f"    {line}")
    return result


def run_antechamber(input_path, output_mol2, charge_method='bcc',
                    net_charge=0, multiplicity=1, gaussian_log=None,
                    atom_type='gaff2', verbose=False, cwd=None):
    """Run antechamber to assign atom types and compute charges.

    For RESP: if gaussian_log is provided, reads ESP from it.
    Otherwise uses AM1-BCC (fast, no QM needed).
    """
    input_path = str(input_path)
    output_mol2 = str(output_mol2)

    if charge_method == 'resp' and gaussian_log:
        # RESP from Gaussian log
        cmd = [
            'antechamber',
            '-i', str(gaussian_log), '-fi', 'gout',
            '-o', output_mol2, '-fo', 'mol2',
            '-c', 'resp', '-at', atom_type,
            '-nc', str(net_charge), '-m', str(multiplicity),
        ]
    else:
        # AM1-BCC or other methods
        ext = Path(input_path).suffix.lower()
        in_fmt = _FORMAT_MAP.get(ext, 'pdb')
        cmd = [
            'antechamber',
            '-i', input_path, '-fi', in_fmt,
            '-o', output_mol2, '-fo', 'mol2',
            '-c', charge_method, '-at', atom_type,
            '-nc', str(net_charge), '-m', str(multiplicity),
        ]

    result = _run_cmd(cmd, verbose=verbose, cwd=cwd)
    if result.returncode != 0:
        print(f"ERROR: antechamber failed:\n{result.stderr}", file=sys.stderr)
        return False
    return True


def run_parmchk2(mol2_path, frcmod_path, ff='gaff2', verbose=False, cwd=None):
    """Run parmchk2 to check for missing parameters."""
    cmd = [
        'parmchk2',
        '-i', str(mol2_path), '-f', 'mol2',
        '-o', str(frcmod_path), '-s', ff,
    ]
    result = _run_cmd(cmd, verbose=verbose, cwd=cwd)
    if result.returncode != 0:
        print(f"ERROR: parmchk2 failed:\n{result.stderr}", file=sys.stderr)
        return False

    # Check for ATTN lines (estimated parameters)
    with open(frcmod_path if os.path.isabs(frcmod_path)
              else os.path.join(cwd or '.', frcmod_path)) as f:
        for line in f:
            if 'ATTN' in line:
                print(f"WARNING: parmchk2 estimated parameters:\n  {line.strip()}",
                      file=sys.stderr)
    return True


def run_tleap(mol2_path, frcmod_path, prmtop_path, inpcrd_path,
              ff='gaff2', verbose=False, cwd=None):
    """Run tleap to generate AMBER topology."""
    script = f"""source leaprc.{ff}
mol = loadmol2 {mol2_path}
loadamberparams {frcmod_path}
check mol
saveamberparm mol {prmtop_path} {inpcrd_path}
quit
"""
    script_path = 'leap.in'
    abs_script = os.path.join(cwd or '.', script_path)
    with open(abs_script, 'w') as f:
        f.write(script)

    result = _run_cmd(['tleap', '-f', script_path], verbose=verbose, cwd=cwd)
    if result.returncode != 0:
        print(f"ERROR: tleap failed:\n{result.stderr}", file=sys.stderr)
        return False

    # Check that output files exist
    prmtop_abs = os.path.join(cwd or '.', prmtop_path)
    if not os.path.exists(prmtop_abs):
        print(f"ERROR: tleap did not produce {prmtop_path}", file=sys.stderr)
        if verbose:
            print(f"  tleap output:\n{result.stdout}", file=sys.stderr)
        return False
    return True


def convert_to_gromacs(prmtop_path, inpcrd_path, output_prefix, mol_name,
                       verbose=False):
    """Convert AMBER topology to GROMACS .itp + .gro using ParmEd."""
    import parmed

    amber = parmed.load_file(str(prmtop_path), str(inpcrd_path))

    # Save full GROMACS topology (includes [ defaults ], [ atomtypes ], etc.)
    top_path = f"{output_prefix}_full.top"
    amber.save(top_path, overwrite=True)
    amber.save(f"{output_prefix}.gro", overwrite=True)

    if verbose:
        total_charge = sum(a.charge for a in amber.atoms)
        print(f"  ParmEd: {len(amber.atoms)} atoms, charge={total_charge:.4f}")

    # Extract [ moleculetype ] section as standalone .itp
    _extract_itp(top_path, f"{output_prefix}.itp", mol_name)

    # Write position restraints
    _write_posre(amber, f"posre_{output_prefix}.itp")

    # Clean up full .top (user only needs .itp)
    os.unlink(top_path)


def _extract_itp(top_path, itp_path, mol_name):
    """Extract moleculetype + atomtypes from ParmEd .top into standalone .itp."""
    with open(top_path) as f:
        lines = f.readlines()

    # Sections to keep (match by normalized name)
    keep_sections = {
        'defaults', 'atomtypes', 'moleculetype',
        'atoms', 'bonds', 'pairs', 'angles', 'dihedrals',
    }

    def _section_name(line):
        """Extract section name from '[ name ]' header."""
        s = line.strip()
        if s.startswith('[') and ']' in s:
            return s.split('[')[1].split(']')[0].strip().lower()
        return None

    with open(itp_path, 'w') as f:
        f.write(f"; Moleculetype: {mol_name}\n")
        f.write(f"; Generated by dvbfixer parametrize (GAFF2)\n\n")

        in_section = False
        rename_next_data = False

        for line in lines:
            sec = _section_name(line)
            if sec is not None:
                if sec in keep_sections:
                    in_section = True
                    if sec == 'moleculetype':
                        rename_next_data = True
                else:
                    in_section = False
                    continue

            if in_section:
                stripped = line.strip()
                # Replace the moleculetype name line
                if rename_next_data and stripped and not stripped.startswith((';', '[')):
                    parts = stripped.split()
                    if len(parts) >= 2:
                        line = f"{mol_name:<16s} {parts[1]}\n"
                    rename_next_data = False
                f.write(line)

        # Add POSRES include
        f.write(f"\n; Include Position restraint file\n")
        f.write(f"#ifdef POSRES\n")
        f.write(f'#include "posre_{mol_name}.itp"\n')
        f.write(f"#endif\n")


def _write_posre(amber_struct, path, fc=1000.0):
    """Write position restraint file for heavy atoms."""
    with open(path, 'w') as f:
        f.write("; Position restraints for heavy atoms\n")
        f.write("; Generated by dvbfixer parametrize\n\n")
        f.write("[ position_restraints ]\n")
        f.write(";  ai  funct    fcx      fcy      fcz\n")
        for i, atom in enumerate(amber_struct.atoms):
            if atom.atomic_number > 1:  # non-hydrogen
                f.write(f"{i + 1:6d}     1  {fc:.1f}  {fc:.1f}  {fc:.1f}\n")


def generate_gaussian_input(input_path, output_path, net_charge=0,
                            multiplicity=1, mem='4GB', nproc=4,
                            method='HF/6-31G*', verbose=False):
    """Generate a Gaussian input file (.com) for RESP charge calculation.

    Writes the .com to the absolute resolved `output_path` (not into a
    tempdir that gets deleted). Antechamber's default gcrt template runs
    HF/6-31G* with Pop=MK ESP grid generation — exactly what RESP needs.

    Post-processes the antechamber output to add:
      - `%mem=<mem>` for memory allocation
      - `%nproc=<nproc>` for parallelism
      - `%chk=<stem>.chk` named after the structure (was 'molecule')
      - Useful remark line with the molecule name and date
      - Optional method override (default keeps antechamber's HF/6-31G* + MK)
    Returns True on success, False otherwise.
    """
    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    stem = output_path.stem

    ext = input_path.suffix.lower()
    in_fmt = _FORMAT_MAP.get(ext, 'pdb')

    with tempfile.TemporaryDirectory(prefix='dvbfixer_gau_') as tmpdir:
        # Copy input into tmpdir so antechamber doesn't pollute the user's cwd
        local_in = Path(tmpdir) / input_path.name
        shutil.copy2(input_path, local_in)
        # antechamber writes the .com to a relative path in cwd; we'll move
        # it to the absolute output_path afterwards.
        local_out = Path(tmpdir) / f'{stem}.com'
        local_gesp = f'{stem}.gesp'  # antechamber needs this as a string

        cmd = [
            'antechamber',
            '-i', local_in.name, '-fi', in_fmt,
            '-o', local_out.name, '-fo', 'gcrt',
            '-at', 'gaff2',
            '-nc', str(net_charge), '-m', str(multiplicity),
            '-gv', '1',
            '-ge', local_gesp,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmpdir)
        if result.returncode != 0 or not local_out.exists():
            print(f"ERROR: antechamber failed to generate Gaussian input:\n"
                  f"{result.stderr}", file=sys.stderr)
            return False

        com_content = local_out.read_text()

    # Post-process: add %mem and %nproc, rename %chk, customise method/remark.
    import datetime
    today = datetime.date.today().isoformat()
    lines = com_content.splitlines()
    new_lines = [
        f'%mem={mem}',
        f'%nproc={nproc}',
        f'%chk={stem}.chk',
    ]
    seen_chk = False
    seen_route = False
    skip_blank = False
    for line in lines:
        s = line.strip()
        # Drop antechamber's default %chk=molecule (we already added one)
        if s.startswith('%chk='):
            seen_chk = True
            continue
        # Override route line with user-chosen method if not the default
        if s.startswith('#') and not seen_route:
            seen_route = True
            if method != 'HF/6-31G*':
                # Substitute method but keep antechamber's MK + IOp directives
                line = line.replace('HF/6-31G*', method, 1)
        # Replace placeholder remark
        if s == 'remark line goes here':
            line = f'RESP charges for {stem} (q={net_charge}, m={multiplicity}) — {today}'
        new_lines.append(line)

    # Write to the FINAL absolute output_path (not in the deleted tempdir)
    output_path.write_text('\n'.join(new_lines) + '\n')
    if verbose:
        print(f"  Wrote {output_path} ({len(new_lines)} lines)")
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog='dvbfixer parametrize',
        description='Parametrize small molecules with GAFF2 for GROMACS MD. '
                    'Uses antechamber + parmchk2 + tleap + ParmEd pipeline.',
    )
    p.add_argument('input',
                   help='Input structure file (.pdb, .mol2, .sdf)')
    p.add_argument('-o', '--output', default=None,
                   help='Output prefix (default: input stem)')
    p.add_argument('-n', '--name', default=None,
                   help='Molecule name for [ moleculetype ] '
                        '(default: from input filename)')
    p.add_argument('-c', '--charge-method', default='bcc',
                   choices=['bcc', 'resp'],
                   help='Charge method: bcc (AM1-BCC, default) or resp')
    p.add_argument('--net-charge', type=int, default=0,
                   help='Net charge of the molecule (default: 0)')
    p.add_argument('--multiplicity', type=int, default=1,
                   help='Spin multiplicity (default: 1)')
    p.add_argument('--gaussian-log', default=None,
                   help='Gaussian log file for RESP charges (output of running '
                        'Gaussian on the .com file from --gen-gaussian)')
    p.add_argument('--gen-gaussian', action='store_true',
                   help='Generate a Gaussian .com input file for RESP charges '
                        'and exit. Run Gaussian on the .com, then re-invoke '
                        'this command with --gaussian-log.')
    p.add_argument('--gaussian-mem', default='4GB',
                   help='Memory for Gaussian %%mem= directive in the generated '
                        '.com (default: 4GB). Examples: 8GB, 16GB.')
    p.add_argument('--gaussian-nproc', type=int, default=4,
                   help='Processors for Gaussian %%nproc= directive in the '
                        'generated .com (default: 4).')
    p.add_argument('--gaussian-method', default='HF/6-31G*',
                   help='QM method for Gaussian ESP calculation (default: '
                        'HF/6-31G*, the AMBER-standard RESP recipe). Override '
                        'only if you know what you are doing — using a '
                        'different method changes the charges in non-obvious '
                        'ways.')
    p.add_argument('--keep-intermediate', action='store_true',
                   help='Keep antechamber/tleap intermediate files')
    p.add_argument('-v', '--verbose', action='store_true',
                   help='Verbose output')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input).resolve()

    if not input_path.exists():
        print(f"Error: {input_path} not found", file=sys.stderr)
        sys.exit(1)

    output_prefix = args.output or input_path.stem
    mol_name = args.name or input_path.stem.upper()

    # RESP without Gaussian log: generate input and exit
    if args.charge_method == 'resp' and not args.gaussian_log:
        if args.gen_gaussian:
            gau_path = Path(f"{output_prefix}.com").resolve()
            print(f"Generating Gaussian input for RESP charges...")
            if generate_gaussian_input(input_path, gau_path,
                                       net_charge=args.net_charge,
                                       multiplicity=args.multiplicity,
                                       mem=args.gaussian_mem,
                                       nproc=args.gaussian_nproc,
                                       method=args.gaussian_method,
                                       verbose=args.verbose):
                print(f"Wrote {gau_path}")
                print()
                print(f"Next steps:")
                print(f"  1. Run Gaussian (this produces {gau_path.stem}.log "
                      f"and {gau_path.stem}.gesp):")
                print(f"       g16 < {gau_path.name} > {gau_path.stem}.log")
                print(f"     (or g09 — both work; computation is "
                      f"HF/6-31G* opt+pop=MK).")
                print(f"  2. Re-run parametrize with the log file:")
                print(f"       dvbfixer parametrize {args.input} "
                      f"-n {mol_name} --net-charge {args.net_charge} "
                      f"-c resp --gaussian-log {gau_path.stem}.log")
            sys.exit(0)
        else:
            print("RESP charges require a Gaussian log file.")
            print("Either provide --gaussian-log FILE, or use --gen-gaussian "
                  "to create a Gaussian input file.")
            print("For quick parametrization, use AM1-BCC (default): "
                  "omit -c resp")
            sys.exit(1)

    # Work in a temp directory to keep things clean
    tmpdir = tempfile.mkdtemp(prefix='dvbfixer_param_')
    try:
        # Copy input to tmpdir
        tmp_input = os.path.join(tmpdir, input_path.name)
        shutil.copy2(input_path, tmp_input)

        if args.gaussian_log:
            glog = Path(args.gaussian_log).resolve()
            shutil.copy2(glog, os.path.join(tmpdir, glog.name))
            gaussian_log_name = glog.name
        else:
            gaussian_log_name = None

        mol2_name = 'mol.mol2'
        frcmod_name = 'mol.frcmod'
        prmtop_name = 'mol.prmtop'
        inpcrd_name = 'mol.inpcrd'

        # Step 1: Antechamber
        print(f"Running antechamber ({args.charge_method} charges, gaff2)...")
        ok = run_antechamber(
            input_path.name, mol2_name,
            charge_method=args.charge_method,
            net_charge=args.net_charge,
            multiplicity=args.multiplicity,
            gaussian_log=gaussian_log_name,
            verbose=args.verbose,
            cwd=tmpdir)
        if not ok:
            sys.exit(1)

        # Step 2: parmchk2
        print(f"Running parmchk2...")
        ok = run_parmchk2(mol2_name, frcmod_name,
                          verbose=args.verbose, cwd=tmpdir)
        if not ok:
            sys.exit(1)

        # Step 3: tleap
        print(f"Running tleap...")
        ok = run_tleap(mol2_name, frcmod_name, prmtop_name, inpcrd_name,
                       verbose=args.verbose, cwd=tmpdir)
        if not ok:
            sys.exit(1)

        # Step 4: Convert to GROMACS
        print(f"Converting to GROMACS format...")
        prmtop_abs = os.path.join(tmpdir, prmtop_name)
        inpcrd_abs = os.path.join(tmpdir, inpcrd_name)
        convert_to_gromacs(prmtop_abs, inpcrd_abs, mol_name,
                           mol_name, verbose=args.verbose)

        # Move output files to current directory / output location
        out_dir = Path(output_prefix).parent if '/' in output_prefix else Path('.')
        out_stem = Path(output_prefix).name if '/' in output_prefix else output_prefix

        for src_name in [f"{mol_name}.itp", f"{mol_name}.gro",
                         f"posre_{mol_name}.itp"]:
            if os.path.exists(src_name):
                dst = out_dir / src_name
                shutil.move(src_name, dst)
                print(f"Wrote {dst}")

        # Keep intermediate files if requested
        if args.keep_intermediate:
            for fname in [mol2_name, frcmod_name, prmtop_name,
                          inpcrd_name, 'leap.in', 'leap.log']:
                src = os.path.join(tmpdir, fname)
                if os.path.exists(src):
                    dst = out_dir / f"{out_stem}_{fname}"
                    shutil.copy2(src, dst)
                    print(f"Wrote {dst} (intermediate)")

        # Print summary
        import parmed
        amber = parmed.load_file(prmtop_abs, inpcrd_abs)
        total_charge = sum(a.charge for a in amber.atoms)
        print(f"\nSummary:")
        print(f"  Molecule:  {mol_name}")
        print(f"  Atoms:     {len(amber.atoms)}")
        print(f"  Charge:    {total_charge:+.4f}")
        print(f"  Method:    GAFF2 + {args.charge_method.upper()}")

    finally:
        # Clean up temp directory
        if not args.keep_intermediate:
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
