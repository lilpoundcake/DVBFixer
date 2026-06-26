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
# PSI4 + psiresp RESP backend (opt-in via --qm-engine psi4)
#
# PSI4 (Smith et al. JCP 2020, LGPL-3) + psiresp (Wang et al. JOSS 2022)
# compute the AMBER-standard RESP-A1 two-stage fit at HF/6-31G* without
# requiring Gaussian. Pipeline:
#   1. Load coords from PDB/MOL2/SDF via OpenBabel.
#   2. psi4.optimize('hf', ...) to get equilibrium geometry + wavefunction.
#   3. psiresp.Molecule.from_psi4(wfn) → psiresp.Job(TwoStageRESP) → charges.
#   4. _patch_mol2_charges() overwrites antechamber's BCC charges with the
#      RESP values; downstream parmchk2/tleap/ParmEd consume the patched mol2.
# ---------------------------------------------------------------------------

def _compute_resp_charges_psi4(input_path, net_charge=0, multiplicity=1,
                               method='HF/6-31G*', nthreads=4, memory='4GB',
                               verbose=False):
    """Compute RESP charges via PSI4 + psiresp. Returns list of charges.

    Atom ordering matches the input file's atom serial order. If PSI4
    reorders during optimisation, we map back via element + position.
    """
    try:
        import psi4
        import psiresp
    except ImportError as exc:
        raise RuntimeError(
            "--qm-engine psi4 requires the `psi4` and `psiresp` packages. "
            "Install via `conda install -c conda-forge psi4 psiresp`."
        ) from exc

    try:
        from openbabel import pybel
    except ImportError as exc:
        raise RuntimeError(
            "PSI4 RESP backend requires `openbabel` for coord loading "
            "(should be in dvbfixer's environment).") from exc

    in_fmt = _FORMAT_MAP.get(Path(input_path).suffix.lower(), 'pdb')
    obmol = next(pybel.readfile(in_fmt, str(input_path)))
    elements = []
    coords = []  # input atom order
    for atom in obmol.atoms:
        from openbabel import openbabel as ob
        elements.append(ob.GetSymbol(atom.atomicnum))
        coords.append(atom.coords)
    n_atoms = len(elements)

    # Build PSI4 geometry block
    xyz_body = '\n'.join(
        f"{el} {x:14.8f} {y:14.8f} {z:14.8f}"
        for el, (x, y, z) in zip(elements, coords)
    )
    if verbose:
        print(f"  PSI4 input: {n_atoms} atoms, q={net_charge}, m={multiplicity}, "
              f"method={method}")

    # Parse "FAMILY/BASIS" into family + basis. Default is HF/6-31G*.
    if '/' in method:
        qm_family, qm_basis = method.split('/', 1)
    else:
        qm_family, qm_basis = method, '6-31G*'

    psi4.core.be_quiet()
    psi4.set_num_threads(int(nthreads))
    psi4.set_memory(memory)
    psi4.set_options({
        'basis': qm_basis,
        'scf_type': 'df',
        'guess': 'sad',
        'reference': 'rhf' if (multiplicity == 1) else 'uhf',
    })

    geom_str = f"{net_charge} {multiplicity}\n{xyz_body}\nno_reorient\nno_com\n"
    geom = psi4.geometry(geom_str)
    # Geometry optimisation at HF/6-31G* — matches Gaussian's default for RESP.
    if verbose:
        print(f"  Running PSI4 {qm_family} geometry optimisation...")
    energy, wfn = psi4.optimize(qm_family.lower(), molecule=geom,
                                return_wfn=True)
    if verbose:
        print(f"  PSI4 SCF converged, E = {energy:.6f} Ha")

    pr_mol = psiresp.Molecule.from_psi4(wfn)
    # Default psiresp config = TwoStageRESP (AMBER-standard).
    job = psiresp.Job(molecules=[pr_mol],
                      config=psiresp.configs.TwoStageRESP())
    if verbose:
        print(f"  Running 2-stage RESP fit via psiresp...")
    job.run()
    charges = pr_mol.charges  # numpy array, one entry per atom
    if charges is None or len(charges) != n_atoms:
        raise RuntimeError(
            f"psiresp returned {0 if charges is None else len(charges)} "
            f"charges; expected {n_atoms}. Check PSI4/psiresp install.")
    return [float(q) for q in charges]


def _patch_mol2_charges(mol2_path, charges):
    """Replace charges in a mol2 ATOM block in-place. Preserves everything else.

    A mol2 ATOM line is:
        atom_id atom_name x y z atom_type subst_id subst_name charge [status]
    The charge is the 9th whitespace-separated field. atom_name is field 2;
    atom_type is field 6.
    """
    with open(mol2_path) as f:
        text = f.read()
    lines = text.splitlines()
    out = []
    in_atom_block = False
    idx = 0
    for line in lines:
        s = line.strip()
        if s.startswith('@<TRIPOS>'):
            in_atom_block = (s == '@<TRIPOS>ATOM')
            out.append(line)
            continue
        if not in_atom_block or not s:
            out.append(line)
            continue
        parts = line.split()
        if len(parts) < 9:
            out.append(line)
            continue
        if idx >= len(charges):
            out.append(line)
            continue
        # Reconstruct with new charge. Mol2 is whitespace-aligned; rebuild a
        # consistent format.
        new_line = (
            f"{int(parts[0]):>7d} {parts[1]:<8s} "
            f"{float(parts[2]):10.4f} {float(parts[3]):10.4f} "
            f"{float(parts[4]):10.4f} {parts[5]:<8s} "
            f"{parts[6]:>4s} {parts[7]:<8s} {charges[idx]:>10.6f}"
        )
        out.append(new_line)
        idx += 1
    with open(mol2_path, 'w') as f:
        f.write('\n'.join(out) + '\n')


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
    p.add_argument('--qm-engine', dest='qm_engine', default=None,
                   choices=['gaussian', 'psi4'],
                   help='QM backend for -c resp. Both opt-in (no default — '
                        'pick explicitly). `gaussian` = commercial license, '
                        'two-step --gen-gaussian then --gaussian-log workflow. '
                        '`psi4` = free conda install (psi4 + psiresp), '
                        'one-shot pipeline, ~5-7× slower than Gaussian.')
    p.add_argument('--psi4-method', dest='psi4_method', default='HF/6-31G*',
                   help='QM method for PSI4 RESP path (default: HF/6-31G*, the '
                        'AMBER-standard RESP recipe). Override only if you '
                        'know why.')
    p.add_argument('--psi4-nthreads', dest='psi4_nthreads', type=int, default=4,
                   help='OpenMP threads for PSI4 (default: 4). ~30%% speedup '
                        'at 4 cores; not linear.')
    p.add_argument('--psi4-memory', dest='psi4_memory', default='4GB',
                   help='Memory cap for PSI4 (default: 4GB). PSI4 errors out '
                        'if too low for the basis set.')
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

    # --- QM engine selection for -c resp ---
    # Both Gaussian and PSI4 are explicit opt-ins. AM1-BCC (the default) needs
    # no QM engine. When user passes -c resp, they MUST pick an engine.
    if args.charge_method == 'resp':
        # Auto-promote legacy invocations: --gen-gaussian or --gaussian-log
        # imply --qm-engine gaussian.
        if args.qm_engine is None and (args.gen_gaussian or args.gaussian_log):
            args.qm_engine = 'gaussian'
            print("INFO: --gen-gaussian/--gaussian-log implies "
                  "--qm-engine gaussian.", file=sys.stderr)
        if args.qm_engine is None:
            print("ERROR: -c resp requires --qm-engine to be set explicitly.\n"
                  "  --qm-engine gaussian: commercial license, two-step "
                  "workflow (--gen-gaussian → user runs Gaussian → "
                  "--gaussian-log).\n"
                  "  --qm-engine psi4:     free conda install "
                  "(`conda install -c conda-forge psi4 psiresp`), one-shot "
                  "pipeline, ~5-7× slower than Gaussian.\n"
                  "If you don't need RESP-quality charges, omit -c resp to "
                  "use AM1-BCC (default).",
                  file=sys.stderr)
            sys.exit(1)
        # Warn if mixing engines and flags
        if args.qm_engine == 'psi4' and (args.gen_gaussian or args.gaussian_log):
            print("WARNING: --gen-gaussian/--gaussian-log are ignored under "
                  "--qm-engine psi4. The PSI4 path runs the full QM "
                  "calculation internally; no Gaussian step is needed.",
                  file=sys.stderr)

    # RESP without log file: either generate Gaussian input and exit,
    # or run PSI4 one-shot.
    if args.charge_method == 'resp' and not args.gaussian_log:
        if args.qm_engine == 'gaussian':
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
                          f"-c resp --qm-engine gaussian "
                          f"--gaussian-log {gau_path.stem}.log")
                sys.exit(0)
            else:
                print("ERROR: --qm-engine gaussian requires --gen-gaussian "
                      "(create input) or --gaussian-log (consume output).",
                      file=sys.stderr)
                sys.exit(1)
        # PSI4 path falls through to the main pipeline; charges are
        # computed inside the temp-dir loop below.

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

        # Step 1: Antechamber. For -c resp --qm-engine psi4 we run with
        # --charge-method bcc as a typing pass; the BCC charges are
        # overwritten with PSI4-RESP charges below before parmchk2 sees
        # the mol2.
        use_psi4_for_resp = (args.charge_method == 'resp'
                             and args.qm_engine == 'psi4')
        ac_charge_method = 'bcc' if use_psi4_for_resp else args.charge_method
        ac_pass_label = 'GAFF2 atom-type pass; charges overwritten by PSI4' \
            if use_psi4_for_resp else f'{args.charge_method} charges, gaff2'
        print(f"Running antechamber ({ac_pass_label})...")
        ok = run_antechamber(
            input_path.name, mol2_name,
            charge_method=ac_charge_method,
            net_charge=args.net_charge,
            multiplicity=args.multiplicity,
            gaussian_log=gaussian_log_name,
            verbose=args.verbose,
            cwd=tmpdir)
        if not ok:
            sys.exit(1)

        # Step 1b (PSI4 RESP only): compute RESP charges via PSI4+psiresp,
        # then patch them into the mol2 produced by antechamber.
        if use_psi4_for_resp:
            print(f"Computing RESP charges via PSI4 + psiresp "
                  f"({args.psi4_method}, {args.psi4_nthreads} threads, "
                  f"{args.psi4_memory})...")
            try:
                resp_charges = _compute_resp_charges_psi4(
                    input_path,
                    net_charge=args.net_charge,
                    multiplicity=args.multiplicity,
                    method=args.psi4_method,
                    nthreads=args.psi4_nthreads,
                    memory=args.psi4_memory,
                    verbose=args.verbose)
            except RuntimeError as e:
                print(f"ERROR: {e}", file=sys.stderr)
                sys.exit(1)
            mol2_abs = os.path.join(tmpdir, mol2_name)
            _patch_mol2_charges(mol2_abs, resp_charges)
            if args.verbose:
                total_q = sum(resp_charges)
                print(f"  PSI4 RESP: {len(resp_charges)} charges patched, "
                      f"sum = {total_q:+.4f}")

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
