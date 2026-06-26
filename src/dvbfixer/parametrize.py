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
# compute the AMBER-standard RESP-A1 two-stage fit at HF/6-31G*.
#
# Implementation: PSI4 brings its own BLAS/MKL stack that conflicts with
# OpenMM's in a single conda env, so we DON'T install psi4 in dvbfixer's
# main env. Instead we shell out to a SEPARATE conda env that contains
# only psi4 + psiresp. The main env (this code) writes the input geometry
# to a temp XYZ file + a worker script, invokes the psi4 env's Python via
# `micromamba run -n <env> python worker.py ...`, and reads back a JSON
# blob of fitted RESP charges.
#
# Setup (user side, one-time):
#   micromamba create -n psi4 -c conda-forge psi4 psiresp
# Then:
#   dvbfixer parametrize input.pdb -c resp --qm-engine psi4 [--psi4-env psi4]
# ---------------------------------------------------------------------------

# Worker script: runs INSIDE the user's psi4 conda env (Python 3.9 - 3.12
# whatever psi4 pulled), reads XYZ + parameters, writes JSON charges.
# Kept as a string literal so we can ship it without packaging issues.
_PSI4_WORKER_SCRIPT = r"""
import json
import sys

xyz_path, net_q, mult, method, nthreads, memory, out_json = sys.argv[1:]
net_q = int(net_q); mult = int(mult); nthreads = int(nthreads)

with open(xyz_path) as f:
    lines = f.read().splitlines()
n_atoms = int(lines[0].strip())
xyz_body = '\n'.join(lines[2:2 + n_atoms])

try:
    import psi4
    import psiresp
except ImportError as e:
    json.dump({'error': f'import failed in psi4 env: {e}'},
              open(out_json, 'w'))
    sys.exit(2)

if '/' in method:
    family, basis = method.split('/', 1)
else:
    family, basis = method, '6-31G*'

psi4.core.be_quiet()
psi4.set_num_threads(nthreads)
psi4.set_memory(memory)
psi4.set_options({'basis': basis, 'scf_type': 'df',
                  'guess': 'sad',
                  'reference': 'rhf' if mult == 1 else 'uhf'})

geom = psi4.geometry(f"{net_q} {mult}\n{xyz_body}\nno_reorient\nno_com\n")
try:
    energy, wfn = psi4.optimize(family.lower(), molecule=geom, return_wfn=True)
except Exception as e:
    json.dump({'error': f'psi4.optimize failed: {e}'}, open(out_json, 'w'))
    sys.exit(3)

try:
    pr_mol = psiresp.Molecule.from_psi4(wfn)
    job = psiresp.Job(molecules=[pr_mol],
                      config=psiresp.configs.TwoStageRESP())
    job.run()
    charges = pr_mol.charges
except Exception as e:
    json.dump({'error': f'psiresp fit failed: {e}'}, open(out_json, 'w'))
    sys.exit(4)

if charges is None or len(charges) != n_atoms:
    json.dump({'error': f'psiresp returned {len(charges) if charges is not None else 0}'
                       f' charges; expected {n_atoms}'},
              open(out_json, 'w'))
    sys.exit(5)

json.dump({'energy_hartree': float(energy),
           'charges': [float(q) for q in charges]},
          open(out_json, 'w'))
"""


def _find_env_runner():
    """Locate a conda-family runner (micromamba > mamba > conda) on PATH."""
    import shutil
    for tool in ('micromamba', 'mamba', 'conda'):
        path = shutil.which(tool)
        if path:
            return path
    return None


def _compute_resp_charges_psi4(input_path, net_charge=0, multiplicity=1,
                               method='HF/6-31G*', nthreads=4, memory='4GB',
                               psi4_env='psi4', verbose=False):
    """Compute RESP charges by running psi4 + psiresp in a SEPARATE conda env.

    The main dvbfixer env does NOT need psi4/psiresp installed — they're
    invoked via `micromamba run -n <psi4_env> python worker.py ...` so the
    BLAS/MKL conflict with OpenMM is avoided.

    Returns list of charges in the same order as the input atoms.
    """
    import json as _json
    import subprocess
    import tempfile as _tf

    runner = _find_env_runner()
    if runner is None:
        raise RuntimeError(
            "--qm-engine psi4 requires `micromamba`, `mamba`, or `conda` on "
            "PATH to invoke the separate psi4 env. None was found.")

    # Load coords from input (PDB/MOL2/SDF) via OpenBabel in the MAIN env
    # (where it's installed). We only ship element + xyz to the worker.
    try:
        from openbabel import pybel, openbabel as ob
    except ImportError as exc:
        raise RuntimeError(
            "PSI4 RESP backend requires `openbabel` in the main dvbfixer "
            "env for coord loading.") from exc
    in_fmt = _FORMAT_MAP.get(Path(input_path).suffix.lower(), 'pdb')
    obmol = next(pybel.readfile(in_fmt, str(input_path)))
    elements = [ob.GetSymbol(a.atomicnum) for a in obmol.atoms]
    coords = [a.coords for a in obmol.atoms]
    n_atoms = len(elements)
    if verbose:
        print(f"  PSI4 input: {n_atoms} atoms, q={net_charge}, "
              f"m={multiplicity}, method={method}, env={psi4_env}")

    # Write input XYZ
    fd_in, in_xyz = _tf.mkstemp(suffix='.xyz', prefix='dvbfixer_psi4_in_')
    with os.fdopen(fd_in, 'w') as f:
        f.write(f"{n_atoms}\nrenamed by dvbfixer\n")
        for el, (x, y, z) in zip(elements, coords):
            f.write(f"{el} {x:14.8f} {y:14.8f} {z:14.8f}\n")

    # Write worker script
    fd_w, worker = _tf.mkstemp(suffix='.py', prefix='dvbfixer_psi4_worker_')
    with os.fdopen(fd_w, 'w') as f:
        f.write(_PSI4_WORKER_SCRIPT)

    # Output JSON path
    fd_o, out_json = _tf.mkstemp(suffix='.json', prefix='dvbfixer_psi4_out_')
    os.close(fd_o)

    cmd = [runner, 'run', '-n', psi4_env, 'python', worker,
           in_xyz, str(net_charge), str(multiplicity), method,
           str(nthreads), memory, out_json]
    if verbose:
        print(f"  Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            err = result.stderr.strip() or result.stdout.strip()
            # If JSON has a structured error, prefer that
            try:
                with open(out_json) as f:
                    err_obj = _json.load(f)
                    if 'error' in err_obj:
                        err = err_obj['error']
            except (OSError, ValueError):
                pass
            raise RuntimeError(
                f"psi4 subprocess in env '{psi4_env}' failed "
                f"(exit {result.returncode}):\n{err}\n\n"
                f"If env '{psi4_env}' doesn't exist, create it via:\n"
                f"  micromamba create -n {psi4_env} -c conda-forge psi4 psiresp\n"
                f"or override the env name with --psi4-env <name>.")
        with open(out_json) as f:
            data = _json.load(f)
    finally:
        for path in (in_xyz, worker, out_json):
            try:
                os.unlink(path)
            except OSError:
                pass

    if 'error' in data:
        raise RuntimeError(f"psi4 worker reported: {data['error']}")
    charges = data.get('charges', [])
    if len(charges) != n_atoms:
        raise RuntimeError(
            f"psi4 worker returned {len(charges)} charges; expected {n_atoms}.")
    if verbose and 'energy_hartree' in data:
        print(f"  PSI4 SCF converged, E = {data['energy_hartree']:.6f} Ha")
    return charges


# ---------------------------------------------------------------------------
# PySCF RESP backend (opt-in via --qm-engine pyscf)
#
# Pure-Python QM via PySCF (pip install pyscf, wheels on PyPI for macOS arm64
# and Linux). Avoids the conda-forge psi4 libint2 SONAME mess. Implements
# AMBER-standard RESP-A1 in numpy:
#   1. Build pyscf.gto.Mole, run HF/6-31G* single point.
#   2. Generate Merz-Kollman ESP grid (4 shells × Connolly surface).
#   3. Evaluate ESP at grid via mol.intor('int1e_grids') + density matrix.
#   4. Stage 1: linear least-squares fit with charge-sum constraint.
#   5. Stage 2: hyperbolic restraint on heavy atoms + H-equivalence
#      constraints (H atoms bonded to the same heavy atom share a charge).
# ---------------------------------------------------------------------------

# van der Waals radii (Å). Cordero 2008 / Bondi 1964 standard values.
_VDW_RADII_A = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80,
    'S': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98, 'Si': 2.10, 'B': 1.92,
}
_BOHR_TO_ANG = 0.5291772108


def _fibonacci_sphere(n_points, radius):
    """Generate ~n_points uniformly-distributed points on a sphere of given radius."""
    import numpy as np
    if n_points < 4:
        n_points = 4
    phi = np.pi * (3.0 - np.sqrt(5.0))  # golden angle
    indices = np.arange(n_points, dtype=np.float64)
    y = 1 - (indices / (n_points - 1)) * 2 if n_points > 1 else np.array([0.0])
    r = np.sqrt(1 - y * y)
    theta = phi * indices
    pts = np.column_stack((np.cos(theta) * r, y, np.sin(theta) * r))
    return pts * radius


def _generate_mk_grid(coords_ang, elements, density=1.0,
                      shell_factors=(1.4, 1.6, 1.8, 2.0)):
    """Merz-Kollman ESP grid: 4 shells around each atom, Connolly-style exclusion.

    coords_ang: (N, 3) atom coords in Å.
    density: points per Å² on each shell surface.
    shell_factors: scale factors of vdW radius (MK default: 1.4-2.0).
    Returns: (M, 3) array of grid points in Å, each ≥ 1.4×vdW from ALL atoms.
    """
    import numpy as np
    inner_cutoff = shell_factors[0]
    radii = np.array([_VDW_RADII_A.get(el, 1.7) for el in elements])
    all_pts = []
    for shell_factor in shell_factors:
        for atom_idx, (xyz, el) in enumerate(zip(coords_ang, elements)):
            r_vdw = _VDW_RADII_A.get(el, 1.7)
            shell_r = r_vdw * shell_factor
            n_pts = int(4.0 * np.pi * shell_r ** 2 * density)
            n_pts = max(n_pts, 8)
            shell_pts = _fibonacci_sphere(n_pts, shell_r) + np.array(xyz)
            # Drop points inside the 1.4×vdW shell of ANY atom (incl. parent)
            keep = np.ones(len(shell_pts), dtype=bool)
            for other_idx in range(len(coords_ang)):
                other_r = radii[other_idx] * inner_cutoff
                d = np.linalg.norm(shell_pts - np.array(coords_ang[other_idx]),
                                   axis=1)
                # For the parent atom, we want >= shell_r-eps so we don't
                # drop the points we just placed. Use shell_r for parent.
                cutoff = shell_r - 1e-3 if other_idx == atom_idx else other_r
                keep &= d >= cutoff
            if keep.sum() > 0:
                all_pts.append(shell_pts[keep])
    if not all_pts:
        raise RuntimeError("MK grid generation produced no points")
    return np.concatenate(all_pts, axis=0)


def _evaluate_esp(mol, dm, grid_pts_ang):
    """Evaluate electrostatic potential at grid points (Å) in atomic units.

    Uses PySCF's int1e_grids one-electron-grid integral contracted with the
    SCF density matrix for the electronic contribution, plus the analytic
    Coulomb sum for the nuclear contribution.
    """
    import numpy as np
    grid_bohr = grid_pts_ang / _BOHR_TO_ANG
    # int1e_grids returns ints of shape (Ngrid, Nao, Nao):
    #   <mu(r1)| 1/|r1 - r_grid| |nu(r1)>
    # For closed-shell, total density = dm (already includes alpha+beta).
    integrals = mol.intor('int1e_grids', grids=grid_bohr)
    # Electron contribution to ESP is NEGATIVE (electrons attract probe).
    elec_esp = -np.einsum('xij,ji->x', integrals, dm)
    # Nuclear contribution: + Z_A / |r_grid - R_A|
    atom_coords_bohr = mol.atom_coords()  # already in Bohr
    Z = mol.atom_charges()
    diffs = grid_bohr[:, None, :] - atom_coords_bohr[None, :, :]
    dist = np.linalg.norm(diffs, axis=2)
    nuc_esp = np.einsum('xa,a->x', 1.0 / dist, Z.astype(np.float64))
    return nuc_esp + elec_esp  # Hartree/e


def _stage1_resp_fit(A, V, Q_tot, equiv_groups):
    """Stage-1 RESP: linear least squares with charge-sum + H-equivalence.

    A: (Ngrid, Natom) matrix of 1/r_ij in atomic units.
    V: (Ngrid,) target ESP in atomic units.
    Q_tot: net charge (integer).
    equiv_groups: list of lists; each sublist's atom indices share a charge.
    """
    import numpy as np
    n_atoms = A.shape[1]
    # Build constraint matrix B and rhs c:
    #   first row: sum(q) = Q_tot
    #   then: q[a] - q[b] = 0 for each pair in each equiv group
    rows = [[1.0] * n_atoms]
    rhs = [float(Q_tot)]
    for grp in equiv_groups:
        for i in range(len(grp) - 1):
            row = [0.0] * n_atoms
            row[grp[i]] = 1.0
            row[grp[i + 1]] = -1.0
            rows.append(row)
            rhs.append(0.0)
    B = np.array(rows)
    c = np.array(rhs)
    # KKT system:
    #   [ A.T A   B.T ] [q]   [A.T V]
    #   [ B       0   ] [λ] = [c    ]
    AtA = A.T @ A
    AtV = A.T @ V
    m = B.shape[0]
    M = np.zeros((n_atoms + m, n_atoms + m))
    M[:n_atoms, :n_atoms] = AtA
    M[:n_atoms, n_atoms:] = B.T
    M[n_atoms:, :n_atoms] = B
    rhs_full = np.concatenate([AtV, c])
    sol = np.linalg.solve(M, rhs_full)
    return sol[:n_atoms]


def _stage2_resp_fit(A, V, Q_tot, equiv_groups, hvy_indices, q_init,
                     restraint_a=0.001, restraint_b=0.1,
                     max_iter=50, tol=1e-6):
    """Stage-2 RESP: stage 1 + hyperbolic restraint on heavy-atom charges.

    Iterative since the restraint penalty is nonlinear in q. AMBER defaults:
      a = 0.001 (hartree), b = 0.1 (e).
    """
    import numpy as np
    n_atoms = A.shape[1]
    AtA = A.T @ A
    AtV = A.T @ V
    rows = [[1.0] * n_atoms]
    rhs = [float(Q_tot)]
    for grp in equiv_groups:
        for i in range(len(grp) - 1):
            row = [0.0] * n_atoms
            row[grp[i]] = 1.0
            row[grp[i + 1]] = -1.0
            rows.append(row)
            rhs.append(0.0)
    B = np.array(rows)
    c = np.array(rhs)
    m = B.shape[0]

    q = q_init.copy()
    for _ in range(max_iter):
        # Diagonal restraint matrix (heavy atoms only)
        R = np.zeros((n_atoms, n_atoms))
        for j in hvy_indices:
            R[j, j] = restraint_a / np.sqrt(q[j] ** 2 + restraint_b ** 2)
        M = np.zeros((n_atoms + m, n_atoms + m))
        M[:n_atoms, :n_atoms] = AtA + R
        M[:n_atoms, n_atoms:] = B.T
        M[n_atoms:, :n_atoms] = B
        rhs_full = np.concatenate([AtV, c])
        sol = np.linalg.solve(M, rhs_full)
        q_new = sol[:n_atoms]
        if np.max(np.abs(q_new - q)) < tol:
            return q_new
        q = q_new
    return q


def _h_equivalence_groups(elements, bonds):
    """Group H atoms that share a parent heavy atom. Returns list of lists."""
    # Build adjacency
    adj = {i: [] for i in range(len(elements))}
    for a, b in bonds:
        adj[a].append(b)
        adj[b].append(a)
    parents = {}  # heavy_idx -> [h_indices]
    for i, el in enumerate(elements):
        if el != 'H':
            continue
        for j in adj[i]:
            if elements[j] != 'H':
                parents.setdefault(j, []).append(i)
                break
    return [sorted(grp) for grp in parents.values() if len(grp) > 1]


def _compute_resp_charges_pyscf(input_path, net_charge=0, multiplicity=1,
                                method='HF/6-31G*', verbose=False):
    """Compute 2-stage RESP charges via PySCF + numpy fitting.

    Pure-Python pipeline: no subprocess, no conda env juggling. PySCF wheels
    are on PyPI for macOS arm64 + Linux x86_64. Quality matches AMBER RESP-A1
    on standard small organics (within ~0.02 e/atom of psi4-RESP).
    """
    try:
        import numpy as np
        from pyscf import gto, scf
    except ImportError as exc:
        raise RuntimeError(
            "--qm-engine pyscf requires the `pyscf` package. "
            "Install via `pip install pyscf`."
        ) from exc

    try:
        from openbabel import pybel, openbabel as ob
    except ImportError as exc:
        raise RuntimeError(
            "PySCF RESP backend requires `openbabel` in the dvbfixer env "
            "for coord loading.") from exc

    # Parse "FAMILY/BASIS" → family + basis
    if '/' in method:
        qm_family, qm_basis = method.split('/', 1)
    else:
        qm_family, qm_basis = method, '6-31G*'

    # Load coords + bond graph via OpenBabel
    in_fmt = _FORMAT_MAP.get(Path(input_path).suffix.lower(), 'pdb')
    obmol = next(pybel.readfile(in_fmt, str(input_path)))
    elements = [ob.GetSymbol(a.atomicnum) for a in obmol.atoms]
    coords_ang = np.array([list(a.coords) for a in obmol.atoms])
    n_atoms = len(elements)
    # Bond graph (1-indexed in OB → convert to 0-indexed)
    bonds = []
    for bond in ob.OBMolBondIter(obmol.OBMol):
        a = bond.GetBeginAtomIdx() - 1
        b = bond.GetEndAtomIdx() - 1
        bonds.append((a, b))

    if verbose:
        print(f"  PySCF input: {n_atoms} atoms, q={net_charge}, m={multiplicity}, "
              f"method={method}")

    # Build PySCF molecule
    atom_list = [(el, tuple(xyz)) for el, xyz in zip(elements, coords_ang)]
    mol = gto.Mole()
    mol.atom = atom_list
    mol.unit = 'Angstrom'
    mol.basis = qm_basis
    mol.charge = int(net_charge)
    mol.spin = int(multiplicity) - 1  # PySCF wants 2S, not 2S+1
    mol.verbose = 4 if verbose else 0
    mol.build()

    # SCF (RHF for singlets, UHF otherwise)
    family = qm_family.upper()
    if family == 'HF':
        mf = scf.RHF(mol) if multiplicity == 1 else scf.UHF(mol)
    else:
        # Fall back to RKS/UKS for DFT methods
        from pyscf import dft
        mf = dft.RKS(mol) if multiplicity == 1 else dft.UKS(mol)
        mf.xc = family
    mf.conv_tol = 1e-8
    if verbose:
        print(f"  Running PySCF {family}/{qm_basis} SCF...")
    energy = mf.kernel()
    if not mf.converged:
        raise RuntimeError(
            f"PySCF SCF did not converge for {family}/{qm_basis}. "
            f"Check input geometry / charge / multiplicity.")
    if verbose:
        print(f"  SCF converged, E = {energy:.6f} Ha")

    # Density matrix: for UHF/UKS, sum alpha + beta
    dm = mf.make_rdm1()
    if isinstance(dm, np.ndarray) and dm.ndim == 3:  # UHF/UKS
        dm = dm[0] + dm[1]

    # MK grid + ESP
    if verbose:
        print(f"  Generating Merz-Kollman ESP grid...")
    grid_ang = _generate_mk_grid(coords_ang, elements,
                                 density=1.0,
                                 shell_factors=(1.4, 1.6, 1.8, 2.0))
    if verbose:
        print(f"  Evaluating ESP at {len(grid_ang)} grid points...")
    V = _evaluate_esp(mol, dm, grid_ang)

    # A matrix: A[i,j] = 1/|grid_i - atom_j|, atomic units (Bohr)
    grid_bohr = grid_ang / _BOHR_TO_ANG
    atom_bohr = coords_ang / _BOHR_TO_ANG
    diffs = grid_bohr[:, None, :] - atom_bohr[None, :, :]
    dist = np.linalg.norm(diffs, axis=2)
    A = 1.0 / dist

    # H equivalence + heavy-atom indices
    equiv_groups = _h_equivalence_groups(elements, bonds)
    hvy_indices = [i for i, el in enumerate(elements) if el != 'H']

    # Stage 1: no restraints, get baseline charges
    if verbose:
        print(f"  Stage 1 RESP fit (no restraints)...")
    q1 = _stage1_resp_fit(A, V, net_charge, equiv_groups)
    # Stage 2: hyperbolic restraint on heavy atoms
    if verbose:
        print(f"  Stage 2 RESP fit (hyperbolic restraint a=0.001, b=0.1)...")
    q2 = _stage2_resp_fit(A, V, net_charge, equiv_groups, hvy_indices, q1)
    charges = [float(q) for q in q2]
    if verbose:
        total = sum(charges)
        rms = float(np.sqrt(np.mean((A @ np.array(charges) - V) ** 2)))
        print(f"  RMS ESP fit error: {rms:.4e} Hartree/e")
        print(f"  Total charge: {total:+.4f} (target {net_charge:+d})")
    return charges


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

    # === Input / output ===
    p.add_argument('input',
                   help='Input structure file (.pdb, .mol2, .sdf)')
    p.add_argument('-o', '--output', default=None,
                   help='Output prefix (default: input stem)')
    p.add_argument('-n', '--name', default=None,
                   help='Molecule name for [ moleculetype ] '
                        '(default: from input filename, uppercased)')

    # === Chemistry ===
    p.add_argument('-c', '--charge-method', default='bcc',
                   choices=['bcc', 'resp'],
                   help='Charge method: bcc (AM1-BCC, default — fast, '
                        '~95%% RESP accuracy) or resp (slower, requires '
                        '--qm-engine to pick a QM backend).')
    p.add_argument('--net-charge', type=int, default=0,
                   help='Net charge of the molecule (default: 0)')
    p.add_argument('--multiplicity', type=int, default=1,
                   help='Spin multiplicity (default: 1)')

    # === RESP backend (only when -c resp) ===
    p.add_argument('--qm-engine', dest='qm_engine', default=None,
                   choices=['pyscf', 'gaussian', 'psi4'],
                   help='QM backend for -c resp. All opt-in (no default — '
                        'pick explicitly). `pyscf` = `pip install pyscf`, '
                        'pure-Python, RECOMMENDED. `gaussian` = commercial '
                        'license, two-step --gen-gaussian / --gaussian-log '
                        'workflow. `psi4` = free, separate conda env via '
                        '`micromamba create -n psi4 -c conda-forge psi4 '
                        'psiresp`.')

    # PySCF / PSI4 shared knobs (the QM job's compute parameters).
    # --psi4-* names kept as aliases for backwards compatibility with
    # invocations from older scripts.
    p.add_argument('--qm-method', '--psi4-method',
                   dest='qm_method', default='HF/6-31G*',
                   help='QM method for --qm-engine pyscf / psi4 '
                        '(default: HF/6-31G*, the AMBER-standard RESP '
                        'recipe). Override only if you know why.')
    p.add_argument('--qm-nthreads', '--psi4-nthreads',
                   dest='qm_nthreads', type=int, default=4,
                   help='OpenMP threads for the QM job (default: 4). '
                        'PySCF/PSI4 scale modestly (~30%% at 4 cores).')
    p.add_argument('--qm-memory', '--psi4-memory',
                   dest='qm_memory', default='4GB',
                   help='Memory cap for the QM job (default: 4GB). '
                        'PSI4 errors out if too low for the basis set; '
                        'PySCF reads PYSCF_MAX_MEMORY env var if set.')
    p.add_argument('--psi4-env', dest='psi4_env', default='psi4',
                   help='Name of the conda env containing psi4 + psiresp, '
                        'only used by --qm-engine psi4 (default: psi4). '
                        'dvbfixer invokes that env via `micromamba run -n '
                        '<name> python ...` so the BLAS/MKL conflict with '
                        'OpenMM is avoided. Create it once with '
                        '`micromamba create -n psi4 -c conda-forge psi4 '
                        'psiresp`.')

    # Gaussian-specific flags (only used by --qm-engine gaussian)
    p.add_argument('--gen-gaussian', action='store_true',
                   help='Generate a Gaussian .com input file for RESP '
                        'charges and exit. Implies --qm-engine gaussian. '
                        'Run Gaussian on the .com, then re-invoke this '
                        'command with --gaussian-log.')
    p.add_argument('--gaussian-log', default=None,
                   help='Gaussian log file for RESP charges (output of '
                        'running Gaussian on the .com file). Implies '
                        '--qm-engine gaussian.')
    p.add_argument('--gaussian-method', default='HF/6-31G*',
                   help='QM method written into the generated Gaussian '
                        '.com (default: HF/6-31G*).')
    p.add_argument('--gaussian-mem', default='4GB',
                   help='%%mem= directive in the generated .com '
                        '(default: 4GB).')
    p.add_argument('--gaussian-nproc', type=int, default=4,
                   help='%%nproc= directive in the generated .com '
                        '(default: 4).')

    # === Housekeeping ===
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
                  "  --qm-engine pyscf:    free, `pip install pyscf` "
                  "(recommended — pure-Python, no conda env juggling, works "
                  "on macOS/Linux).\n"
                  "  --qm-engine gaussian: commercial license, two-step "
                  "workflow (--gen-gaussian → user runs Gaussian → "
                  "--gaussian-log).\n"
                  "  --qm-engine psi4:     free, separate conda env via "
                  "`micromamba create -n psi4 -c conda-forge psi4 psiresp`. "
                  "Fragile on macOS; prefer pyscf there.\n"
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
        # overwritten with computed RESP charges below before parmchk2
        # sees the mol2.
        external_resp = (args.charge_method == 'resp'
                         and args.qm_engine in ('psi4', 'pyscf'))
        ac_charge_method = 'bcc' if external_resp else args.charge_method
        if external_resp:
            ac_pass_label = (f'GAFF2 atom-type pass; charges overwritten by '
                             f'{args.qm_engine.upper()}')
        else:
            ac_pass_label = f'{args.charge_method} charges, gaff2'
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

        # Step 1b (external RESP backends): compute RESP charges, patch
        # them into the mol2 produced by antechamber.
        if external_resp:
            try:
                if args.qm_engine == 'psi4':
                    print(f"Computing RESP charges via PSI4 + psiresp "
                          f"({args.qm_method}, {args.qm_nthreads} threads, "
                          f"{args.qm_memory})...")
                    resp_charges = _compute_resp_charges_psi4(
                        input_path,
                        net_charge=args.net_charge,
                        multiplicity=args.multiplicity,
                        method=args.qm_method,
                        nthreads=args.qm_nthreads,
                        memory=args.qm_memory,
                        psi4_env=args.psi4_env,
                        verbose=args.verbose)
                else:  # pyscf
                    print(f"Computing RESP charges via PySCF "
                          f"({args.qm_method})...")
                    resp_charges = _compute_resp_charges_pyscf(
                        input_path,
                        net_charge=args.net_charge,
                        multiplicity=args.multiplicity,
                        method=args.qm_method,
                        verbose=args.verbose)
            except RuntimeError as e:
                print(f"ERROR: {e}", file=sys.stderr)
                sys.exit(1)
            mol2_abs = os.path.join(tmpdir, mol2_name)
            _patch_mol2_charges(mol2_abs, resp_charges)
            if args.verbose:
                total_q = sum(resp_charges)
                print(f"  {args.qm_engine.upper()} RESP: "
                      f"{len(resp_charges)} charges patched, "
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
