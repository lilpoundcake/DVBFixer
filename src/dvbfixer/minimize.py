"""Energy-minimize a PDB structure with OpenMM using selective restraints.

Reads a .dat file (from 'dvbfixer prepare') to determine which atoms were
added by PDBFixer. Original atoms get strong positional restraints while
newly added atoms are free to relax into physically reasonable positions.

Three-tier restraint system:
  - Original heavy atoms: strong restraint (default 100 kcal/mol/A^2)
  - New backbone (N, CA, C, O, CB): weak restraint (default 5 kcal/mol/A^2)
  - New sidechain + all hydrogens: free (no restraint)

Minimization runs in two phases:
  1. Full restraints (default 1000 iterations)
  2. Restraints reduced 10x (default 1000 iterations)
"""

import argparse
import json
import sys
from pathlib import Path

from openmm import CustomExternalForce, LangevinMiddleIntegrator
from openmm.app import ForceField, Modeller, PDBFile, PME, Simulation
from openmm.unit import kelvin, nanometer, picosecond


DEFAULT_FF = ['amber19/protein.ff19SB.xml', 'amber19/tip3p.xml']
DEFAULT_PH = 7.0
DEFAULT_PADDING = 1.0  # nm
DEFAULT_RESTRAINT_K = 100.0  # kcal/mol/A^2 for original atoms
DEFAULT_WEAK_K = 5.0  # kcal/mol/A^2 for added atoms (backbone only)
DEFAULT_MAX_ITER = 1000

BACKBONE_NAMES = {'N', 'CA', 'C', 'O', 'CB'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer minimize",
        description="Energy-minimize a PDB structure with OpenMM. Uses selective "
        "restraints: original atoms are restrained, newly added atoms (from "
        "PDBFixer .dat file) are free to relax.",
    )
    p.add_argument("input", help="Input PDB file")
    p.add_argument("-o", "--output", help="Output minimized PDB (default: <input>_minimized.pdb)")
    p.add_argument("--dat", help="Restraint data file from 'dvbfixer prepare' (default: <input>.dat)")
    p.add_argument("--ph", type=float, default=DEFAULT_PH,
                   help=f"pH for hydrogen addition if needed (default: {DEFAULT_PH})")
    p.add_argument("--ff", nargs='+', default=DEFAULT_FF,
                   help=f"Force field XML files (default: {' '.join(DEFAULT_FF)})")
    p.add_argument("--padding", type=float, default=DEFAULT_PADDING,
                   help=f"Solvent padding in nm (default: {DEFAULT_PADDING})")
    p.add_argument("--restraint-k", type=float, default=DEFAULT_RESTRAINT_K,
                   help=f"Restraint force constant for original atoms in kcal/mol/A^2 (default: {DEFAULT_RESTRAINT_K})")
    p.add_argument("--weak-k", type=float, default=DEFAULT_WEAK_K,
                   help=f"Restraint force constant for added backbone atoms in kcal/mol/A^2 (default: {DEFAULT_WEAK_K})")
    p.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER,
                   help=f"Max minimization iterations per phase (default: {DEFAULT_MAX_ITER})")
    p.add_argument("--rebuild-h", action="store_true",
                   help="Strip and re-add hydrogens via OpenMM (default: keep existing)")
    p.add_argument("--no-solvent", action="store_true",
                   help="Minimize in vacuum (no solvent box)")
    p.add_argument("--platform", choices=["CPU", "CUDA", "OpenCL", "Reference"],
                   help="OpenMM platform (default: auto-select fastest)")
    p.add_argument("--rename", action="store_true",
                   help="Rename non-canonical residues (AMBER/CHARMM) to standard names before processing")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print detailed progress")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# .dat file loading
# ---------------------------------------------------------------------------

def load_dat(path):
    """Load .dat file, return set of (chain, resid, icode, atom_name) for added atoms."""
    with open(path) as f:
        dat = json.load(f)

    added_keys = set()
    for entry in dat["added_atoms"]:
        added_keys.add((entry["chain"], entry["resid"],
                        entry["icode"], entry["atom"]))

    print(f"Loaded restraint data: {path} ({len(added_keys)} added atoms)")
    return added_keys


def resolve_new_atom_indices(topology, added_keys, verbose=False):
    """Match (chain, resid, icode, atom_name) keys from .dat against a topology.
    Returns set of atom indices that are 'new' (added by PDBFixer).
    """
    normalized_keys = {(c, r, ic.strip(), a) for c, r, ic, a in added_keys}

    indices = set()
    matched_keys = set()
    for atom in topology.atoms():
        res = atom.residue
        key = (res.chain.id, res.id, res.insertionCode.strip(), atom.name)
        if key in normalized_keys:
            indices.add(atom.index)
            matched_keys.add(key)

    unmatched = normalized_keys - matched_keys
    if unmatched and verbose:
        print(f"  {len(unmatched)} atoms from .dat not in topology "
              f"(hydrogens are re-added during minimization)")

    return indices


# ---------------------------------------------------------------------------
# Restraints & minimization
# ---------------------------------------------------------------------------

def build_restraint_force(topology, positions, new_atom_indices, strong_k, weak_k):
    """Build a CustomExternalForce with positional restraints.

    - Original heavy atoms: strong restraint to initial position
    - Newly added backbone atoms: weak restraint (keeps loop shape reasonable)
    - Newly added sidechain atoms & all hydrogens: no restraint (free to move)
    """
    # 1 kcal/mol/A^2 = 418.4 kJ/mol/nm^2
    conv = 418.4

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    n_strong = 0
    n_weak = 0
    n_free = 0

    for atom in topology.atoms():
        pos = positions[atom.index]
        x0 = pos[0].value_in_unit(nanometer)
        y0 = pos[1].value_in_unit(nanometer)
        z0 = pos[2].value_in_unit(nanometer)

        is_hydrogen = atom.element.symbol == 'H'
        is_new = atom.index in new_atom_indices

        if is_hydrogen:
            n_free += 1
            continue
        elif not is_new:
            k = strong_k * conv
            force.addParticle(atom.index, [k, x0, y0, z0])
            n_strong += 1
        elif atom.name in BACKBONE_NAMES:
            k = weak_k * conv
            force.addParticle(atom.index, [k, x0, y0, z0])
            n_weak += 1
        else:
            n_free += 1

    return force, n_strong, n_weak, n_free


def _has_hydrogens(topology):
    """Check if topology already contains hydrogen atoms."""
    return any(a.element.symbol == 'H' for a in topology.atoms())


# AMBER protonation variant names that need explicit variants in addHydrogens
AMBER_VARIANTS = {'HIE', 'HID', 'HIP', 'ASH', 'GLH', 'CYM', 'CYX', 'LYN'}


def _read_amber_renames(pdb_path):
    """Read PDB text to find AMBER protonation names before OpenMM normalizes them.

    OpenMM's PDBFile converts GLH->GLU, HIE->HIS, CYX->CYS on read.
    This reads the raw PDB to capture the original names.

    Returns dict: (chain_id, resid_str) -> amber_name
    """
    renames = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip()
            if resname in AMBER_VARIANTS:
                chain = line[21]
                resid = line[22:26].strip()
                renames[(chain, resid)] = resname
    return renames


def _build_variants(topology, amber_renames):
    """Build variants list from AMBER renames dict.

    Returns list of variant names (or None) for each residue, suitable for
    Modeller.addHydrogens(variants=...).
    """
    if not amber_renames:
        return None
    variants = []
    has_any = False
    for res in topology.residues():
        key = (res.chain.id, res.id)
        if key in amber_renames:
            variants.append(amber_renames[key])
            has_any = True
        else:
            variants.append(None)
    return variants if has_any else None


def _get_known_residues():
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    return PROTEIN_RESIDUES | SOLVENT_IONS


def _strip_hetatm(topology, positions):
    """Remove non-protein/solvent/ion residues. Returns (new_top, new_pos)."""
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    to_delete = [res for res in topology.residues() if res.name not in known]

    if not to_delete:
        return topology, positions

    modeller = Modeller(topology, positions)
    modeller.delete(to_delete)
    return modeller.topology, modeller.positions


def minimize(topology, positions, new_atom_indices, args, amber_renames=None):
    """Set up system with solvent, apply selective restraints, minimize.

    Non-protein residues (glycans, ligands) are stripped before parametrization
    and restored afterward with original coordinates.
    """
    # Strip non-protein HETATM
    stripped_top, stripped_pos = _strip_hetatm(topology, positions)
    n_stripped = sum(1 for _ in topology.residues()) - sum(1 for _ in stripped_top.residues())
    if n_stripped > 0:
        print(f"Stripped {n_stripped} non-protein residues for minimization")

    # Remap new_atom_indices to stripped topology
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    old_to_new = {}
    new_idx = 0
    for atom in topology.atoms():
        if atom.residue.name in known:
            old_to_new[atom.index] = new_idx
            new_idx += 1
    stripped_new_indices = {old_to_new[i] for i in new_atom_indices if i in old_to_new}

    print("Loading force field...")
    forcefield = ForceField(*args.ff)

    modeller = Modeller(stripped_top, stripped_pos)

    # Rename all HIS to explicit variants before any OpenMM operation.
    # OpenMM's template matching and auto-detection crash on ambiguous HIS.
    for res in modeller.topology.residues():
        if res.name == 'HIS':
            atom_names = {a.name for a in res.atoms()}
            if 'HD1' in atom_names and 'HE2' in atom_names:
                res.name = 'HIP'
            elif 'HD1' in atom_names:
                res.name = 'HID'
            elif 'HE2' in atom_names:
                res.name = 'HIE'
            elif 'ND1' in atom_names and 'NE2' not in atom_names:
                res.name = 'HID'
            elif 'NE2' in atom_names and 'ND1' not in atom_names:
                res.name = 'HIE'
            else:
                res.name = 'HIE'  # default

    # Default: keep existing hydrogens. --rebuild-h strips and re-adds via OpenMM.
    keep_h = not args.rebuild_h
    if keep_h:
        if not _has_hydrogens(modeller.topology):
            print("No hydrogens found, adding them...")
            modeller.addHydrogens(forcefield, pH=args.ph)
        else:
            # Check for residues missing H (e.g. mutated HIS with no H added yet).
            # OpenMM's addHydrogens can't add H to residues that have none —
            # it needs a valid template match first. Use PDBFixer instead,
            # which handles missing atoms properly via addMissingAtoms.
            has_incomplete = False
            for res in modeller.topology.residues():
                res_has_h = any(a.element.symbol == 'H' for a in res.atoms())
                if not res_has_h and res.name not in ('HOH', 'WAT', 'NA', 'CL',
                                                       'K', 'MG', 'CA', 'ZN'):
                    has_incomplete = True
                    break
            if has_incomplete:
                print("Found residues missing H, using PDBFixer to add them...")
                import tempfile as _tf
                from pdbfixer import PDBFixer as _PDBFixer
                with _tf.NamedTemporaryFile(mode='w', suffix='.pdb',
                                            delete=False) as _tmp:
                    PDBFile.writeFile(modeller.topology, modeller.positions,
                                     _tmp, keepIds=True)
                    _tmp_path = _tmp.name
                try:
                    with open(_tmp_path) as _f:
                        _fixer = _PDBFixer(pdbfile=_f)
                finally:
                    Path(_tmp_path).unlink()
                _fixer.findMissingResidues()
                _fixer.findMissingAtoms()
                _fixer.addMissingAtoms()
                _fixer.addMissingHydrogens(args.ph)
                modeller = Modeller(_fixer.topology, _fixer.positions)
            else:
                print("Keeping existing hydrogens from input")
    else:
        # Strip H, fix any missing heavy atoms via PDBFixer, then re-add correct H.
        # addHydrogens adds correct H for every residue, so stripping is safe.
        h_to_delete = [a for a in modeller.topology.atoms() if a.element.symbol == 'H']
        if h_to_delete:
            modeller.delete(h_to_delete)

        # Use PDBFixer to detect and add any missing heavy atoms
        import tempfile as _tf
        from pdbfixer import PDBFixer as _PDBFixer
        with _tf.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as _tmp:
            PDBFile.writeFile(modeller.topology, modeller.positions, _tmp, keepIds=True)
            _tmp_path = _tmp.name
        try:
            with open(_tmp_path) as _f:
                _fixer = _PDBFixer(pdbfile=_f)
        finally:
            Path(_tmp_path).unlink()
        _fixer.findMissingResidues()
        _fixer.findMissingAtoms()
        n_missing = sum(len(v) for v in _fixer.missingAtoms.values())
        n_terminals = sum(len(v) for v in _fixer.missingTerminals.values())
        if n_missing or n_terminals:
            print(f"Fixing {n_missing} missing atom(s), {n_terminals} terminal(s)")
        _fixer.addMissingAtoms()
        modeller = Modeller(_fixer.topology, _fixer.positions)

        print(f"Adding hydrogens (pH {args.ph})...")
        variants = _build_variants(modeller.topology, amber_renames)
        if variants:
            n_var = sum(1 for v in variants if v is not None)
            print(f"  {n_var} AMBER protonation variants detected")
            modeller.addHydrogens(forcefield, pH=args.ph, variants=variants)
        else:
            modeller.addHydrogens(forcefield, pH=args.ph)

    if not args.no_solvent:
        print("Adding solvent...")
        modeller.addSolvent(forcefield, model='tip3p',
                            padding=args.padding * nanometer)

    print("Creating system...")
    if args.no_solvent:
        from openmm.app import NoCutoff
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff)
    else:
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)

    restraint, n_strong, n_weak, n_free = build_restraint_force(
        modeller.topology, modeller.positions,
        stripped_new_indices, args.restraint_k, args.weak_k
    )
    system.addForce(restraint)
    print(f"Restraints: {n_strong} strong (original), {n_weak} weak (new backbone), {n_free} free")

    integrator = LangevinMiddleIntegrator(300 * kelvin, 1.0 / picosecond, 0.002 * picosecond)

    if args.platform:
        from openmm import Platform
        platform = Platform.getPlatformByName(args.platform)
        simulation = Simulation(modeller.topology, system, integrator, platform)
    else:
        simulation = Simulation(modeller.topology, system, integrator)

    simulation.context.setPositions(modeller.positions)

    state = simulation.context.getState(getEnergy=True)
    print(f"Energy before minimization: {state.getPotentialEnergy()}")

    # Phase 1: full restraints
    print(f"Minimizing (phase 1: {args.max_iter} iterations, original atoms restrained)...")
    simulation.minimizeEnergy(maxIterations=args.max_iter)

    state = simulation.context.getState(getEnergy=True)
    print(f"Energy after phase 1: {state.getPotentialEnergy()}")

    # Phase 2: 10x reduced restraints
    print("Minimizing (phase 2: reduced restraints)...")
    for i in range(restraint.getNumParticles()):
        idx, params = restraint.getParticleParameters(i)
        restraint.setParticleParameters(i, idx, [params[0] / 10.0, params[1], params[2], params[3]])
    restraint.updateParametersInContext(simulation.context)

    simulation.minimizeEnergy(maxIterations=args.max_iter)

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"Energy after phase 2: {state.getPotentialEnergy()}")

    min_topology = simulation.topology
    min_positions = state.getPositions()

    # Restore non-protein residues with original coordinates
    if n_stripped > 0:
        # Strip solvent from minimized result first
        if not args.no_solvent:
            min_topology, min_positions = strip_solvent(min_topology, min_positions)

        # Merge: use minimized positions for protein atoms, original for HETATM
        # Match by (chain, resid, atomname) since addHydrogens may have changed indices
        from openmm.unit import nanometer as nm_unit
        import numpy as np

        # Build position lookup from minimized protein
        min_pos_map = {}
        for atom in min_topology.atoms():
            res = atom.residue
            key = (res.chain.id, res.id, atom.name)
            min_pos_map[key] = min_positions[atom.index].value_in_unit(nm_unit)

        # Build result: original positions, overwritten by minimized where available
        n_atoms = len(positions)
        result = np.zeros((n_atoms, 3))
        for i in range(n_atoms):
            result[i] = positions[i].value_in_unit(nm_unit)

        n_updated = 0
        for atom in topology.atoms():
            res = atom.residue
            key = (res.chain.id, res.id, atom.name)
            if key in min_pos_map:
                result[atom.index] = min_pos_map[key]
                n_updated += 1

        if args.verbose:
            print(f"Restored: {n_updated} protein atoms updated, "
                  f"{n_atoms - n_updated} HETATM atoms kept original")

        return topology, result * nm_unit

    return min_topology, min_positions


def strip_solvent(topology, positions):
    """Remove water and ions from the final structure."""
    modeller = Modeller(topology, positions)
    modeller.deleteWater()

    ions_to_delete = []
    for res in modeller.topology.residues():
        if res.name.upper() in ('NA', 'CL', 'NA+', 'CL-', 'K', 'K+', 'MG', 'MG2+', 'CA2+'):
            ions_to_delete.append(res)
    if ions_to_delete:
        modeller.delete(ions_to_delete)

    return modeller.topology, modeller.positions


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_minimized")

    # .dat is optional: if --dat given, use it; otherwise auto-detect; otherwise skip
    if args.dat:
        dat_path = Path(args.dat)
    else:
        dat_path = input_path.with_suffix(".dat")

    if args.rename:
        from dvbfixer.rename import canonicalize_pdb
        import tempfile as _tf
        _tmp = Path(_tf.mktemp(suffix='.pdb'))
        n = canonicalize_pdb(input_path, _tmp, args.verbose)
        if n > 0:
            print(f"Canonicalized {n} non-canonical residue(s)")
            input_path = _tmp
        elif _tmp.exists():
            _tmp.unlink()

    print(f"=== Minimization: {input_path} ===")

    # Read AMBER protonation names before OpenMM normalizes them
    amber_renames = _read_amber_renames(str(input_path))
    if amber_renames:
        print(f"Detected {len(amber_renames)} AMBER protonation variants in input")

    pdb = PDBFile(str(input_path))
    topology = pdb.topology
    positions = pdb.positions

    # Load restraint data from .dat if available
    if args.dat and not dat_path.exists():
        print(f"Error: .dat file not found: {dat_path}", file=sys.stderr)
        sys.exit(1)

    if dat_path.exists():
        added_keys = load_dat(dat_path)
        new_atom_indices = resolve_new_atom_indices(topology, added_keys, args.verbose)

        # Also load variant overrides from .dat (saved by prepare)
        import json
        with open(dat_path) as _df:
            _dat = json.load(_df)
        for key_str, var_name in _dat.get('variant_overrides', {}).items():
            ch, rn = key_str.split(':', 1)
            if (ch, rn) not in amber_renames:
                amber_renames[(ch, rn)] = var_name
        if amber_renames:
            print(f"  Total protonation variants (PDB + .dat): {len(amber_renames)}")
    else:
        print("No .dat file — all atoms get uniform strong restraints")
        new_atom_indices = set()

    final_topology, final_positions = minimize(
        topology, positions, new_atom_indices, args,
        amber_renames=amber_renames,
    )

    # Strip solvent from output (if HETATM wasn't stripped, solvent stripping
    # wasn't done inside minimize)
    has_hetatm = any(r.name not in _get_known_residues() for r in topology.residues())
    if not args.no_solvent and not has_hetatm:
        final_topology, final_positions = strip_solvent(final_topology, final_positions)

    with open(output_path, 'w') as f:
        PDBFile.writeFile(final_topology, final_positions, f, keepIds=True)
    print(f"\nSaved minimized structure: {output_path}")
