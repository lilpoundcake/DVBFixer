"""Pull two atoms together to form a bond using OpenMM partial minimization.

Selected atoms near the bond endpoints are free to move, while all other atoms
are frozen (mass=0) but still provide forces. This is equivalent to Schrödinger's
"minimize selected atoms" workspace operation.

Force-field selection flows through `ffutils.resolve_ff` (short-name aliases,
auto-detection). For arbitrary unknown ligands, use `minimize
--parametrize-ligands` beforehand to register GAFF2 + AM1-BCC templates.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from openmm import CustomBondForce, CustomExternalForce, LangevinMiddleIntegrator
from openmm.app import Modeller, NoCutoff, PDBFile, Simulation
from openmm.unit import angstrom, kelvin, nanometer, picosecond

from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

# Default force field (protein only)
DEFAULT_FF = 'auto'

# Default target distances by bond type (angstroms)
TARGET_DISTANCES = {
    'SS': 2.05,         # CYS SG - CYS SG disulfide
    'glycosidic': 1.45, # ASN ND2/OD1 - sugar C1
    'generic': 1.50,
}

SUGAR_RESNAMES = {'NAG', 'NDG', 'BMA', 'MAN', 'FUC', 'FUL', 'GAL', 'BGC', 'GLC', 'SIA'}

# Max bonds per element in biological context
MAX_BONDS = {'C': 4, 'N': 4, 'O': 2, 'S': 2, 'H': 1, 'P': 5, 'SE': 2}

# Expected bond length ranges (angstroms) by element pair
BOND_LENGTH_RANGE = {
    ('S', 'S'): (1.8, 2.2),
    ('C', 'N'): (1.3, 1.6),
    ('C', 'O'): (1.2, 1.6),
    ('C', 'C'): (1.2, 1.6),
    ('C', 'S'): (1.7, 1.9),
    ('N', 'C'): (1.3, 1.6),
}

# Pull force constant (kcal/mol/A^2 -> kJ/mol/nm^2)
PULL_K = 500.0 * 418.4  # 500 kcal/mol/A^2 in kJ/mol/nm^2


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer pull",
        description="Pull two atoms together to form a bond using OpenMM partial "
        "minimization. Selected atoms near the bond are free, rest are frozen "
        "(mass=0). Works with any residue type (glycans auto-parametrized via OpenFF).",
    )
    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB file")
    io.add_argument("-o", "--output", help="Output PDB (default: <input>_pulled.pdb)")

    bonds = p.add_argument_group("Bond specification")
    bonds.add_argument("--bond", nargs=2, required=True, action="append", metavar="SPEC",
                       help="Two atom specs: chain:resnum:atomname (e.g. --bond H:239:SG K:239:SG). "
                            "Can be specified multiple times for multiple bonds.")
    bonds.add_argument("--target-distance", type=float,
                       help="Target bond distance in angstroms (default: auto by bond type)")
    bonds.add_argument("--anchor", metavar="SPEC",
                       help="Anchor one endpoint (freeze its side, only move the other)")

    physics = p.add_argument_group("Physics / restraints")
    physics.add_argument("--radius", type=float, default=10.0,
                         help="Radius of free region around bond endpoints in angstroms (default: 10.0)")
    physics.add_argument("--max-iter", type=int, default=1000,
                         help="Max minimization iterations (default: 1000)")

    ff = p.add_argument_group("Force field")
    ff.add_argument("--ff", nargs='+', default=[DEFAULT_FF],
                    help="Force field selection. Accepts a short name "
                         "(auto, amber, amber+glycam, charmm, ...) or an "
                         "explicit list of OpenMM XML paths. Default: 'auto' — "
                         "detect from residue names in the input. "
                         "See docs/force-fields.md.")

    content = p.add_argument_group("Content selection")
    content.add_argument("--rename", action="store_true",
                         help="Rename non-canonical residues (AMBER/CHARMM) to standard names before processing")

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument("-v", "--verbose", action="store_true",
                      help="Print detailed progress")

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(p, batch=True)
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_atom_spec(spec):
    """Parse 'chain:resnum:atomname' -> (chain_id, resSeq, atom_name)."""
    parts = spec.split(':')
    if len(parts) != 3:
        print(f"Invalid atom spec '{spec}': expected chain:resnum:atomname", file=sys.stderr)
        sys.exit(1)
    return parts[0], int(parts[1]), parts[2]


def find_atom_index(topology, chain_id, resSeq, atomName):
    """Find atom index in OpenMM topology."""
    for atom in topology.atoms():
        res = atom.residue
        if (res.chain.id == chain_id and
                res.id == str(resSeq) and
                atom.name == atomName):
            return atom.index
    return -1


def get_target_distance(topology, idx_a, idx_b, user_override=None):
    """Determine target bond distance from atom/residue types."""
    if user_override is not None:
        return user_override

    atoms = list(topology.atoms())
    a, b = atoms[idx_a], atoms[idx_b]

    if (a.residue.name == 'CYS' and a.name == 'SG' and
            b.residue.name == 'CYS' and b.name == 'SG'):
        return TARGET_DISTANCES['SS']

    if a.residue.name in SUGAR_RESNAMES or b.residue.name in SUGAR_RESNAMES:
        return TARGET_DISTANCES['glycosidic']

    return TARGET_DISTANCES['generic']


def identify_free_atoms(topology, positions, idx_a, idx_b, radius, anchor_idx=None):
    """Find atom indices within radius of either bond endpoint.

    If anchor_idx is set, only atoms near the OTHER endpoint are free.
    Returns set of free atom indices.
    """
    coords = np.array([[p.value_in_unit(nanometer) for p in pos] for pos in positions])
    pos_a = coords[idx_a]
    pos_b = coords[idx_b]

    d_a = np.linalg.norm(coords - pos_a, axis=1) * 10.0  # nm -> A
    d_b = np.linalg.norm(coords - pos_b, axis=1) * 10.0

    if anchor_idx == idx_a:
        mask = d_b <= radius
    elif anchor_idx == idx_b:
        mask = d_a <= radius
    else:
        mask = (d_a <= radius) | (d_b <= radius)

    return set(np.where(mask)[0])


# ---------------------------------------------------------------------------
# Pre-pull hydrogen cleanup
# ---------------------------------------------------------------------------

def remove_bond_hydrogens(topology, positions, bond_info, verbose=False):
    """Remove hydrogens that conflict with new bonds.

    - CYS SG: remove HG (thiol hydrogen) before disulfide formation
    - ASN ND2: remove HD22 (extra amide H) before glycosidic bond
    """
    atoms_to_delete = []
    atoms_list = list(topology.atoms())

    for idx_a, idx_b in bond_info:
        for idx in (idx_a, idx_b):
            atom = atoms_list[idx]
            res = atom.residue

            # CYS SG -> remove HG
            if res.name == 'CYS' and atom.name == 'SG':
                for bond in topology.bonds():
                    a1, a2 = bond
                    partner = a2 if a1.index == idx else (a1 if a2.index == idx else None)
                    if partner and partner.element.symbol == 'H':
                        atoms_to_delete.append(partner)
                        if verbose:
                            print(f"  Removing {res.chain.id}:{res.name}{res.id}:{partner.name} "
                                  f"(thiol H before disulfide)")

            # ASN ND2 -> remove last H (HD22)
            if res.name == 'ASN' and atom.name == 'ND2':
                h_atoms = []
                for bond in topology.bonds():
                    a1, a2 = bond
                    partner = a2 if a1.index == idx else (a1 if a2.index == idx else None)
                    if partner and partner.element.symbol == 'H':
                        h_atoms.append(partner)
                if len(h_atoms) > 1:
                    to_remove = sorted(h_atoms, key=lambda a: a.name)[-1]
                    atoms_to_delete.append(to_remove)
                    if verbose:
                        print(f"  Removing {res.chain.id}:{res.name}{res.id}:{to_remove.name} "
                              f"(extra H before glycosidic bond)")

    if not atoms_to_delete:
        return topology, positions

    # Deduplicate
    seen = set()
    unique = []
    for a in atoms_to_delete:
        if a.index not in seen:
            seen.add(a.index)
            unique.append(a)

    modeller = Modeller(topology, positions)
    modeller.delete(unique)
    print(f"Removed {len(unique)} hydrogen(s) conflicting with new bonds")
    return modeller.topology, modeller.positions


# ---------------------------------------------------------------------------
# Valence validation
# ---------------------------------------------------------------------------

def count_bonds(topology, atom_index):
    """Count bonds for a given atom in the topology."""
    n = 0
    for bond in topology.bonds():
        if bond[0].index == atom_index or bond[1].index == atom_index:
            n += 1
    return n


def validate_pre_pull(topology, resolved_bonds, verbose=False):
    """Check that proposed bonds are chemically reasonable.

    Warns (does not error) if valence would be exceeded or bond type is unusual.
    Only examines the two atoms at each bond endpoint.
    """
    atoms_list = list(topology.atoms())

    for idx_a, idx_b, target_dist, spec_a, spec_b in resolved_bonds:
        a, b = atoms_list[idx_a], atoms_list[idx_b]
        elem_a = a.element.symbol
        elem_b = b.element.symbol

        # Check bond type
        pair = (elem_a, elem_b)
        pair_rev = (elem_b, elem_a)
        if pair not in BOND_LENGTH_RANGE and pair_rev not in BOND_LENGTH_RANGE:
            print(f"  Warning: unusual bond type {elem_a}-{elem_b} "
                  f"({spec_a} <-> {spec_b})")

        # Check expected bond length range
        brange = BOND_LENGTH_RANGE.get(pair) or BOND_LENGTH_RANGE.get(pair_rev)
        if brange and not (brange[0] <= target_dist <= brange[1]):
            print(f"  Warning: target distance {target_dist:.2f} A outside expected "
                  f"range {brange[0]:.1f}-{brange[1]:.1f} A for {elem_a}-{elem_b} bond")

        # Check valence
        for idx, atom, spec in [(idx_a, a, spec_a), (idx_b, b, spec_b)]:
            elem = atom.element.symbol
            max_b = MAX_BONDS.get(elem)
            if max_b is None:
                continue
            current = count_bonds(topology, idx)
            if current >= max_b:
                print(f"  Warning: {spec} ({atom.residue.name}:{atom.name}) "
                      f"already has {current} bond(s) (max for {elem} = {max_b})")
            elif verbose:
                print(f"  {spec}: {current}/{max_b} bonds ({elem})")


def validate_post_pull(topology, positions, resolved_bonds, verbose=False):
    """Check pulled bonds converged and no steric clashes in pulling residues."""
    atoms_list = list(topology.atoms())

    for idx_a, idx_b, target_dist, spec_a, spec_b in resolved_bonds:
        final_dist = _measure_dist(positions, idx_a, idx_b)
        a, b = atoms_list[idx_a], atoms_list[idx_b]
        elem_a, elem_b = a.element.symbol, b.element.symbol

        # Distance convergence
        print(f"Final: {spec_a} <-> {spec_b}: {final_dist:.2f} A (target: {target_dist:.2f} A)")
        if abs(final_dist - target_dist) > 0.5:
            print(f"  Warning: pull did not converge — {abs(final_dist - target_dist):.2f} A off target")

        # Bond length range
        pair = (elem_a, elem_b)
        pair_rev = (elem_b, elem_a)
        brange = BOND_LENGTH_RANGE.get(pair) or BOND_LENGTH_RANGE.get(pair_rev)
        if brange and not (brange[0] - 0.3 <= final_dist <= brange[1] + 0.3):
            print(f"  Warning: final distance {final_dist:.2f} A outside expected range "
                  f"for {elem_a}-{elem_b} bond ({brange[0]:.1f}-{brange[1]:.1f} A)")

        # Steric clashes within pulling residues
        res_a, res_b = a.residue, b.residue
        res_atoms = [at for at in atoms_list
                     if at.element.symbol != 'H' and at.residue in (res_a, res_b)]
        bonded_pairs = set()
        for bond in topology.bonds():
            bonded_pairs.add((bond[0].index, bond[1].index))
            bonded_pairs.add((bond[1].index, bond[0].index))

        for i, at1 in enumerate(res_atoms):
            for at2 in res_atoms[i+1:]:
                if (at1.index, at2.index) in bonded_pairs:
                    continue
                d = _measure_dist(positions, at1.index, at2.index)
                if d < 1.0:
                    print(f"  Warning: steric clash {at1.residue.name}:{at1.name} - "
                          f"{at2.residue.name}:{at2.name} = {d:.2f} A")


# ---------------------------------------------------------------------------
# OpenMM minimization
# ---------------------------------------------------------------------------

def strip_hetatm(topology, positions):
    """Remove non-protein residues from topology. Returns (new_top, new_pos, removed_info).

    removed_info is a list of (chain_id, res_id, res_name, atom_data) for restoration.
    """
    from openmm.app import Modeller

    known = PROTEIN_RESIDUES | SOLVENT_IONS
    to_delete = []
    for res in topology.residues():
        if res.name not in known:
            to_delete.append(res)

    if not to_delete:
        return topology, positions, []

    modeller = Modeller(topology, positions)
    modeller.delete(to_delete)
    return modeller.topology, modeller.positions, to_delete


def run_pull(topology, positions, free_indices, bond_pairs, ff_files, max_iter, verbose):
    """Pull atoms together using OpenMM with frozen/free atoms.

    Args:
        bond_pairs: list of (idx_a, idx_b, target_dist) tuples

    Strips HETATM residues before minimization (they can't be parametrized),
    minimizes protein with mass=0 freezing, then restores HETATM.
    For bonds where one endpoint is HETATM, uses CustomExternalForce to pull
    the protein atom toward the fixed HETATM position.
    """
    from openmm.app import ForceField
    from openmm.unit import nanometer as nm_unit

    # Strip non-protein residues
    stripped_top, stripped_pos, removed = strip_hetatm(topology, positions)

    if removed and verbose:
        print(f"Stripped {len(removed)} non-protein residues for minimization")

    # Remap atom indices after stripping
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    old_to_new = {}
    new_idx = 0
    for atom in topology.atoms():
        if atom.residue.name in known:
            old_to_new[atom.index] = new_idx
            new_idx += 1

    # Classify bond pairs: protein-protein vs protein-HETATM
    pp_bonds = []   # both in protein: (new_a, new_b, target_dist)
    ph_bonds = []   # one in HETATM: (new_protein_idx, hetatm_position_nm, target_dist)

    for idx_a, idx_b, target_dist in bond_pairs:
        new_a = old_to_new.get(idx_a)
        new_b = old_to_new.get(idx_b)

        if new_a is not None and new_b is not None:
            pp_bonds.append((new_a, new_b, target_dist))
        elif new_a is not None and new_b is None:
            # b is HETATM — pull a toward b's fixed position
            pos_b = [positions[idx_b][i].value_in_unit(nm_unit) for i in range(3)]
            ph_bonds.append((new_a, pos_b, target_dist))
        elif new_a is None and new_b is not None:
            # a is HETATM — pull b toward a's fixed position
            pos_a = [positions[idx_a][i].value_in_unit(nm_unit) for i in range(3)]
            ph_bonds.append((new_b, pos_a, target_dist))
        else:
            print("Error: both bond endpoints are in non-protein residues. "
                  "Pull between two HETATM atoms is not supported.",
                  file=sys.stderr)
            sys.exit(1)

    if verbose:
        if pp_bonds:
            print(f"  {len(pp_bonds)} protein-protein bond(s)")
        if ph_bonds:
            print(f"  {len(ph_bonds)} protein-HETATM bond(s) (pulling toward fixed position)")

    new_free = {old_to_new[i] for i in free_indices if i in old_to_new}

    if verbose:
        print("Loading force field...")

    forcefield = ForceField(*ff_files)

    if verbose:
        print("Creating system (vacuum, no cutoff)...")

    system = forcefield.createSystem(stripped_top, nonbondedMethod=NoCutoff)

    # Freeze all atoms except free ones (mass=0)
    n_frozen = 0
    for i in range(system.getNumParticles()):
        if i not in new_free:
            system.setParticleMass(i, 0.0)
            n_frozen += 1

    if verbose:
        print(f"Atoms: {len(new_free)} free, {n_frozen} frozen")

    # Add harmonic pull forces for protein-protein bonds
    if pp_bonds:
        bond_force = CustomBondForce("0.5*k*(r-r0)^2")
        bond_force.addPerBondParameter("k")
        bond_force.addPerBondParameter("r0")
        for new_a, new_b, target_dist in pp_bonds:
            target_nm = target_dist / 10.0
            bond_force.addBond(new_a, new_b, [PULL_K, target_nm])
        system.addForce(bond_force)

    # Add pull-to-point forces for protein-HETATM bonds
    # E = 0.5*k*(dist - r0)^2 where dist = distance from atom to fixed HETATM position
    if ph_bonds:
        ext_force = CustomExternalForce(
            "0.5*k*(sqrt((x-tx)^2+(y-ty)^2+(z-tz)^2)-r0)^2"
        )
        ext_force.addPerParticleParameter("k")
        ext_force.addPerParticleParameter("r0")
        ext_force.addPerParticleParameter("tx")
        ext_force.addPerParticleParameter("ty")
        ext_force.addPerParticleParameter("tz")
        for prot_idx, target_pos, target_dist in ph_bonds:
            target_nm = target_dist / 10.0
            ext_force.addParticle(prot_idx,
                                  [PULL_K, target_nm, target_pos[0], target_pos[1], target_pos[2]])
        system.addForce(ext_force)

    integrator = LangevinMiddleIntegrator(300 * kelvin, 1.0 / picosecond, 0.002 * picosecond)
    simulation = Simulation(stripped_top, system, integrator)
    simulation.context.setPositions(stripped_pos)

    state = simulation.context.getState(getEnergy=True)
    if verbose:
        print(f"Energy before: {state.getPotentialEnergy()}")

    print(f"Minimizing ({max_iter} iterations, {len(bond_pairs)} bond(s))...")
    simulation.minimizeEnergy(maxIterations=max_iter)

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    if verbose:
        print(f"Energy after: {state.getPotentialEnergy()}")

    # Restore original topology with HETATM — merge minimized protein coords
    # with original HETATM coords
    from openmm.unit import nanometer as nm_unit
    minimized_pos = state.getPositions()
    new_to_old = {v: k for k, v in old_to_new.items()}

    n_atoms = len(positions)
    result = np.zeros((n_atoms, 3))
    for i in range(n_atoms):
        result[i] = positions[i].value_in_unit(nm_unit)
    for new_i in range(len(minimized_pos)):
        old_i = new_to_old[new_i]
        result[old_i] = minimized_pos[new_i].value_in_unit(nm_unit)

    return result * nm_unit


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(topology, positions, output_path, bond_info_list, input_path=None):
    """Write PDB with SSBOND/CONECT for new bonds.

    bond_info_list: list of (idx_a, idx_b) tuples for new bonds.
    """
    from dvbfixer.pdbutils import append_before_end, build_serial_map

    with open(output_path, 'w') as f:
        PDBFile.writeFile(topology, positions, f, keepIds=True)

    serial_map = build_serial_map(output_path)
    atoms = list(topology.atoms())

    extra_lines = []
    ss_count = 0

    for idx_a, idx_b in bond_info_list:
        a, b = atoms[idx_a], atoms[idx_b]

        key_a = (a.residue.chain.id, a.residue.id, a.name)
        key_b = (b.residue.chain.id, b.residue.id, b.name)
        serial_a = serial_map.get(key_a)
        serial_b = serial_map.get(key_b)

        # SSBOND for CYS-CYS
        if (a.residue.name == 'CYS' and b.residue.name == 'CYS' and
                a.name == 'SG' and b.name == 'SG'):
            ss_count += 1
            ch_a, res_a = a.residue.chain.id, int(a.residue.id)
            ch_b, res_b = b.residue.chain.id, int(b.residue.id)
            ssbond = f"SSBOND {ss_count:3d} CYS {ch_a} {res_a:4d}    CYS {ch_b} {res_b:4d}"
            extra_lines.append(ssbond.ljust(80) + "\n")

        # CONECT for the new bond
        if serial_a and serial_b:
            extra_lines.append(f"CONECT{serial_a:5d}{serial_b:5d}".ljust(80) + "\n")
            extra_lines.append(f"CONECT{serial_b:5d}{serial_a:5d}".ljust(80) + "\n")

    append_before_end(output_path, extra_lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _measure_dist(positions, idx_a, idx_b):
    """Measure distance between two atoms in angstroms."""
    pa = positions[idx_a]
    pb = positions[idx_b]
    return np.sqrt(sum((a.value_in_unit(angstrom) - b.value_in_unit(angstrom))**2
                       for a, b in zip(pa, pb)))


def _resolve_bonds(topology, positions, bond_specs, target_distance_override=None):
    """Resolve bond specs to atom indices. Returns list of (idx_a, idx_b, target_dist, spec_a, spec_b)."""
    resolved = []
    for spec_a, spec_b in bond_specs:
        chain_a, res_a, name_a = parse_atom_spec(spec_a)
        chain_b, res_b, name_b = parse_atom_spec(spec_b)

        idx_a = find_atom_index(topology, chain_a, res_a, name_a)
        idx_b = find_atom_index(topology, chain_b, res_b, name_b)
        if idx_a < 0:
            print(f"Atom not found: {spec_a}", file=sys.stderr)
            sys.exit(1)
        if idx_b < 0:
            print(f"Atom not found: {spec_b}", file=sys.stderr)
            sys.exit(1)

        target_dist = get_target_distance(topology, idx_a, idx_b, target_distance_override)
        resolved.append((idx_a, idx_b, target_dist, spec_a, spec_b))
    return resolved


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_pulled")

    from dvbfixer.ffutils import print_ff_selection, resolve_ff
    args.ff, _ff_alias, _ff_reason = resolve_ff(
        args.ff, input_path, verbose=args.verbose)
    print_ff_selection(_ff_alias, _ff_reason, args.ff)

    if args.rename:
        from dvbfixer.rename import canonicalize_pdb
        from dvbfixer.tempfiles import make_temp_path
        _tmp = make_temp_path(suffix='.pdb')
        n = canonicalize_pdb(input_path, _tmp, args.verbose)
        if n > 0:
            print(f"Canonicalized {n} non-canonical residue(s)")
            input_path = _tmp
        elif _tmp.exists():
            _tmp.unlink()

    pdb = PDBFile(str(input_path))
    topology = pdb.topology
    positions = pdb.positions

    # Resolve bond specs to atom indices
    resolved = _resolve_bonds(topology, positions, args.bond, args.target_distance)

    # Print initial distances and validate
    for idx_a, idx_b, target_dist, spec_a, spec_b in resolved:
        initial_dist = _measure_dist(positions, idx_a, idx_b)
        print(f"=== Bond: {spec_a} <-> {spec_b} ===")
        print(f"  Initial distance: {initial_dist:.2f} A")
        print(f"  Target distance:  {target_dist:.2f} A")

    validate_pre_pull(topology, resolved, args.verbose)

    # Handle anchor
    anchor_idx = None
    if args.anchor:
        ac, ar, an = parse_atom_spec(args.anchor)
        anchor_idx = find_atom_index(topology, ac, ar, an)
        if anchor_idx < 0:
            print(f"Anchor atom not found: {args.anchor}", file=sys.stderr)
            sys.exit(1)

    # Build bond pairs and free atom regions
    bond_pairs = []  # (idx_a, idx_b, target_dist)
    bond_info = []   # (idx_a, idx_b) for write_output
    free_indices = set()
    needs_pull = False

    for idx_a, idx_b, target_dist, spec_a, spec_b in resolved:
        bond_info.append((idx_a, idx_b))
        initial_dist = _measure_dist(positions, idx_a, idx_b)

        if initial_dist <= target_dist + 0.5:
            print(f"  {spec_a} <-> {spec_b}: already close enough, skipping.")
        else:
            needs_pull = True
            bond_pairs.append((idx_a, idx_b, target_dist))

            fi = identify_free_atoms(
                topology, positions, idx_a, idx_b, args.radius, anchor_idx
            )
            if idx_a not in fi or idx_b not in fi:
                print("Error: bond endpoints outside free region. Increase --radius.",
                      file=sys.stderr)
                sys.exit(1)
            free_indices |= fi

    if needs_pull:
        print(f"\nFree atoms: {len(free_indices)} (within {args.radius} A radius)")
        new_positions = run_pull(
            topology, positions, free_indices, bond_pairs,
            args.ff, args.max_iter, args.verbose
        )
    else:
        new_positions = positions

    # Remove conflicting hydrogens after minimization (HG on CYS SG, HD22 on glycosylated ASN ND2)
    topology, new_positions = remove_bond_hydrogens(topology, new_positions, bond_info, args.verbose)

    # Re-resolve indices after hydrogen removal (topology changed)
    resolved = _resolve_bonds(topology, new_positions, args.bond, args.target_distance)
    bond_info = [(idx_a, idx_b) for idx_a, idx_b, _, _, _ in resolved]

    # Validate results
    print()
    validate_post_pull(topology, new_positions, resolved, args.verbose)

    write_output(topology, new_positions, output_path,
                 bond_info, input_path=input_path)
    print(f"Wrote {output_path}")
