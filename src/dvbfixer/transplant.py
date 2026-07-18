"""dvbfixer transplant — Transplant molecules from a graft PDB into an acceptor PDB.

Workflow for glycoprotein preparation with GLYCAM:
1. Extract glycosylation site residues from protein → donor.pdb
2. Submit to GLYCAM-Web → get default.pdb with renamed residues + glycans
3. Transplant: align donor to acceptor, apply transform to graft,
   replace donor residues in acceptor with graft content

Also works with CHARMM-GUI output: use simple transplant mode (--donor + --select)
to copy glycan chains or other molecules from CHARMM-GUI PDB into your structure.

Also works as a general molecule transplant tool (--donor only, no --graft).
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from dvbfixer.acpype_export import (
    add_glycam_bonds,
    export_gromacs,
    prepare_for_openmm,
)
from dvbfixer.ffutils import PROTEIN_RESIDUES


def _parse_pdb(path):
    """Parse PDB into list of dicts with all fields preserved."""
    atoms = []
    other_lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                atoms.append({
                    'record': line[:6].strip(),
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'altloc': line[16],
                    'resname': line[17:20].strip(),
                    'chain': line[21],
                    'resseq': int(line[22:26]),
                    'icode': line[26] if len(line) > 26 else ' ',
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'occupancy': line[54:60].strip() if len(line) > 54 else '1.00',
                    'bfactor': line[60:66].strip() if len(line) > 60 else '0.00',
                    'element': line[76:78].strip() if len(line) > 76 else '',
                    'raw': line,
                })
            elif line.startswith(('CONECT', 'SSBOND', 'LINK')):
                other_lines.append(line)
    return atoms, other_lines


def _get_chain_atoms(atoms, chain_id):
    return [a for a in atoms if a['chain'] == chain_id]


def _get_residue_key(atom):
    return (atom['chain'], atom['resseq'], atom['icode'])


def _get_chains(atoms):
    seen = set()
    chains = []
    for a in atoms:
        if a['chain'] not in seen:
            seen.add(a['chain'])
            chains.append(a['chain'])
    return chains


# ---------------------------------------------------------------------------
# Kabsch superposition
# ---------------------------------------------------------------------------

def _kabsch(P, Q):
    """Find R, t such that Q @ R + t aligns Q onto P."""
    centroid_P = P.mean(axis=0)
    centroid_Q = Q.mean(axis=0)
    P_c = P - centroid_P
    Q_c = Q - centroid_Q
    H = Q_c.T @ P_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    R = U @ sign_matrix @ Vt
    t = centroid_P - centroid_Q @ R
    return R, t


def _find_common_ca(atoms_a, atoms_b, align_chains=None):
    """Find matching CA atoms between two structures.

    Returns (coords_a, coords_b) arrays for Kabsch alignment.
    """
    ca_a = {}
    for a in atoms_a:
        if a['name'] == 'CA':
            ca_a[(a['chain'], a['resseq'], a['icode'])] = np.array([a['x'], a['y'], a['z']])

    ca_b = {}
    for a in atoms_b:
        if a['name'] == 'CA':
            ca_b[(a['chain'], a['resseq'], a['icode'])] = np.array([a['x'], a['y'], a['z']])

    if align_chains:
        chain_map = {}
        for pair in align_chains:
            ch_a, ch_b = pair.split(':')
            chain_map[ch_a] = ch_b
    else:
        chains_a = set(a['chain'] for a in atoms_a if a['name'] == 'CA')
        chains_b = set(a['chain'] for a in atoms_b if a['name'] == 'CA')
        common = chains_a & chains_b
        chain_map = {c: c for c in common}

    if not chain_map:
        print("ERROR: No common protein chains for alignment.", file=sys.stderr)
        print(f"  Structure A chains: {sorted(set(a['chain'] for a in atoms_a))}", file=sys.stderr)
        print(f"  Structure B chains: {sorted(set(a['chain'] for a in atoms_b))}", file=sys.stderr)
        print("  Use --align A:B to specify chain mapping.", file=sys.stderr)
        sys.exit(1)

    coords_a, coords_b = [], []
    for ch_a, ch_b in chain_map.items():
        for (ch, resseq, icode), coord in ca_a.items():
            if ch != ch_a:
                continue
            key_b = (ch_b, resseq, icode)
            if key_b in ca_b:
                coords_a.append(coord)
                coords_b.append(ca_b[key_b])

    if len(coords_a) < 3:
        print(f"ERROR: Only {len(coords_a)} matching CA atoms. Need >= 3.", file=sys.stderr)
        sys.exit(1)

    return np.array(coords_a), np.array(coords_b)


def _transform_atoms(atoms, R, t):
    for a in atoms:
        coord = np.array([a['x'], a['y'], a['z']])
        new = coord @ R + t
        a['x'], a['y'], a['z'] = new


# ---------------------------------------------------------------------------
# PDB writer
# ---------------------------------------------------------------------------

def _write_pdb(atoms_groups, path, conect_maps=None):
    """Write PDB. Each group gets a TER record.

    conect_maps: list of (conect_lines, source_serial_map) tuples.
    source_serial_map maps source serial -> (chain, resseq, atomname) for remapping.
    If None, no CONECT written.
    """
    serial = 0
    # Build (chain, resseq, atomname) -> new_serial for CONECT remapping
    atom_key_to_serial = {}

    with open(path, 'w') as f:
        f.write("REMARK    Generated by dvbfixer transplant\n")

        for group in atoms_groups:
            for atom in group:
                serial += 1
                atom_key = (atom['chain'], atom['resseq'], atom['name'])
                atom_key_to_serial[atom_key] = serial

                name = atom['name']
                if len(name) < 4:
                    name_field = f" {name:<3s}"
                else:
                    name_field = f"{name:<4s}"

                rec = 'HETATM' if atom['record'] == 'HETATM' else 'ATOM  '
                elem = atom.get('element', '')
                if not elem:
                    elem = name[0] if name[0].isalpha() else name[1]

                f.write(
                    f"{rec}{serial:5d} {name_field}"
                    f" {atom['resname']:>3s} "
                    f"{atom['chain']}"
                    f"{atom['resseq']:4d}"
                    f"{atom['icode'] if atom['icode'].strip() else ' '}"
                    f"   "
                    f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                    f"{float(atom['occupancy']):6.2f}{float(atom['bfactor']):6.2f}"
                    f"          "
                    f"{elem:>2s}\n"
                )

            if group:
                last = group[-1]
                serial += 1
                f.write(
                    f"TER   {serial:5d}      "
                    f"{last['resname']:>3s} "
                    f"{last['chain']}"
                    f"{last['resseq']:4d}\n"
                )

        # Write CONECT using atom identity (chain, resseq, name) matching
        if conect_maps:
            written = set()
            for conect_lines, src_serial_map in conect_maps:
                for line in conect_lines:
                    serials = []
                    s = line[6:]
                    while len(s) >= 5:
                        chunk = s[:5].strip()
                        if chunk:
                            serials.append(int(chunk))
                        s = s[5:]
                    if len(serials) < 2:
                        continue

                    new_serials = []
                    ok = True
                    for old_s in serials:
                        key = src_serial_map.get(old_s)
                        if key and key in atom_key_to_serial:
                            new_serials.append(atom_key_to_serial[key])
                        else:
                            ok = False
                            break

                    if ok:
                        pair = (new_serials[0], tuple(new_serials[1:]))
                        if pair not in written:
                            written.add(pair)
                            conect = f"CONECT{new_serials[0]:5d}"
                            for ns in new_serials[1:]:
                                conect += f"{ns:5d}"
                            f.write(conect + "\n")

        f.write("END\n")


# ---------------------------------------------------------------------------
# Selection parsing
# ---------------------------------------------------------------------------

def _parse_selection(spec):
    """Parse 'A,B' or 'A:NAG' or 'A:301-310'."""
    selections = []
    for part in spec.split(','):
        part = part.strip()
        if ':' in part:
            chain, filt = part.split(':', 1)
            if '-' in filt and filt.replace('-', '').isdigit():
                start, end = filt.split('-')
                selections.append((chain, None, (int(start), int(end))))
            elif filt.isdigit():
                selections.append((chain, None, (int(filt), int(filt))))
            else:
                selections.append((chain, filt, None))
        else:
            selections.append((part, None, None))
    return selections


def _select_atoms(atoms, selections):
    selected = []
    for atom in atoms:
        for chain, resname_filter, resseq_range in selections:
            if atom['chain'] != chain:
                continue
            if resname_filter and atom['resname'] != resname_filter:
                continue
            if resseq_range:
                start, end = resseq_range
                if not (start <= atom['resseq'] <= end):
                    continue
            selected.append(atom)
            break
    return selected


def _select_conect(conect_lines, selected_serials):
    selected = []
    for line in conect_lines:
        serials = []
        s = line[6:]
        while len(s) >= 5:
            chunk = s[:5].strip()
            if chunk:
                serials.append(int(chunk))
            s = s[5:]
        if serials and all(s in selected_serials for s in serials):
            selected.append(line)
    return selected


def _build_serial_map(atoms):
    """Build serial -> (chain, resseq, atomname) map for CONECT remapping."""
    smap = {}
    for a in atoms:
        smap[a['serial']] = (a['chain'], a['resseq'], a['name'])
    return smap


# ---------------------------------------------------------------------------
# Residue replacement: remove donor residues from acceptor, insert graft
# ---------------------------------------------------------------------------

def _identify_donor_residues(donor_atoms):
    """Get set of (chain, resseq, icode) from donor — these get removed from acceptor."""
    residues = set()
    for a in donor_atoms:
        residues.add((a['chain'], a['resseq'], a['icode']))
    return residues


def _build_graft_residue_map(graft_atoms, donor_residues):
    """Map graft protein residues to donor residue numbers.

    The graft has renumbered residues (1,2,3...) but they correspond to the
    donor residues. Protein residues in graft are matched by order per chain
    to donor residues sorted by resseq.

    Returns dict mapping (graft_chain, graft_resseq, graft_icode) ->
                         (acceptor_chain, acceptor_resseq, acceptor_icode)
    """
    # Group donor residues by chain, sorted
    donor_by_chain = {}
    for ch, resseq, icode in sorted(donor_residues):
        donor_by_chain.setdefault(ch, []).append((resseq, icode))

    # Group graft protein residues by chain, in order
    graft_protein_by_chain = {}
    for a in graft_atoms:
        if a['resname'] in PROTEIN_RESIDUES:
            key = (a['chain'], a['resseq'], a['icode'])
            if a['chain'] not in graft_protein_by_chain:
                graft_protein_by_chain[a['chain']] = []
            if not graft_protein_by_chain[a['chain']] or \
               graft_protein_by_chain[a['chain']][-1] != key:
                graft_protein_by_chain[a['chain']].append(key)

    # Deduplicate
    for ch in graft_protein_by_chain:
        seen = set()
        deduped = []
        for k in graft_protein_by_chain[ch]:
            if k not in seen:
                seen.add(k)
                deduped.append(k)
        graft_protein_by_chain[ch] = deduped

    mapping = {}
    for ch in graft_protein_by_chain:
        graft_residues = graft_protein_by_chain[ch]
        donor_residues_ch = donor_by_chain.get(ch, [])
        if len(graft_residues) != len(donor_residues_ch):
            print(f"  WARNING: Chain {ch}: {len(graft_residues)} graft protein residues "
                  f"vs {len(donor_residues_ch)} donor residues", file=sys.stderr)
        for i, graft_key in enumerate(graft_residues):
            if i < len(donor_residues_ch):
                d_resseq, d_icode = donor_residues_ch[i]
                mapping[graft_key] = (ch, d_resseq, d_icode)

    return mapping


def _renumber_graft(graft_atoms, residue_map, donor_residues):
    """Renumber graft atoms: protein residues get donor numbering,
    non-protein residues (glycans) get new numbers after the last
    acceptor residue in their chain.

    Detects resseq backward jumps in non-protein residues to split
    duplicate glycan trees into separate chains (e.g. two glycan trees
    from two different donor sites that share the same graft chain ID).

    Also removes graft backbone atoms that overlap with acceptor
    (N-terminal N,H and C-terminal C,O of the replaced fragment).
    """
    # Find max resseq per chain from donor residues (for glycan numbering)
    max_resseq = {}
    for ch, resseq, icode in donor_residues:
        max_resseq[ch] = max(max_resseq.get(ch, 0), resseq)

    # Collect non-protein residue keys in order, per chain
    nonprot_keys_by_chain = {}  # chain -> [keys in order]
    seen_keys = set()
    for a in graft_atoms:
        key = (a['chain'], a['resseq'], a['icode'])
        if key in residue_map or a['resname'] in PROTEIN_RESIDUES:
            continue
        if key not in seen_keys:
            seen_keys.add(key)
            ch = a['chain']
            if ch not in nonprot_keys_by_chain:
                nonprot_keys_by_chain[ch] = []
            nonprot_keys_by_chain[ch].append(key)

    # Detect resseq backward jumps within each chain's non-protein residues
    # and assign separate target chains for each segment
    used_chains = set(ch for ch, _, _ in donor_residues)
    all_ids = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')

    glycan_remap = {}
    for ch, keys in nonprot_keys_by_chain.items():
        # Split on resseq backward jumps
        segments = [[keys[0]]]
        for i in range(1, len(keys)):
            if keys[i][1] < keys[i-1][1]:
                segments.append([])
            segments[-1].append(keys[i])

        for seg_idx, seg_keys in enumerate(segments):
            # First segment keeps the original target chain
            if seg_idx == 0:
                target_ch = ch
            else:
                # Assign a new unique chain ID
                target_ch = None
                for cid in all_ids:
                    if cid not in used_chains:
                        target_ch = cid
                        break
                if target_ch is None:
                    target_ch = ch  # fallback
            used_chains.add(target_ch)

            counter = max_resseq.get(target_ch, 0) + 1000
            for key in seg_keys:
                counter += 1
                glycan_remap[key] = (target_ch, counter, ' ')

    # Apply renumbering
    for a in graft_atoms:
        key = (a['chain'], a['resseq'], a['icode'])
        if key in residue_map:
            new_ch, new_resseq, new_icode = residue_map[key]
            a['chain'] = new_ch
            a['resseq'] = new_resseq
            a['icode'] = new_icode
        elif key in glycan_remap:
            new_ch, new_resseq, new_icode = glycan_remap[key]
            a['chain'] = new_ch
            a['resseq'] = new_resseq
            a['icode'] = new_icode


# ---------------------------------------------------------------------------
# OpenMM relaxation with AMBER + GLYCAM
# ---------------------------------------------------------------------------

def _relax_structure(path, output_path, stages, verbose=False):
    """Multi-stage energy minimization with AMBER + GLYCAM.

    Protein heavy atoms are restrained; glycan atoms move freely.
    Stages define progressively reducing restraint strength.
    """
    from openmm import CustomExternalForce, LangevinMiddleIntegrator, unit
    from openmm.app import ForceField, HBonds, Modeller, NoCutoff, PDBFile, Simulation

    print("\nRelaxing structure with AMBER + GLYCAM...")

    # Preprocess: CYS→CYX for SS bonds
    temp_dir = path.parent
    temp_pdb = temp_dir / '_relax_temp.pdb'
    _, _ = prepare_for_openmm(path, temp_pdb)

    pdb = PDBFile(str(temp_pdb))
    topology = pdb.topology
    positions = pdb.positions

    # Clean up temp file
    temp_pdb.unlink(missing_ok=True)

    forcefield = ForceField('amber14-all.xml', 'amber14/GLYCAM_06j-1.xml')

    # Add missing bonds for GLYCAM residues
    # OpenMM PDBFile only infers bonds for standard residues. GLYCAM residues
    # (NLN, OLS, OLT, sugar residues) need bonds added from FF templates.
    add_glycam_bonds(topology, forcefield, verbose)

    # Add missing hydrogens (e.g., backbone H on NLN from GLYCAM)
    Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')
    modeller = Modeller(topology, positions)
    modeller.addHydrogens(forcefield)
    topology = modeller.topology
    positions = modeller.positions
    n_after = sum(1 for _ in topology.atoms())
    print(f"  After addHydrogens: {n_after} atoms")

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
    )

    # Add position restraints on protein heavy atoms
    restraint = CustomExternalForce(
        'k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    restraint.addGlobalParameter('k', 0.0)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')

    restrained_count = 0
    for atom in topology.atoms():
        if (atom.residue.name in PROTEIN_RESIDUES
                and atom.element.symbol != 'H'):
            pos = positions[atom.index]
            restraint.addParticle(
                atom.index,
                [pos.x, pos.y, pos.z],
            )
            restrained_count += 1

    system.addForce(restraint)
    print(f"  Restrained {restrained_count} protein heavy atoms")

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)

    # Get initial energy
    state = simulation.context.getState(getEnergy=True)
    e0 = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Initial energy: {e0:.1f} kJ/mol")

    KJ_PER_NM2 = unit.kilojoules_per_mole / unit.nanometer**2

    for i, (k_val, max_iter) in enumerate(stages):
        k_actual = k_val * KJ_PER_NM2
        simulation.context.setParameter('k', k_actual.value_in_unit(KJ_PER_NM2))
        simulation.minimizeEnergy(maxIterations=max_iter)
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        if k_val > 0:
            print(f"  Stage {i+1}: k={k_val:.0f} kJ/mol/nm2, {max_iter} iters -> {energy:.1f} kJ/mol")
        else:
            print(f"  Stage {i+1}: unrestrained, {max_iter} iters -> {energy:.1f} kJ/mol")

    # Write output
    state = simulation.context.getState(getPositions=True)
    final_positions = state.getPositions()

    with open(output_path, 'w') as f:
        PDBFile.writeFile(topology, final_positions, f, keepIds=True)

    print(f"  Final energy: {energy:.1f} kJ/mol")
    return output_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer transplant",
        description="Transplant molecules from graft PDB into acceptor PDB. "
                    "Supports GLYCAM-Web and CHARMM-GUI output. "
                    "Aligns via Kabsch superposition on donor CA atoms.",
    )
    io = p.add_argument_group("Input / output")
    io.add_argument("acceptor", help="Acceptor PDB file (receives molecules)")
    io.add_argument("--donor", required=True,
                    help="Donor PDB: original residues extracted from acceptor "
                         "(used for alignment and identifying replacement sites)")
    io.add_argument("--graft",
                    help="Graft PDB: modified donor + added molecules "
                         "(e.g. GLYCAM output). If omitted, donor is used as graft.")
    io.add_argument("-o", "--output",
                    help="Output PDB (default: <acceptor>_transplant.pdb)")

    selection = p.add_argument_group("Molecule selection")
    selection.add_argument("--select",
                           help="What to transplant (if no --graft): chain IDs, "
                                "'A,B' or 'A:NAG' or 'A:301-310'")

    alignment = p.add_argument_group("Alignment")
    alignment.add_argument("--align", action='append', default=[],
                           help="Enable Kabsch superposition and specify chain mapping: "
                                "DONOR:ACCEPTOR (e.g. H:H). Repeatable. "
                                "If given without value, auto-detects matching chains.")
    alignment.add_argument("--superpose", action='store_true',
                           help="Enable Kabsch superposition (auto-detect chain mapping)")

    relax = p.add_argument_group("Relaxation (OpenMM)")
    relax.add_argument("--relax", action='store_true',
                       help="Run OpenMM minimization with AMBER+GLYCAM after transplant")
    relax.add_argument("--relax-stages", default='1000:5000,100:5000,10:5000,0:5000',
                       help="Relaxation stages as k1:iter1,k2:iter2,... "
                            "k in kJ/mol/nm2 (default: 1000:5000,100:5000,10:5000,0:5000)")

    export = p.add_argument_group("GROMACS export")
    export.add_argument("--gromacs", metavar="DIR",
                        help="Export GROMACS topology via ACPYPE to DIR. "
                             "Uses AMBER+GLYCAM with per-pair 1-4 scaling ([ pairs_nb ]).")

    content = p.add_argument_group("Content selection")
    content.add_argument("--no-infer-conect", dest="no_infer_conect",
                         action="store_true",
                         help="Skip automatic CONECT inference on acceptor/donor/graft "
                              "(default: infer missing bonds before transplant).")

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument("-v", "--verbose", action='store_true')

    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    acceptor_path = Path(args.acceptor)
    donor_path = Path(args.donor)
    graft_path = Path(args.graft) if args.graft else None

    for p in [acceptor_path, donor_path] + ([graft_path] if graft_path else []):
        if not p.exists():
            print(f"File not found: {p}", file=sys.stderr)
            sys.exit(1)

    # Compute output path from ORIGINAL acceptor BEFORE we replace the path
    # with a materialised temp copy.
    output_path = Path(args.output) if args.output else \
        acceptor_path.with_stem(acceptor_path.stem + "_transplant")

    # Auto-infer CONECT on each input separately (different serial spaces).
    if not args.no_infer_conect:
        from dvbfixer.pdbutils import _materialise_inferred_pdb
        acceptor_path = Path(_materialise_inferred_pdb(
            acceptor_path, verbose=args.verbose))
        donor_path = Path(_materialise_inferred_pdb(
            donor_path, verbose=args.verbose))
        if graft_path:
            graft_path = Path(_materialise_inferred_pdb(
                graft_path, verbose=args.verbose))

    # Parse structures
    donor_atoms, donor_other = _parse_pdb(donor_path)
    acceptor_atoms, acceptor_other = _parse_pdb(acceptor_path)

    print(f"Acceptor: {acceptor_path.name} ({len(acceptor_atoms)} atoms, "
          f"chains {','.join(_get_chains(acceptor_atoms))})")
    print(f"Donor:    {donor_path.name} ({len(donor_atoms)} atoms, "
          f"chains {','.join(_get_chains(donor_atoms))})")

    if graft_path:
        graft_atoms, graft_other = _parse_pdb(graft_path)
        print(f"Graft:    {graft_path.name} ({len(graft_atoms)} atoms, "
              f"chains {','.join(_get_chains(graft_atoms))})")

        # --- Replace + Add workflow ---

        # 1. Identify donor residues to remove from acceptor
        donor_residues = _identify_donor_residues(donor_atoms)
        print(f"\nReplacing {len(donor_residues)} donor residues in acceptor:")
        for ch, resseq, icode in sorted(donor_residues):
            # Find resname from acceptor
            resname = next(
                (a['resname'] for a in acceptor_atoms
                 if a['chain'] == ch and a['resseq'] == resseq and a['icode'] == icode),
                '???')
            print(f"  {ch}:{resname}{resseq}")

        # 2. Align donor to acceptor using CA atoms
        if args.superpose or args.align:
            donor_coords, acc_coords = _find_common_ca(
                donor_atoms, acceptor_atoms, args.align or None)
            pre_rmsd = np.sqrt(np.mean(np.sum((donor_coords - acc_coords) ** 2, axis=1)))
            R, t = _kabsch(acc_coords, donor_coords)

            aligned = donor_coords @ R + t
            post_rmsd = np.sqrt(np.mean(np.sum((aligned - acc_coords) ** 2, axis=1)))
            print(f"\nAligned on {len(donor_coords)} CA atoms: "
                  f"RMSD {pre_rmsd:.2f} -> {post_rmsd:.3f} A")

            # Apply same transform to graft
            _transform_atoms(graft_atoms, R, t)
        else:
            print("\nNo superposition (use --superpose or --align to enable)")

        # 3. Renumber graft residues to match acceptor numbering
        residue_map = _build_graft_residue_map(graft_atoms, donor_residues)
        _renumber_graft(graft_atoms, residue_map, donor_residues)

        # Count graft residues
        graft_residue_keys = set(_get_residue_key(a) for a in graft_atoms)
        protein_graft = set(k for k in graft_residue_keys
                           if any(a['resname'] in PROTEIN_RESIDUES
                                  for a in graft_atoms if _get_residue_key(a) == k))
        glycan_graft = graft_residue_keys - protein_graft
        print(f"\nGraft contains: {len(protein_graft)} protein residues, "
              f"{len(glycan_graft)} non-protein residues")

        # 4. Remove donor residues from acceptor
        removed = 0
        filtered_acceptor = []
        for a in acceptor_atoms:
            key = (a['chain'], a['resseq'], a['icode'])
            if key in donor_residues:
                removed += 1
            else:
                filtered_acceptor.append(a)
        print(f"Removed {removed} atoms from acceptor at donor sites")

        # 5. Insert graft into acceptor at the right positions
        # Split graft into protein part (replaces donor) and non-protein (appended)
        graft_protein = [a for a in graft_atoms if a['resname'] in PROTEIN_RESIDUES]
        graft_nonprotein = [a for a in graft_atoms if a['resname'] not in PROTEIN_RESIDUES]

        # Build output: for each chain, insert graft protein residues at correct position
        output_atoms = []
        for chain_id in _get_chains(filtered_acceptor):
            chain_atoms = _get_chain_atoms(filtered_acceptor, chain_id)
            graft_chain_protein = [a for a in graft_protein if a['chain'] == chain_id]

            if not graft_chain_protein:
                output_atoms.extend(chain_atoms)
                continue

            # Group graft by resseq for insertion
            graft_by_resseq = {}
            for a in graft_chain_protein:
                graft_by_resseq.setdefault(a['resseq'], []).append(a)

            # Insert graft residues at correct position in chain
            prev_resseq = -9999
            inserted_resseqs = set()
            for a in chain_atoms:
                # Check if any graft residues should be inserted before this atom
                for g_resseq in sorted(graft_by_resseq.keys()):
                    if g_resseq not in inserted_resseqs and \
                       prev_resseq < g_resseq <= a['resseq']:
                        output_atoms.extend(graft_by_resseq[g_resseq])
                        inserted_resseqs.add(g_resseq)
                output_atoms.append(a)
                prev_resseq = a['resseq']

            # Append any remaining graft residues for this chain
            for g_resseq in sorted(graft_by_resseq.keys()):
                if g_resseq not in inserted_resseqs:
                    output_atoms.extend(graft_by_resseq[g_resseq])

        # Append non-protein graft atoms
        output_atoms.extend(graft_nonprotein)

        # Group for writing (by chain, keeping insertion order)
        groups = []
        current_chain = None
        current_group = []
        for a in output_atoms:
            if a['chain'] != current_chain:
                if current_group:
                    groups.append(current_group)
                current_group = []
                current_chain = a['chain']
            current_group.append(a)
        if current_group:
            groups.append(current_group)

        # Build CONECT maps: acceptor uses acceptor serials, graft uses graft serials
        # Graft atoms have been renumbered, so we need the serial map AFTER renumbering
        acc_serial_map = _build_serial_map(acceptor_atoms)
        graft_serial_map = _build_serial_map(graft_atoms)

        # Filter acceptor CONECT: remove any referencing donor residues (being replaced)
        acc_conect = [l for l in acceptor_other if l.startswith('CONECT')]
        graft_conect = [l for l in graft_other if l.startswith('CONECT')]

        conect_maps = [
            (acc_conect, acc_serial_map),
            (graft_conect, graft_serial_map),
        ]

    else:
        # --- Simple transplant (no --graft) ---
        if not args.select:
            print("ERROR: --select required when not using --graft", file=sys.stderr)
            sys.exit(1)

        selections = _parse_selection(args.select)
        transplant_atoms = _select_atoms(donor_atoms, selections)

        if not transplant_atoms:
            print(f"ERROR: No atoms matched '{args.select}'", file=sys.stderr)
            sys.exit(1)

        transplant_residues = set(_get_residue_key(a) for a in transplant_atoms)
        print(f"Selected: {len(transplant_atoms)} atoms, "
              f"{len(transplant_residues)} residues")

        if args.superpose or args.align:
            donor_coords, acc_coords = _find_common_ca(
                donor_atoms, acceptor_atoms, args.align or None)
            pre_rmsd = np.sqrt(np.mean(np.sum((donor_coords - acc_coords) ** 2, axis=1)))
            R, t = _kabsch(acc_coords, donor_coords)
            _transform_atoms(transplant_atoms, R, t)
            aligned = donor_coords @ R + t
            post_rmsd = np.sqrt(np.mean(np.sum((aligned - acc_coords) ** 2, axis=1)))
            print(f"Aligned on {len(donor_coords)} CA atoms: "
                  f"RMSD {pre_rmsd:.2f} -> {post_rmsd:.3f} A")

        # Reassign chains if collision
        used_chains = set(a['chain'] for a in acceptor_atoms)
        all_ids = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')
        for a in transplant_atoms:
            if a['chain'] in used_chains:
                for cid in all_ids:
                    if cid not in used_chains:
                        old = a['chain']
                        for a2 in transplant_atoms:
                            if a2['chain'] == old:
                                a2['chain'] = cid
                        used_chains.add(cid)
                        print(f"  Remapped chain {old} -> {cid}")
                        break

        groups = []
        for chain_id in _get_chains(acceptor_atoms):
            groups.append(_get_chain_atoms(acceptor_atoms, chain_id))
        for chain_id in _get_chains(transplant_atoms):
            groups.append(_get_chain_atoms(transplant_atoms, chain_id))

        acc_serial_map = _build_serial_map(acceptor_atoms)
        donor_serial_map = _build_serial_map(donor_atoms)
        transplant_serials = set(a['serial'] for a in transplant_atoms)
        donor_conect = _select_conect(donor_other, transplant_serials)
        acc_conect = [l for l in acceptor_other if l.startswith('CONECT')]

        conect_maps = [
            (acc_conect, acc_serial_map),
            (donor_conect, donor_serial_map),
        ]

    # Write output
    _write_pdb(groups, output_path, conect_maps=conect_maps)
    total = sum(len(g) for g in groups)
    print(f"\nWrote {output_path.name} ({total} atoms)")

    # Relax if requested
    if args.relax:
        stages = []
        for part in args.relax_stages.split(','):
            k_str, iter_str = part.split(':')
            stages.append((float(k_str), int(iter_str)))

        relax_output = output_path.with_stem(output_path.stem + '_relaxed')
        _relax_structure(output_path, relax_output, stages, verbose=args.verbose)
        print(f"Wrote {relax_output.name}")

        # Use relaxed structure for GROMACS export
        gmx_source = relax_output
    else:
        gmx_source = output_path

    # Export GROMACS topology if requested
    if args.gromacs:
        export_gromacs(gmx_source, args.gromacs, verbose=args.verbose)


if __name__ == "__main__":
    main()
