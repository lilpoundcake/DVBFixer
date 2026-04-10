"""Convert PDB glycan nomenclature to GLYCAM force field naming.

GLYCAM uses 3-character residue codes encoding [linkage][sugar][anomer]:
  - Linkage: 0=terminal, 2-9=single position, V/W/U/Z/X/Y=multi-position
  - Sugar: G=Glc, L=Gal, M=Man, Y=GlcNAc, V=GalNAc, f=Fuc, S=Neu5Ac, ...
  - Anomer: A=alpha, B=beta (lowercase sugar code = L-sugar)

Detects glycosidic bonds from CONECT records (or distance-based fallback),
determines linkage patterns, and renames residues and atoms accordingly.
Optionally adds ROH cap at the reducing end.
"""

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path


# PDB residue name -> (GLYCAM sugar code, GLYCAM anomer code)
# Lowercase sugar code = L-sugar
PDB_TO_GLYCAM = {
    'BGC': ('G', 'B'),   # beta-D-glucopyranose
    'GLC': ('G', 'A'),   # alpha-D-glucopyranose
    'GAL': ('L', 'B'),   # beta-D-galactopyranose
    'BGA': ('L', 'A'),   # alpha-D-galactopyranose
    'MAN': ('M', 'A'),   # alpha-D-mannopyranose
    'BMA': ('M', 'B'),   # beta-D-mannopyranose
    'AMA': ('M', 'A'),   # alpha-D-mannopyranose (alt PDB code)
    'NAG': ('Y', 'B'),   # N-acetyl-beta-D-glucosamine
    'NDG': ('Y', 'A'),   # N-acetyl-alpha-D-glucosamine
    'BGL': ('Y', 'B'),   # beta-GlcNAc (alt PDB code)
    'NGA': ('V', 'B'),   # N-acetyl-beta-D-galactosamine
    'A2G': ('V', 'A'),   # N-acetyl-alpha-D-galactosamine
    'FUC': ('f', 'A'),   # alpha-L-fucopyranose
    'FUL': ('f', 'B'),   # beta-L-fucopyranose
    'AFU': ('f', 'A'),   # alpha-L-fucopyranose (alt PDB code)
    'SIA': ('S', 'A'),   # alpha-Neu5Ac (sialic acid)
    'XYS': ('X', 'B'),   # beta-D-xylopyranose
    'XYP': ('X', 'A'),   # alpha-D-xylopyranose
    'RIB': ('R', 'B'),   # beta-D-ribose
    'GCU': ('Z', 'B'),   # beta-D-glucuronic acid
    'IDS': ('U', 'A'),   # alpha-L-iduronic acid
    'RAM': ('h', 'A'),   # alpha-L-rhamnose
}

# Multi-linkage position sets -> GLYCAM linkage code
MULTI_LINKAGE = {
    frozenset({2, 3}): 'Z',
    frozenset({2, 4}): 'Y',
    frozenset({2, 6}): 'X',
    frozenset({3, 4}): 'W',
    frozenset({3, 6}): 'V',
    frozenset({4, 6}): 'U',
    frozenset({2, 3, 4}): 'T',
    frozenset({2, 3, 6}): 'S',
    frozenset({2, 4, 6}): 'R',
    frozenset({3, 4, 6}): 'Q',
    frozenset({2, 3, 4, 6}): 'P',
}

# PDB atom name -> GLYCAM atom name (per PDB residue type)
# Universal hydroxyl H rename: PDB HOx → GLYCAM HxO (applies to ALL sugars)
_HYDROXYL_H_RENAME = {
    'HO1': 'H1O', 'HO2': 'H2O', 'HO3': 'H3O', 'HO4': 'H4O',
    'HO6': 'H6O', 'HO7': 'H7O', 'HO8': 'H8O', 'HO9': 'H9O',
}

# N-acetyl group rename: PDB standard → GLYCAM
_NACETYL_RENAME_PDB = {
    'C7': 'C2N', 'O7': 'O2N', 'C8': 'CME',           # PDB standard (NAG/NGA)
}
# N-acetyl rename for CHARMM-GUI style names
_NACETYL_RENAME_CHARMM = {
    'N': 'N2', 'HN': 'H2N',                            # amide N + H
    'C': 'C2N', 'O': 'O2N',                            # carbonyl C + O
    'CT': 'CME', 'HT1': 'H1M', 'HT2': 'H2M', 'HT3': 'H3M',  # methyl
}

# N-acetyl sugar PDB names (these get both PDB and CHARMM N-acetyl renames)
_NACETYL_SUGARS = {'NAG', 'NDG', 'BGL', 'NGA', 'A2G'}

# Per-residue specific renames (overrides)
GLYCAM_ATOM_MAP = {
    # Sialic acid renames (different from standard sugars)
    'SIA': {'C10': 'C5N', 'C11': 'CME', 'O10': 'O5N'},
}

# Protein residues that can be glycosylated
PROTEIN_TO_GLYCAM = {
    'ASN': ('NLN', 'ND2'),   # N-linked glycosylation
    'SER': ('OLS', 'OG'),    # O-linked
    'THR': ('OLT', 'OG1'),   # O-linked
}

# Sialic acid residues (anomeric carbon is C2, not C1)
SIALIC_ACID_RESIDUES = {'SIA'}


def _parse_pdb(path):
    """Parse PDB file into atoms, residues, and bond graph."""
    with open(path) as f:
        lines = f.readlines()

    atoms = []          # list of atom dicts (preserving order)
    serial_to_idx = {}  # serial -> index in atoms list
    bond_graph = defaultdict(set)  # serial -> set of bonded serials
    link_records = []   # parsed LINK records

    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            serial = int(line[6:11])
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resseq = int(line[22:26])
            icode = line[26] if len(line) > 26 else ' '
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) > 76 else name[0]

            atom = {
                'serial': serial, 'name': name, 'resname': resname,
                'chain': chain, 'resseq': resseq, 'icode': icode,
                'x': x, 'y': y, 'z': z, 'element': element,
                'line': line,
            }
            serial_to_idx[serial] = len(atoms)
            atoms.append(atom)

        elif line.startswith('CONECT'):
            serials = []
            s = line[6:]
            while len(s) >= 5:
                chunk = s[:5].strip()
                if chunk:
                    serials.append(int(chunk))
                s = s[5:]
            if len(serials) >= 2:
                src = serials[0]
                for dst in serials[1:]:
                    bond_graph[src].add(dst)
                    bond_graph[dst].add(src)

        elif line.startswith('LINK  '):
            # LINK records: atom1 (13-26) -- atom2 (43-56)
            try:
                name1 = line[12:16].strip()
                resname1 = line[17:20].strip()
                chain1 = line[21]
                resseq1 = int(line[22:26])
                name2 = line[42:46].strip()
                resname2 = line[47:50].strip()
                chain2 = line[51]
                resseq2 = int(line[52:56])
                link_records.append({
                    'name1': name1, 'resname1': resname1,
                    'chain1': chain1, 'resseq1': resseq1,
                    'name2': name2, 'resname2': resname2,
                    'chain2': chain2, 'resseq2': resseq2,
                })
            except (ValueError, IndexError):
                pass

    # Build residue dict: (chain, resseq) -> list of atom indices
    residues = defaultdict(list)
    for i, atom in enumerate(atoms):
        residues[(atom['chain'], atom['resseq'])].append(i)

    return atoms, residues, bond_graph, link_records


def _detect_glycosidic_bonds(atoms, residues, bond_graph):
    """Detect glycosidic bonds from CONECT records.

    Returns:
        glyco_bonds: list of (parent_reskey, child_reskey, linkage_position)
        protein_links: list of (protein_reskey, sugar_reskey)
    """
    glyco_bonds = []
    protein_links = []

    # Build serial -> atom lookup
    serial_to_atom = {a['serial']: a for a in atoms}

    seen_bonds = set()
    for src_serial, neighbors in bond_graph.items():
        if src_serial not in serial_to_atom:
            continue
        src = serial_to_atom[src_serial]
        src_reskey = (src['chain'], src['resseq'])

        for dst_serial in neighbors:
            if dst_serial not in serial_to_atom:
                continue
            dst = serial_to_atom[dst_serial]
            dst_reskey = (dst['chain'], dst['resseq'])

            if src_reskey == dst_reskey:
                continue  # intra-residue

            bond_key = tuple(sorted([src_serial, dst_serial]))
            if bond_key in seen_bonds:
                continue
            seen_bonds.add(bond_key)

            # Check for protein-sugar bond
            if src['resname'] in PROTEIN_TO_GLYCAM and dst['resname'] in PDB_TO_GLYCAM:
                glycam_name, link_atom = PROTEIN_TO_GLYCAM[src['resname']]
                if src['name'] == link_atom:
                    protein_links.append((src_reskey, dst_reskey))
                    continue
            if dst['resname'] in PROTEIN_TO_GLYCAM and src['resname'] in PDB_TO_GLYCAM:
                glycam_name, link_atom = PROTEIN_TO_GLYCAM[dst['resname']]
                if dst['name'] == link_atom:
                    protein_links.append((dst_reskey, src_reskey))
                    continue

            # Check for sugar-sugar glycosidic bond: C1-O (or C2-O for sialic acid)
            anomeric = None
            oxygen = None
            if _is_anomeric_carbon(src) and dst['name'].startswith('O'):
                anomeric, oxygen = src, dst
            elif _is_anomeric_carbon(dst) and src['name'].startswith('O'):
                anomeric, oxygen = dst, src

            if anomeric and oxygen:
                child_reskey = (anomeric['chain'], anomeric['resseq'])
                parent_reskey = (oxygen['chain'], oxygen['resseq'])
                # Extract linkage position from oxygen name (O4 -> 4)
                o_name = oxygen['name']
                try:
                    pos = int(o_name[1:])
                    glyco_bonds.append((parent_reskey, child_reskey, pos))
                except ValueError:
                    pass

    return glyco_bonds, protein_links


def _detect_glycosidic_bonds_by_distance(atoms, residues):
    """Fallback: detect glycosidic bonds by C1-O distance < 2.0 A."""
    glyco_bonds = []
    protein_links = []

    # Collect anomeric carbons and oxygens per residue
    anomeric_atoms = {}  # reskey -> atom dict for C1 (or C2 for SIA)
    oxygen_atoms = defaultdict(list)  # reskey -> list of O atoms

    for atom in atoms:
        reskey = (atom['chain'], atom['resseq'])
        if _is_anomeric_carbon(atom):
            anomeric_atoms[reskey] = atom
        if atom['name'].startswith('O') and atom['resname'] in PDB_TO_GLYCAM:
            oxygen_atoms[reskey].append(atom)

    # Check each anomeric C against all O atoms on other residues
    for child_key, c_atom in anomeric_atoms.items():
        for parent_key, o_list in oxygen_atoms.items():
            if child_key == parent_key:
                continue
            for o_atom in o_list:
                dx = c_atom['x'] - o_atom['x']
                dy = c_atom['y'] - o_atom['y']
                dz = c_atom['z'] - o_atom['z']
                dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                if dist < 2.0:
                    try:
                        pos = int(o_atom['name'][1:])
                        glyco_bonds.append((parent_key, child_key, pos))
                    except ValueError:
                        pass

    # Check protein-sugar links
    for atom in atoms:
        if atom['resname'] not in PROTEIN_TO_GLYCAM:
            continue
        _, link_atom = PROTEIN_TO_GLYCAM[atom['resname']]
        if atom['name'] != link_atom:
            continue
        prot_key = (atom['chain'], atom['resseq'])
        for sugar_key, c_atom in anomeric_atoms.items():
            dx = atom['x'] - c_atom['x']
            dy = atom['y'] - c_atom['y']
            dz = atom['z'] - c_atom['z']
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < 2.0:
                protein_links.append((prot_key, sugar_key))

    return glyco_bonds, protein_links


def _is_anomeric_carbon(atom):
    """Check if atom is an anomeric carbon (C1, or C2 for sialic acid)."""
    if atom['resname'] in SIALIC_ACID_RESIDUES:
        return atom['name'] == 'C2'
    return atom['name'] == 'C1'


def _determine_linkage_code(positions):
    """Convert set of child linkage positions to GLYCAM linkage code."""
    if not positions:
        return '0'
    if len(positions) == 1:
        return str(next(iter(positions)))
    key = frozenset(positions)
    if key in MULTI_LINKAGE:
        return MULTI_LINKAGE[key]
    # Unknown multi-linkage — use lowest position as fallback
    print(f"  WARNING: Unknown multi-linkage combination {sorted(positions)}, "
          f"using '{min(positions)}'", file=sys.stderr)
    return str(min(positions))


def _rename_atom(pdb_resname, atom_name):
    """Rename a PDB atom name to GLYCAM atom name.

    Applies in order:
    1. Residue-specific overrides (GLYCAM_ATOM_MAP)
    2. N-acetyl renames (PDB standard + CHARMM-GUI style) for NAG/NGA/BGL/NDG/A2G
    3. Universal hydroxyl H rename (HOx → HxO) for all sugars
    """
    # 1. Residue-specific
    specific = GLYCAM_ATOM_MAP.get(pdb_resname, {})
    if atom_name in specific:
        return specific[atom_name]

    # 2. N-acetyl renames
    if pdb_resname in _NACETYL_SUGARS:
        if atom_name in _NACETYL_RENAME_PDB:
            return _NACETYL_RENAME_PDB[atom_name]
        if atom_name in _NACETYL_RENAME_CHARMM:
            return _NACETYL_RENAME_CHARMM[atom_name]

    # 3. Universal hydroxyl H rename
    if atom_name in _HYDROXYL_H_RENAME:
        return _HYDROXYL_H_RENAME[atom_name]

    return atom_name


def _format_atom_line(atom, new_resname, new_serial, new_resseq=None):
    """Format an ATOM/HETATM PDB line with updated names."""
    name = atom['name']
    # PDB atom name formatting: 4-char names start at col 13, shorter at col 14
    if len(name) < 4:
        name_field = f" {name:<3s}"
    else:
        name_field = f"{name:<4s}"

    resseq = new_resseq if new_resseq is not None else atom['resseq']
    icode = atom.get('icode', ' ')
    if icode.strip() == '':
        icode = ' '

    element = atom.get('element', name[0])

    return (
        f"HETATM{new_serial:5d} {name_field} "
        f"{new_resname:>3s} {atom['chain']}{resseq:4d}{icode}"
        f"   {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
        f"                      {element:>2s}  \n"
    )


def convert_to_glycam(input_path, output_path, add_roh=True, verbose=False):
    """Convert PDB glycan nomenclature to GLYCAM naming."""
    atoms, residues, bond_graph, link_records = _parse_pdb(input_path)

    # Detect glycosidic bonds
    if bond_graph:
        glyco_bonds, protein_links = _detect_glycosidic_bonds(
            atoms, residues, bond_graph)
        if verbose:
            print(f"  Detected {len(glyco_bonds)} glycosidic bonds from CONECT")
    else:
        glyco_bonds, protein_links = _detect_glycosidic_bonds_by_distance(
            atoms, residues)
        if verbose:
            print(f"  Detected {len(glyco_bonds)} glycosidic bonds by distance")

    if verbose and protein_links:
        print(f"  Detected {len(protein_links)} protein-sugar links")

    # Build parent->children map and child->parent map
    children_positions = defaultdict(set)  # reskey -> set of linkage positions
    child_to_parent = {}                   # child_reskey -> parent_reskey
    for parent_key, child_key, pos in glyco_bonds:
        children_positions[parent_key].add(pos)
        child_to_parent[child_key] = parent_key

    # Identify sugar residues and protein-linked sugars
    sugar_reskeys = set()
    for reskey, atom_indices in residues.items():
        resname = atoms[atom_indices[0]]['resname']
        if resname in PDB_TO_GLYCAM:
            sugar_reskeys.add(reskey)

    protein_linked_sugars = {sugar_key for _, sugar_key in protein_links}

    # Find reducing-end sugars (roots): sugars with no parent sugar and not protein-linked
    root_sugars = set()
    for reskey in sugar_reskeys:
        if reskey not in child_to_parent and reskey not in protein_linked_sugars:
            root_sugars.add(reskey)

    # Build GLYCAM residue names
    glycam_names = {}  # reskey -> GLYCAM 3-char name
    for reskey in sugar_reskeys:
        atom_idx = residues[reskey][0]
        pdb_resname = atoms[atom_idx]['resname']

        if pdb_resname not in PDB_TO_GLYCAM:
            if verbose:
                print(f"  WARNING: Unknown sugar {pdb_resname} at "
                      f"{reskey[0]}:{reskey[1]}, skipping")
            continue

        sugar_code, anomer_code = PDB_TO_GLYCAM[pdb_resname]
        linkage_code = _determine_linkage_code(children_positions.get(reskey, set()))
        glycam_name = linkage_code + sugar_code + anomer_code
        glycam_names[reskey] = glycam_name

        if verbose:
            positions = children_positions.get(reskey, set())
            pos_str = ','.join(f'O{p}' for p in sorted(positions)) if positions else 'terminal'
            print(f"  {reskey[0]}:{pdb_resname}{reskey[1]} -> {glycam_name} ({pos_str})")

    # Rename atoms in sugar residues
    for reskey in sugar_reskeys:
        for atom_idx in residues[reskey]:
            atom = atoms[atom_idx]
            atom['name'] = _rename_atom(atom['resname'], atom['name'])

    # Handle protein-linked residues
    protein_renames = {}
    for prot_key, sugar_key in protein_links:
        prot_idx = residues[prot_key][0]
        prot_resname = atoms[prot_idx]['resname']
        if prot_resname in PROTEIN_TO_GLYCAM:
            glycam_prot_name, _ = PROTEIN_TO_GLYCAM[prot_resname]
            protein_renames[prot_key] = glycam_prot_name
            if verbose:
                print(f"  {prot_key[0]}:{prot_resname}{prot_key[1]} -> {glycam_prot_name}")

    # Write output PDB
    serial = 0
    out_lines = []
    roh_atoms = []  # (atom_dict, resseq) for ROH cap

    # Determine residue order from input
    reskey_order = []
    seen_reskeys = set()
    for atom in atoms:
        rk = (atom['chain'], atom['resseq'])
        if rk not in seen_reskeys:
            seen_reskeys.add(rk)
            reskey_order.append(rk)

    # Track old serial -> new serial for CONECT remapping
    serial_map = {}
    roh_serial_map = {}  # old serial -> new serial for ROH atoms

    for reskey in reskey_order:
        atom_indices = residues[reskey]
        first_atom = atoms[atom_indices[0]]
        pdb_resname = first_atom['resname']

        # Determine output residue name
        if reskey in glycam_names:
            out_resname = glycam_names[reskey]
        elif reskey in protein_renames:
            out_resname = protein_renames[reskey]
        else:
            out_resname = pdb_resname

        for atom_idx in atom_indices:
            atom = atoms[atom_idx]

            # For reducing-end sugars: extract O1 for ROH cap
            if add_roh and reskey in root_sugars and atom['name'] == 'O1':
                roh_atoms.append((dict(atom), reskey[1]))
                # Don't write O1 in the sugar itself
                continue

            # Skip HD22 on NLN (glycosylated ASN)
            if reskey in protein_renames and protein_renames[reskey] == 'NLN':
                if atom['name'] == 'HD22':
                    continue

            serial += 1
            serial_map[atom['serial']] = serial
            line = _format_atom_line(atom, out_resname, serial)
            out_lines.append(line)

    # Add ROH cap residues at end
    for roh_atom, orig_resseq in roh_atoms:
        serial += 1
        roh_serial_map[roh_atom['serial']] = serial
        serial_map[roh_atom['serial']] = serial
        roh_atom['name'] = 'O1'
        line = _format_atom_line(roh_atom, 'ROH', serial, new_resseq=orig_resseq)
        out_lines.append(line)
        if verbose:
            print(f"  Added ROH cap at {roh_atom['chain']}:{orig_resseq}")

    # Write TER
    out_lines.append('TER\n')

    # Remap and write CONECT records
    for src_serial, neighbors in sorted(bond_graph.items()):
        if src_serial not in serial_map:
            continue
        new_src = serial_map[src_serial]
        new_neighbors = []
        for dst in sorted(neighbors):
            if dst in serial_map:
                new_neighbors.append(serial_map[dst])
        if new_neighbors:
            conect = f"CONECT{new_src:5d}"
            for n in new_neighbors:
                conect += f"{n:5d}"
            out_lines.append(conect + '\n')

    out_lines.append('END\n')

    with open(output_path, 'w') as f:
        f.writelines(out_lines)

    return len(glycam_names), len(protein_renames)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer glycam",
        description="Convert PDB glycan nomenclature to GLYCAM force field naming. "
        "Renames sugar residues to GLYCAM 3-character codes encoding "
        "[linkage][sugar][anomer], renames atoms, and optionally adds ROH cap "
        "at the reducing end.",
    )
    p.add_argument("input", help="Input PDB file with glycan structure")
    p.add_argument("-o", "--output",
                   help="Output PDB file (default: <input>_glycam.pdb)")
    p.add_argument("--no-roh", action="store_true",
                   help="Do not add ROH cap at the reducing end")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print conversion details")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = (Path(args.output) if args.output
                   else input_path.with_stem(input_path.stem + "_glycam"))

    n_sugars, n_protein = convert_to_glycam(
        input_path, output_path,
        add_roh=not args.no_roh,
        verbose=args.verbose,
    )

    print(f"Converted {n_sugars} sugar residue(s)" +
          (f", {n_protein} protein residue(s)" if n_protein else ""))
    print(f"Wrote {output_path}")
