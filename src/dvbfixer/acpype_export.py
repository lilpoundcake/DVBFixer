"""ACPYPE-based GROMACS topology export for AMBER+GLYCAM systems.

Pipeline: PDB -> OpenMM (AMBER14+GLYCAM parametrize) -> ParmEd (prmtop/inpcrd)
-> ACPYPE (GROMACS .top/.gro with per-pair 1-4 parameters via [ pairs_nb ]).

Handles the mixed 1-4 scaling problem: AMBER uses fudgeLJ=0.5/fudgeQQ=0.8333,
GLYCAM uses 1.0/1.0. GROMACS only supports one global value. ACPYPE solves this
using [ pairs_nb ] directive with per-pair LJ/Coulomb parameters.
"""

from pathlib import Path

from dvbfixer.ffutils import PROTEIN_RESIDUES


def detect_ss_bonds(pdb_path):
    """Detect disulfide bonds from CONECT records between SG atoms.

    Returns set of (chain, resseq) for CYS residues involved in SS bonds.
    """
    with open(pdb_path) as f:
        lines = f.readlines()

    # Build serial -> atom info
    serial_to_atom = {}
    sg_serials = set()
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            serial = int(line[6:11])
            chain = line[21]
            resseq = int(line[22:26])
            resname = line[17:20].strip()
            atomname = line[12:16].strip()
            serial_to_atom[serial] = {
                'chain': chain, 'resseq': resseq,
                'resname': resname, 'name': atomname,
            }
            if atomname == 'SG' and resname == 'CYS':
                sg_serials.add(serial)

    # Check CONECT for SG-SG bonds
    ss_residues = set()
    for line in lines:
        if not line.startswith('CONECT'):
            continue
        serials = []
        s = line[6:]
        while len(s) >= 5:
            chunk = s[:5].strip()
            if chunk:
                serials.append(int(chunk))
            s = s[5:]
        if len(serials) >= 2 and serials[0] in sg_serials:
            for s in serials[1:]:
                if s in sg_serials:
                    a1 = serial_to_atom[serials[0]]
                    a2 = serial_to_atom[s]
                    ss_residues.add((a1['chain'], a1['resseq']))
                    ss_residues.add((a2['chain'], a2['resseq']))

    return ss_residues


# GLYCAM sugar residue detection
_GLYCAM_LINKAGE = set('0123456789VWUZXYTSRQPvwuzxytsr')
_GLYCAM_ANOMER = {'A', 'B'}
_PDB_SUGARS = {'NAG', 'NDG', 'BGL', 'BMA', 'MAN', 'GAL', 'BGC', 'GLC',
               'FUC', 'FUL', 'AFU', 'SIA', 'NGA', 'A2G', 'AMA', 'BGA'}


def _is_glycam_sugar(resname):
    if resname in _PDB_SUGARS:
        return True
    if len(resname) == 3 and resname[0] in _GLYCAM_LINKAGE and resname[2] in _GLYCAM_ANOMER:
        return True
    return False


def prepare_for_openmm(pdb_path, temp_path, extra_ss=None, strip_glycam_h=True):
    """Preprocess PDB for OpenMM:
    - CYS->CYX for disulfide bonds (from CONECT + extra_ss), strip HG
    - Strip H from GLYCAM protein residues (NLN/OLS/OLT), re-added by addHydrogens
    - Remove terminal atoms (OXT, H2, H3) from mid-chain residues

    extra_ss: optional set of (chain, resseq) to force CYX renaming on,
              in addition to CONECT-detected ones.
    """
    ss_residues = detect_ss_bonds(pdb_path)
    if extra_ss:
        ss_residues |= extra_ss

    with open(pdb_path) as f:
        lines = f.readlines()

    # Detect which residues have neighbors on both sides (mid-chain)
    atom_lines = [l for l in lines if l.startswith(('ATOM  ', 'HETATM'))]
    residue_order = []
    seen = set()
    for l in atom_lines:
        key = (l[21], int(l[22:26]))
        if key not in seen:
            seen.add(key)
            residue_order.append(key)

    # Build set of (chain, resseq) that have both a predecessor and successor
    midchain = set()
    for i in range(1, len(residue_order) - 1):
        if residue_order[i-1][0] == residue_order[i][0] == residue_order[i+1][0]:
            midchain.add(residue_order[i])

    # AMBER protonation variants → standard names (OpenMM needs standard names
    # in topology, variants passed separately to addHydrogens)
    _AMBER_TO_STD = {
        'ASH': 'ASP', 'GLH': 'GLU',
        'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',
        'LYN': 'LYS',
    }
    # OpenMM variant names for addHydrogens
    _OPENMM_VARIANTS = {'ASH', 'GLH', 'HIE', 'HID', 'HIP', 'CYX', 'LYN'}

    amber_variants = {}  # (chain, resseq) -> variant name
    nln_fix = 0
    terminal_fix = 0
    # Collect lines, separating protein from glycan for reordering
    protein_out = []
    glycan_out = []
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            chain = line[21]
            resseq = int(line[22:26])
            resname = line[17:20].strip()
            atomname = line[12:16].strip()

            if resname == 'CYS' and (chain, resseq) in ss_residues:
                amber_variants[(chain, resseq)] = 'CYX'
                if atomname == 'HG':
                    continue
                line = line[:17] + 'CYX' + line[20:]

            # Capture and rename AMBER protonation variants to standard.
            # Capture AMBER protonation variant names (for later use).
            # Do NOT rename — AMBER14 has templates for ASH/GLH/HIE etc.
            # Renaming would mismatch the H atoms.
            if resname in _AMBER_TO_STD:
                amber_variants[(chain, resseq)] = resname

            # GLYCAM protein residues: strip all H (will be re-added by addHydrogens)
            if strip_glycam_h and resname in ('NLN', 'OLS', 'OLT') and atomname[0] == 'H':
                nln_fix += 1
                continue

            # Remove terminal atoms from mid-chain residues
            if (chain, resseq) in midchain and atomname in ('OXT', 'H2', 'H3'):
                terminal_fix += 1
                continue

            # Rename PDB-standard glycan atom names to GLYCAM convention
            _GLYCAM_ATOM_RENAME = {
                'N': 'N2', 'HN': 'H2N', 'C': 'C2N', 'O': 'O2N',
                'CT': 'CME', 'HT1': 'H1M', 'HT2': 'H2M', 'HT3': 'H3M',
                'HO1': 'H1O', 'HO2': 'H2O', 'HO3': 'H3O', 'HO4': 'H4O',
                'HO6': 'H6O', 'HO7': 'H7O', 'HO8': 'H8O', 'HO9': 'H9O',
            }

            # Sort into protein vs glycan (glycan goes after all protein)
            if _is_glycam_sugar(resname):
                # Rename atom if needed
                new_aname = _GLYCAM_ATOM_RENAME.get(atomname)
                if new_aname:
                    if len(new_aname) < 4:
                        name_field = f' {new_aname:<3s}'
                    else:
                        name_field = f'{new_aname:<4s}'
                    line = line[:12] + name_field + line[16:]
                glycan_out.append(line)
            else:
                protein_out.append(line)
        elif line.startswith('CONECT'):
            # Keep CONECT — needed for glycam intra-residue bonds and ND2→C1
            glycan_out.append(line)
        elif line.startswith('TER'):
            # Skip TER records from glycan chains — they break protein continuity
            if len(line) > 20:
                ter_rn = line[17:20].strip()
                if _is_glycam_sugar(ter_rn):
                    continue
            protein_out.append(line)
        else:
            protein_out.append(line)

    # Add missing HD21 to NLN residues (GLYCAM glycosylated ASN needs HD21 on ND2)
    import numpy as np
    nln_atoms = {}  # (chain, resseq) -> {atomname: (coords, serial)}
    for line in protein_out:
        if line.startswith(('ATOM', 'HETATM')):
            rn = line[17:20].strip()
            if rn == 'NLN':
                an = line[12:16].strip()
                ch = line[21]
                rs = int(line[22:26])
                key = (ch, rs)
                if key not in nln_atoms:
                    nln_atoms[key] = {}
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    nln_atoms[key][an] = (np.array([x, y, z]), int(line[6:11]))
                except ValueError:
                    pass

    hd21_lines = {}
    for key, atoms in nln_atoms.items():
        if 'HD21' not in atoms and 'ND2' in atoms and 'CG' in atoms:
            nd2, nd2_ser = atoms['ND2']
            cg, _ = atoms['CG']
            vec = nd2 - cg
            norm = np.linalg.norm(vec)
            if norm > 0.1:
                hd21_pos = nd2 + (vec / norm) * 1.01
                ch, rs = key
                hd21_lines[key] = (
                    f"ATOM  {(nd2_ser+1) % 100000:5d} HD21 NLN {ch}{rs:4d}    "
                    f"{hd21_pos[0]:8.3f}{hd21_pos[1]:8.3f}{hd21_pos[2]:8.3f}"
                    f"  1.00  0.00           H  \n")

    # Write protein chains first (contiguous), then glycan chains.
    with open(temp_path, 'w') as f:
        for line in protein_out:
            f.write(line)
            # Insert HD21 after ND2 line
            if line.startswith(('ATOM', 'HETATM')):
                rn = line[17:20].strip()
                an = line[12:16].strip()
                ch = line[21]
                rs = int(line[22:26]) if len(line) > 25 else 0
                if rn == 'NLN' and an == 'ND2' and (ch, rs) in hd21_lines:
                    f.write(hd21_lines[(ch, rs)])
        for line in glycan_out:
            f.write(line)

    if ss_residues:
        print(f"  Renamed {len(ss_residues)} CYS -> CYX (disulfide)")
    if amber_variants:
        print(f"  Renamed {len(amber_variants)} AMBER protonation variants to standard")
    if nln_fix:
        print(f"  Stripped {nln_fix} H from GLYCAM protein residues (will re-add)")
    if terminal_fix:
        print(f"  Removed {terminal_fix} spurious terminal atoms (OXT/H2/H3) from mid-chain")

    return temp_path, amber_variants


def add_glycam_bonds(topology, forcefield, verbose=False):
    """Add intra-residue and inter-residue bonds for GLYCAM residues.

    OpenMM PDBFile only infers bonds for standard amino acids. GLYCAM residues
    (NLN, OLS, OLT, sugars) get no intra-residue bonds. This function uses
    the force field templates to add the missing bonds.
    """
    # Standard residues that PDBFile already handles
    standard_res = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLN', 'GLU', 'GLY',
        'HIS', 'HIE', 'HID', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
        'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ACE', 'NME',
    }

    # GLYCAM protein residues that should form peptide bonds
    glycam_protein = {'NLN', 'OLS', 'OLT'}

    # Collect existing bonds for fast lookup
    existing_bonds = set()
    for b in topology.bonds():
        existing_bonds.add((b[0].index, b[1].index))
        existing_bonds.add((b[1].index, b[0].index))

    added_intra = 0
    added_inter = 0

    for chain in topology.chains():
        residues = list(chain.residues())
        for i, res in enumerate(residues):
            if res.name in standard_res:
                continue

            # This is a non-standard residue -- add intra-residue bonds from template
            atom_map = {a.name: a for a in res.atoms()}

            # Try to find matching template
            matched = False
            for tname in [res.name, 'N' + res.name, 'C' + res.name]:
                if tname in forcefield._templates:
                    template = forcefield._templates[tname]
                    for b in template.bonds:
                        a1_name = template.atoms[b[0]].name
                        a2_name = template.atoms[b[1]].name
                        if a1_name in atom_map and a2_name in atom_map:
                            a1 = atom_map[a1_name]
                            a2 = atom_map[a2_name]
                            if (a1.index, a2.index) not in existing_bonds:
                                topology.addBond(a1, a2)
                                existing_bonds.add((a1.index, a2.index))
                                existing_bonds.add((a2.index, a1.index))
                                added_intra += 1
                    matched = True
                    break

            if not matched and verbose:
                print(f"    WARNING: No FF template for {res.name}:{res.id}")

            # Peptide bonds: connect to previous and next residue
            if res.name in glycam_protein:
                # Bond to previous: prev C -> this N
                if i > 0:
                    prev_atoms = {a.name: a for a in residues[i-1].atoms()}
                    if 'C' in prev_atoms and 'N' in atom_map:
                        c = prev_atoms['C']
                        n = atom_map['N']
                        if (c.index, n.index) not in existing_bonds:
                            topology.addBond(c, n)
                            existing_bonds.add((c.index, n.index))
                            existing_bonds.add((n.index, c.index))
                            added_inter += 1

                # Bond to next: this C -> next N
                if i < len(residues) - 1:
                    next_atoms = {a.name: a for a in residues[i+1].atoms()}
                    if 'C' in atom_map and 'N' in next_atoms:
                        c = atom_map['C']
                        n = next_atoms['N']
                        if (c.index, n.index) not in existing_bonds:
                            topology.addBond(c, n)
                            existing_bonds.add((c.index, n.index))
                            existing_bonds.add((n.index, c.index))
                            added_inter += 1

    if added_intra or added_inter:
        print(f"  Added {added_intra} intra-residue + {added_inter} inter-residue bonds for GLYCAM")


_SOLVENT_IONS_BLOCK = """\
; TIP3P water model (AMBER)
[ moleculetype ]
; name  nrexcl
SOL     2

[ atoms ]
;  nr  type  resi  res  atom  cgnr    charge      mass
    1   OW    1    SOL   OW    1    -0.834000    16.00000
    2   HW    1    SOL  HW1    2     0.417000     1.00800
    3   HW    1    SOL  HW2    3     0.417000     1.00800

[ bonds ]
;  ai   aj  funct   r      k
    1    2    1    0.09572  462750.4
    1    3    1    0.09572  462750.4

[ angles ]
;  ai   aj   ak  funct  theta    cth
    2    1    3    1    104.52   836.800

; Ion moleculetypes
[ moleculetype ]
; name  nrexcl
NA      1

[ atoms ]
;  nr  type  resi  res  atom  cgnr    charge      mass
    1   Na+   1    NA    NA    1     1.000000    22.99000

[ moleculetype ]
; name  nrexcl
CL      1

[ atoms ]
;  nr  type  resi  res  atom  cgnr    charge      mass
    1   Cl-   1    CL    CL    1    -1.000000    35.45000

[ moleculetype ]
; name  nrexcl
K       1

[ atoms ]
;  nr  type  resi  res  atom  cgnr    charge      mass
    1   K+    1     K     K    1     1.000000    39.10000

"""

# Ion atomtypes to add to [ atomtypes ] section
_ION_ATOMTYPES = """\
 OW      OW          0.00000  0.00000   A     3.15061e-01   6.36386e-01 ; 1.77  0.1521
 HW      HW          0.00000  0.00000   A     0.00000e+00   0.00000e+00 ; 0.00  0.0000
 Na+     Na+         0.00000  0.00000   A     2.43928e-01   3.65846e-01 ; 1.37  0.0874
 Cl-     Cl-         0.00000  0.00000   A     4.47766e-01   1.48913e-01 ; 2.51  0.0356
 K+      K+          0.00000  0.00000   A     3.03796e-01   8.10369e-01 ; 1.71  0.1937
"""


def _insert_posres_include(top_path, stem):
    """Insert #ifdef POSRES / #include posre / #endif into the .top file.

    The include goes at the end of the main (first) moleculetype section,
    just before the solvent/ion moleculetypes or [ system ].
    """
    with open(top_path) as f:
        lines = f.readlines()

    posres_block = (
        f'\n; Include position restraint file\n'
        f'#ifdef POSRES\n'
        f'#include "posre_{stem}.itp"\n'
        f'#endif\n\n'
    )

    # Find the end of the first moleculetype: look for second [ moleculetype ]
    # or [ system ] — whichever comes first
    mt_count = 0
    insert_at = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == '[ moleculetype ]':
            mt_count += 1
            if mt_count == 2:
                # Insert before the line preceding 2nd moleculetype
                # (skip back over blank/comment lines to insert cleanly)
                insert_at = i
                while insert_at > 0 and lines[insert_at - 1].strip() in ('', ';'):
                    insert_at -= 1
                break
        elif stripped == '[ system ]':
            insert_at = i
            while insert_at > 0 and lines[insert_at - 1].strip() == '':
                insert_at -= 1
            break

    if insert_at is not None:
        lines.insert(insert_at, posres_block)

    with open(top_path, 'w') as f:
        f.writelines(lines)


def _append_solvent_ions(top_path):
    """Append water (TIP3P) and ion moleculetypes to ACPYPE .top file.

    ACPYPE output only contains the solute. gmx solvate/genion add SOL/NA/CL
    to [ molecules ] but need matching moleculetype definitions.
    Inserts atomtypes and moleculetypes before [ system ].
    """
    with open(top_path) as f:
        content = f.read()

    # Add ion/water atomtypes to the [ atomtypes ] section
    # Find end of [ atomtypes ] (next section header)
    lines = content.split('\n')
    in_atomtypes = False
    atomtypes_end = None
    for i, line in enumerate(lines):
        if line.strip() == '[ atomtypes ]':
            in_atomtypes = True
        elif in_atomtypes and line.strip().startswith('['):
            atomtypes_end = i
            break

    if atomtypes_end is not None:
        # Check which atomtypes already exist
        existing = set()
        for line in lines[:atomtypes_end]:
            parts = line.split()
            if parts and not parts[0].startswith(';') and not parts[0].startswith('['):
                existing.add(parts[0])

        new_types = []
        for line in _ION_ATOMTYPES.strip().split('\n'):
            parts = line.split()
            if parts and parts[0] not in existing:
                new_types.append(line)

        if new_types:
            insert = '\n'.join(new_types) + '\n'
            lines.insert(atomtypes_end, insert)

    content = '\n'.join(lines)

    # Insert moleculetypes before [ system ]
    system_idx = content.find('[ system ]')
    if system_idx != -1:
        content = content[:system_idx] + _SOLVENT_IONS_BLOCK + content[system_idx:]

    with open(top_path, 'w') as f:
        f.write(content)


def export_gromacs(pdb_path, output_dir, basename=None, extra_ss=None, verbose=False):
    """Export GROMACS topology files using ACPYPE.

    Pipeline: PDB -> OpenMM (AMBER+GLYCAM parametrize) -> ParmEd (prmtop/inpcrd)
    -> ACPYPE (GROMACS .top/.gro with per-pair 1-4 parameters via [ pairs_nb ]).

    Args:
        pdb_path: Input PDB file
        output_dir: Directory for output files
        basename: Stem for output filenames (default: pdb_path.stem)
        extra_ss: Optional set of (chain, resseq) to force CYX renaming
        verbose: Print detailed output
    """
    from openmm.app import ForceField, Modeller, PDBFile, NoCutoff
    import parmed
    from acpype.topol import MolTopol
    import shutil

    pdb_path = Path(pdb_path)
    output_dir = Path(output_dir)
    stem = basename or pdb_path.stem

    print("\nExporting GROMACS topology via ACPYPE...")

    # Prepare PDB for OpenMM (CYX, GLYCAM bonds, H)
    temp_pdb = pdb_path.parent / '_gmx_temp.pdb'
    _, amber_variants = prepare_for_openmm(pdb_path, temp_pdb, extra_ss=extra_ss,
                                           strip_glycam_h=False)

    pdb = PDBFile(str(temp_pdb))
    topology = pdb.topology
    positions = pdb.positions

    # Restore AMBER variant names that PDBFile normalized (ASH→ASP, GLH→GLU, etc.)
    # For N/C-terminal protonation variants: no NASH/NGLH templates in AMBER14,
    # so strip the protonation H (HD2 for ASH, HE2 for GLH) and use standard name.
    n_terminal_res = set()
    c_terminal_res = set()
    for chain in topology.chains():
        res_list = list(chain.residues())
        if res_list:
            n_terminal_res.add(res_list[0])
            c_terminal_res.add(res_list[-1])

    terminal_h_to_delete = []
    for res in topology.residues():
        key = (res.chain.id, int(res.id))
        orig = amber_variants.get(key)
        if orig:
            if res in n_terminal_res or res in c_terminal_res:
                # Terminal: strip protonation H, keep standard name
                h_name = {'ASH': 'HD2', 'GLH': 'HE2'}.get(orig)
                if h_name:
                    for atom in res.atoms():
                        if atom.name == h_name:
                            terminal_h_to_delete.append(atom)
                    if verbose:
                        print(f"  Terminal {orig} {res.chain.id}:{res.id} → "
                              f"standard {res.name} (stripped {h_name})")
            else:
                # Internal: restore AMBER variant name
                res.name = orig

    if terminal_h_to_delete:
        modeller_tmp = Modeller(topology, positions)
        modeller_tmp.delete(terminal_h_to_delete)
        topology = modeller_tmp.topology
        positions = modeller_tmp.positions

    # Fix missing peptide bonds (C→N distance > 1.9 Å from GLYCAM transplant)
    existing_bonds = set()
    for bond in topology.bonds():
        existing_bonds.add((bond[0].index, bond[1].index))
        existing_bonds.add((bond[1].index, bond[0].index))
    for chain in topology.chains():
        res_list = list(chain.residues())
        for i in range(len(res_list) - 1):
            c_atom = n_atom = None
            for a in res_list[i].atoms():
                if a.name == 'C': c_atom = a
            for a in res_list[i+1].atoms():
                if a.name == 'N': n_atom = a
            if c_atom and n_atom and (c_atom.index, n_atom.index) not in existing_bonds:
                topology.addBond(c_atom, n_atom)
                if verbose:
                    print(f"  Added peptide bond: {res_list[i].name}{res_list[i].id}:C → "
                          f"{res_list[i+1].name}{res_list[i+1].id}:N")

    forcefield = ForceField('amber14-all.xml', 'amber14/GLYCAM_06j-1.xml')
    add_glycam_bonds(topology, forcefield, verbose)

    # Build variants list from captured AMBER protonation names.
    # Skip variants for N/C-terminal residues (no NASH/CGLH in AMBER14).
    _OPENMM_VARIANTS = {'ASH', 'GLH', 'HIE', 'HID', 'HIP', 'CYX', 'LYN'}
    n_terminal = set()
    c_terminal = set()
    for chain in topology.chains():
        res_list = list(chain.residues())
        if res_list:
            n_terminal.add((chain.id, int(res_list[0].id)))
            c_terminal.add((chain.id, int(res_list[-1].id)))

    variants = None
    if amber_variants:
        variants = []
        for res in topology.residues():
            key = (res.chain.id, int(res.id))
            var = amber_variants.get(key)
            if (key in n_terminal or key in c_terminal) and var in _OPENMM_VARIANTS:
                variants.append(None)  # no terminal variant templates
            elif var and var in _OPENMM_VARIANTS:
                variants.append(var)
            else:
                variants.append(None)

    # Always call addHydrogens — adds missing H without disturbing existing ones.
    # glycam-hydrogens.xml provides H definitions for GLYCAM sugar residues.
    Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')
    modeller = Modeller(topology, positions)
    n_before = sum(1 for _ in topology.atoms())
    modeller.addHydrogens(forcefield, variants=variants)
    topology = modeller.topology
    positions = modeller.positions
    n_after = sum(1 for _ in topology.atoms())
    print(f"  After addHydrogens: {n_after} atoms ({n_after - n_before} added)")
    positions = modeller.positions
    print(f"  Parametrized: {sum(1 for _ in topology.atoms())} atoms")

    # Create system WITHOUT constraints (ParmEd needs all bond types)
    # Restore AMBER variant names and build residueTemplates for createSystem.
    n_term_keys = set()
    c_term_keys = set()
    for chain in topology.chains():
        rl = list(chain.residues())
        if rl:
            n_term_keys.add((chain.id, int(rl[0].id)))
            c_term_keys.add((chain.id, int(rl[-1].id)))

    res_templates = {}
    for res in topology.residues():
        key = (res.chain.id, int(res.id))
        orig = amber_variants.get(key)
        if orig and orig in forcefield._templates:
            if key in n_term_keys or key in c_term_keys:
                continue  # keep standard name for terminals
            res.name = orig
            res_templates[res] = orig
        # Force template for GLYCAM sugars (ignoreExternalBonds can cause mismatches)
        elif res.name in forcefield._templates and _is_glycam_sugar(res.name):
            res_templates[res] = res.name

    # ignoreExternalBonds=True: N-linked glycan UYB has C1→ND2 (N atom)
    # but GLYCAM template expects C1→O (glycosidic O). The bond element
    # mismatch is correct chemistry — just not what the template expects.
    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff,
                                     ignoreExternalBonds=True,
                                     residueTemplates=res_templates)

    # ParmEd: OpenMM -> AMBER prmtop/inpcrd
    structure = parmed.openmm.load_topology(topology, system, positions)

    prmtop = pdb_path.parent / f'_{stem}.prmtop'
    inpcrd = pdb_path.parent / f'_{stem}.inpcrd'
    structure.save(str(prmtop), overwrite=True)
    structure.save(str(inpcrd), overwrite=True)
    print(f"  Saved AMBER files: {prmtop.name}, {inpcrd.name}")

    # ACPYPE: AMBER -> GROMACS (handles mixed 1-4 scaling via [ pairs_nb ])
    # ACPYPE creates .amb2gmx dir in CWD, so chdir to pdb_path.parent
    old_cwd = Path.cwd()
    try:
        import os
        os.chdir(pdb_path.parent)
        mol = MolTopol(
            acFileXyz=str(inpcrd),
            acFileTop=str(prmtop),
            amb2gmx=True,
            basename=stem,
        )
        mol.writeGromacsTopolFiles()
    finally:
        os.chdir(old_cwd)

    # Move ACPYPE output to target directory
    acpype_dir = pdb_path.parent / f'{stem}.amb2gmx'
    output_dir.mkdir(parents=True, exist_ok=True)

    gmx_top = acpype_dir / f'{stem}_GMX.top'
    gmx_gro = acpype_dir / f'{stem}_GMX.gro'
    posre = acpype_dir / f'posre_{stem}.itp'

    copied = []
    for src, dst_name in [
        (gmx_top, 'topol.top'),
        (gmx_gro, f'{stem}.gro'),
        (posre, f'posre_{stem}.itp'),
    ]:
        if src.exists():
            dst = output_dir / dst_name
            shutil.copy2(src, dst)
            copied.append(dst_name)

    # Append water and ion moleculetypes before [ system ] so gmx solvate/genion work
    # Insert position restraint include in the main moleculetype
    top_path = output_dir / 'topol.top'
    if top_path.exists():
        _insert_posres_include(top_path, stem)
        _append_solvent_ions(top_path)

    # Cleanup temp files
    temp_pdb.unlink(missing_ok=True)
    prmtop.unlink(missing_ok=True)
    inpcrd.unlink(missing_ok=True)
    shutil.rmtree(acpype_dir, ignore_errors=True)

    print(f"  GROMACS files: {', '.join(copied)} -> {output_dir}/")
    return output_dir
