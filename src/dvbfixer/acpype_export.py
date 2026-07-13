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
    """Detect disulfide bonds from CONECT records between SG atoms, with
    distance-based fallback (SG-SG within 2.5 Å) for inputs missing CONECTs.

    Recognizes SG on both CYS and CYX residues — input may already use either
    name (e.g. crystal structures use CYS, dvbfixer prot output uses CYX).

    Returns set of (chain, resseq) for cysteines involved in SS bonds.
    """
    with open(pdb_path) as f:
        lines = f.readlines()

    # Build serial -> atom info. Capture SG atoms of all cysteine variants
    # (CYS / CYX / CYM) for downstream SS detection.
    serial_to_atom = {}
    sg_atoms = []  # list of dicts with chain/resseq/serial/x/y/z
    _CYS_NAMES = {'CYS', 'CYX', 'CYM'}
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            try:
                serial = int(line[6:11])
            except (ValueError, IndexError):
                continue
            chain = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            resname = line[17:20].strip()
            atomname = line[12:16].strip()
            info = {
                'chain': chain, 'resseq': resseq,
                'resname': resname, 'name': atomname,
                'serial': serial,
            }
            serial_to_atom[serial] = info
            if atomname == 'SG' and resname in _CYS_NAMES:
                try:
                    info['x'] = float(line[30:38])
                    info['y'] = float(line[38:46])
                    info['z'] = float(line[46:54])
                except ValueError:
                    continue
                sg_atoms.append(info)

    sg_serial_set = {a['serial'] for a in sg_atoms}

    # Pass 1: CONECT-based detection
    ss_residues = set()
    for line in lines:
        if not line.startswith('CONECT'):
            continue
        serials = []
        s = line[6:]
        while len(s) >= 5:
            chunk = s[:5].strip()
            if chunk:
                try:
                    serials.append(int(chunk))
                except ValueError:
                    pass
            s = s[5:]
        if len(serials) >= 2 and serials[0] in sg_serial_set:
            for s2 in serials[1:]:
                if s2 in sg_serial_set:
                    a1 = serial_to_atom[serials[0]]
                    a2 = serial_to_atom[s2]
                    ss_residues.add((a1['chain'], a1['resseq']))
                    ss_residues.add((a2['chain'], a2['resseq']))

    # Pass 2: distance-based fallback. Catches inputs without CONECT records.
    # SG-SG bond length is ~2.05 Å; use a generous 2.5 Å cutoff.
    _SS_CUTOFF_A2 = 2.5 ** 2
    for i, a1 in enumerate(sg_atoms):
        if 'x' not in a1:
            continue
        for a2 in sg_atoms[i + 1:]:
            if 'x' not in a2:
                continue
            d2 = ((a1['x'] - a2['x']) ** 2
                  + (a1['y'] - a2['y']) ** 2
                  + (a1['z'] - a2['z']) ** 2)
            if d2 < _SS_CUTOFF_A2:
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


def prepare_for_openmm(pdb_path, temp_path, extra_ss=None, strip_glycam_h=True,
                       prot_overrides=None, keep_all_hydrogens=False):
    """Preprocess PDB for OpenMM:
    - CYS->CYX for disulfide bonds (from CONECT + extra_ss), strip HG
    - Strip H from GLYCAM protein residues (NLN/OLS/OLT), re-added by addHydrogens
    - Remove terminal atoms (OXT, H2, H3) from mid-chain residues
    - Apply protonation overrides: rename residues to ASH/GLH/HIE/HID/HIP/CYX/LYN

    extra_ss: optional set of (chain, resseq) to force CYX renaming on,
              in addition to CONECT-detected ones.
    prot_overrides: optional dict {(chain, resseq): variant_name} to rename
                    residues to AMBER protonation names before OpenMM loads them.
    """
    ss_residues = detect_ss_bonds(pdb_path)
    if extra_ss:
        ss_residues |= extra_ss
    prot_overrides = prot_overrides or {}

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

            # Apply --protonate override: rename to AMBER variant
            override = prot_overrides.get((chain, resseq))
            if override and override in _OPENMM_VARIANTS:
                amber_variants[(chain, resseq)] = override
                # Strip wrong-protonation H (e.g. HE2 when overriding to HID,
                # HD1 when overriding to HIE, HD2 when overriding HIS→HIP).
                # addHydrogens will add the correct H.
                _STRIP_FOR_VARIANT = {
                    'ASH': set(), 'GLH': set(),
                    'HIE': {'HD1'}, 'HID': {'HE2'}, 'HIP': set(),
                    'CYX': {'HG'}, 'LYN': {'HZ3'},
                }
                if atomname in _STRIP_FOR_VARIANT.get(override, set()):
                    continue
                # Keep standard residue name in topology (HIS/ASP/GLU/CYS/LYS);
                # the AMBER variant name is restored later via amber_variants.
                std = _AMBER_TO_STD.get(resname, resname)
                line = line[:17] + f"{std:>3s}" + line[20:]
                # Skip the auto-capture below so override isn't overwritten.

            elif resname == 'CYS' and (chain, resseq) in ss_residues:
                amber_variants[(chain, resseq)] = 'CYX'
                if atomname == 'HG':
                    continue
                line = line[:17] + 'CYX' + line[20:]

            elif resname == 'CYX':
                # Already named CYX (e.g. from `dvbfixer protonate` output).
                # Capture in amber_variants so res_templates can disambiguate
                # against CYM after PDBFile normalizes CYX→CYS internally.
                amber_variants[(chain, resseq)] = 'CYX'
                # Strip HG defensively — CYX must not have it.
                if atomname == 'HG':
                    continue

            elif resname == 'CYM':
                # Explicitly deprotonated cysteine; preserve like CYX.
                amber_variants[(chain, resseq)] = 'CYM'
                if atomname in ('HG',):
                    continue

            elif resname in _AMBER_TO_STD:
                # Capture AMBER protonation variant from input PDB (no override).
                # Rename to standard so PDBFile loads cleanly; restore via amber_variants.
                amber_variants[(chain, resseq)] = resname
                line = line[:17] + f"{_AMBER_TO_STD[resname]:>3s}" + line[20:]

            # GLYCAM protein residues: strip all H (will be re-added by addHydrogens)
            # Only strip H atoms whose names don't match GLYCAM convention.
            # Previously stripped EVERY H on NLN/OLS/OLT (including correctly-
            # placed HD21), which forced addHydrogens to re-place it at a
            # generic position — destroying the canonical cis-OD1 geometry
            # from minimize's _rigid_track_glycan_trees.
            # Keep: H, HA, HB2, HB3, HD21 (NLN); H, HA, HB2, HB3 (OLS);
            #       H, HA, HB, HG21/HG22/HG23 (OLT).
            # Strip: CHARMM-style (HN, HT*); HD22 on NLN (only HD21 in
            #        glycosylated ASN); HG/HG1 on OLS/OLT (replaced by sugar).
            _GLYCAM_KEEP_H = {
                'NLN': {'H', 'HA', 'HB2', 'HB3', 'HD21'},
                'OLS': {'H', 'HA', 'HB2', 'HB3'},
                'OLT': {'H', 'HA', 'HB', 'HG21', 'HG22', 'HG23'},
            }
            # --keep-all-hydrogens skips the GLYCAM H strip: every input H
            # passes through even if it's not on the KEEP allowlist.
            if (not keep_all_hydrogens
                    and strip_glycam_h and resname in _GLYCAM_KEEP_H
                    and atomname[0] == 'H'
                    and atomname not in _GLYCAM_KEEP_H[resname]):
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


def add_glycam_bonds(topology, forcefield, verbose=False, positions=None):
    """Add intra-residue and inter-residue bonds for GLYCAM residues.

    OpenMM PDBFile only infers bonds for standard amino acids. GLYCAM residues
    (NLN, OLS, OLT, sugars) get no intra-residue bonds. This function uses
    the force field templates to add the missing bonds.

    When `positions` is provided, ALSO detects sugar-sugar glycosidic bonds
    by distance (C1/C2 anomeric atom of one sugar within 2.0 Å of an
    O2/O3/O4/O6 atom of another sugar) and adds them. Without these inter-
    sugar bonds, OpenMM's template matching fails on linkage-position sugars
    like 6LB ("missing 1 externally bonded O atom").
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

    # Sugar-sugar glycosidic bond detection (distance-based). Anomeric carbon
    # (C1 for most sugars, C2 for sialic acid) of one sugar within 2.0 Å of
    # a linkage O atom (O2/O3/O4/O6) of another sugar.
    # Required so GLYCAM templates like 6LB / 4YB / VMB which declare an
    # external O bond can match. Without these bonds, addHydrogens and
    # createSystem fail with "missing externally bonded O atom".
    added_sugar = 0
    if positions is not None:
        try:
            from openmm.unit import nanometer
            from dvbfixer.ffutils import is_glycam_sugar, KNOWN_GLYCAN_SMILES
        except Exception:
            is_glycam_sugar = lambda n: False
            KNOWN_GLYCAN_SMILES = {}
            nanometer = None

        def _is_sugar(name):
            return is_glycam_sugar(name) or name in KNOWN_GLYCAN_SMILES

        def _is_sialic(name):
            # PDB sialic acid + GLYCAM sialic (3-char with sugar code S/s)
            if name in ('SIA', 'NAN'):
                return True
            return len(name) == 3 and name[1] in ('S', 's')

        anomeric_names_default = {'C1'}
        linkage_o_names = {'O2', 'O3', 'O4', 'O6'}
        sugar_residues = [r for r in topology.residues() if _is_sugar(r.name)]
        cutoff_nm = 0.20  # 2.0 Å

        # Pre-extract anomeric C and linkage O lists per residue
        anomeric_atoms = []  # list of (atom, position_nm)
        linkage_atoms = []   # list of (atom, position_nm)
        for r in sugar_residues:
            allowed_c = {'C2'} if _is_sialic(r.name) else anomeric_names_default
            for a in r.atoms():
                p = positions[a.index]
                pv = p.value_in_unit(nanometer) if nanometer is not None else p
                if a.name in allowed_c:
                    anomeric_atoms.append((a, pv))
                elif a.name in linkage_o_names:
                    linkage_atoms.append((a, pv))

        for atom_c, c_pos in anomeric_atoms:
            cx, cy, cz = float(c_pos[0]), float(c_pos[1]), float(c_pos[2])
            best = None
            best_d2 = cutoff_nm * cutoff_nm
            for atom_o, o_pos in linkage_atoms:
                if atom_o.residue.index == atom_c.residue.index:
                    continue
                ox, oy, oz = float(o_pos[0]), float(o_pos[1]), float(o_pos[2])
                d2 = (cx - ox) ** 2 + (cy - oy) ** 2 + (cz - oz) ** 2
                if d2 < best_d2:
                    best_d2 = d2
                    best = atom_o
            if best is not None:
                pair = (atom_c.index, best.index)
                rpair = (best.index, atom_c.index)
                if pair not in existing_bonds and rpair not in existing_bonds:
                    topology.addBond(atom_c, best)
                    existing_bonds.add(pair)
                    existing_bonds.add(rpair)
                    added_sugar += 1

    if added_intra or added_inter or added_sugar:
        print(f"  Added {added_intra} intra-residue + {added_inter} inter-residue "
              f"+ {added_sugar} sugar-sugar bonds for GLYCAM")


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


def export_gromacs(pdb_path, output_dir, basename=None, extra_ss=None,
                   prot_overrides=None, verbose=False,
                   keep_all_hydrogens=False):
    """Export GROMACS topology files using ACPYPE.

    Pipeline: PDB -> OpenMM (AMBER+GLYCAM parametrize) -> ParmEd (prmtop/inpcrd)
    -> ACPYPE (GROMACS .top/.gro with per-pair 1-4 parameters via [ pairs_nb ]).

    Args:
        pdb_path: Input PDB file
        output_dir: Directory for output files
        basename: Stem for output filenames (default: pdb_path.stem)
        extra_ss: Optional set of (chain, resseq) to force CYX renaming
        prot_overrides: Optional dict {(chain, resseq): variant} for protonation
                        renames (ASH/GLH/HIE/HID/HIP/CYX/LYN)
        verbose: Print detailed output
        keep_all_hydrogens: If True, skip the _GLYCAM_KEEP_H allowlist filter
                            in prepare_for_openmm — every input H atom passes
                            through even if not on the KEEP list.
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
    # strip_glycam_h=True: NLN/OLS/OLT often have CHARMM-style atom names
    # (HN instead of H) that don't match GLYCAM templates. Strip and let
    # addHydrogens regenerate with correct names.
    temp_pdb = pdb_path.parent / '_gmx_temp.pdb'
    _, amber_variants = prepare_for_openmm(pdb_path, temp_pdb, extra_ss=extra_ss,
                                           strip_glycam_h=True,
                                           prot_overrides=prot_overrides,
                                           keep_all_hydrogens=keep_all_hydrogens)

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
                    import warnings
                    warnings.warn(
                        f"Terminal {orig} {res.chain.id}:{res.id} → {res.name}: "
                        f"AMBER14 has no N/C-terminal protonated template "
                        f"(NASH/NGLH). Using standard {res.name} (stripped {h_name})."
                    )
            # NOTE: Do NOT rename to AMBER variant here — addHydrogens uses
            # hydrogens.xml which only has standard residue names (HIS, ASP,
            # GLU, CYS, LYS) with variant attributes. Renaming to HIE/HID/etc.
            # makes addHydrogens unable to find H definitions. Variant info is
            # passed via the variants list. Renaming happens later for
            # createSystem (line ~637).

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

    # tip3pfb.xml included for ion templates (Ca2+, Mg2+, Zn2+, Na+, Cl-,
    # K+, ...). AMBER ships ion params inside the water-model XML, so
    # without it any structure containing ions fails template matching.
    forcefield = ForceField('amber14-all.xml', 'amber14/GLYCAM_06j-1.xml',
                             'amber14/tip3pfb.xml')
    # Pass positions so sugar-sugar glycosidic bonds are detected by distance.
    # Without these bonds, GLYCAM templates for linkage-position sugars (e.g.
    # 6LB declares O6 must be externally bonded) fail to match and addHydrogens
    # raises "No template found for residue X. The atoms and bonds match X,
    # but the set of externally bonded atoms is missing 1 O atom."
    add_glycam_bonds(topology, forcefield, verbose, positions=positions)

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

    # Variants without terminal templates in AMBER14 (must drop at termini)
    _NO_TERMINAL_VARIANT = {'ASH', 'GLH'}

    variants = None
    if amber_variants:
        variants = []
        for res in topology.residues():
            key = (res.chain.id, int(res.id))
            var = amber_variants.get(key)
            if (key in n_terminal or key in c_terminal) and var in _NO_TERMINAL_VARIANT:
                variants.append(None)  # no NASH/NGLH/CASH/CGLH templates
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
        if orig:
            is_terminal = key in n_term_keys or key in c_term_keys
            # ASH/GLH have no terminal templates — keep standard name
            if is_terminal and orig in _NO_TERMINAL_VARIANT:
                continue
            # CYX/HIE/HID/HIP/LYN: AMBER14 has terminal versions (NCYX, CCYX,
            # NHIE, CHIE, etc.). Use them at terminals.
            if is_terminal:
                term_prefix = 'N' if key in n_term_keys else 'C'
                term_name = term_prefix + orig
                if term_name in forcefield._templates:
                    res.name = orig
                    res_templates[res] = term_name
                    continue
            if orig in forcefield._templates:
                res.name = orig
                res_templates[res] = orig
        # Defensive disambiguation for any remaining CYS not yet templated.
        # After addHydrogens, an isolated CYS might still match BOTH CYM and
        # CYX templates (e.g. when the protein-glycan workflow stripped HG on
        # an SS bond we missed). Force the right template:
        #   - HG present → protonated CYS (template CYS)
        #   - HG absent → CYX (assume disulfide; safer than CYM which is
        #     charged and rare). For CYM, use --protonate flag explicitly.
        elif res.name == 'CYS' and res not in res_templates:
            atom_names = {a.name for a in res.atoms()}
            is_terminal = key in n_term_keys or key in c_term_keys
            if 'HG' in atom_names:
                term_name = ('N' if key in n_term_keys else 'C') + 'CYS'
                if is_terminal and term_name in forcefield._templates:
                    res_templates[res] = term_name
                else:
                    res_templates[res] = 'CYS'
            else:
                # Missing HG and not yet flagged as CYX — assume disulfide.
                if is_terminal:
                    term_name = ('N' if key in n_term_keys else 'C') + 'CYX'
                    if term_name in forcefield._templates:
                        res.name = 'CYX'
                        res_templates[res] = term_name
                        continue
                if 'CYX' in forcefield._templates:
                    res.name = 'CYX'
                    res_templates[res] = 'CYX'
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
