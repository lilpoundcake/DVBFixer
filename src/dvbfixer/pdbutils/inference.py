"""Automatic CONECT-record inference.

Split out of the flat ``pdbutils.py`` in Phase 1.4. Contains the
OpenBabel-based bond perception, distance-based fallback, filtering
against standard-AA templates, and the three "domain override" rules
(SS bonds, glycosidic linkages, N/O glycosylation) that give dvbfixer
its CONECT idempotency.

Line-level PDB I/O helpers (``build_serial_map`` etc.) live in
:mod:`dvbfixer.pdbutils.io`. Both are re-exported from the
:mod:`dvbfixer.pdbutils` package init so ``from dvbfixer.pdbutils
import X`` keeps working across the split.
"""

import atexit
import os
import tempfile

# Sugar / glycoprotein-attachment knowledge used by the inference domain
# overrides. Imported lazily where needed to avoid circular imports.
#
# The PDB (3-char) and CHARMM-GUI (4-char) sugar name sets used to be
# maintained here as separate, independently-drifted copies of the ones in
# `dvbfixer.ffutils` — e.g. `BNE5AC` (beta-anomer sialic acid) was missing
# from this file's own CHARMM-GUI set even though `ffutils` already had it.
# Both name sets now live solely in `ffutils`; `_is_sugar_resname` below
# delegates to `ffutils.is_pdb_sugar_resname`, the canonical check.

# Anomeric carbon atom names. Sialic acid uses C2 (it's a ketose);
# all aldose hexoses/pentoses use C1.
_ANOMERIC_C_HEXOSE = 'C1'
_ANOMERIC_C_SIALIC = 'C2'
_SIALIC_RESNAMES = {'SIA', 'NAN', 'ANE5', 'BNE5', 'ANE5AC', 'BNE5AC'}

# Ring oxygens that participate in glycosidic linkages (acceptor side).
_GLYCOSIDIC_LINKAGE_OXYGENS = {'O2', 'O3', 'O4', 'O6'}

# Protein attachment atoms for N- and O-glycosylation (donor side).
_GLYCOSYLATION_DONOR_ATOMS = {
    'ASN': 'ND2', 'NLN': 'ND2',
    'SER': 'OG',  'OLS': 'OG',
    'THR': 'OG1', 'OLT': 'OG1',
}

# CYS family for SS bond detection.
_CYS_RESNAMES = {'CYS', 'CYX', 'CYM'}

# Cutoffs used by the domain overrides (mirror existing scattered fallbacks).
_SS_CUTOFF_A = 2.5
_GLYCOSIDIC_CUTOFF_A = 2.0
_GLYCOSYLATION_CUTOFF_A = 2.5


def _is_sugar_resname(resname):
    """True if resname is a recognised sugar (PDB, CHARMM-GUI, or GLYCAM)."""
    from dvbfixer.ffutils import is_pdb_sugar_resname
    return is_pdb_sugar_resname(resname)


def _anomeric_carbon(resname):
    """Return the anomeric carbon atom name for a sugar resname."""
    if resname in _SIALIC_RESNAMES:
        return _ANOMERIC_C_SIALIC
    # GLYCAM sialic codes match `?S?` (e.g. 0SA, 4SB)
    if len(resname) == 3 and resname[1] in ('S', 's'):
        return _ANOMERIC_C_SIALIC
    return _ANOMERIC_C_HEXOSE


def _has_any_conect(pdb_path):
    """Quick scan: does the file contain at least one CONECT record?"""
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('CONECT'):
                return True
    return False


def _parse_existing_conect(pdb_path):
    """Parse CONECT records into a set of canonical (min, max) serial pairs."""
    pairs = set()
    with open(pdb_path) as f:
        for line in f:
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
            if len(serials) < 2:
                continue
            src = serials[0]
            for tgt in serials[1:]:
                if src == tgt:
                    continue
                a, b = (src, tgt) if src < tgt else (tgt, src)
                pairs.add((a, b))
    return pairs


def _load_atom_table(pdb_path):
    """Parse all ATOM/HETATM lines once. Returns:

      atoms: list of dicts, in file order, with keys
        serial, name, altloc, resname, chain, resid, icode, x, y, z, element, is_het
      by_serial: {serial: dict}
    """
    atoms = []
    by_serial = {}
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            try:
                serial = int(line[6:11])
            except ValueError:
                continue
            name = line[12:16].strip()
            altloc = line[16:17].strip()
            resname = line[17:21].strip()  # 4-char tolerant (CHARMM-GUI / AMBER variants)
            if len(resname) > 3 and resname[3] == ' ':
                resname = resname[:3]
            chain = line[21:22].strip()
            try:
                resid = int(line[22:26])
            except ValueError:
                resid = 0
            icode = line[26:27].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            element = line[76:78].strip()
            if not element:
                # Fall back to atom-name heuristic: first letter, skip leading digit.
                en = name[0].upper() if name and not name[0].isdigit() else (
                    name[1].upper() if len(name) > 1 else 'X')
                element = en
            element = element[0].upper() + element[1:].lower()
            atom = {
                'serial': serial, 'name': name, 'altloc': altloc,
                'resname': resname, 'chain': chain, 'resid': resid,
                'icode': icode, 'x': x, 'y': y, 'z': z,
                'element': element, 'is_het': line.startswith('HETATM'),
            }
            atoms.append(atom)
            by_serial[serial] = atom
    return atoms, by_serial


def _residue_key(atom):
    """Composite key uniquely identifying a residue."""
    return (atom['chain'], atom['resid'], atom['icode'])


def _detect_coarse_grained(atoms):
    """Heuristic: ≥80% C-only beads with no H/N/O → likely coarse-grained."""
    if not atoms:
        return False
    n = len(atoms)
    n_c = sum(1 for a in atoms if a['element'] == 'C')
    n_hno = sum(1 for a in atoms if a['element'] in ('H', 'N', 'O'))
    return (n_c / n >= 0.8) and (n_hno / n < 0.05)


def _openbabel_bonds(pdb_path):
    """Run OpenBabel ConnectTheDots and return raw bond serials.

    Returns set of canonical (min, max) tuples. Returns None on import or
    load failure (caller falls back to element-aware cutoffs).
    """
    try:
        from openbabel import openbabel as ob
    except ImportError:
        return None
    # OpenBabel emits a "Failed to kekulize" warning while PerceiveBondOrders
    # runs on protein/glycan ring systems. We don't use bond orders so this is
    # noise — silence the global error log level for the duration of the call.
    ob.obErrorLog.SetOutputLevel(0)
    obconv = ob.OBConversion()
    obconv.SetInFormat('pdb')
    mol = ob.OBMol()
    if not obconv.ReadFile(mol, str(pdb_path)):
        return None
    # NOTE: ReadFile already runs ConnectTheDots+PerceiveBondOrders for PDB
    # input. Do NOT call mol.ConnectTheDots() again — it wipes the PDB
    # serial mapping (GetSerialNum returns 0 for all atoms afterwards).
    bonds = set()
    for bond in ob.OBMolBondIter(mol):
        # OpenBabel atom indices are 1-based and match PDB serials when the
        # input PDB has no serial gaps. Map via atom iter to be safe.
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        # ob.OBAtom.GetResidue().GetIdx() gives 0-based residue idx; we need
        # the original PDB serials from the atom's PDBSerial.
        try:
            sa = a.GetResidue().GetSerialNum(a)
        except Exception:
            sa = a.GetIdx()
        try:
            sb = b.GetResidue().GetSerialNum(b)
        except Exception:
            sb = b.GetIdx()
        if sa == sb:
            continue
        lo, hi = (sa, sb) if sa < sb else (sb, sa)
        bonds.add((lo, hi))
    return bonds


def _fallback_element_cutoffs():
    """Element-aware bond cutoffs (Å). Used if OpenBabel unavailable."""
    return {
        ('C', 'C'): 1.7, ('C', 'N'): 1.55, ('C', 'O'): 1.55, ('C', 'S'): 1.85,
        ('N', 'N'): 1.5, ('N', 'O'): 1.45, ('O', 'P'): 1.7, ('S', 'S'): 2.5,
        ('C', 'H'): 1.2, ('N', 'H'): 1.15, ('O', 'H'): 1.1, ('S', 'H'): 1.4,
        ('C', 'F'): 1.45, ('C', 'Cl'): 1.85, ('P', 'O'): 1.7,
        ('Fe', 'N'): 2.5, ('Fe', 'O'): 2.5, ('Fe', 'S'): 2.7,
        ('Zn', 'N'): 2.3, ('Zn', 'O'): 2.3, ('Zn', 'S'): 2.5,
    }


def _fallback_distance_bonds(atoms, by_serial):
    """Distance-based bond inference fallback (no OpenBabel)."""
    try:
        import numpy as np
        from scipy.spatial import cKDTree
    except ImportError:
        np = None
        cKDTree = None

    cutoffs = _fallback_element_cutoffs()
    bonds = set()
    n = len(atoms)

    if cKDTree is not None and n >= 2000:
        coords = np.array([[a['x'], a['y'], a['z']] for a in atoms],
                          dtype=np.float32)
        tree = cKDTree(coords)
        # Largest cutoff is around 2.7 Å (Fe-S); query at 3.0 to be safe.
        pairs = tree.query_pairs(r=3.0)
        for i, j in pairs:
            ai, aj = atoms[i], atoms[j]
            ek = tuple(sorted((ai['element'], aj['element'])))
            cutoff = cutoffs.get(ek, 1.7)
            d = ((ai['x']-aj['x'])**2 + (ai['y']-aj['y'])**2 +
                 (ai['z']-aj['z'])**2) ** 0.5
            if 0.4 <= d <= cutoff:
                lo, hi = (ai['serial'], aj['serial']) if ai['serial'] < aj['serial'] else (aj['serial'], ai['serial'])
                bonds.add((lo, hi))
    else:
        # O(n²) — acceptable for small files
        for i in range(n):
            ai = atoms[i]
            for j in range(i + 1, n):
                aj = atoms[j]
                dx = ai['x'] - aj['x']
                if abs(dx) > 3.0:
                    continue
                dy = ai['y'] - aj['y']
                if abs(dy) > 3.0:
                    continue
                dz = ai['z'] - aj['z']
                if abs(dz) > 3.0:
                    continue
                d2 = dx*dx + dy*dy + dz*dz
                if d2 > 9.0 or d2 < 0.16:  # 3.0² and 0.4²
                    continue
                ek = tuple(sorted((ai['element'], aj['element'])))
                cutoff = cutoffs.get(ek, 1.7)
                if d2 <= cutoff * cutoff:
                    lo, hi = (ai['serial'], aj['serial']) if ai['serial'] < aj['serial'] else (aj['serial'], ai['serial'])
                    bonds.add((lo, hi))
    return bonds


def _apply_filter(bonds, atoms, by_serial):
    """Drop bonds we don't want to emit.

    Drops:
      - bonds where BOTH atoms are in the same standard amino-acid residue
        (FF templates own intra-AA chemistry)
      - bonds involving any water atom
      - bonds between two DIFFERENT standard-AA residues that AREN'T a
        canonical peptide backbone C-N or disulfide SG-SG. OpenBabel's
        ConnectTheDots occasionally infers proximity-based false-positive
        bonds (e.g. ARG NH1 → some Ser OG, at ~1.7 Å H-bond distance,
        looks like a bond). These slip past OpenMM's template matcher
        as "extra external X atom" errors, which are confusing because
        the residue itself is fine.
    Keeps:
      - any bond with at least one HETATM in a non-AA residue
      - SS bonds (CYS family SG-SG) regardless of residue classification
      - real peptide backbone C-N bonds between adjacent protein residues
      - intra-residue bonds when the residue name is an AMBER protonation
        variant (LYN / ASH / GLH / CYX / CYM / HID / HIE / HIP). OpenMM's
        ``PDBFile._standardResidues`` set only covers the 20 canonical
        AAs, so PDBFile does NOT infer intra-residue bonds for these
        variant names. Without CONECT records, the loaded topology has
        zero bonds for the variant residues and every downstream
        template-match step fails ("residue X has no bonds").
      - same-residue bonds on a sugar (GLYCAM or plain PDB name) ONLY if
        within a real covalent distance (~1.7 Å). OpenBabel's
        ``ConnectTheDots`` has no template for GLYCAM/PDB sugar residues
        (unlike the 20 canonical AAs) and can propose same-residue
        "bonds" at 2.5-7+ Å (ring/branch atoms that are merely close in
        3D, not bonded) — these were previously kept unconditionally,
        producing chemically-impossible CONECT entries in the output.
    """
    from dvbfixer.ffutils import PROTEIN_RESIDUES, is_pdb_sugar_resname
    from dvbfixer.ffutils.variants import ALL_VARIANTS
    _WATER = {'HOH', 'WAT', 'TIP3', 'TIP4', 'TIP5', 'SOL', 'SPC', 'SPCE'}
    _SUGAR_COVALENT_CUTOFF2 = 1.7 * 1.7  # Å², same default as the fallback distance-bonder
    # OpenMM's PDBFile parser recognises these — no need to emit CONECT
    # for their intra-residue bonds (FF templates own the chemistry AND
    # the parser knows the names).
    _OPENMM_PARSER_STANDARD_AA = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
        'TYR', 'VAL',
    }
    out = set()
    # File-order neighbours define the only ordinary peptide links we allow.
    # A distance bonder can otherwise connect the N of one residue to the C
    # of any spatially close residue, producing OpenMM's misleading "1 C atom
    # too many" template error during addHydrogens.
    protein_order = []
    seen_residues = set()
    for atom in atoms:
        key = _residue_key(atom)
        if atom['resname'] in PROTEIN_RESIDUES and key not in seen_residues:
            seen_residues.add(key)
            protein_order.append(key)
    peptide_pairs = {
        (left, right)
        for left, right in zip(protein_order, protein_order[1:])
        if left[0] == right[0]
    }
    for lo, hi in bonds:
        a = by_serial.get(lo)
        b = by_serial.get(hi)
        if a is None or b is None:
            continue
        if a['resname'] in _WATER or b['resname'] in _WATER:
            continue
        same_res = _residue_key(a) == _residue_key(b)
        is_ss = (a['name'] == 'SG' and b['name'] == 'SG'
                 and a['resname'] in _CYS_RESNAMES
                 and b['resname'] in _CYS_RESNAMES)
        if is_ss:
            out.add((lo, hi))
            continue
        # Intra-residue bond: keep it if the residue name is an AMBER
        # variant OpenMM's PDBFile can't infer bonds for; otherwise drop
        # (standard AAs are template-owned).
        if same_res and a['resname'] in _OPENMM_PARSER_STANDARD_AA:
            continue
        if same_res and a['resname'] in ALL_VARIANTS:
            out.add((lo, hi))
            continue
        if same_res and a['resname'] in PROTEIN_RESIDUES:
            continue
        if same_res and is_pdb_sugar_resname(a['resname']):
            dx = a['x'] - b['x']
            dy = a['y'] - b['y']
            dz = a['z'] - b['z']
            if dx * dx + dy * dy + dz * dz > _SUGAR_COVALENT_CUTOFF2:
                continue
            out.add((lo, hi))
            continue
        # Inter-residue where BOTH are standard-AA: keep ONLY the
        # canonical peptide-backbone C-N bond. Everything else (N-O,
        # N-N, S-N, sidechain-sidechain contacts) is spurious.
        if (not same_res
                and a['resname'] in PROTEIN_RESIDUES
                and b['resname'] in PROTEIN_RESIDUES):
            names = {a['name'], b['name']}
            if names != {'C', 'N'}:
                continue
            c_atom, n_atom = (a, b) if a['name'] == 'C' else (b, a)
            if (_residue_key(c_atom), _residue_key(n_atom)) not in peptide_pairs:
                continue
        out.add((lo, hi))
    return out


def _domain_overrides(atoms):
    """Force-add the three known glycoprotein bond patterns.

    Returns set of canonical (min, max) serial pairs.
    """
    extras = set()
    # Index by residue and by element for fast lookup.
    by_res = {}   # res_key -> list of atom dicts
    for a in atoms:
        by_res.setdefault(_residue_key(a), []).append(a)

    # --- SS bonds ---
    # Proper 1:1 nearest-neighbor matching, not naive all-pairs-within-
    # cutoff: when several CYS SG atoms cluster close together (e.g. two
    # chain copies whose N-termini pack near each other), a plain cutoff
    # scan force-adds every pair within range, giving one SG atom two or
    # three "partners". (This is belt-and-suspenders — the final dedup
    # in infer_conect_records() also enforces this across every bond
    # source, including OpenBabel's own SG-SG detection.)
    sg_atoms = [a for a in atoms
                if a['name'] == 'SG' and a['resname'] in _CYS_RESNAMES]
    candidates = []
    for i, a in enumerate(sg_atoms):
        for b in sg_atoms[i + 1:]:
            if _residue_key(a) == _residue_key(b):
                continue
            d2 = ((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2)
            if d2 <= _SS_CUTOFF_A ** 2:
                candidates.append((d2, a['serial'], b['serial']))
    candidates.sort()
    matched_sg = set()
    for _d2, s1, s2 in candidates:
        if s1 in matched_sg or s2 in matched_sg:
            continue
        matched_sg.add(s1)
        matched_sg.add(s2)
        lo, hi = (s1, s2) if s1 < s2 else (s2, s1)
        extras.add((lo, hi))

    # --- Glycosidic linkages: anomeric C of sugar i → O on sugar j ---
    sugar_residues = {}  # res_key -> (resname, atom_dict_by_name)
    for rk, ratoms in by_res.items():
        rname = ratoms[0]['resname']
        if _is_sugar_resname(rname):
            sugar_residues[rk] = {a['name']: a for a in ratoms}
    sugar_keys = list(sugar_residues.keys())
    for i, rk_i in enumerate(sugar_keys):
        rname_i = by_res[rk_i][0]['resname']
        anomeric_name = _anomeric_carbon(rname_i)
        ac = sugar_residues[rk_i].get(anomeric_name)
        if ac is None:
            continue
        for rk_j in sugar_keys:
            if rk_j == rk_i:
                continue
            for o_name in _GLYCOSIDIC_LINKAGE_OXYGENS:
                ot = sugar_residues[rk_j].get(o_name)
                if ot is None:
                    continue
                d2 = ((ac['x']-ot['x'])**2 + (ac['y']-ot['y'])**2 +
                      (ac['z']-ot['z'])**2)
                if d2 <= _GLYCOSIDIC_CUTOFF_A ** 2:
                    lo, hi = (ac['serial'], ot['serial']) if ac['serial'] < ot['serial'] else (ot['serial'], ac['serial'])
                    extras.add((lo, hi))

    # --- Glycosylation: protein ND2/OG/OG1 → sugar anomeric C ---
    for rk, ratoms in by_res.items():
        rname = ratoms[0]['resname']
        donor_atom_name = _GLYCOSYLATION_DONOR_ATOMS.get(rname)
        if not donor_atom_name:
            continue
        donor = next((a for a in ratoms if a['name'] == donor_atom_name), None)
        if donor is None:
            continue
        for rk_s, sugar_atoms in sugar_residues.items():
            anomeric_name = _anomeric_carbon(by_res[rk_s][0]['resname'])
            ac = sugar_atoms.get(anomeric_name)
            if ac is None:
                continue
            d2 = ((donor['x']-ac['x'])**2 + (donor['y']-ac['y'])**2 +
                  (donor['z']-ac['z'])**2)
            if d2 <= _GLYCOSYLATION_CUTOFF_A ** 2:
                lo, hi = (donor['serial'], ac['serial']) if donor['serial'] < ac['serial'] else (ac['serial'], donor['serial'])
                extras.add((lo, hi))

    return extras


def infer_conect_records(pdb_path, *, preserve_existing=True,
                         include_protein_backbone=False, verbose=False):
    """Infer connectivity for a PDB and return canonical sorted (s1, s2) pairs.

    Algorithm:
      1. Parse existing CONECT (kept if preserve_existing=True).
      2. Run OpenBabel ConnectTheDots; fall back to element-aware cutoffs
         + scipy.spatial.cKDTree if OpenBabel is unavailable.
      3. Filter: drop intra-AA bonds (FF-template owned), water bonds.
      4. Domain overrides: force-add SS, glycosidic, glycosylation bonds.
      5. Drop bonds referencing serials that no longer exist.
      6. Return union, deduplicated, sorted.
    """
    atoms, by_serial = _load_atom_table(pdb_path)
    if not atoms:
        return []

    if _detect_coarse_grained(atoms):
        if verbose:
            print(f"  [conect] {pdb_path}: looks coarse-grained; "
                  f"skipping inference, passing existing CONECT through.")
        existing = _parse_existing_conect(pdb_path) if preserve_existing else set()
        return sorted(existing)

    existing = _parse_existing_conect(pdb_path) if preserve_existing else set()
    # Sanitize inherited CONECT too: a false bond inferred by an earlier
    # pipeline stage must not become permanently trusted on the next pass.
    existing = _apply_filter(existing, atoms, by_serial)

    inferred = _openbabel_bonds(pdb_path)
    used_fallback = False
    if inferred is None:
        if verbose:
            print("  [conect] OpenBabel unavailable or load failed; "
                  "falling back to element-aware distance cutoffs.")
        inferred = _fallback_distance_bonds(atoms, by_serial)
        used_fallback = True

    inferred = _apply_filter(inferred, atoms, by_serial)
    inferred |= _domain_overrides(atoms)

    # Union with existing, drop stale serials.
    union = (existing | inferred)
    valid = {(a, b) for (a, b) in union if a in by_serial and b in by_serial}
    valid, _n_dropped_ss = _dedupe_ss_bonds(valid, existing, by_serial)

    if verbose:
        print(f"  [conect] {pdb_path}: "
              f"existing={len(existing)}, "
              f"inferred={len(inferred)}, "
              f"total={len(valid)}"
              f"{' (fallback)' if used_fallback else ''}"
              f"{f', dropped {_n_dropped_ss} spurious SG-SG' if _n_dropped_ss else ''}")

    return sorted(valid)


def _dedupe_ss_bonds(pairs, existing, by_serial):
    """Enforce at most one SG-SG partner per CYS-family SG atom.

    Every one of OpenBabel's ConnectTheDots, the element-aware distance
    fallback, and ``_domain_overrides`` can independently flag an SG-SG
    contact — none of them is aware of what the *other* two already
    decided. When several CYS SG atoms sit close together (e.g. two
    chain copies whose N-termini pack near each other, or a genuinely
    clashing test structure), the union ends up with one SG bonded to
    two or three different partners, and OpenMM's forcefield matcher
    then fails with "N S atom(s) too many" for the CYX template — a
    disulfide SG must have exactly one external S bond.

    Precedence: pairs already present in the file's own (preserved)
    CONECT records win outright (trust whatever upstream step —
    typically ``prepare``'s PROPKA+Reduce pass — already established
    the correct pairing for). Remaining candidates are resolved by
    greedy nearest-distance 1:1 matching.

    Returns ``(deduped_pairs, n_dropped)``.
    """
    ss_pairs = []
    other_pairs = []
    for a, b in pairs:
        aa = by_serial.get(a)
        bb = by_serial.get(b)
        if (aa is not None and bb is not None
                and aa['name'] == 'SG' and bb['name'] == 'SG'
                and aa['resname'] in _CYS_RESNAMES and bb['resname'] in _CYS_RESNAMES):
            ss_pairs.append((a, b))
        else:
            other_pairs.append((a, b))
    if len(ss_pairs) < 2:
        return pairs, 0

    locked_serial = set()
    kept_ss = []
    # Existing (preserved) pairs win outright.
    for a, b in ss_pairs:
        if (a, b) in existing or (b, a) in existing:
            kept_ss.append((a, b))
            locked_serial.add(a)
            locked_serial.add(b)
    # Remaining candidates: greedy nearest-distance 1:1 matching among
    # SG atoms not already locked by an existing pair.
    remaining = [(a, b) for (a, b) in ss_pairs if (a, b) not in kept_ss]
    scored = sorted(
        (((by_serial[a]['x']-by_serial[b]['x'])**2
          + (by_serial[a]['y']-by_serial[b]['y'])**2
          + (by_serial[a]['z']-by_serial[b]['z'])**2), a, b)
        for a, b in remaining
    )
    for _d2, a, b in scored:
        if a in locked_serial or b in locked_serial:
            continue
        kept_ss.append((a, b))
        locked_serial.add(a)
        locked_serial.add(b)
    return set(other_pairs) | set(kept_ss), len(ss_pairs) - len(kept_ss)


def write_conect_block(bonds):
    """Format a list of (s1, s2) pairs as PDB CONECT lines (≤4 partners each)."""
    if not bonds:
        return []
    # Group partners per source serial.
    partners = {}
    for a, b in bonds:
        partners.setdefault(a, set()).add(b)
        partners.setdefault(b, set()).add(a)
    lines = []
    for src in sorted(partners):
        ps = sorted(partners[src])
        for i in range(0, len(ps), 4):
            chunk = ps[i:i + 4]
            line = f"CONECT{src:5d}"
            for p in chunk:
                line += f"{p:5d}"
            lines.append(line.ljust(80) + "\n")
    return lines


def _strip_existing_conect(lines):
    return [ln for ln in lines if not ln.startswith('CONECT')]


def _materialise_inferred_pdb(pdb_path, *, verbose=False):
    """Write a temp PDB copy with merged CONECT records and return its path.

    Registers atexit cleanup of the temp file. Never modifies pdb_path.
    """
    bonds = infer_conect_records(pdb_path, verbose=verbose)
    with open(pdb_path) as f:
        lines = f.readlines()
    stripped = _strip_existing_conect(lines)
    new_conect = write_conect_block(bonds)
    # Insert new CONECT block just before END/ENDMDL, or at end if neither.
    out_lines = []
    inserted = False
    for ln in stripped:
        if not inserted and ln.startswith(('END', 'ENDMDL')):
            out_lines.extend(new_conect)
            inserted = True
        out_lines.append(ln)
    if not inserted:
        out_lines.extend(new_conect)

    fd, tmp_path = tempfile.mkstemp(suffix='.pdb', prefix='dvbfixer_conect_')
    with os.fdopen(fd, 'w') as f:
        f.writelines(out_lines)

    def _cleanup(p=tmp_path):
        try:
            os.unlink(p)
        except OSError:
            pass
    atexit.register(_cleanup)
    return tmp_path

