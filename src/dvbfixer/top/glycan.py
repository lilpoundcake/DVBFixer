"""Glycan/glycolipid link detection for `dvbfixer top`.

Split out of `top/pipeline.py`: `detect_glycan_links` + `build_glycan_trees`
are plain top-level functions called once each from `main()`, not part of
`TopologyBuilder` — they run BEFORE any topology building starts, to figure
out which residues form glycan trees so the chain classifier can route them
away from the generic per-chain builder.
"""
import math
from collections import defaultdict

from dvbfixer.top.ff_data import CERAMIDE_RTP, PDB_TO_CARB, PDB_TO_LIPID


def _is_ceramide(resname):
    """Check if a residue name is a ceramide (PDB or RTP name)."""
    return resname in PDB_TO_LIPID or resname in CERAMIDE_RTP


def _parse_conect_bonds(pdb_path):
    """Parse CONECT and LINK records from PDB to get explicit bond pairs.

    Returns set of (serial1, serial2) tuples from CONECT (unordered),
    serial→(chain, resseq, atomname) mapping,
    and list of LINK bonds as (chain1, resseq1, atom1, chain2, resseq2, atom2).
    """
    bonds = set()
    serial_map = {}
    link_bonds = []  # (ch1, rs1, aname1, ch2, rs2, aname2)

    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                try:
                    serial = int(line[6:11])
                    aname = line[12:16].strip()
                    chain = line[21]
                    resseq = int(line[22:26])
                    serial_map[serial] = (chain, resseq, aname)
                except (ValueError, IndexError):
                    pass
            elif line.startswith('CONECT'):
                parts = line[6:].split()
                if len(parts) >= 2:
                    try:
                        s1 = int(parts[0])
                        for p in parts[1:]:
                            s2 = int(p)
                            if s1 != s2:
                                bonds.add((min(s1, s2), max(s1, s2)))
                    except ValueError:
                        pass
            elif line.startswith('LINK'):
                try:
                    a1 = line[12:16].strip()
                    ch1 = line[21]
                    rs1 = int(line[22:26])
                    a2 = line[42:46].strip()
                    ch2 = line[51]
                    rs2 = int(line[52:56])
                    link_bonds.append((ch1, rs1, a1, ch2, rs2, a2))
                except (ValueError, IndexError):
                    pass
    return bonds, serial_map, link_bonds


def detect_glycan_links(chains, pdb_path=None):
    """Auto-detect glycosidic bonds from CONECT records and C1-O distances.

    Uses CONECT records as primary source (reliable), distance-based
    detection as fallback.

    Detects:
    - Sugar-sugar links: C1 of one sugar bonded to Ox of another
    - Protein-sugar links: C1 bonded to ASN ND2 / SER OG / THR OG1
    - Sugar-ceramide links: sugar O bonded to ceramide C1S

    Returns list of GlycanLink tuples.
    """
    sugar_names = set(PDB_TO_CARB.keys())
    links = []

    # Also detect sugars already named with CHARMM RTP names
    charmm_sugar_names = set(PDB_TO_CARB.values())

    # Collect C1 atoms from sugars and ND2/O* from sugars+ASN
    c1_atoms = []
    target_atoms = []
    # Ceramide C1S atoms (for sugar O1-ceramide C1S detection)
    cer_c1s_atoms = []
    # Sugar O atoms (for reverse check against ceramide)
    sugar_o_atoms = []

    for chain in chains:
        ch = chain.chain_id
        for res in chain.residues:
            is_sugar = res.resname in sugar_names or res.resname in charmm_sugar_names
            # Sialic acid variants: anomeric carbon is C2 (not C1)
            is_sialic = res.resname in ('SIA', 'ANE5', 'ANE5AC', 'BNE5AC')
            for aname, x, y, z in res.atoms:
                if is_sugar and aname == 'C1' and not is_sialic:
                    c1_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))
                # Sialic acid links via C2
                if is_sugar and is_sialic and aname == 'C2':
                    c1_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))

                if is_sugar and aname.startswith('O') and aname != 'O5':
                    target_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))
                    sugar_o_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))

                # Protein residues that can be glycosylated
                if res.resname in ('ASN', 'NLN') and aname == 'ND2':
                    target_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))
                elif res.resname in ('SER', 'OLS') and aname == 'OG':
                    target_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))
                elif res.resname in ('THR', 'OLT') and aname == 'OG1':
                    target_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))

                # Collect ceramide C1S for sugar-ceramide bond detection
                if _is_ceramide(res.resname) and aname == 'C1S':
                    cer_c1s_atoms.append((ch, res.resname, res.resseq, x, y, z))

    # --- Primary: CONECT/LINK-based detection (explicit bonds, most reliable) ---
    conect_found = set()  # (ch1, rs1, ch2, rs2) pairs found via CONECT/LINK
    if pdb_path is not None:
        conect_bonds, serial_map, link_bonds = _parse_conect_bonds(str(pdb_path))

        # Process LINK records first (inter-residue bonds, always reliable)
        # Build sets of known C1/C2 and target atoms for fast lookup
        c1_set = {(ch, rs, aname) for ch, rn, rs, aname, x, y, z in c1_atoms}
        tgt_set = {(ch, rs, aname) for ch, rn, rs, aname, x, y, z in target_atoms}

        for ch1, rs1, a1, ch2, rs2, a2 in link_bonds:
            # Check if this is a sugar C1/C2 → target bond
            if (ch1, rs1, a1) in c1_set and (ch2, rs2, a2) in tgt_set:
                links.append((ch2, rs2, a2, ch1, rs1, a1))
                conect_found.add((ch1, rs1, ch2, rs2))
                conect_found.add((ch2, rs2, ch1, rs1))
            elif (ch2, rs2, a2) in c1_set and (ch1, rs1, a1) in tgt_set:
                links.append((ch1, rs1, a1, ch2, rs2, a2))
                conect_found.add((ch1, rs1, ch2, rs2))
                conect_found.add((ch2, rs2, ch1, rs1))

        # Process CONECT records (atom serial-based bonds)
        if conect_bonds:
            id_to_serial = {}
            for serial, (ch, rs, aname) in serial_map.items():
                id_to_serial[(ch, rs, aname)] = serial

            for ch1, rn1, rs1, aname1, x1, y1, z1 in c1_atoms:
                s1 = id_to_serial.get((ch1, rs1, aname1))
                if s1 is None:
                    continue
                for ch2, rn2, rs2, aname2, x2, y2, z2 in target_atoms:
                    if ch1 == ch2 and rs1 == rs2:
                        continue
                    if (ch1, rs1, ch2, rs2) in conect_found:
                        continue
                    s2 = id_to_serial.get((ch2, rs2, aname2))
                    if s2 is None:
                        continue
                    bond_key = (min(s1, s2), max(s1, s2))
                    if bond_key in conect_bonds:
                        links.append((ch2, rs2, aname2, ch1, rs1, aname1))
                        conect_found.add((ch1, rs1, ch2, rs2))
                        conect_found.add((ch2, rs2, ch1, rs1))

    # --- Fallback: distance-based detection (for atoms not found via CONECT) ---
    _PROTEIN_ATOMS = {'ND2', 'OG', 'OG1'}
    for ch1, rn1, rs1, aname1, x1, y1, z1 in c1_atoms:
        for ch2, rn2, rs2, aname2, x2, y2, z2 in target_atoms:
            if ch1 == ch2 and rs1 == rs2:
                continue
            # Skip if already found via CONECT
            if (ch1, rs1, ch2, rs2) in conect_found:
                continue
            d = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            cutoff = 2.5 if aname2 in _PROTEIN_ATOMS else 2.0
            if d < cutoff:
                links.append((ch2, rs2, aname2, ch1, rs1, aname1))

    # Find sugar O - ceramide C1S bonds within 2.0 A
    # In glycolipids, sugar O1 bridges to ceramide C1S (ceramide's own O1/HO1
    # are already removed in CHARMM-GUI output)
    for ch1, rn1, rs1, aname1, x1, y1, z1 in sugar_o_atoms:
        for ch2, rn2, rs2, x2, y2, z2 in cer_c1s_atoms:
            d = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            if d < 2.0:
                # Link: ceramide (donor, has C1S) <- sugar (has Ox)
                links.append((ch2, rs2, 'C1S', ch1, rs1, aname1))

    return links


def build_glycan_trees(glycan_links, chains):
    """Build glycan trees from linkage list.

    Returns list of glycan trees, each a list of (chain_id, resseq, resname)
    in topological order (root first).
    Also detects ceramide-sugar links for glycolipids.
    """
    # Build adjacency: sugar -> list of child sugars
    sugar_residues = set()
    charmm_sugar_names = set(PDB_TO_CARB.values())
    for chain in chains:
        for res in chain.residues:
            if res.resname in PDB_TO_CARB or res.resname in charmm_sugar_names:
                sugar_residues.add((chain.chain_id, res.resseq))

    # Identify ceramide residues
    ceramide_residues = set()
    for chain in chains:
        for res in chain.residues:
            if _is_ceramide(res.resname):
                ceramide_residues.add((chain.chain_id, res.resseq))

    # Find root sugars (connected to protein ASN, not to another sugar)
    children = defaultdict(list)  # parent -> [child]
    parent_of = {}  # child -> parent
    protein_links = []  # [(protein_chain, protein_resseq, protein_atom, sugar_chain, sugar_resseq)]
    ceramide_links = []  # [(cer_chain, cer_resseq, cer_atom, sugar_chain, sugar_resseq, sugar_atom)]

    for don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom in glycan_links:
        acc_key = (acc_ch, acc_rs)
        don_key = (don_ch, don_rs)
        if don_key in sugar_residues and acc_key in sugar_residues:
            # sugar-sugar link: donor is parent, acceptor is child
            children[don_key].append((acc_key, don_atom, acc_atom))
            parent_of[acc_key] = (don_key, don_atom)
        elif don_key in ceramide_residues and acc_key in sugar_residues:
            # ceramide-sugar link — do NOT add to parent_of so the sugar
            # becomes a root of its glycan tree (ceramide handled separately)
            ceramide_links.append((don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom))
        elif acc_key in sugar_residues:
            # protein-sugar link
            protein_links.append((don_ch, don_rs, don_atom, acc_ch, acc_rs))

    # Find roots (sugars with no sugar parent)
    roots = []
    for ch, rs in sugar_residues:
        if (ch, rs) not in parent_of:
            roots.append((ch, rs))

    # BFS from each root to build tree
    trees = []
    for root in roots:
        tree = []
        # Find linkage info for this root
        link_atoms = {}  # (ch, rs) -> set of linked O positions
        queue = [root]
        visited = set()
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            tree.append(node)
            for child, donor_atom, child_atom in children.get(node, []):
                # Record which O on parent is linked
                if node not in link_atoms:
                    link_atoms[node] = set()
                link_atoms[node].add(donor_atom)
                # Child's anomeric O is linked: O1 for C1 linkage, O2 for C2 linkage
                if child_atom in ('C1', 'C2'):
                    linked_o = 'O1' if child_atom == 'C1' else 'O2'
                    if child not in link_atoms:
                        link_atoms[child] = set()
                    link_atoms[child].add(linked_o)
                queue.append(child)

        # Mark root's O1 as linked if bonded to ceramide
        for cl in ceramide_links:
            if (cl[3], cl[4]) == root:
                if root not in link_atoms:
                    link_atoms[root] = set()
                link_atoms[root].add('O1')

        # Also mark root's O1 as linked (bonds to protein or is reducing end)
        for pl in protein_links:
            if (pl[3], pl[4]) == root:
                if root not in link_atoms:
                    link_atoms[root] = set()
                link_atoms[root].add('O1')

        # Filter ceramide_links relevant to this tree
        tree_set = set(tree)
        tree_cer_links = [cl for cl in ceramide_links
                          if (cl[3], cl[4]) in tree_set]

        trees.append((tree, link_atoms, protein_links, tree_cer_links))

    return trees
