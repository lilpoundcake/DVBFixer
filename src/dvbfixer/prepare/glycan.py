"""Glycan-aware helpers for ``dvbfixer prepare``.

Split out of the flat ``prepare.py`` in the Phase 2.3 follow-up work.

- :func:`find_glycosylated_atoms` / :func:`find_glycosylated_atoms_with_sugar`
  — locate ASN/SER/THR residues bonded to a sugar via CONECT records plus a
  distance fallback. The `_with_sugar` variant also returns the bonded
  sugar's residue name, needed by
  :func:`rename_glycosylated_protein_residues` to decide when to
  rename ASN → NLN, etc.
- :func:`remove_extra_glycan_hydrogens` — after PDBFixer adds hydrogens,
  drop the extra HD22 / HG / HG1 from glycosylated ASN / SER / THR
  respectively (the sidechain N/O is bonded to the sugar, not to H).
- :func:`rename_glycosylated_protein_residues` — rename ASN/SER/THR →
  NLN/OLS/OLT when the bonded sugar is GLYCAM-named. PDB-named sugars
  leave the protein residue standard-named.
- :func:`_build_glycan_trees` — BFS the glycan connectivity from each
  protein anchor, needed by the heterogen-H hydrogen adders below.
- :func:`add_heterogen_h_via_rdkit` — the preferred H-addition path
  (SDF round-trip through RDKit). Falls back to OpenBabel via
  :func:`add_heterogen_h_via_openbabel` when RDKit can't handle the
  input (unknown ligand, bad valence, etc.).
"""

from __future__ import annotations

from openmm.app import Modeller

from dvbfixer.prepare.cli import GLYCOSYLATED_RESIDUES, SUGAR_RESNAMES


def find_glycosylated_atoms(input_path):
    """Find protein atoms bonded to sugars (CONECT + distance fallback).

    Returns a dict-like set of `(chain_id, resid, atom_name)` tuples
    (e.g. ASN ND2 bonded to NAG C1). Use `find_glycosylated_atoms_with_sugar`
    for the additional bonded-sugar resname info needed by GLYCAM rename
    logic.
    """
    return set(find_glycosylated_atoms_with_sugar(input_path).keys())


def find_glycosylated_atoms_with_sugar(input_path):
    """Find protein atoms bonded to sugars (CONECT records + distance
    fallback). Returns dict `{(chain_id, resid, atom_name): sugar_resname}`.

    Distance fallback is needed when the input is missing CONECT records
    for some glycosidic bonds. ND2/OG/OG1 within 2.0 Å of a sugar anomeric
    C (C1, or C2 for sialic acid) → glycosidic bond.

    The sugar resname lets downstream rename logic distinguish GLYCAM-named
    sugars (UYB/4YB/VMB/...) — which need the protein anchor renamed to
    NLN/OLS/OLT for GLYCAM template matching — from PDB-named sugars
    (NAG/NDG/BMA/MAN/...) — which are handled by whatever FF the user
    routes them through and don't need the rename (extra HD22 is still
    removed regardless).
    """
    from dvbfixer.ffutils import is_glycam_sugar

    with open(input_path) as f:
        lines = f.readlines()

    serials = {}
    # Also collect atom coords for distance fallback
    by_id = {}  # (chain, resname, resid, atomname) → (x, y, z)
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                serial = int(line[6:11])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            chain = line[21]
            resname = line[17:20].strip()
            resid = line[22:26].strip()
            atomname = line[12:16].strip()
            serials[serial] = (chain, resname, resid, atomname)
            by_id[(chain, resname, resid, atomname)] = (x, y, z)

    glycosylated = {}  # (chain, resid, atom) → sugar_resname

    def _sugar_class(name):
        if name in SUGAR_RESNAMES:
            return name
        if is_glycam_sugar(name):
            return name
        return None

    # Pass 1: CONECT-based detection.
    for line in lines:
        if not line.startswith('CONECT'):
            continue
        parts = []
        s = line[6:]
        while len(s) >= 5:
            chunk = s[:5].strip()
            if chunk:
                try:
                    parts.append(int(chunk))
                except ValueError:
                    pass
            s = s[5:]
        if len(parts) < 2:
            continue
        src = parts[0]
        for dst in parts[1:]:
            if src not in serials or dst not in serials:
                continue
            s_info, d_info = serials[src], serials[dst]
            d_sugar = _sugar_class(d_info[1])
            s_sugar = _sugar_class(s_info[1])
            if s_info[1] in GLYCOSYLATED_RESIDUES and d_sugar:
                glycosylated[(s_info[0], s_info[2], s_info[3])] = d_sugar
            elif d_info[1] in GLYCOSYLATED_RESIDUES and s_sugar:
                glycosylated[(d_info[0], d_info[2], d_info[3])] = s_sugar

    # Pass 2: distance-based fallback. Iterate over protein anchor atoms
    # (ASN/SER/THR ND2/OG/OG1) and find nearby sugar anomeric carbons.
    # 2.0 Å cutoff catches typical glycosidic bonds (~1.45 Å) with margin.
    _ANCHOR_ATOMS = {'ASN': 'ND2', 'SER': 'OG', 'THR': 'OG1',
                     'NLN': 'ND2', 'OLS': 'OG', 'OLT': 'OG1'}
    _ANOMERIC_NAMES = {'C1', 'C2'}  # C2 for sialic acid
    SIALIC_RESNAMES = {'SIA', 'NAN'}
    sugar_anomeric = []  # list of (chain, resname, resid, atom, coord)
    for (ch, rn, rs, an), pos in by_id.items():
        if _sugar_class(rn) is None:
            continue
        if an not in _ANOMERIC_NAMES:
            continue
        # C2 anomeric is only for sialic acid (or GLYCAM *S* codes)
        if an == 'C2':
            is_sialic = rn in SIALIC_RESNAMES or (
                len(rn) == 3 and rn[1] in ('S', 's')
            )
            if not is_sialic:
                continue
        sugar_anomeric.append((ch, rn, rs, an, pos))

    cutoff2 = 2.0 ** 2
    for (ch, rn, rs, an), pos in by_id.items():
        if rn not in _ANCHOR_ATOMS:
            continue
        if an != _ANCHOR_ATOMS[rn]:
            continue
        key = (ch, rs, an)
        if key in glycosylated:
            continue  # already detected via CONECT
        # Find nearest sugar anomeric C within cutoff
        best = None
        best_d2 = cutoff2
        for (s_ch, s_rn, s_rs, s_an, s_pos) in sugar_anomeric:
            d2 = ((pos[0] - s_pos[0]) ** 2
                  + (pos[1] - s_pos[1]) ** 2
                  + (pos[2] - s_pos[2]) ** 2)
            if d2 < best_d2:
                best_d2 = d2
                best = s_rn
        if best is not None:
            glycosylated[key] = best

    return glycosylated


def remove_extra_glycan_hydrogens(fixer, glycosylated_atoms, verbose=False):
    """Remove extra hydrogen from glycosylated atoms (e.g. HD22 from ASN ND2).

    PDBFixer adds hydrogens assuming standard templates, but glycosylated
    ASN ND2 has an external bond to the sugar, so it should have only 1 H
    instead of 2.
    """
    if not glycosylated_atoms:
        return

    # Find atoms to delete: for each glycosylated atom, find the last H bonded to it
    atoms_to_delete = []
    for atom in fixer.topology.atoms():
        res = atom.residue
        key = (res.chain.id, res.id, atom.name)
        if key not in glycosylated_atoms:
            continue

        # Find hydrogens bonded to this atom
        h_atoms = []
        for bond in fixer.topology.bonds():
            a1, a2 = bond
            if a1.index == atom.index and a2.element.symbol == 'H':
                h_atoms.append(a2)
            elif a2.index == atom.index and a1.element.symbol == 'H':
                h_atoms.append(a1)

        if len(h_atoms) > 1:
            # Remove the last hydrogen (HD22 for ASN ND2)
            to_remove = sorted(h_atoms, key=lambda a: a.name)[-1]
            atoms_to_delete.append(to_remove)
            if verbose:
                print(f"  Removing extra H: {res.chain.id}:{res.name}{res.id}:{to_remove.name} "
                      f"(glycosylated {atom.name})")

    if atoms_to_delete:
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.delete(atoms_to_delete)
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions
        print(f"Removed {len(atoms_to_delete)} extra hydrogen(s) from glycosylated residues")


# Map standard protein residue → GLYCAM glycosylated variant
_GLYCAM_RENAME = {'ASN': 'NLN', 'SER': 'OLS', 'THR': 'OLT'}
# H atoms that GLYCAM glycoprotein templates do NOT expect (one bond replaced
# by sugar). NLN keeps HD21, drops HD22. OLS/OLT drop the hydroxyl H.
_GLYCAM_DROP_H = {
    'NLN': {'HD22'},
    'OLS': {'HG', 'HG1'},
    'OLT': {'HG1'},
}


def _write_heterogen_conects(pdb_path, topology):
    """Emit CONECT records for all topology bonds where either atom is a
    heterogen (not in PROTEIN_RESIDUES | SOLVENT_IONS). Appends before END.

    Without this, PyMOL/VMD show no protein-glycan or sugar-sugar bonds, and
    downstream minimize loses cross-residue connectivity.
    """
    from collections import defaultdict

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    from dvbfixer.pdbutils import append_before_end, build_serial_map

    known = PROTEIN_RESIDUES | SOLVENT_IONS
    serial_map = build_serial_map(pdb_path)

    adj = defaultdict(set)
    for a1, a2 in topology.bonds():
        r1, r2 = a1.residue, a2.residue
        # Skip protein-only bonds — those are inferred from templates
        if r1.name in known and r2.name in known:
            continue
        k1 = (r1.chain.id, str(r1.id), a1.name)
        k2 = (r2.chain.id, str(r2.id), a2.name)
        s1 = serial_map.get(k1)
        s2 = serial_map.get(k2)
        if s1 is None or s2 is None:
            continue
        adj[s1].add(s2)
        adj[s2].add(s1)

    if not adj:
        return

    lines = []
    for src in sorted(adj):
        partners = sorted(adj[src])
        # CONECT supports up to 4 partners per line
        for i in range(0, len(partners), 4):
            chunk = partners[i:i + 4]
            lines.append(
                "CONECT" + f"{src:5d}"
                + "".join(f"{p:5d}" for p in chunk)
                + "\n"
            )
    append_before_end(pdb_path, lines)


def _build_glycan_trees(topology, positions, known_set):
    """Group heterogen residues into connected trees.

    Edges come from topology.bonds() PLUS distance-based sugar-sugar links
    (C1/C2 of one sugar within ≤ 2.0 Å of any O atom of another) since input
    PDBs often omit sugar-sugar CONECTs.

    Returns:
      trees: list of sets of residue.index
      protein_anchors: dict tree_idx → list of protein atoms bonded to the tree
      extra_inter_bonds: list of (atom1, atom2) — sugar-sugar bonds found by
        distance that are NOT already in topology.bonds()
    """
    # Adjacency between heterogen residues
    adj = {}  # res.index → set of res.index
    het_residues = [r for r in topology.residues() if r.name not in known_set]
    het_idx = {r.index for r in het_residues}
    for r in het_residues:
        adj.setdefault(r.index, set())

    # Track protein atoms bonded to each heterogen residue
    protein_links = {}  # het_res.index → list of protein atoms

    # Existing bonds from topology
    existing_edge = set()  # frozenset({res1.index, res2.index}) for inter-residue
    for a1, a2 in topology.bonds():
        r1, r2 = a1.residue, a2.residue
        if r1.index == r2.index:
            continue
        in1, in2 = r1.index in het_idx, r2.index in het_idx
        if in1 and in2:
            adj[r1.index].add(r2.index)
            adj[r2.index].add(r1.index)
            existing_edge.add(frozenset({r1.index, r2.index}))
        elif in1 and r2.name in known_set and r2.name not in {'HOH', 'WAT', 'TIP3', 'SOL'}:
            protein_links.setdefault(r1.index, []).append(a2)
        elif in2 and r1.name in known_set and r1.name not in {'HOH', 'WAT', 'TIP3', 'SOL'}:
            protein_links.setdefault(r2.index, []).append(a1)

    # Distance-based sugar-sugar bond perception: C1/C2 of one sugar within
    # ≤ 2.0 Å of any O atom of another sugar residue.
    extra_inter_bonds = []
    # Collect anomeric C atoms and O atoms grouped by residue
    anomeric = {}  # res.index → list of (atom, pos_a)
    oxygens = {}   # res.index → list of (atom, pos_a)
    for r in het_residues:
        for a in r.atoms():
            if a.element.symbol == 'C' and a.name in ('C1', 'C2'):
                p = positions[a.index]
                anomeric.setdefault(r.index, []).append(
                    (a, (p.x * 10.0, p.y * 10.0, p.z * 10.0))
                )
            elif a.element.symbol == 'O':
                p = positions[a.index]
                oxygens.setdefault(r.index, []).append(
                    (a, (p.x * 10.0, p.y * 10.0, p.z * 10.0))
                )

    for r1 in het_residues:
        for atom_c, p_c in anomeric.get(r1.index, ()):
            for r2 in het_residues:
                if r1.index == r2.index:
                    continue
                if frozenset({r1.index, r2.index}) in existing_edge:
                    continue
                for atom_o, p_o in oxygens.get(r2.index, ()):
                    dx = p_c[0] - p_o[0]
                    dy = p_c[1] - p_o[1]
                    dz = p_c[2] - p_o[2]
                    d2 = dx * dx + dy * dy + dz * dz
                    if d2 <= 4.0:  # 2.0 Å
                        adj[r1.index].add(r2.index)
                        adj[r2.index].add(r1.index)
                        existing_edge.add(frozenset({r1.index, r2.index}))
                        extra_inter_bonds.append((atom_c, atom_o))
                        break

    # BFS to build trees
    seen = set()
    trees = []
    tree_anchors = []
    for r in het_residues:
        if r.index in seen:
            continue
        tree = set()
        stack = [r.index]
        while stack:
            ri = stack.pop()
            if ri in seen:
                continue
            seen.add(ri)
            tree.add(ri)
            for ni in adj[ri]:
                if ni not in seen:
                    stack.append(ni)
        trees.append(tree)
        # Collect protein anchors for this tree
        anchors = []
        for ri in tree:
            anchors.extend(protein_links.get(ri, []))
        tree_anchors.append(anchors)

    return trees, tree_anchors, extra_inter_bonds


def add_heterogen_h_via_rdkit(topology, positions, output_pdb_path,
                               verbose=False):
    """BioLuminate-style H addition for all heterogens via RDKit.

    Pipes the prepared PDB (with CONECT records for inter-residue bonds) through
    RDKit:  MolFromPDBBlock (uses CONECT + distance perception for full-molecule
    graph) → AddHs(addCoords=True) (valence-based H placement with proper 3D
    coords, respecting all bonds). Then re-parses the augmented PDB to update
    OpenMM topology.

    Returns new (topology, positions).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        if verbose:
            print("  RDKit not available — skipping heterogen H polish")
        return topology, positions

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    known = PROTEIN_RESIDUES | SOLVENT_IONS

    # Need any heterogens without H?
    needs_h = False
    for r in topology.residues():
        if r.name in known:
            continue
        if not any(a.element.symbol == 'H' for a in r.atoms()):
            needs_h = True
            break
    if not needs_h:
        return topology, positions

    # Read current PDB (which has CONECTs we just wrote)
    pdb_text = open(output_pdb_path).read()

    # Suppress RDKit warnings
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')

    mol = Chem.MolFromPDBBlock(pdb_text, sanitize=False, removeHs=False,
                                proximityBonding=True)
    if mol is None:
        if verbose:
            print("  RDKit failed to parse PDB — skipping heterogen H polish")
        return topology, positions

    # Atoms in known PDB sugars that should NOT carry hydrogens. These are
    # carbonyl atoms (C=O) and atoms forming external bonds we know about.
    # RDKit's AddHs after MolFromPDBBlock can't perceive double bonds from
    # PDB CONECT records alone, so we post-filter the added H by parent.
    _NO_H_ATOM_NAMES = {
        'NAG': {'C7', 'O7'},       # acetamide C=O
        'NDG': {'C7', 'O7'},
        'NGA': {'C7', 'O7'},
        'A2G': {'C7', 'O7'},
        'SIA': {'C1', 'O1A', 'O1B'},  # carboxylate
    }

    n_before = mol.GetNumAtoms()
    try:
        molh = Chem.AddHs(mol, addCoords=True)
    except Exception as e:
        if verbose:
            print(f"  RDKit AddHs failed: {e}")
        return topology, positions

    # MMFF94 geometry refinement of heterogen atoms (BioLuminate-style).
    # Freeze protein atoms; let glycan H atoms + glycosidic geometry relax.
    # This gives proper bond lengths/angles for glycans before any FF minimize.
    # RDKit's MMFF setup emits multi-line C++ stack traces to stderr when
    # SSSR wasn't precomputed on an unusual ligand (harmless — we catch the
    # exception and skip the pass). Silence its log for the duration.
    _RDKIT_LOG_QUIETED = False
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        _RDKIT_LOG_QUIETED = True
    except Exception:
        pass

    import contextlib as _cl
    import io as _io
    _rd_stderr = _io.StringIO()
    try:
        with _cl.redirect_stderr(_rd_stderr):
            from rdkit.Chem import AllChem
            # Get atom indices of heterogen heavy atoms + their newly-added H
            het_indices = set()
            for atom in molh.GetAtoms():
                info = atom.GetPDBResidueInfo()
                if info and info.GetResidueName().strip() not in known:
                    het_indices.add(atom.GetIdx())
                elif atom.GetIdx() >= n_before:
                    # Added H — keep its position frozen if parent is protein,
                    # let it move if parent is heterogen
                    for bond in atom.GetBonds():
                        other = bond.GetOtherAtom(atom)
                        other_info = other.GetPDBResidueInfo()
                        if (other_info and
                                other_info.GetResidueName().strip() not in known):
                            het_indices.add(atom.GetIdx())
                        break

            if het_indices and len(het_indices) < molh.GetNumAtoms():
                # NOTE: do NOT pre-compute SSSR on `molh` here. `molh` is
                # the whole PDB (protein + heterogens, thousands of atoms)
                # and RDKit's SSSR can run unboundedly long on large mixed
                # systems. If MMFFGetMoleculeProperties needs ring info
                # and it isn't there, it will raise — the except below
                # takes the "MMFF skipped" branch, which is the correct
                # behaviour for ligands MMFF can't type.
                props = AllChem.MMFFGetMoleculeProperties(
                    molh, mmffVariant='MMFF94s')
                if props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(
                        molh, props, ignoreInterfragInteractions=False)
                    if ff is not None:
                        for i in range(molh.GetNumAtoms()):
                            if i not in het_indices:
                                ff.AddFixedPoint(i)
                        ff.Minimize(maxIts=400)
                        if verbose:
                            print(f"  MMFF94 refined {len(het_indices)} "
                                  f"heterogen atoms")
    except Exception:
        if verbose:
            print("  MMFF94 refinement skipped (ligand not MMFF-typeable — "
                  "harmless; coords unchanged)")
    finally:
        if _RDKIT_LOG_QUIETED:
            try:
                RDLogger.EnableLog('rdApp.*')
            except Exception:
                pass

    n_after = molh.GetNumAtoms()
    added = n_after - n_before
    if added <= 0:
        return topology, positions

    # RDKit's PDB writer loses residue info for atoms added by AddHs (they
    # all get tagged as UNL/Unknown Ligand). Instead, extract each new H atom
    # with its parent (bonded heavy atom) and use the parent's residue info
    # to assign H to the correct OpenMM residue.
    from openmm import Vec3
    from openmm.app import Topology, element
    from openmm.unit import Quantity, nanometer

    new_h_specs = []  # (parent_omm, h_name, (x_nm, y_nm, z_nm))
    omm_atoms = list(topology.atoms())

    # Map RDKit heavy atom index → OpenMM atom (positional, since both built
    # from same PDB in order).
    rdkit_heavy_to_omm = {}
    rdkit_idx = 0
    for omm_atom in omm_atoms:
        if rdkit_idx >= mol.GetNumAtoms():
            break
        rdkit_heavy_to_omm[rdkit_idx] = omm_atom
        rdkit_idx += 1

    h_count_by_parent = {}
    conf = molh.GetConformer()
    for ai in range(n_before, n_after):
        new_atom = molh.GetAtomWithIdx(ai)
        if new_atom.GetAtomicNum() != 1:
            continue
        # Find parent heavy atom (must be in original)
        parent_idx = None
        for bond in new_atom.GetBonds():
            other = bond.GetOtherAtom(new_atom)
            if other.GetIdx() < n_before:
                parent_idx = other.GetIdx()
                break
        if parent_idx is None:
            continue
        parent_omm = rdkit_heavy_to_omm.get(parent_idx)
        if parent_omm is None or parent_omm.residue.name in known:
            continue
        # Skip H if parent atom is on the no-H list for its residue type
        # (carbonyl C/O atoms that RDKit can't perceive as double-bonded).
        if parent_omm.name in _NO_H_ATOM_NAMES.get(parent_omm.residue.name, ()):
            continue
        # Get H position from conformer (in Å)
        pos = conf.GetAtomPosition(ai)
        x_nm, y_nm, z_nm = pos.x / 10.0, pos.y / 10.0, pos.z / 10.0
        # RDKit AddHs(addCoords=True) sometimes leaves an H at origin if it
        # couldn't compute a 3D position. Replace with a small offset from
        # the parent so downstream tools (xtb etc.) don't choke on (0,0,0).
        if abs(pos.x) < 1e-6 and abs(pos.y) < 1e-6 and abs(pos.z) < 1e-6:
            p_parent = positions[parent_omm.index]
            v_par = p_parent.value_in_unit(nanometer)
            # Offset 1.0 Å in arbitrary direction; minimization will fix
            x_nm = float(v_par[0]) + 0.1
            y_nm = float(v_par[1]) + 0.05
            z_nm = float(v_par[2]) + 0.05
        # Name H based on parent
        pname = parent_omm.name
        if pname.startswith(('O', 'N', 'S')):
            h_base = f"H{pname}"
        elif len(pname) > 1:
            h_base = f"H{pname[1:]}"
        else:
            h_base = "H"
        n_existing = h_count_by_parent.get(parent_omm.index, 0)
        h_name = h_base if n_existing == 0 else f"{h_base}{n_existing + 1}"
        h_count_by_parent[parent_omm.index] = n_existing + 1
        new_h_specs.append((parent_omm, h_name, (x_nm, y_nm, z_nm)))

    if not new_h_specs:
        return topology, positions

    # Build new topology with H inserted after each parent atom
    new_top = Topology()
    new_pos_list = []
    omm_to_new = {}
    h_by_parent = {}
    for parent, hname, hpos in new_h_specs:
        h_by_parent.setdefault(parent.index, []).append((hname, hpos))

    for chain in topology.chains():
        new_chain = new_top.addChain(chain.id)
        for res in chain.residues():
            new_res = new_top.addResidue(res.name, new_chain, res.id,
                                          res.insertionCode)
            for atom in res.atoms():
                new_atom = new_top.addAtom(atom.name, atom.element, new_res)
                omm_to_new[atom.index] = new_atom
                p = positions[atom.index]
                new_pos_list.append(Vec3(p.x, p.y, p.z))
                for (h_name, h_pos) in h_by_parent.get(atom.index, []):
                    h_atom = new_top.addAtom(h_name, element.hydrogen, new_res)
                    # Bond new H to its parent — without this, the H has no
                    # CONECT record in output PDB; OpenBabel's proximityBonding
                    # may miss the bond if RDKit placed H >1.9 Å away, and the
                    # H drifts free during obminimize/xtb refinement.
                    new_top.addBond(new_atom, h_atom)
                    new_pos_list.append(Vec3(*h_pos))

    # Carry over original OpenMM bonds (protein, protein-glycan from CONECT)
    added_pairs = set()
    def _add_bond(a1, a2):
        if a1 is None or a2 is None:
            return
        pair = (min(a1.index, a2.index), max(a1.index, a2.index))
        if pair in added_pairs:
            return
        added_pairs.add(pair)
        new_top.addBond(a1, a2)

    for b in topology.bonds():
        _add_bond(omm_to_new.get(b[0].index), omm_to_new.get(b[1].index))

    # Carry over intra-heterogen bonds that RDKit perceived (sugar-sugar,
    # intra-sugar) — these were missing from PDBFixer topology.
    for bond in molh.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        # Only consider bonds between original heavy atoms (not added H)
        if a >= n_before or b >= n_before:
            continue
        omm_a = rdkit_heavy_to_omm.get(a)
        omm_b = rdkit_heavy_to_omm.get(b)
        if omm_a is None or omm_b is None:
            continue
        # Skip protein-only bonds
        if omm_a.residue.name in known and omm_b.residue.name in known:
            continue
        _add_bond(omm_to_new.get(omm_a.index), omm_to_new.get(omm_b.index))

    new_positions = Quantity(new_pos_list, nanometer)
    print(f"RDKit added {len(new_h_specs)} H atoms to heterogens")
    return new_top, new_positions


def _process_single_residue(res, external_bond_counts, positions, known, ob,
                             new_h_specs, perceived_intra_bonds, verbose):
    """Per-residue H placement via OpenBabel PDB roundtrip.

    external_bond_counts: dict atom.index → number of inter-residue bonds.
    Used to suppress extra H atoms at linkage positions.
    """
    from openbabel import pybel

    res_atoms = list(res.atoms())
    # Build PDB string for just this residue
    lines = []
    for i, a in enumerate(res_atoms):
        p = positions[a.index]
        x_a, y_a, z_a = p.x * 10.0, p.y * 10.0, p.z * 10.0
        name = a.name
        if len(name) >= 4:
            name_field = name[:4]
        else:
            name_field = f" {name:<3s}"
        lines.append(
            f"HETATM{i+1:5d} {name_field} "
            f"{res.name:>3s} A{1:4d}    "
            f"{x_a:8.3f}{y_a:8.3f}{z_a:8.3f}  1.00  0.00          "
            f"{a.element.symbol:>2s}\n"
        )
    lines.append("END\n")
    pdb_str = "".join(lines)

    try:
        mol = pybel.readstring('pdb', pdb_str)
    except Exception as e:
        if verbose:
            print(f"  OpenBabel read failed for {res.name}{res.id}: {e}")
        return

    n_before = mol.OBMol.NumAtoms()
    if n_before == 0:
        return

    # Harvest perceived intra-residue bonds BEFORE addh (read-only access)
    for bond in ob.OBMolBondIter(mol.OBMol):
        b_idx = bond.GetBeginAtomIdx()
        e_idx = bond.GetEndAtomIdx()
        if 1 <= b_idx <= n_before and 1 <= e_idx <= n_before:
            a1 = res_atoms[b_idx - 1]
            a2 = res_atoms[e_idx - 1]
            perceived_intra_bonds.append((a1, a2))

    mol.addh()
    pdb_text = mol.write('pdb')

    # Parse PDB text — H atoms and CONECT bonds
    serial_to_info = {}  # serial → (elem, x, y, z)
    for line in pdb_text.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            try:
                serial = int(line[6:11])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                elem = line[76:78].strip()
                if not elem:
                    name = line[12:16].strip()
                    elem = name[0] if name and name[0].isalpha() else 'H'
                serial_to_info[serial] = (elem, x, y, z)
            except (ValueError, IndexError):
                continue

    # CONECT records: for each H, find its parent
    h_to_parent = {}
    for line in pdb_text.split('\n'):
        if not line.startswith('CONECT'):
            continue
        try:
            parts = [line[i:i + 5] for i in range(6, len(line.rstrip()), 5)]
            serials = [int(p) for p in parts if p.strip()]
        except ValueError:
            continue
        if len(serials) < 2:
            continue
        src = serials[0]
        src_info = serial_to_info.get(src)
        if not src_info:
            continue
        for dst in serials[1:]:
            dst_info = serial_to_info.get(dst)
            if not dst_info:
                continue
            if src_info[0] == 'H' and dst_info[0] != 'H':
                h_to_parent[src] = dst
            elif dst_info[0] == 'H' and src_info[0] != 'H':
                h_to_parent[dst] = src

    # Group H atoms by parent. Only H atoms with serial > n_before are new.
    new_h_by_parent = {}  # parent_serial → list of (x, y, z)
    for h_serial, parent_serial in h_to_parent.items():
        if h_serial <= n_before:
            continue
        if parent_serial < 1 or parent_serial > n_before:
            continue
        info = serial_to_info.get(h_serial)
        if not info:
            continue
        new_h_by_parent.setdefault(parent_serial, []).append(info[1:])

    # Add H to topology, BUT skip extras at atoms with external bonds.
    h_count_by_parent = {}
    for parent_serial, h_positions in new_h_by_parent.items():
        parent_omm = res_atoms[parent_serial - 1]
        n_external = external_bond_counts.get(parent_omm.index, 0)
        # n_h_to_keep = n_h_openbabel_wanted - n_external_bonds_already_present.
        # OpenBabel doesn't know about external bonds, so it adds H to satisfy
        # full valence. We drop n_external H atoms to account for the actual
        # external linkages.
        n_keep = max(0, len(h_positions) - n_external)
        for (x_a, y_a, z_a) in h_positions[:n_keep]:
            pname = parent_omm.name
            if pname.startswith(('O', 'N', 'S')):
                h_base = f"H{pname}"
            elif len(pname) > 1:
                h_base = f"H{pname[1:]}"
            else:
                h_base = "H"
            n_existing = h_count_by_parent.get(parent_omm.index, 0)
            h_name = h_base if n_existing == 0 else f"{h_base}{n_existing + 1}"
            h_count_by_parent[parent_omm.index] = n_existing + 1
            new_h_specs.append((parent_omm, parent_omm.residue, h_name,
                                (x_a / 10.0, y_a / 10.0, z_a / 10.0)))


def add_heterogen_h_via_openbabel(topology, positions, verbose=False):
    """Use OpenBabel to add H atoms to heterogens (sugars, ligands), respecting
    inter-residue bonds (protein-glycan, sugar-sugar).

    Groups heterogens into connected trees, builds ONE OBMol per tree (including
    protein-anchor atoms with their real elements so OpenBabel sees the linkage
    and doesn't add H at linkage carbons). Distance-based sugar-sugar bonds are
    discovered too. Perceived intra-residue and sugar-sugar bonds are added to
    the returned OpenMM topology so downstream FF parametrization sees them.

    Returns new (topology, positions).
    """
    try:
        from openbabel import openbabel as ob
    except ImportError:
        if verbose:
            print("  OpenBabel not available — skipping heterogen H polish")
        return topology, positions

    from openmm import Vec3
    from openmm.app import Topology, element
    from openmm.unit import Quantity, nanometer

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    known = PROTEIN_RESIDUES | SOLVENT_IONS

    # Find heterogen residues that need H. A residue needs H if:
    # (a) it has no H at all, or
    # (b) its H count is < heavy_count - external_bonds (insufficient — e.g.
    #     a NAG with only methyl H but no ring/hydroxyl/amide H).
    # We strip existing H from (b) residues before OpenBabel regenerates them,
    # so the placement is consistent (not a mix of stale + new positions).
    res_by_index = {r.index: r for r in topology.residues()}

    # Find distance-perceived inter-residue bonds (sugar-sugar) so they get
    # added to the output topology and counted as external bonds for H-filtering.
    _, _, extra_inter_bonds = _build_glycan_trees(topology, positions, known)

    # Tentative external bond count (topology + distance-perceived), used for
    # the "insufficient H" heuristic. Per-residue external bonds in the strip
    # block ignore intra-residue bonds.
    tentative_ext = {}
    for a1, a2 in topology.bonds():
        if a1.residue.index != a2.residue.index:
            tentative_ext[a1.index] = tentative_ext.get(a1.index, 0) + 1
            tentative_ext[a2.index] = tentative_ext.get(a2.index, 0) + 1
    for a1, a2 in extra_inter_bonds:
        tentative_ext[a1.index] = tentative_ext.get(a1.index, 0) + 1
        tentative_ext[a2.index] = tentative_ext.get(a2.index, 0) + 1

    needs_h = set()
    h_to_strip = []
    for r in res_by_index.values():
        if r.name in known:
            continue
        heavy = [a for a in r.atoms() if a.element.symbol != 'H']
        h_atoms = [a for a in r.atoms() if a.element.symbol == 'H']
        # Crude expected H count: sum (default_valence - existing_bonds) over
        # heavy atoms. With incomplete bond info this overcounts but is fine
        # for an "insufficient H" trigger.
        n_ext = sum(tentative_ext.get(a.index, 0) for a in heavy)
        # Empirical minimum: each heavy atom (excluding O/N involved in
        # external bonds) typically carries ≥1 H. So expect ≥ len(heavy) - n_ext.
        min_expected = max(0, len(heavy) - n_ext)
        if len(h_atoms) == 0:
            needs_h.add(r.index)
        elif len(h_atoms) < min_expected:
            # Strip stale H so OpenBabel regenerates a consistent set.
            needs_h.add(r.index)
            h_to_strip.extend(h_atoms)

    if not needs_h:
        return topology, positions

    # Strip stale H from partially-H'd residues before regeneration.
    if h_to_strip:
        from openmm.app import Modeller as _Modeller
        m = _Modeller(topology, positions)
        m.delete(h_to_strip)
        topology = m.topology
        positions = m.positions
        # Rebuild res_by_index after deletion (atom indices shift)
        res_by_index = {r.index: r for r in topology.residues()}
        # Re-run glycan tree detection on new topology
        _, _, extra_inter_bonds = _build_glycan_trees(topology, positions, known)
        # Re-identify needs_h by (chain, id) since indices shifted
        needs_h_keys = set()
        for r in topology.residues():
            if r.name in known:
                continue
            if not any(a.element.symbol == 'H' for a in r.atoms()):
                needs_h_keys.add((r.chain.id, r.id))
        needs_h = {r.index for r in res_by_index.values()
                   if (r.chain.id, r.id) in needs_h_keys}
        if verbose:
            print(f"  Stripped {len(h_to_strip)} stale heterogen H atoms before regeneration")

    # Count external bonds per atom (inter-residue bonds from topology + distance).
    # Used to suppress extra H at linkage positions in per-residue OpenBabel.
    external_bond_counts = {}
    for a1, a2 in topology.bonds():
        if a1.residue.index != a2.residue.index:
            external_bond_counts[a1.index] = external_bond_counts.get(a1.index, 0) + 1
            external_bond_counts[a2.index] = external_bond_counts.get(a2.index, 0) + 1
    for a1, a2 in extra_inter_bonds:
        external_bond_counts[a1.index] = external_bond_counts.get(a1.index, 0) + 1
        external_bond_counts[a2.index] = external_bond_counts.get(a2.index, 0) + 1

    new_h_specs = []
    perceived_intra_bonds = []

    for res in res_by_index.values():
        if res.index not in needs_h:
            continue
        try:
            _process_single_residue(
                res, external_bond_counts, positions, known, ob,
                new_h_specs, perceived_intra_bonds, verbose,
            )
        except Exception as e:
            if verbose:
                print(f"  OpenBabel residue {res.name}{res.id} failed: {e}")

    if not new_h_specs and not extra_inter_bonds and not perceived_intra_bonds:
        return topology, positions

    # Build a fresh topology preserving order, inserting new H after parents.
    new_top = Topology()
    new_pos_list = []
    omm_to_new = {}

    h_by_parent = {}
    for parent, res, hname, hpos in new_h_specs:
        h_by_parent.setdefault(parent.index, []).append((res, hname, hpos))

    for chain in topology.chains():
        new_chain = new_top.addChain(chain.id)
        for res in chain.residues():
            new_res = new_top.addResidue(res.name, new_chain, res.id,
                                          res.insertionCode)
            for atom in res.atoms():
                new_atom = new_top.addAtom(atom.name, atom.element, new_res)
                omm_to_new[atom.index] = new_atom
                p = positions[atom.index]
                new_pos_list.append(Vec3(p.x, p.y, p.z))
                for (h_res, h_name, h_pos) in h_by_parent.get(atom.index, []):
                    h_atom = new_top.addAtom(h_name, element.hydrogen, new_res)
                    # Bond new H to its parent — without this, downstream tools
                    # (OpenBabel proximityBonding, xtb distance perception) may
                    # miss the bond if RDKit placed H >1.9 Å from parent, and
                    # the H drifts free during refinement.
                    new_top.addBond(new_atom, h_atom)
                    new_pos_list.append(Vec3(*h_pos))

    # Re-add bonds from original topology + perceived intra + sugar-sugar
    added_pairs = set()
    def _add_bond(a1, a2):
        n1 = omm_to_new.get(a1.index)
        n2 = omm_to_new.get(a2.index)
        if n1 is None or n2 is None:
            return
        pair = (min(n1.index, n2.index), max(n1.index, n2.index))
        if pair in added_pairs:
            return
        added_pairs.add(pair)
        new_top.addBond(n1, n2)

    for b in topology.bonds():
        _add_bond(b[0], b[1])
    for (a1, a2) in perceived_intra_bonds:
        _add_bond(a1, a2)
    for (a1, a2) in extra_inter_bonds:
        _add_bond(a1, a2)

    new_positions = Quantity(new_pos_list, nanometer)

    added = len(new_h_specs)
    if added:
        print(f"OpenBabel added {added} H atoms to heterogens")
    if extra_inter_bonds and verbose:
        print(f"  Detected {len(extra_inter_bonds)} sugar-sugar bonds by distance")
    return new_top, new_positions


def rename_glycosylated_protein_residues(topology, positions, glycosylated_atoms,
                                          verbose=False, sugar_by_anchor=None):
    """Rename ASN/SER/THR with glycosidic bonds to NLN/OLS/OLT — but ONLY
    when the bonded sugar is GLYCAM-named (UYB/4YB/VMB/...).

    For PDB-named sugars (NAG/NDG/BMA/MAN/...) the renames are NOT applied:
    NLN/OLS/OLT are GLYCAM-specific templates and the user is expected to
    either (a) convert PDB sugars to GLYCAM naming via `dvbfixer convert
    --to-amber` first, or (b) rely on downstream tools that don't need the
    glycoprotein rename. Extra HD22 (etc.) is still removed in both cases
    — that's handled by the caller via `remove_extra_glycan_hydrogens`.

    `sugar_by_anchor` (optional): dict {(chain, resid, atom): sugar_resname}
    from `find_glycosylated_atoms_with_sugar`. If None, every entry in
    `glycosylated_atoms` is treated as if bonded to a PDB-named sugar
    (no rename).

    Returns (new_topology, new_positions, renamed_keys).
    """
    from dvbfixer.ffutils import is_glycam_sugar

    by_res = {}
    for ch, rid, atom in glycosylated_atoms:
        by_res.setdefault((ch, rid), set()).add(atom)

    # Determine which anchor residues are bonded to GLYCAM-named sugars.
    glycam_anchors = set()
    if sugar_by_anchor:
        for (ch, rid, atom), sugar_rn in sugar_by_anchor.items():
            if is_glycam_sugar(sugar_rn):
                glycam_anchors.add((ch, rid))

    renamed = set()
    h_to_drop = []
    for res in topology.residues():
        key = (res.chain.id, res.id)
        if key not in by_res:
            continue
        if key not in glycam_anchors:
            # PDB-named sugar → keep ASN/SER/THR name. HD22/HG removal is
            # handled by remove_extra_glycan_hydrogens downstream.
            continue
        new_name = _GLYCAM_RENAME.get(res.name)
        if new_name is None:
            continue
        expected = {'NLN': 'ND2', 'OLS': 'OG', 'OLT': 'OG1'}[new_name]
        atom_names = {a.name for a in res.atoms()}
        if expected not in atom_names:
            if verbose:
                print(f"  Skip glycam-rename: {res.chain.id}:{res.name}{res.id} "
                      f"missing {expected}")
            continue
        res.name = new_name
        renamed.add(key)
        for a in res.atoms():
            if a.name in _GLYCAM_DROP_H[new_name]:
                h_to_drop.append(a)
        if verbose:
            print(f"  GLYCAM rename: {res.chain.id}:{res.id} → {new_name}")

    if h_to_drop:
        m = Modeller(topology, positions)
        m.delete(h_to_drop)
        return m.topology, m.positions, renamed
    return topology, positions, renamed

