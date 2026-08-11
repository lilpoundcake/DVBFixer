"""GROMACS topology construction for `dvbfixer top`.

Split out of `top/pipeline.py`: `TopologyBuilder` (protein/glycan/
glycolipid/small-molecule chain builder) plus its two private helpers
(`_match_atom_names`, `_resolve_sugar_rtp`) that are used only by it.
`main()` in `top/pipeline.py` still owns orchestration — reading the
PDB, classifying chains, calling `detect_glycan_links`/
`build_glycan_trees` (top/glycan.py), and writing output — this module
only builds the in-memory `ChainTopology` per chain.
"""
import sys
from collections import defaultdict
from pathlib import Path

from dvbfixer.rtp_parser import (
    parse_arn,
    parse_atomtypes,
    parse_r2b,
    parse_rtp,
    parse_tdb,
)
from dvbfixer.top.ff_data import (
    _EXPLICIT_RENAMES,
    CARB_ATOM_MAP,
    PDB_TO_CARB,
    PDB_TO_GMX,
    PDB_TO_LIPID,
    STANDARD_AA,
)
from dvbfixer.top.types import AtomEntry, ChainTopology


def _match_atom_names(rtp_names, pdb_names, arn_rtp_to_pdb=None):
    """Match RTP atom names to PDB atom names, handling naming conventions.

    PDB/OpenMM uses IUPAC naming (HB2/HB3), AMBER RTP uses old naming
    (HB1/HB2). For prochiral methylene hydrogens the numbering is shifted:
      AMBER HB1 = IUPAC HB2,  AMBER HB2 = IUPAC HB3

    Also handles:
    - ILE CD/HD1-3 → CD1/HD11-13
    - C-terminal OC1/OC2 → OXT/O
    - N-terminal H1 ← H (when only H exists in PDB)
    - ARN-based renames (e.g. CHARMM HN → PDB H)

    arn_rtp_to_pdb: dict[rtp_name -> pdb_name] from ARN file (reverse mapping).

    Returns dict[rtp_name -> pdb_name] for matched atoms.
    """
    rtp_set = set(rtp_names)
    pdb_set = set(pdb_names)
    mapping = {}
    used_pdb = set()

    # Pass 0a: ARN-based renames (e.g. CHARMM HN→H)
    if arn_rtp_to_pdb:
        for rtp_name in rtp_names:
            if rtp_name in arn_rtp_to_pdb:
                pdb_name = arn_rtp_to_pdb[rtp_name]
                if pdb_name in pdb_set and rtp_name not in pdb_set:
                    mapping[rtp_name] = pdb_name
                    used_pdb.add(pdb_name)

    # Pass 0b: explicit renames (ILE CD→CD1, C-terminal OC1→OXT, etc.)
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if rtp_name in _EXPLICIT_RENAMES:
            pdb_name = _EXPLICIT_RENAMES[rtp_name]
            if pdb_name in pdb_set and rtp_name not in pdb_set:
                mapping[rtp_name] = pdb_name
                used_pdb.add(pdb_name)

    # Detect if shift mapping applies for H-atom prefixes:
    # RTP has HB1,HB2 and PDB has HB2,HB3 → shift applies for HB prefix
    # Only shift when: lowest RTP number is NOT in PDB, but lowest+1 IS
    shift_prefixes = set()
    # Group RTP H-atoms by prefix
    prefix_groups = defaultdict(list)
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if len(rtp_name) >= 2 and rtp_name[-1].isdigit() and rtp_name[0] == 'H':
            prefix = rtp_name[:-1]
            num = int(rtp_name[-1])
            prefix_groups[prefix].append(num)

    for prefix, nums in prefix_groups.items():
        min_num = min(nums)
        max_num = max(nums)
        min_name = f"{prefix}{min_num}"
        max_shifted = f"{prefix}{max_num + 1}"
        # Shift applies if: lowest RTP name not in PDB AND highest+1 IS in PDB
        # This means the PDB numbering is offset by +1 from RTP for this prefix
        # e.g. RTP {HB1,HB2}, PDB has HB3 but not HB1 → shift
        if min_name not in pdb_set and max_shifted in pdb_set:
            shift_prefixes.add(prefix)

    # N-terminal special case: RTP H1 → PDB H (when PDB has H but not H1)
    # N-terminal residues have H1/H2/H3 in RTP but PDB may have H/H2/H3
    if 'H1' in rtp_set and 'H1' not in pdb_set and 'H' in pdb_set and 'H' not in rtp_set:
        mapping['H1'] = 'H'
        used_pdb.add('H')

    # Canonical PDB cap aliases.  The bundled RTPs retain force-field-native
    # HH31/HH32/HH33 and CH3 names, while OpenMM writes H1/H2/H3 and writes
    # the NME methyl carbon as C.
    for rtp_name, pdb_name in (("HH31", "H1"), ("HH32", "H2"),
                               ("HH33", "H3")):
        if (rtp_name in rtp_set and rtp_name not in mapping
                and pdb_name in pdb_set and pdb_name not in used_pdb):
            mapping[rtp_name] = pdb_name
            used_pdb.add(pdb_name)
    if ('CH3' in rtp_set and 'CH3' not in pdb_set and 'C' in pdb_set
            and 'C' not in used_pdb):
        mapping['CH3'] = 'C'
        used_pdb.add('C')

    # Pass 1: shift matching for detected prefixes
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if len(rtp_name) >= 2 and rtp_name[-1].isdigit() and rtp_name[0] == 'H':
            prefix = rtp_name[:-1]
            if prefix in shift_prefixes:
                num = int(rtp_name[-1])
                shifted = f"{prefix}{num + 1}"
                if shifted in pdb_set and shifted not in used_pdb:
                    mapping[rtp_name] = shifted
                    used_pdb.add(shifted)
                elif num == 1:
                    # H1 → H (N-terminal special case)
                    base = prefix
                    if base in pdb_set and base not in used_pdb:
                        mapping[rtp_name] = base
                        used_pdb.add(base)

    # Pass 2: exact matches for remaining
    for name in rtp_names:
        if name not in mapping and name in pdb_set and name not in used_pdb:
            mapping[name] = name
            used_pdb.add(name)

    # Pass 3: shift for any remaining unmatched
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if len(rtp_name) >= 2 and rtp_name[-1].isdigit():
            prefix = rtp_name[:-1]
            num = int(rtp_name[-1])
            shifted = f"{prefix}{num + 1}"
            if shifted in pdb_set and shifted not in used_pdb:
                mapping[rtp_name] = shifted
                used_pdb.add(shifted)

    # Pass 4: singleton numbered atom → base name (e.g. CHARMM HG1 → PDB HG)
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if len(rtp_name) >= 2 and rtp_name[-1].isdigit():
            prefix = rtp_name[:-1]
            # Only if this is the sole RTP atom with this prefix
            count = sum(1 for n in rtp_names if n.startswith(prefix) and n != prefix
                        and len(n) > len(prefix) and n[len(prefix):].isdigit())
            if count == 1 and prefix in pdb_set and prefix not in used_pdb:
                mapping[rtp_name] = prefix
                used_pdb.add(prefix)

    return mapping


def _resolve_sugar_rtp(resname, atoms, residues_dict):
    """Resolve a sugar PDB/CHARMM name to its correct RTP entry.

    Handles special cases like BGAL with N-acetyl atoms -> BGALNA.
    """
    # First try direct RTP match (CHARMM-GUI native names)
    if resname in residues_dict:
        rtp_name = resname
    else:
        rtp_name = PDB_TO_CARB.get(resname)

    if rtp_name is None:
        return None

    # Auto-detect: BGAL/AGAL with N-acetyl atoms -> BGALNA/AGALNA
    if rtp_name in ('BGAL', 'AGAL'):
        pdb_atom_names = {a[0] for a in atoms} if isinstance(atoms, list) else atoms
        has_nacetyl = bool(pdb_atom_names & {'N', 'HN', 'CT'})
        if has_nacetyl:
            alt = 'BGALNA' if rtp_name == 'BGAL' else 'AGALNA'
            if alt in residues_dict:
                rtp_name = alt

    if rtp_name not in residues_dict:
        return None

    return rtp_name


class TopologyBuilder:
    def __init__(self, ff_dir, ff_type='amber', verbose=False,
                 keep_all_hydrogens=False):
        self.ff_dir = Path(ff_dir)
        self.ff_type = ff_type
        self.verbose = verbose
        # When True, skip stripping HO1/HO2/HO3/HO4/HO6 at glycosidic
        # linkage sites and skip the associated charge redistribution.
        # Used for free reducing-end sugars where the H is real.
        self.keep_all_hydrogens = keep_all_hydrogens

        # Parse all FF files
        rtp_file = self.ff_dir / ('aminoacids.rtp' if ff_type == 'amber' else 'aminoacids.rtp')
        self.bonded_types, self.residues = parse_rtp(rtp_file)

        # For CHARMM, also parse merged.rtp if it exists (has more residues)
        merged_rtp = self.ff_dir / 'merged.rtp'
        if merged_rtp.exists():
            bt2, res2 = parse_rtp(merged_rtp)
            # merged.rtp is the main file for CHARMM, aminoacids.rtp may be subset
            self.residues.update(res2)
            self.bonded_types = bt2

        # Also load carb.rtp for sugar residues
        carb_rtp = self.ff_dir / 'carb.rtp'
        if carb_rtp.exists():
            carb_bt, carb_res = parse_rtp(carb_rtp)
            self.residues.update(carb_res)
            self.carb_bonded_types = carb_bt
        else:
            self.carb_bonded_types = None

        # Load all other molecule-type RTP files (CHARMM has lipids, NA, etc.)
        for rtp_name in ['lipid.rtp', 'na.rtp', 'cgenff.rtp', 'ethers.rtp',
                         'metals.rtp', 'silicates.rtp', 'solvent.rtp']:
            rtp_path = self.ff_dir / rtp_name
            if rtp_path.exists():
                _, extra_res = parse_rtp(rtp_path)
                self.residues.update(extra_res)

        self.r2b = parse_r2b(self.ff_dir / 'aminoacids.r2b')

        # Load all R2B files
        for r2b_name in ['carb.r2b', 'lipid.r2b', 'na.r2b', 'cgenff.r2b',
                         'ethers.r2b', 'metals.r2b', 'silicates.r2b', 'solvent.r2b']:
            r2b_path = self.ff_dir / r2b_name
            if r2b_path.exists():
                extra_r2b = parse_r2b(r2b_path)
                self.r2b.update(extra_r2b)

        self.arn = parse_arn(self.ff_dir / 'aminoacids.arn')
        self.atom_masses = parse_atomtypes(self.ff_dir / 'atomtypes.atp')

        # Terminal patches (CHARMM uses these, AMBER has empty TDB)
        n_tdb = self.ff_dir / 'aminoacids.n.tdb'
        c_tdb = self.ff_dir / 'aminoacids.c.tdb'
        self.n_patches = parse_tdb(n_tdb) if n_tdb.exists() else {}
        self.c_patches = parse_tdb(c_tdb) if c_tdb.exists() else {}

        # Build reverse ARN: (resname, ff_name) -> gromacs_name
        # For PDB->FF mapping we need (resname, pdb_name) -> ff_name
        # ARN file format is: resname gromacs_name ff_name
        # gromacs_name = what PDB uses, ff_name = what RTP uses
        self.arn_reverse = {}
        for (resname, gmx_name), ff_name in self.arn.items():
            self.arn_reverse[(resname, ff_name)] = gmx_name

    def _resolve_resname(self, pdb_resname, position, chain_ss_residues):
        """Map PDB residue name to RTP building block name.

        position: 'nter', 'cter', 'mid', 'twter' (both terminals = single residue chain)
        chain_ss_residues: set of resseq numbers involved in SS bonds
        """
        # Normalize non-canonical names
        gmx_name = PDB_TO_GMX.get(pdb_resname, pdb_resname)

        # Check if in r2b mapping
        if gmx_name not in self.r2b:
            # Try as-is (some residues use their PDB name directly)
            if pdb_resname in self.r2b:
                gmx_name = pdb_resname
            elif pdb_resname in self.residues:
                return pdb_resname  # Direct RTP match, no terminal variant needed
            else:
                return None

        main, nter, cter, twter = self.r2b[gmx_name]

        if position == 'twter' and twter != '-':
            return twter
        elif position == 'nter' and nter != '-':
            return nter
        elif position == 'cter' and cter != '-':
            return cter
        else:
            return main

    def _map_atom_name(self, rtp_resname, pdb_atom_name):
        """Map a PDB atom name to the name used in the RTP entry.

        ARN maps gromacs_name -> ff_name. PDB names match gromacs_names.
        We need to convert PDB name to RTP name (which is the ff_name).
        """
        # Check specific residue mapping first
        key = (rtp_resname, pdb_atom_name)
        if key in self.arn:
            return self.arn[key]

        # Check wildcard mapping
        key = ('*', pdb_atom_name)
        if key in self.arn:
            return self.arn[key]

        # No mapping needed — name is the same
        return pdb_atom_name

    def _get_pdb_name(self, rtp_resname, rtp_atom_name):
        """Get the PDB/GROMACS name for an RTP atom (reverse of _map_atom_name)."""
        key = (rtp_resname, rtp_atom_name)
        if key in self.arn_reverse:
            return self.arn_reverse[key]
        # Check wildcard
        for (resn, ff_name), gmx_name in self.arn_reverse.items():
            if resn == '*' and ff_name == rtp_atom_name:
                return gmx_name
        return rtp_atom_name

    def build_chain(self, chain, ss_residues=None, ss_pairs=None):
        """Build topology for a single chain.

        ss_residues: set of resseq involved in SS bonds in this chain.
        ss_pairs: list of (resseq1, resseq2) intra-chain SS bond pairs.
        """
        if ss_residues is None:
            ss_residues = set()
        if ss_pairs is None:
            ss_pairs = []

        n_res = len(chain.residues)
        if n_res == 0:
            return None

        # Step 1: Resolve RTP names for each residue
        rtp_names = []
        for i, res in enumerate(chain.residues):
            if n_res == 1:
                pos = 'twter'
            elif i == 0:
                pos = 'nter'
            elif i == n_res - 1:
                pos = 'cter'
            else:
                pos = 'mid'

            # Check SS bond
            pdb_name = res.resname
            if pdb_name == 'CYS' and res.resseq in ss_residues:
                pdb_name = 'CYX'

            rtp_name = self._resolve_resname(pdb_name, pos, ss_residues)
            if rtp_name is None:
                print(f"WARNING: Residue {res.resname} {res.chain_id}:{res.resseq} "
                      f"not found in force field", file=sys.stderr)
                return None

            if rtp_name not in self.residues:
                print(f"WARNING: RTP entry '{rtp_name}' not found for "
                      f"{res.resname} {res.chain_id}:{res.resseq}", file=sys.stderr)
                return None

            rtp_names.append(rtp_name)
            if self.verbose:
                label = f"  {res.chain_id}:{res.resname}{res.resseq}"
                if rtp_name != res.resname:
                    label += f" -> {rtp_name}"
                print(label)

        # Step 2: Build atom list from RTP entries, matched to PDB
        chain_top = ChainTopology(
            name=f"Protein_chain_{chain.chain_id.strip() or 'X'}",
            nrexcl=self.bonded_types.nrexcl,
        )

        # Map: (residue_index, rtp_atom_name) -> global atom index (1-based)
        atom_index_map = {}
        global_idx = 0
        cgnr_offset = 0

        for res_i, (res, rtp_name) in enumerate(zip(chain.residues, rtp_names)):
            rtp_res = self.residues[rtp_name]
            pdb_atom_names = {a[0] for a in res.atoms}
            rtp_atom_names = [a[0] for a in rtp_res.atoms]

            # Build per-residue ARN reverse mapping (rtp_name → pdb_name)
            arn_rtp_to_pdb = {}
            for rtp_aname in rtp_atom_names:
                # Check residue-specific ARN
                key = (rtp_name, rtp_aname)
                if key in self.arn_reverse:
                    arn_rtp_to_pdb[rtp_aname] = self.arn_reverse[key]
                # Check wildcard ARN
                for (resn, ff_name), gmx_name in self.arn_reverse.items():
                    if resn == '*' and ff_name == rtp_aname:
                        arn_rtp_to_pdb[rtp_aname] = gmx_name
                        break

            # Build RTP→PDB atom name mapping
            rtp_to_pdb = _match_atom_names(rtp_atom_names, pdb_atom_names,
                                           arn_rtp_to_pdb)

            # Build PDB atom coordinate lookup
            pdb_coords = {a[0]: (a[1], a[2], a[3]) for a in res.atoms}

            # Protonation H atoms (HD2/HE2/HD1) are added to res.atoms
            # by _add_protonation_hydrogens() before build_chain is called.
            # They should already be in pdb_coords via the res.atoms loop above.

            # Pre-scan: find skipped atoms and redistribute their charge
            # to bonded neighbors. Handles glycosylated ASN (HD21 or HD22
            # removed when ND2 bonds to sugar) and other missing atoms.
            rtp_bonds_local = {}
            for a1, a2 in rtp_res.bonds:
                if not a1.startswith(('+', '-')) and not a2.startswith(('+', '-')):
                    rtp_bonds_local.setdefault(a1, []).append(a2)
                    rtp_bonds_local.setdefault(a2, []).append(a1)
            skip_charge = {}
            for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
                pdb_name = rtp_to_pdb.get(atom_name)
                if pdb_name is None:
                    for bonded in rtp_bonds_local.get(atom_name, []):
                        if rtp_to_pdb.get(bonded) is not None:
                            skip_charge[bonded] = \
                                skip_charge.get(bonded, 0.0) + charge
                            break
                    if self.verbose:
                        dest = [b for b in rtp_bonds_local.get(atom_name, [])
                                if rtp_to_pdb.get(b) is not None]
                        print(f"    Skipping {rtp_name}:{atom_name} "
                              f"(not in PDB {res.chain_id}:{res.resseq},"
                              f" charge {charge:+.4f} → "
                              f"{dest[0] if dest else 'LOST'})")

            for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
                pdb_name = rtp_to_pdb.get(atom_name)
                if pdb_name is None:
                    continue

                charge = charge + skip_charge.get(atom_name, 0.0)

                global_idx += 1
                mass = self.atom_masses.get(atom_type, 0.0)
                x, y, z = pdb_coords.get(pdb_name, (0.0, 0.0, 0.0))
                chain_top.atoms.append(AtomEntry(
                    index=global_idx,
                    atom_type=atom_type,
                    resnr=res_i + 1,
                    resname=rtp_name,
                    atomname=pdb_name,
                    cgnr=cgnr + cgnr_offset,
                    charge=charge,
                    mass=mass,
                    x=x, y=y, z=z,
                    chain_id=res.chain_id,
                    orig_resseq=res.resseq,
                    orig_resname=res.resname,
                ))
                atom_index_map[(res_i, atom_name)] = global_idx

            cgnr_offset += max((cgnr for _, _, _, cgnr in rtp_res.atoms), default=0)

        # Step 3: Build bonds, resolving inter-residue references
        for res_i, rtp_name in enumerate(rtp_names):
            rtp_res = self.residues[rtp_name]
            for a1, a2 in rtp_res.bonds:
                idx1 = self._resolve_atom_ref(a1, res_i, atom_index_map, n_res)
                idx2 = self._resolve_atom_ref(a2, res_i, atom_index_map, n_res)
                if idx1 is not None and idx2 is not None:
                    bond = (min(idx1, idx2), max(idx1, idx2))
                    chain_top.bonds.append(bond)

        # Deduplicate bonds
        chain_top.bonds = sorted(set(chain_top.bonds))

        # Step 3b: Add intra-chain SS bonds (SG-SG)
        # RTP CYS2 has CB-SG but not the inter-residue SG-SG bond
        resseq_to_resi = {res.resseq: i for i, res in enumerate(chain.residues)}
        for res1_seq, res2_seq in ss_pairs:
            resi1 = resseq_to_resi.get(res1_seq)
            resi2 = resseq_to_resi.get(res2_seq)
            if resi1 is None or resi2 is None:
                continue
            sg1 = atom_index_map.get((resi1, 'SG'))
            sg2 = atom_index_map.get((resi2, 'SG'))
            if sg1 is not None and sg2 is not None:
                bond = (min(sg1, sg2), max(sg1, sg2))
                chain_top.bonds.append(bond)
                if self.verbose:
                    print(f"  Added intra-chain SS bond: "
                          f"CYS2 {chain.chain_id}:{res1_seq}:SG - "
                          f"CYS2 {chain.chain_id}:{res2_seq}:SG")
        chain_top.bonds = sorted(set(chain_top.bonds))

        # Step 4: Apply terminal patches (CHARMM)
        if self.ff_type == 'charmm' and n_res > 0:
            self._apply_terminal_patches(chain, chain_top, rtp_names, atom_index_map)

        # Step 5: Build bond graph for angle/dihedral generation
        adj = defaultdict(set)
        for i, j in chain_top.bonds:
            adj[i].add(j)
            adj[j].add(i)

        # Step 6: Generate angles
        chain_top.angles = self._generate_angles(adj)

        # Step 7: Generate proper dihedrals
        chain_top.dihedrals = self._generate_dihedrals(
            adj, rtp_names, atom_index_map, n_res
        )

        # Step 8: Generate 1-4 pairs from dihedrals
        chain_top.pairs = self._generate_pairs(chain_top.dihedrals)

        # Step 9: Resolve impropers from RTP
        chain_top.impropers = self._resolve_impropers(rtp_names, atom_index_map, n_res)

        # Step 10: Resolve CMAP (CHARMM)
        chain_top.cmap = self._resolve_cmap(rtp_names, atom_index_map, n_res)

        # Step 11: Renumber atoms to be contiguous (patches may create gaps)
        self._renumber_atoms(chain_top)

        return chain_top

    def _resolve_atom_ref(self, ref, res_i, atom_index_map, n_res):
        """Resolve atom reference like '-C', '+N', or 'CA' to global index."""
        if ref.startswith('-'):
            target_res = res_i - 1
            atom_name = ref[1:]
        elif ref.startswith('+'):
            target_res = res_i + 1
            atom_name = ref[1:]
        else:
            target_res = res_i
            atom_name = ref

        if target_res < 0 or target_res >= n_res:
            return None

        return atom_index_map.get((target_res, atom_name))

    def _renumber_atoms(self, chain_top):
        """Renumber all atom indices to be contiguous 1..N based on list order."""
        if not chain_top.atoms:
            return
        # Build remap based on current position in the atoms list
        remap = {}
        for new_idx, atom in enumerate(chain_top.atoms, 1):
            remap[atom.index] = new_idx

        for atom in chain_top.atoms:
            atom.index = remap[atom.index]

        def remap_tuple(t):
            return tuple(remap.get(x, x) if isinstance(x, int) else x for x in t)

        chain_top.bonds = [remap_tuple(b) for b in chain_top.bonds]
        chain_top.pairs = [remap_tuple(p) for p in chain_top.pairs]
        chain_top.angles = [remap_tuple(a) for a in chain_top.angles]
        chain_top.dihedrals = [remap_tuple(d) for d in chain_top.dihedrals]
        chain_top.impropers = [remap_tuple(i) for i in chain_top.impropers]
        chain_top.cmap = [remap_tuple(c) for c in chain_top.cmap]

    def _apply_terminal_patches(self, chain, chain_top, rtp_names, atom_index_map):
        """Apply CHARMM terminal patches (TDB). Only for protein chains."""
        n_res = len(chain.residues)

        first_res = chain.residues[0].resname
        last_res = chain.residues[-1].resname

        # Treat the two ends independently.  A previous all-or-nothing guard
        # skipped both patches whenever ACE was first, which also left an
        # otherwise uncapped C terminus unpatched (and vice versa for NME).
        if first_res != 'ACE':
            if first_res == 'GLY':
                n_patch_name = 'GLY-NH3+'
            elif first_res == 'PRO':
                n_patch_name = 'PRO-NH2+'
            elif first_res in STANDARD_AA or first_res in PDB_TO_GMX:
                n_patch_name = 'NH3+'
            else:
                n_patch_name = None
            if n_patch_name in self.n_patches:
                self._apply_patch(self.n_patches[n_patch_name], 0,
                                  chain_top, atom_index_map)

        if last_res != 'NME' and (last_res in STANDARD_AA or last_res in PDB_TO_GMX):
            if 'COO-' in self.c_patches:
                self._apply_patch(self.c_patches['COO-'], n_res - 1,
                                  chain_top, atom_index_map)

    def _apply_patch(self, patch, res_i, chain_top, atom_index_map):
        """Apply a single terminal patch to the chain topology."""
        # Save coordinates of atoms that will be deleted (for reuse by added atoms)
        deleted_coords = {}
        for atom_name in patch.delete:
            idx = atom_index_map.get((res_i, atom_name))
            if idx is not None:
                for atom in chain_top.atoms:
                    if atom.index == idx:
                        deleted_coords[atom_name] = (atom.x, atom.y, atom.z)
                        break

        # Delete atoms
        delete_indices = set()
        for atom_name in patch.delete:
            idx = atom_index_map.get((res_i, atom_name))
            if idx is not None:
                delete_indices.add(idx)

        if delete_indices:
            chain_top.atoms = [a for a in chain_top.atoms if a.index not in delete_indices]
            chain_top.bonds = [(i, j) for i, j in chain_top.bonds
                               if i not in delete_indices and j not in delete_indices]

        # Replace atoms (change type, mass, charge)
        for name, new_type, mass, charge in patch.replace:
            idx = atom_index_map.get((res_i, name))
            if idx is not None:
                for atom in chain_top.atoms:
                    if atom.index == idx:
                        atom.atom_type = new_type
                        atom.mass = mass
                        atom.charge = charge
                        break

        # Add atoms — insert after reference atom in the residue
        if patch.add:
            max_idx = max(a.index for a in chain_top.atoms) if chain_top.atoms else 0
            resnr = None
            resname = 'UNK'
            for atom in chain_top.atoms:
                if atom_index_map.get((res_i, atom.atomname)) == atom.index:
                    resnr = atom.resnr
                    resname = atom.resname
                    break

            for add_entry in patch.add:
                count = add_entry['count']
                name_base = add_entry['name']
                atype = add_entry['type']
                mass = add_entry['mass']
                charge = add_entry['charge']
                cgnr_val = add_entry['cgnr']

                # Find insertion position: after the first reference atom
                ref_atoms = add_entry['ref_atoms']
                insert_after_idx = None
                if ref_atoms:
                    insert_after_idx = atom_index_map.get((res_i, ref_atoms[0]))

                # Find position in atoms list to insert + get reference coords
                insert_pos = len(chain_top.atoms)
                ref_x, ref_y, ref_z = 0.0, 0.0, 0.0
                ref_chain_id = ' '
                if insert_after_idx is not None:
                    for pos, atom in enumerate(chain_top.atoms):
                        if atom.index == insert_after_idx:
                            insert_pos = pos + 1
                            ref_x, ref_y, ref_z = atom.x, atom.y, atom.z
                            ref_chain_id = atom.chain_id
                            break

                # Build coordinate list for added atoms:
                # - Use deleted atom coords when available (e.g. OT1 gets O's coords)
                # - Otherwise use reference atom coords with small offsets to avoid overlap
                add_coords = []
                for k in range(count):
                    atom_name = f"{name_base}{k + 1}" if count > 1 else name_base
                    # Try to find a matching deleted atom's coordinates
                    # COO-: deleted O → use for OT1; offset for OT2
                    # NH3+: deleted HN → use for H1; offset for H2, H3
                    coord_found = False
                    if k == 0:
                        # First added atom: try deleted atoms that share the same
                        # element (O→OT1, HN→H1)
                        for dname, dcoord in deleted_coords.items():
                            # Match by element: H* deleted → first H added
                            if dname[0] == name_base[0]:
                                add_coords.append(dcoord)
                                coord_found = True
                                break
                    if not coord_found:
                        # Offset from reference atom to avoid overlap
                        offset = 0.1 * (k + 1)  # 1 Angstrom increments
                        add_coords.append((ref_x + offset, ref_y + offset, ref_z))

                new_atoms = []
                for k in range(count):
                    max_idx += 1
                    atom_name = f"{name_base}{k + 1}" if count > 1 else name_base
                    actual_cgnr = cgnr_val if cgnr_val > 0 else (
                        chain_top.atoms[insert_pos - 1].cgnr if insert_pos > 0 else 1
                    )
                    ax, ay, az = add_coords[k]
                    new_atom = AtomEntry(
                        index=max_idx,
                        atom_type=atype,
                        resnr=resnr or (res_i + 1),
                        resname=resname,
                        atomname=atom_name,
                        cgnr=actual_cgnr,
                        charge=charge,
                        mass=mass,
                        x=ax, y=ay, z=az,
                        chain_id=ref_chain_id,
                        orig_resseq=chain_top.atoms[insert_pos - 1].orig_resseq if insert_pos > 0 else 0,
                        orig_resname=chain_top.atoms[insert_pos - 1].orig_resname if insert_pos > 0 else '',
                    )
                    new_atoms.append(new_atom)
                    atom_index_map[(res_i, atom_name)] = max_idx

                    # Add bond to reference atom
                    if insert_after_idx is not None:
                        bond = (min(insert_after_idx, max_idx), max(insert_after_idx, max_idx))
                        chain_top.bonds.append(bond)

                # Insert at the right position
                for offset, new_atom in enumerate(new_atoms):
                    chain_top.atoms.insert(insert_pos + offset, new_atom)

        # Add impropers from patch
        for imp in patch.impropers:
            indices = []
            for atom_name in imp:
                idx = atom_index_map.get((res_i, atom_name))
                if idx is not None:
                    indices.append(idx)
            if len(indices) == 4:
                chain_top.impropers.append(tuple(indices))

    def _generate_angles(self, adj):
        """Generate all angles from bond graph."""
        angles = set()
        for j in sorted(adj.keys()):
            neighbors = sorted(adj[j])
            for idx_a, i in enumerate(neighbors):
                for k in neighbors[idx_a + 1:]:
                    angles.add((i, j, k))
        return sorted(angles)

    def _generate_dihedrals(self, adj, rtp_names, atom_index_map, n_res):
        """Generate proper dihedrals from bond graph + RTP explicit dihedrals."""
        # Generate all possible dihedrals from connectivity
        generated = set()
        for j in sorted(adj.keys()):
            for k in sorted(adj[j]):
                if k <= j:
                    continue
                for i in sorted(adj[j]):
                    if i == k:
                        continue
                    for l in sorted(adj[k]):
                        if l == j:
                            continue
                        generated.add((i, j, k, l))

        dihedrals = sorted(generated)

        # Add explicit RTP dihedrals (AMBER ILDN corrections etc.)
        explicit = []
        for res_i, rtp_name in enumerate(rtp_names):
            rtp_res = self.residues[rtp_name]
            for dih in rtp_res.dihedrals:
                indices = []
                for ref in dih[:4]:
                    idx = self._resolve_atom_ref(ref, res_i, atom_index_map, n_res)
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 4:
                    if len(dih) == 5:
                        explicit.append((*indices, dih[4]))
                    else:
                        explicit.append(tuple(indices))

        # Merge: explicit dihedrals with type names go at the end
        result = [(i, j, k, l) for i, j, k, l in dihedrals]
        result.extend(explicit)
        return result

    def _generate_pairs(self, dihedrals):
        """Generate 1-4 pairs from proper dihedrals."""
        pairs = set()
        for dih in dihedrals:
            i, l = dih[0], dih[3]
            pair = (min(i, l), max(i, l))
            pairs.add(pair)
        return sorted(pairs)

    def _resolve_impropers(self, rtp_names, atom_index_map, n_res):
        """Resolve all improper dihedrals from RTP entries."""
        impropers = []
        for res_i, rtp_name in enumerate(rtp_names):
            rtp_res = self.residues[rtp_name]
            for imp in rtp_res.impropers:
                indices = []
                for ref in imp:
                    idx = self._resolve_atom_ref(ref, res_i, atom_index_map, n_res)
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 4:
                    impropers.append(tuple(indices))
        return impropers

    def _resolve_cmap(self, rtp_names, atom_index_map, n_res):
        """Resolve CMAP entries from RTP."""
        cmaps = []
        for res_i, rtp_name in enumerate(rtp_names):
            rtp_res = self.residues[rtp_name]
            for cm in rtp_res.cmap:
                indices = []
                for ref in cm:
                    idx = self._resolve_atom_ref(ref, res_i, atom_index_map, n_res)
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 5:
                    cmaps.append(tuple(indices))
        return cmaps

    def build_glycan_chain(self, tree, link_atoms, all_chains, glycan_links):
        """Build topology for a glycan tree (one moleculetype).

        tree: list of (chain_id, resseq) in topological order
        link_atoms: dict (chain_id, resseq) -> set of O atoms that are linked
        all_chains: list of PDBChain objects
        glycan_links: list of (don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom)
        """
        # Build residue lookup from PDB chains
        res_lookup = {}  # (chain, resseq) -> PDBResidue
        for chain in all_chains:
            for res in chain.residues:
                res_lookup[(chain.chain_id, res.resseq)] = res

        bt = self.carb_bonded_types or self.bonded_types
        first_ch, first_rs = tree[0]
        first_res = res_lookup.get((first_ch, first_rs))
        chain_name = f"Glycan_{first_ch.strip() or 'X'}_{first_rs}"

        chain_top = ChainTopology(
            name=chain_name,
            nrexcl=bt.nrexcl,
        )

        atom_index_map = {}  # (tree_idx, rtp_atom_name) -> global atom index
        global_idx = 0

        for tree_idx, (ch, rs) in enumerate(tree):
            res = res_lookup.get((ch, rs))
            if res is None:
                print(f"WARNING: Sugar {ch}:{rs} not found in PDB", file=sys.stderr)
                continue

            # Map PDB name -> CHARMM name (with auto-detect for BGALNA etc.)
            rtp_name = _resolve_sugar_rtp(res.resname, res.atoms, self.residues)
            if rtp_name is None:
                print(f"WARNING: No CHARMM RTP for {res.resname} ({ch}:{rs})",
                      file=sys.stderr)
                continue

            rtp_res = self.residues[rtp_name]
            pdb_atom_names = {a[0] for a in res.atoms}
            rtp_atom_names = [a[0] for a in rtp_res.atoms]

            # Determine which HO atoms to remove at linked positions
            # and redistribute their charge to the bonded O atom.
            # --keep-all-hydrogens skips this stripping (user opts to keep
            # every input H, e.g. for a free reducing end).
            linked_os = link_atoms.get((ch, rs), set())
            remove_ho = {}  # ho_name -> o_name (for charge transfer)
            rtp_charges = {a[0]: a[2] for a in rtp_res.atoms}
            if not self.keep_all_hydrogens:
                for o_name in linked_os:
                    ho_name = 'HO' + o_name[1:]
                    if ho_name in rtp_charges:
                        remove_ho[ho_name] = o_name

            # Build RTP->PDB mapping using sugar-specific atom map
            # Invert the PDB->CHARMM map to get CHARMM->PDB
            carb_rtp_to_pdb = {}
            pdb_resname = res.resname
            if pdb_resname in CARB_ATOM_MAP:
                for pdb_aname, charmm_aname in CARB_ATOM_MAP[pdb_resname].items():
                    carb_rtp_to_pdb[charmm_aname] = pdb_aname

            # Also add ARN-based renames
            for rtp_aname in rtp_atom_names:
                if rtp_aname in carb_rtp_to_pdb:
                    continue
                for (resn, ff_name), gmx_name in self.arn_reverse.items():
                    if resn == '*' and ff_name == rtp_aname:
                        carb_rtp_to_pdb[rtp_aname] = gmx_name
                        break

            rtp_to_pdb = _match_atom_names(rtp_atom_names, pdb_atom_names,
                                           carb_rtp_to_pdb)
            pdb_coords = {a[0]: (a[1], a[2], a[3]) for a in res.atoms}

            # Build charge adjustments: O gets HO's charge when HO is removed
            charge_adjust = {}  # o_name -> extra charge
            # When glycosidic bond forms, O type changes from hydroxyl to ether
            type_change = {}  # o_name -> new_type
            for ho_name, o_name in remove_ho.items():
                pdb_o_name = rtp_to_pdb.get(o_name, o_name)
                if pdb_o_name not in pdb_coords and o_name in linked_os:
                    # O not in PDB (removed by CHARMM-GUI at linkage):
                    # redistribute O + HO combined charge to anomeric carbon
                    c_name = 'C2' if rtp_name in ('ANE5AC', 'BNE5AC') else 'C1'
                    charge_adjust[c_name] = charge_adjust.get(c_name, 0.0) + \
                        rtp_charges.get(o_name, 0.0) + rtp_charges[ho_name]
                else:
                    charge_adjust[o_name] = charge_adjust.get(o_name, 0.0) + rtp_charges[ho_name]
                    # OC311 (hydroxyl) -> OC3C61 (ether) for linked O
                    type_change[o_name] = 'OC3C61'

            # Defensive: if O1/O2 is not in PDB but wasn't detected as linked,
            # still redistribute its charge to anomeric C. Also skip its HO if
            # the HO is also not in PDB (both removed at glycosidic bond site).
            # --keep-all-hydrogens skips this defensive redistribution too.
            if not self.keep_all_hydrogens:
                for o_name in ('O1', 'O2'):
                    if o_name in linked_os:
                        continue  # already handled above
                    pdb_o = rtp_to_pdb.get(o_name, o_name)
                    if pdb_o not in pdb_coords and o_name in rtp_charges:
                        c_name = 'C2' if rtp_name in ('ANE5AC', 'BNE5AC') else 'C1'
                        charge_adjust[c_name] = charge_adjust.get(c_name, 0.0) + \
                            rtp_charges[o_name]
                        # Also handle corresponding HO if it's also not in PDB
                        ho_name = 'HO' + o_name[1:]
                        if ho_name in rtp_charges:
                            ho_pdb = rtp_to_pdb.get(ho_name, ho_name)
                            if ho_pdb not in pdb_coords:
                                charge_adjust[c_name] += rtp_charges[ho_name]
                                remove_ho[ho_name] = o_name

            for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
                # Skip HO atoms at linked positions
                if atom_name in remove_ho:
                    continue

                # Apply charge redistribution and type change for linked O atoms
                adj_charge = charge + charge_adjust.get(atom_name, 0.0)
                adj_type = type_change.get(atom_name, atom_type)

                pdb_name = rtp_to_pdb.get(atom_name, atom_name)
                # For carbs, skip atoms not in PDB: H atoms, linked O atoms,
                # and any O1/O2 not present (glycosidic bond sites where
                # CHARMM-GUI removes the bridging O)
                if pdb_name not in pdb_coords:
                    if (atom_name.startswith('H') or atom_name in linked_os
                            or atom_name in ('O1', 'O2')):
                        if self.verbose:
                            print(f"    Skipping {rtp_name}:{atom_name} "
                                  f"(not in PDB {ch}:{rs})")
                        continue

                global_idx += 1
                mass = self.atom_masses.get(adj_type, 0.0)
                x, y, z = pdb_coords.get(pdb_name, (0.0, 0.0, 0.0))
                chain_top.atoms.append(AtomEntry(
                    index=global_idx,
                    atom_type=adj_type,
                    resnr=tree_idx + 1,
                    resname=rtp_name,
                    atomname=pdb_name if pdb_name in pdb_coords else atom_name,
                    cgnr=cgnr,
                    charge=adj_charge,
                    mass=mass,
                    x=x, y=y, z=z,
                    chain_id=ch,
                    orig_resseq=rs,
                    orig_resname=res.resname,
                ))
                atom_index_map[(tree_idx, atom_name)] = global_idx

            # Intra-residue bonds from RTP (skip bonds involving removed atoms)
            for a1, a2 in rtp_res.bonds:
                if a1 in remove_ho or a2 in remove_ho:
                    continue
                idx1 = atom_index_map.get((tree_idx, a1))
                idx2 = atom_index_map.get((tree_idx, a2))
                if idx1 is not None and idx2 is not None:
                    chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # Add inter-residue glycosidic bonds
        tree_pos = {(ch, rs): i for i, (ch, rs) in enumerate(tree)}
        for don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom in glycan_links:
            don_key = (don_ch, don_rs)
            acc_key = (acc_ch, acc_rs)
            if don_key in tree_pos and acc_key in tree_pos:
                don_idx = tree_pos[don_key]
                acc_idx = tree_pos[acc_key]
                idx1 = atom_index_map.get((don_idx, don_atom))
                idx2 = atom_index_map.get((acc_idx, acc_atom))
                if idx1 is not None and idx2 is not None:
                    chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # Deduplicate bonds
        chain_top.bonds = sorted(set(chain_top.bonds))

        # Build bond graph and generate angles/dihedrals/pairs
        adj = defaultdict(set)
        for i, j in chain_top.bonds:
            adj[i].add(j)
            adj[j].add(i)

        # Generate angles — all use default ftype from [ bondedtypes ]
        chain_top.angles = list(self._generate_angles(adj))

        # For carbs, generate all dihedrals from connectivity
        generated = set()
        for j in sorted(adj.keys()):
            for k in sorted(adj[j]):
                if k <= j:
                    continue
                for i in sorted(adj[j]):
                    if i == k:
                        continue
                    for l in sorted(adj[k]):
                        if l == j:
                            continue
                        generated.add((i, j, k, l))
        chain_top.dihedrals = sorted(generated)

        chain_top.pairs = self._generate_pairs(chain_top.dihedrals)

        # Resolve impropers from RTP
        for tree_idx, (ch, rs) in enumerate(tree):
            res = res_lookup.get((ch, rs))
            if res is None:
                continue
            rtp_name = _resolve_sugar_rtp(res.resname, res.atoms, self.residues)
            if rtp_name is None:
                continue
            rtp_res = self.residues[rtp_name]
            for imp in rtp_res.impropers:
                indices = []
                for ref in imp:
                    idx = atom_index_map.get((tree_idx, ref))
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 4:
                    chain_top.impropers.append(tuple(indices))

        # Renumber
        self._renumber_atoms(chain_top)

        return chain_top

    def build_glycolipid_chain(self, ceramide_res, tree, link_atoms,
                               all_chains, glycan_links, ceramide_link):
        """Build topology for a glycolipid (ceramide + sugar tree as one moleculetype).

        ceramide_res: PDBResidue for the ceramide
        tree: list of (chain_id, resseq) for sugars in topological order
        link_atoms: dict (chain_id, resseq) -> set of linked O atoms
        all_chains: list of PDBChain objects
        glycan_links: list of sugar-sugar links
        ceramide_link: (cer_chain, cer_resseq, cer_atom, sugar_chain, sugar_resseq, sugar_atom)
        """
        # Build residue lookup
        res_lookup = {}
        for chain in all_chains:
            for res in chain.residues:
                res_lookup[(chain.chain_id, res.resseq)] = res

        bt = self.carb_bonded_types or self.bonded_types
        cer_ch = ceramide_link[0]
        cer_rs = ceramide_link[1]
        chain_name = f"Glycolipid_{cer_ch.strip() or 'X'}_{cer_rs}"

        chain_top = ChainTopology(
            name=chain_name,
            nrexcl=bt.nrexcl,
        )

        atom_index_map = {}  # (resnr, rtp_atom_name) -> global atom index
        global_idx = 0

        # --- Step 1: Build ceramide residue ---
        cer_rtp_name = PDB_TO_LIPID.get(ceramide_res.resname, ceramide_res.resname)
        if cer_rtp_name not in self.residues:
            print(f"WARNING: No RTP entry for ceramide {ceramide_res.resname} "
                  f"(tried {cer_rtp_name})", file=sys.stderr)
            return None

        rtp_res = self.residues[cer_rtp_name]
        pdb_atom_names = {a[0] for a in ceramide_res.atoms}
        pdb_coords = {a[0]: (a[1], a[2], a[3]) for a in ceramide_res.atoms}
        rtp_charges = {a[0]: a[2] for a in rtp_res.atoms}

        # At linkage: remove ceramide HO1 and O1 (not in PDB — CHARMM-GUI
        # already removed them). Redistribute their combined charge to C1S.
        # --keep-all-hydrogens skips this to preserve the input H verbatim.
        cer_remove_ho = {}
        cer_charge_adjust = {}
        cer_type_change = {}
        if 'HO1' in rtp_charges and not self.keep_all_hydrogens:
            cer_remove_ho['HO1'] = 'O1'
            # O1 not in PDB: redistribute O1+HO1 combined charge to C1S
            if 'O1' not in pdb_atom_names:
                o1_charge = rtp_charges.get('O1', 0.0)
                ho1_charge = rtp_charges['HO1']
                cer_charge_adjust['C1S'] = o1_charge + ho1_charge
            else:
                # O1 is in PDB: standard redistribution
                cer_charge_adjust['O1'] = rtp_charges['HO1']
                cer_type_change['O1'] = 'OC301'

        resnr = 1  # ceramide is residue 1
        for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
            if atom_name in cer_remove_ho:
                continue

            adj_charge = charge + cer_charge_adjust.get(atom_name, 0.0)
            adj_type = cer_type_change.get(atom_name, atom_type)

            # Skip atoms not in PDB (H atoms, and O1 if already removed)
            if atom_name not in pdb_coords:
                if atom_name.startswith('H') or atom_name in ('O1', 'HO1'):
                    if self.verbose:
                        print(f"    Skipping {cer_rtp_name}:{atom_name} (not in PDB)")
                    continue

            global_idx += 1
            mass = self.atom_masses.get(adj_type, 0.0)
            x, y, z = pdb_coords.get(atom_name, (0.0, 0.0, 0.0))
            chain_top.atoms.append(AtomEntry(
                index=global_idx,
                atom_type=adj_type,
                resnr=resnr,
                resname=cer_rtp_name,
                atomname=atom_name,
                cgnr=cgnr,
                charge=adj_charge,
                mass=mass,
                x=x, y=y, z=z,
                chain_id=cer_ch,
                orig_resseq=cer_rs,
                orig_resname=ceramide_res.resname,
            ))
            atom_index_map[(0, atom_name)] = global_idx

        # Intra-residue bonds for ceramide
        for a1, a2 in rtp_res.bonds:
            if a1 in cer_remove_ho or a2 in cer_remove_ho:
                continue
            idx1 = atom_index_map.get((0, a1))
            idx2 = atom_index_map.get((0, a2))
            if idx1 is not None and idx2 is not None:
                chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # Ceramide impropers
        for imp in rtp_res.impropers:
            indices = []
            for ref in imp:
                idx = atom_index_map.get((0, ref))
                if idx is not None:
                    indices.append(idx)
            if len(indices) == 4:
                chain_top.impropers.append(tuple(indices))

        # --- Step 2: Build sugar tree (reusing glycan chain logic) ---
        for tree_idx, (ch, rs) in enumerate(tree):
            res = res_lookup.get((ch, rs))
            if res is None:
                print(f"WARNING: Sugar {ch}:{rs} not found in PDB", file=sys.stderr)
                continue

            # Map PDB name -> CHARMM name (with auto-detect for BGALNA etc.)
            rtp_name = _resolve_sugar_rtp(res.resname, res.atoms, self.residues)
            if rtp_name is None:
                print(f"WARNING: No CHARMM RTP for {res.resname} ({ch}:{rs})",
                      file=sys.stderr)
                continue

            rtp_res = self.residues[rtp_name]
            pdb_atom_names_set = {a[0] for a in res.atoms}
            rtp_atom_names = [a[0] for a in rtp_res.atoms]

            # Determine linked O atoms; skip HO-strip if --keep-all-hydrogens
            linked_os = link_atoms.get((ch, rs), set())
            remove_ho = {}
            rtp_charges = {a[0]: a[2] for a in rtp_res.atoms}
            if not self.keep_all_hydrogens:
                for o_name in linked_os:
                    ho_name = 'HO' + o_name[1:]
                    if ho_name in rtp_charges:
                        remove_ho[ho_name] = o_name

            # Build RTP->PDB atom name mapping
            carb_rtp_to_pdb = {}
            pdb_resname = res.resname
            if pdb_resname in CARB_ATOM_MAP:
                for pdb_aname, charmm_aname in CARB_ATOM_MAP[pdb_resname].items():
                    carb_rtp_to_pdb[charmm_aname] = pdb_aname

            for rtp_aname in rtp_atom_names:
                if rtp_aname in carb_rtp_to_pdb:
                    continue
                for (resn, ff_name), gmx_name in self.arn_reverse.items():
                    if resn == '*' and ff_name == rtp_aname:
                        carb_rtp_to_pdb[rtp_aname] = gmx_name
                        break

            rtp_to_pdb = _match_atom_names(rtp_atom_names, pdb_atom_names_set,
                                           carb_rtp_to_pdb)
            pdb_coords = {a[0]: (a[1], a[2], a[3]) for a in res.atoms}

            charge_adjust = {}
            type_change = {}
            # Determine which O atom links to ceramide (if any)
            cer_linked_o = ceramide_link[5] if (ch, rs) == (ceramide_link[3], ceramide_link[4]) else None
            for ho_name, o_name in remove_ho.items():
                pdb_o_name = rtp_to_pdb.get(o_name, o_name)
                if pdb_o_name not in pdb_coords and o_name in linked_os:
                    # O not in PDB (removed at linkage): redistribute to anomeric C
                    c_name = 'C2' if rtp_name in ('ANE5AC', 'BNE5AC') else 'C1'
                    charge_adjust[c_name] = charge_adjust.get(c_name, 0.0) + \
                        rtp_charges.get(o_name, 0.0) + rtp_charges[ho_name]
                else:
                    charge_adjust[o_name] = charge_adjust.get(o_name, 0.0) + rtp_charges[ho_name]
                    # Ceramide-linked O becomes OC301 (linear ether),
                    # sugar-sugar linked O becomes OC3C61 (cyclic ether)
                    if o_name == cer_linked_o:
                        type_change[o_name] = 'OC301'
                    else:
                        type_change[o_name] = 'OC3C61'

            # Defensive: redistribute charge of O1/O2 not in PDB even if not
            # in linked_os. --keep-all-hydrogens skips this too.
            if not self.keep_all_hydrogens:
                for o_name in ('O1', 'O2'):
                    if o_name in linked_os:
                        continue
                    pdb_o = rtp_to_pdb.get(o_name, o_name)
                    if pdb_o not in pdb_coords and o_name in rtp_charges:
                        c_name = 'C2' if rtp_name in ('ANE5AC', 'BNE5AC') else 'C1'
                        charge_adjust[c_name] = charge_adjust.get(c_name, 0.0) + \
                            rtp_charges[o_name]
                        ho_name = 'HO' + o_name[1:]
                        if ho_name in rtp_charges:
                            ho_pdb = rtp_to_pdb.get(ho_name, ho_name)
                            if ho_pdb not in pdb_coords:
                                charge_adjust[c_name] += rtp_charges[ho_name]
                                remove_ho[ho_name] = o_name

            resnr_sugar = tree_idx + 2  # ceramide is resnr 1
            for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
                if atom_name in remove_ho:
                    continue

                adj_charge = charge + charge_adjust.get(atom_name, 0.0)
                adj_type = type_change.get(atom_name, atom_type)

                pdb_name = rtp_to_pdb.get(atom_name, atom_name)
                if pdb_name not in pdb_coords:
                    if (atom_name.startswith('H') or atom_name in linked_os
                            or atom_name in ('O1', 'O2')):
                        if self.verbose:
                            print(f"    Skipping {rtp_name}:{atom_name} "
                                  f"(not in PDB {ch}:{rs})")
                        continue

                global_idx += 1
                mass = self.atom_masses.get(adj_type, 0.0)
                x, y, z = pdb_coords.get(pdb_name, (0.0, 0.0, 0.0))
                chain_top.atoms.append(AtomEntry(
                    index=global_idx,
                    atom_type=adj_type,
                    resnr=resnr_sugar,
                    resname=rtp_name,
                    atomname=pdb_name if pdb_name in pdb_coords else atom_name,
                    cgnr=cgnr,
                    charge=adj_charge,
                    mass=mass,
                    x=x, y=y, z=z,
                    chain_id=ch,
                    orig_resseq=rs,
                    orig_resname=res.resname,
                ))
                atom_index_map[(tree_idx + 1, atom_name)] = global_idx

            # Intra-residue bonds
            for a1, a2 in rtp_res.bonds:
                if a1 in remove_ho or a2 in remove_ho:
                    continue
                idx1 = atom_index_map.get((tree_idx + 1, a1))
                idx2 = atom_index_map.get((tree_idx + 1, a2))
                if idx1 is not None and idx2 is not None:
                    chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # --- Step 3: Add inter-residue bonds ---
        # Ceramide C1S — sugar O1 bond (sugar O1 bridges ceramide C1S to sugar C1)
        cer_atom = ceramide_link[2]  # e.g. 'C1S'
        sugar_atom = ceramide_link[5]  # e.g. 'O1'
        cer_idx = atom_index_map.get((0, cer_atom))
        root_sugar_idx = atom_index_map.get((1, sugar_atom))
        if cer_idx is not None and root_sugar_idx is not None:
            chain_top.bonds.append((min(cer_idx, root_sugar_idx),
                                    max(cer_idx, root_sugar_idx)))

        # Sugar-sugar glycosidic bonds
        tree_pos = {(ch, rs): i + 1 for i, (ch, rs) in enumerate(tree)}
        for don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom in glycan_links:
            don_key = (don_ch, don_rs)
            acc_key = (acc_ch, acc_rs)
            if don_key in tree_pos and acc_key in tree_pos:
                don_tidx = tree_pos[don_key]
                acc_tidx = tree_pos[acc_key]
                idx1 = atom_index_map.get((don_tidx, don_atom))
                idx2 = atom_index_map.get((acc_tidx, acc_atom))
                if idx1 is not None and idx2 is not None:
                    chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # Deduplicate bonds
        chain_top.bonds = sorted(set(chain_top.bonds))

        # --- Step 4: Build bond graph and enumerate angles/dihedrals/pairs ---
        adj = defaultdict(set)
        for i, j in chain_top.bonds:
            adj[i].add(j)
            adj[j].add(i)

        chain_top.angles = list(self._generate_angles(adj))

        generated = set()
        for j in sorted(adj.keys()):
            for k in sorted(adj[j]):
                if k <= j:
                    continue
                for i in sorted(adj[j]):
                    if i == k:
                        continue
                    for l in sorted(adj[k]):
                        if l == j:
                            continue
                        generated.add((i, j, k, l))
        chain_top.dihedrals = sorted(generated)

        chain_top.pairs = self._generate_pairs(chain_top.dihedrals)

        # Sugar impropers
        for tree_idx, (ch, rs) in enumerate(tree):
            res = res_lookup.get((ch, rs))
            if res is None:
                continue
            rtp_name = _resolve_sugar_rtp(res.resname, res.atoms, self.residues)
            if rtp_name is None:
                continue
            rtp_res = self.residues[rtp_name]
            for imp in rtp_res.impropers:
                indices = []
                for ref in imp:
                    idx = atom_index_map.get((tree_idx + 1, ref))
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 4:
                    chain_top.impropers.append(tuple(indices))

        # Renumber
        self._renumber_atoms(chain_top)

        return chain_top
