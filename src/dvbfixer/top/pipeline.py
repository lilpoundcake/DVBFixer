"""dvbfixer top — Generate GROMACS .itp/.top topology files from PDB.

Parses GROMACS force field RTP/ARN/R2B/TDB files directly and builds
correct topology with proper atom types, charges, bonds, angles,
dihedrals, impropers, and CMAP (CHARMM).
"""

import sys
from collections import defaultdict
from pathlib import Path

from dvbfixer.top.cli import FF_CHOICES, FF_DIR, parse_args
from dvbfixer.top.ff_data import (
    _KNOWN_4CHAR_RESNAMES,
    _WATER_ATOMS_PER_MOL,
    _WATER_DEFAULT_ION_SET,
    _WATER_ION_ALIAS,
    _WATER_RESNAMES,
    CERAMIDE_RTP,
    PDB_TO_CARB,
    PDB_TO_GMX,
    PDB_TO_LIPID,
    STANDARD_AA,
)
from dvbfixer.top.glycan import _is_ceramide, build_glycan_trees, detect_glycan_links
from dvbfixer.top.topology_builder import TopologyBuilder
from dvbfixer.top.types import AtomEntry, ChainTopology, PDBChain, PDBResidue
from dvbfixer.top.writers import (
    _write_moleculetype,
    write_pdb,
    write_posre,
    write_top,
)


# ---------------------------------------------------------------------------
# PDB reader
# ---------------------------------------------------------------------------
def read_pdb_chains(path):
    """Read PDB file and extract chains with residues and atoms.

    Detects resseq backward jumps within a chain (e.g. two glycan trees
    with same chain ID and overlapping residue numbers) and splits them
    into separate sub-chains with generated chain IDs.
    """
    # First pass: collect lines per original chain ID, preserving order
    # Handle 4-char resnames (CHARMM-GUI style): if col 21 is not a space
    # and col 17-20 is not blank, the resname extends to col 21 and chain ID
    # is effectively blank.
    chain_lines = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            # Detect 4-char resnames: col 17-20 is the resname field,
            # col 21 is normally chain ID. If col 20 is not space and col 21
            # is not space either, it's likely a 4-char resname with no chain.
            resname_3 = line[17:20].strip()
            chain_id = line[21]
            if resname_3 and chain_id != ' ' and not chain_id.isalpha():
                # Could be 4-char resname (e.g. CER1, BGAL, ANE5, AGLC)
                resname_4 = line[17:21].strip()
                if (resname_4 in PDB_TO_LIPID or resname_4 in PDB_TO_CARB
                        or resname_4 in CERAMIDE_RTP
                        or (len(resname_4) == 4 and resname_3 not in STANDARD_AA)):
                    chain_id = ' '
            chain_lines[chain_id].append(line)

    # Second pass: split chains on resseq backward jumps
    chains = []
    used_ids = set(chain_lines.keys())
    # Pool of available chain IDs for sub-chains
    all_ids = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')

    for orig_id, lines in chain_lines.items():
        # Detect breaks (resseq goes backwards)
        segments = [[]]  # list of line groups
        prev_resseq = -999999
        for line in lines:
            resseq = int(line[22:26])
            if resseq < prev_resseq:
                segments.append([])
            segments[-1].append(line)
            prev_resseq = resseq

        for seg_idx, seg_lines in enumerate(segments):
            if seg_idx == 0:
                cid = orig_id
            else:
                # Assign a new chain ID
                cid = None
                for candidate in all_ids:
                    if candidate not in used_ids:
                        cid = candidate
                        break
                if cid is None:
                    cid = f'{orig_id}{seg_idx}'
                used_ids.add(cid)

            chain = PDBChain(chain_id=cid)
            for line in seg_lines:
                resname = line[17:20].strip()
                # Check for 4-char resname (CHARMM-GUI style or GROMACS output).
                # Use 4-char if it's in known sets, OR if the 3-char truncation
                # isn't a standard amino acid (catches all CHARMM carb/lipid names).
                resname_4 = line[17:21].strip()
                if len(resname_4) == 4:
                    if (resname_4 in PDB_TO_LIPID or resname_4 in PDB_TO_CARB
                            or resname_4 in CERAMIDE_RTP
                            or resname_4 in _KNOWN_4CHAR_RESNAMES
                            or resname not in STANDARD_AA):
                        resname = resname_4
                resseq = int(line[22:26])
                icode = line[26].strip()
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                key = (resseq, icode)
                if not chain.residues or (chain.residues[-1].resseq, chain.residues[-1].icode) != key:
                    chain.residues.append(PDBResidue(
                        chain_id=cid, resname=resname,
                        resseq=resseq, icode=icode,
                    ))
                chain.residues[-1].atoms.append((atom_name, x, y, z))

            chains.append(chain)

    return chains


def _parse_ion_names(ions_itp_path):
    """Parse moleculetype names from ions.itp file."""
    names = set()
    with open(ions_itp_path) as f:
        in_moltype = False
        for line in f:
            stripped = line.strip()
            if stripped == '[ moleculetype ]':
                in_moltype = True
                continue
            if in_moltype and stripped and not stripped.startswith(';'):
                names.add(stripped.split()[0])
                in_moltype = False
            if stripped.startswith('[') and 'moleculetype' not in stripped:
                in_moltype = False
    return names


def _count_molecules(pdb_path, mol_names):
    """Count ion/buffer molecules in PDB by residue name.

    Returns list of (name, count) in order of first appearance,
    preserving the PDB residue ordering for [ molecules ].
    """
    from collections import OrderedDict
    counts = OrderedDict()
    seen_residues = set()  # (chain, resseq, resname) to avoid double-counting atoms
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            resname = line[17:20].strip()
            # Also check 4-char resname (cols 17-20)
            resname4 = line[17:21].strip()
            name = None
            if resname in mol_names:
                name = resname
            elif resname4 in mol_names:
                name = resname4
            if name is None:
                continue
            chain = line[21]
            resseq = int(line[22:26])
            key = (chain, resseq, name)
            if key not in seen_residues:
                seen_residues.add(key)
                counts[name] = counts.get(name, 0) + 1
    return list(counts.items())


def _gro_to_pdb(gro_path):
    """Convert a GROMACS .gro file to a temporary PDB via MDAnalysis."""
    import tempfile
    import warnings

    import MDAnalysis as mda
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(str(gro_path))
        tmp = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
        tmp.close()
        u.atoms.write(tmp.name)
    return Path(tmp.name)


def _split_chain_by_distance(chain, gap_cutoff=4.0):
    """Split a chain into sub-chains where consecutive residues are > gap_cutoff apart.

    Uses nearest-atom distance between consecutive residues (same approach as
    split_chains.py criterion 3). Molecules that are physically separate
    (e.g. ACET, ACEH in buffer) get split into individual chains.
    """
    import numpy as np

    if len(chain.residues) <= 1:
        return [chain]

    # Build coordinate arrays per residue
    res_coords = []
    for r in chain.residues:
        coords = np.array([(x, y, z) for _, x, y, z in r.atoms])
        res_coords.append(coords)

    # Find breaks: where nearest-atom distance > gap_cutoff
    breaks = [0]
    for i in range(1, len(chain.residues)):
        prev = res_coords[i - 1]
        cur = res_coords[i]
        # Nearest-atom distance
        diff = prev[:, None, :] - cur[None, :, :]
        min_dist = np.sqrt((diff ** 2).sum(axis=2)).min()
        if min_dist > gap_cutoff:
            breaks.append(i)

    if len(breaks) == 1:
        return [chain]  # no splits needed

    # Split into sub-chains
    result = []
    for bi in range(len(breaks)):
        start = breaks[bi]
        end = breaks[bi + 1] if bi + 1 < len(breaks) else len(chain.residues)
        sub = PDBChain(chain_id=chain.chain_id)
        sub.residues = chain.residues[start:end]
        result.append(sub)

    return result


def _count_water(pdb_path):
    """Count water molecules (SOL/HOH/WAT/TIP3/TIP4/TIP5/SPC/SPCE) in PDB.

    Uses atom count / atoms-per-molecule (per resname, via
    `_WATER_ATOMS_PER_MOL`) to handle resseq overflow in large systems
    where PDB wraps at 9999 — per-residue counting would undercount once
    multiple distinct residues collide onto the same wrapped resSeq. A
    single hardcoded "//3" divisor is wrong for TIP4/TIP5 (4/5-site
    models) and HOH (often deposited O-only, 1 atom); atoms are tallied
    separately per matched resname so each uses its own divisor.
    """
    water_atom_counts: dict[str, int] = {}
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            resname = line[17:20].strip()
            resname4 = line[17:21].strip()
            matched = resname if resname in _WATER_RESNAMES else (
                resname4 if resname4 in _WATER_RESNAMES else None)
            if matched:
                water_atom_counts[matched] = water_atom_counts.get(matched, 0) + 1
    return sum(count // _WATER_ATOMS_PER_MOL.get(name, 3)
               for name, count in water_atom_counts.items())


def _add_protonation_hydrogens(protein_chains, pdb_path, ff_type, verbose=False):
    """Add missing protonation H atoms using OpenMM Modeller with variants.

    Same approach as protonate.py: load PDB → strip H → strip non-protein →
    build variants list → PDBFixer fix missing atoms → Modeller.addHydrogens
    with CHARMM FF and variants → extract only needed protonation H coords.

    OpenMM variant names (ASH, GLH, HIP, HID, HIE) work with both
    charmm36.xml and amber14-all.xml force fields.
    """
    import os
    import tempfile

    # Map protonated names to OpenMM variant names
    _PROT_TO_VARIANT = {
        'ASPP': 'ASH', 'ASH': 'ASH', 'ASPH': 'ASH',
        'GLUP': 'GLH', 'GLH': 'GLH', 'GLUH': 'GLH',
        'HSP': 'HIP', 'HIP': 'HIP', 'HISH': 'HIP',
        'HSD': 'HID', 'HID': 'HID', 'HISD': 'HID',
        'HSE': 'HIE', 'HIE': 'HIE', 'HISE': 'HIE',
    }
    # Which H atoms we want for each protonated form
    _PROT_H_ATOMS = {
        'ASPP': {'HD2'}, 'ASH': {'HD2'}, 'ASPH': {'HD2'},
        'GLUP': {'HE2'}, 'GLH': {'HE2'}, 'GLUH': {'HE2'},
        'HSP': {'HD1', 'HE2'}, 'HIP': {'HD1', 'HE2'}, 'HISH': {'HD1', 'HE2'},
        'HSD': {'HD1'}, 'HID': {'HD1'}, 'HISD': {'HD1'},
        'HSE': {'HE2'}, 'HIE': {'HE2'}, 'HISE': {'HE2'},
    }

    # Collect residues that need protonation H
    need_h = {}  # (chain_id, resseq) -> (chain_ref, res_ref, prot_name, missing_h)
    for chain in protein_chains:
        for res in chain.residues:
            prot_name = res.resname.upper()
            if prot_name in _PROT_H_ATOMS:
                existing = {a[0] for a in res.atoms}
                missing = _PROT_H_ATOMS[prot_name] - existing
                if missing:
                    need_h[(chain.chain_id, res.resseq)] = (
                        chain, res, prot_name, missing)

    if not need_h:
        return

    if verbose:
        for (cid, rseq), (ch, res, pn, mh) in need_h.items():
            print(f"  Need H for {pn} {cid}:{rseq}: "
                  f"{', '.join(sorted(mh))}")

    from openmm import unit
    from openmm.app import ForceField, Modeller, PDBFile
    from pdbfixer import PDBFixer

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    from dvbfixer.protonate import _strip_hydrogens

    # Build variant lookup: (chain_id, resseq) -> OpenMM variant name
    variant_lookup = {}
    for (cid, rseq), (ch, res, pn, mh) in need_h.items():
        variant = _PROT_TO_VARIANT.get(pn)
        if variant:
            variant_lookup[(cid, rseq)] = variant

    # Load PDB with OpenMM (same approach as protonate.py)
    pdb = PDBFile(str(pdb_path))

    # Strip existing hydrogens
    topology, positions = _strip_hydrogens(pdb.topology, pdb.positions)

    # Strip non-protein residues (glycans, ligands, etc.) AND GLYCAM
    # glycosylated residues (NLN/OLS/OLT) — they don't need protonation H,
    # and their templates live in GLYCAM_06j-1.xml not amber14-all.xml.
    # Their coords are preserved in protein_chains, so this stripping is
    # only for OpenMM's hydrogen-addition pass.
    _GLYCAM_PROT = {'NLN', 'OLS', 'OLT'}
    known = (PROTEIN_RESIDUES | SOLVENT_IONS) - _GLYCAM_PROT
    to_delete = [res for res in topology.residues() if res.name not in known]
    if to_delete:
        modeller = Modeller(topology, positions)
        modeller.delete(to_delete)
        topology, positions = modeller.topology, modeller.positions

    # Build variants list. At N/C-terminals, drop ASH/GLH variants if using
    # AMBER FF since AMBER14 has no NASH/NGLH/CASH/CGLH templates (no RESP
    # charges were ever computed for terminal protonated ASP/GLU). HID/HIE/HIP
    # have terminal templates (NHIE/CHIE etc.) so they work fine.
    _AMBER_NO_TERMINAL = {'ASH', 'GLH'}
    _VARIANT_TO_STD = {'ASH': 'ASP', 'GLH': 'GLU'}

    def _build_variants(topo):
        terminals = set()
        for chain in topo.chains():
            res_list = list(chain.residues())
            if res_list:
                terminals.add(res_list[0].index)
                terminals.add(res_list[-1].index)
        vlist = []
        skipped_terminals = []
        for res in topo.residues():
            cid = res.chain.id
            try:
                rseq = int(res.id)
            except ValueError:
                vlist.append(None)
                continue
            key = (cid, rseq)
            var = variant_lookup.get(key)
            if (var in _AMBER_NO_TERMINAL and ff_type != 'charmm'
                    and res.index in terminals):
                vlist.append(None)
                skipped_terminals.append((var, cid, rseq))
                # Drop from need_h so HD2/HE2 isn't expected later
                need_h.pop(key, None)
                # Revert protein_chains rename (ASH→ASP, GLH→GLU) so the
                # output topology uses the standard terminal RTP entry
                std = _VARIANT_TO_STD[var]
                for pc in protein_chains:
                    if pc.chain_id == cid:
                        for r in pc.residues:
                            if r.resseq == rseq:
                                r.resname = std
            else:
                from dvbfixer.ffutils.variants import openmm_hydrogen_variant
                vlist.append(openmm_hydrogen_variant(var))
        if skipped_terminals:
            import warnings
            for var, cid, rseq in skipped_terminals:
                std = _VARIANT_TO_STD[var]
                warnings.warn(
                    f"Terminal {var} {cid}:{rseq} → {std}: AMBER14 has no "
                    f"terminal protonated template (NASH/NGLH/CASH/CGLH). "
                    f"Using standard {std} (no HD2/HE2 added).", stacklevel=2
                )
        return vlist

    variants = _build_variants(topology)

    # Use PDBFixer to fix any missing heavy atoms
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
                mode='w', suffix='.pdb', delete=False) as f:
            PDBFile.writeFile(topology, positions, f, keepIds=True)
            tmp_path = f.name

        fixer = PDBFixer(filename=tmp_path)
        fixer.findMissingResidues()
        fixer.missingResidues = {}  # Only fix atoms, not residues
        fixer.findMissingAtoms()
        from dvbfixer.ffutils.geometry import rebuild_missing_atoms_with_retry
        rebuild_missing_atoms_with_retry(fixer, verbose=verbose, log_prefix="[top] ")

        modeller = Modeller(fixer.topology, fixer.positions)

        # Rebuild variants for potentially reordered topology
        variants = _build_variants(modeller.topology)

        # Add H with proper geometry — use CHARMM or AMBER FF
        if ff_type == 'charmm':
            ff = ForceField('charmm36.xml', 'charmm36/water.xml')
        else:
            ff = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller.addHydrogens(ff, variants=variants)
        # Post-addHydrogens guards: PDBFixer's addMissingAtoms just
        # above can rebuild sidechains with D-Cα; reflect them before
        # extracting the protonation H set.
        from dvbfixer.ffutils.geometry import (
            fix_ca_chirality,
            repair_misplaced_hydrogens,
        )
        repair_misplaced_hydrogens(modeller.topology, modeller.positions,
                                    verbose=verbose)
        fix_ca_chirality(modeller.topology, modeller.positions,
                          verbose=verbose)

        # Extract only the specific protonation H atoms we need
        added = 0
        for atom in modeller.topology.atoms():
            res = atom.residue
            cid = res.chain.id
            try:
                rseq = int(res.id)
            except ValueError:
                continue
            key = (cid, rseq)
            if key not in need_h:
                continue
            chain_ref, res_ref, prot_name, missing_h = need_h[key]
            if atom.name in missing_h:
                pos = modeller.positions[atom.index]
                xyz = pos.value_in_unit(unit.angstrom)
                res_ref.atoms.append((atom.name, xyz[0], xyz[1], xyz[2]))
                added += 1
                if verbose:
                    print(f"    Placed {atom.name} at {res_ref.resname} "
                          f"{cid}:{rseq} "
                          f"({xyz[0]:.3f}, {xyz[1]:.3f}, {xyz[2]:.3f})")

    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.unlink(tmp_path)

    if verbose:
        print(f"  Added {added} protonation H atom(s) via OpenMM Modeller")


def _extract_molecule_lines(pdb_path, mol_names):
    """Extract ATOM/HETATM lines for ion/buffer molecules from PDB."""
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            resname = line[17:20].strip()
            resname4 = line[17:21].strip()
            if resname in mol_names or resname4 in mol_names:
                lines.append(line)
    return lines


def read_ssbonds(path):
    """Read SS bonds from SSBOND records or auto-detect from CYS SG-SG distances.

    Returns [(chain1, resseq1, chain2, resseq2)].
    """
    ssbonds = []

    # Try SSBOND records first
    with open(path) as f:
        for line in f:
            if line.startswith('SSBOND'):
                ch1 = line[15]
                res1 = int(line[17:21])
                ch2 = line[29]
                res2 = int(line[31:35])
                ssbonds.append((ch1, res1, ch2, res2))

    if ssbonds:
        return ssbonds

    # Auto-detect from CYS SG atom distances (< 2.5 A)
    import math
    sg_atoms = []  # [(chain, resseq, x, y, z)]
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            resname = line[17:20].strip()
            atom_name = line[12:16].strip()
            if resname == 'CYS' and atom_name == 'SG':
                chain = line[21]
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                sg_atoms.append((chain, resseq, x, y, z))

    # Find SG pairs within 2.5 A
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            ch1, r1, x1, y1, z1 = sg_atoms[i]
            ch2, r2, x2, y2, z2 = sg_atoms[j]
            dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            if dist < 2.5:
                ssbonds.append((ch1, r1, ch2, r2))

    return ssbonds


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)

    if not input_path.exists():
        print(f"Error: {input_path} not found", file=sys.stderr)
        sys.exit(1)

    # Validate (--ff, --water) compatibility and resolve --ion-set auto.
    # CHARMM ions (SOD/CLA/POT/CAL/MGA) are fitted to CHARMM-TIP3P; mixing them
    # with OPC/TIP4P/TIP4P-Ew is not supported by CHARMM developers.
    if args.ff == 'charmm':
        if args.water in {'tip4p', 'tip4pew', 'opc'}:
            print(f"ERROR: --water {args.water} is not parametrized for CHARMM36. "
                  f"Use --ff amber to combine those waters with matched ions, or "
                  f"pick --water tip3p|spc|spce for CHARMM.",
                  file=sys.stderr)
            sys.exit(1)
        if args.ion_set != 'auto':
            print("INFO: --ion-set is ignored with --ff charmm "
                  "(CHARMM ions come from the bundled ions.itp).",
                  file=sys.stderr)
    else:  # AMBER
        if args.ion_set == 'auto':
            args.ion_set = _WATER_DEFAULT_ION_SET[args.water]
        # Warn about water-model substitutions for non-JC waters
        if args.water in _WATER_ION_ALIAS and args.ion_set.startswith('jc-'):
            alias = _WATER_ION_ALIAS[args.water]
            print(f"WARNING: plain {args.water.upper()} was not parametrized by "
                  f"Joung-Cheatham; using {alias.upper()} ions ({args.ion_set}).",
                  file=sys.stderr)

    # Convert GRO to temp PDB if needed
    tmp_pdb = None
    orig_input_path = input_path  # preserve for output path defaults
    if input_path.suffix.lower() == '.gro':
        print("Converting GRO to PDB via MDAnalysis...")
        tmp_pdb = _gro_to_pdb(input_path)
        input_path = tmp_pdb

    # Auto-infer CONECT records so SS detection, glycosidic-bond detection,
    # and glycosylation-site detection work on inputs without CONECT.
    if not args.no_infer_conect:
        from dvbfixer.pdbutils import _materialise_inferred_pdb
        input_path = Path(_materialise_inferred_pdb(
            input_path, verbose=args.verbose))

    # ACPYPE mode: OpenMM + ParmEd + ACPYPE pipeline (delegated to
    # top.acpype so the RTP-based `--ff amber|charmm` path stays clean).
    if args.acpype:
        from dvbfixer.top.acpype import run_acpype_mode
        run_acpype_mode(input_path, args)
        return

    # Determine FF directory
    if args.ff_dir:
        ff_dir = Path(args.ff_dir)
    else:
        ff_name = FF_CHOICES[args.ff]
        ff_dir = FF_DIR / ff_name
        if not ff_dir.exists():
            print(f"Error: Force field directory not found: {ff_dir}", file=sys.stderr)
            print("Use --ff-dir to specify the path", file=sys.stderr)
            sys.exit(1)

    ff_name = ff_dir.name
    print(f"Using force field: {ff_name}")

    # Read PDB
    chains = read_pdb_chains(input_path)
    if not chains:
        print("Error: No chains found in PDB", file=sys.stderr)
        sys.exit(1)

    # Build topology builder first (need its residue dict for chain filtering)
    builder = TopologyBuilder(ff_dir, args.ff, args.verbose,
                              keep_all_hydrogens=args.keep_all_hydrogens)

    # Strip water and ion residues (handled separately via counting)
    ion_names_for_filter = set()
    ions_path_check = ff_dir / 'ions.itp'
    if ions_path_check.exists():
        ion_names_for_filter = _parse_ion_names(ions_path_check)
    skip_resnames = _WATER_RESNAMES | ion_names_for_filter

    for chain in chains:
        chain.residues = [r for r in chain.residues
                          if r.resname not in skip_resnames]

    # Remove empty chains (were all water/ions)
    chains = [c for c in chains if c.residues]

    # Split chains by nearest-atom distance (detect separate molecules in
    # GROMACS PDB output where multiple molecules share a chain ID)
    split_chains_list = []
    for chain in chains:
        split_chains_list.extend(_split_chain_by_distance(chain, gap_cutoff=4.0))
    chains = split_chains_list

    if chains:
        resnames = set()
        for c in chains:
            for r in c.residues:
                resnames.add(r.resname)
        print(f"Found {len(chains)} chain(s), residue types: {', '.join(sorted(resnames))}")

    protein_chains = []
    other_chains = []
    small_mol_counts = {}  # resname -> count (for single-residue small molecules)
    sugar_names = set(PDB_TO_CARB.keys())
    charmm_sugar_names_filter = set(PDB_TO_CARB.values())
    for chain in chains:
        has_protein = any(
            r.resname in STANDARD_AA or r.resname in PDB_TO_GMX
            for r in chain.residues
        )
        if has_protein:
            # Keep only protein residues in protein chains
            chain.residues = [
                r for r in chain.residues
                if r.resname in STANDARD_AA or r.resname in PDB_TO_GMX
            ]
            protein_chains.append(chain)
        elif len(chain.residues) == 1 and chain.residues[0].resname in builder.residues:
            # Single-residue molecule (detected by distance split) — count it
            rn = chain.residues[0].resname
            small_mol_counts[rn] = small_mol_counts.get(rn, 0) + 1
        else:
            # Check if chain has any FF-recognized residues
            has_known = any(
                r.resname in builder.residues or r.resname in sugar_names
                or _is_ceramide(r.resname)
                for r in chain.residues
            )
            if has_known:
                other_chains.append(chain)
            else:
                unknown_resnames = sorted({r.resname for r in chain.residues})
                print(f"WARNING: chain {chain.chain_id} has no residue(s) "
                      f"recognized by the {args.ff} force field or GLYCAM/"
                      f"CHARMM sugar/ceramide tables — dropping it entirely "
                      f"from the topology (resname(s): "
                      f"{', '.join(unknown_resnames)}). No .itp will be "
                      f"generated for this chain; add it to PDB_TO_LIPID or "
                      f"parametrize it separately (e.g. GAFF2/antechamber) "
                      f"if it needs a topology.", file=sys.stderr)

    if not protein_chains and not other_chains and not small_mol_counts:
        print("Error: No recognized chains found", file=sys.stderr)
        sys.exit(1)

    # Detect SS bonds
    ss_bonds = read_ssbonds(input_path)

    # Parse explicit --ss flags
    for ss_spec in args.ss:
        parts = ss_spec.split(':')
        if len(parts) == 4:
            ss_bonds.append((parts[0], int(parts[1]), parts[2], int(parts[3])))

    # Build per-chain SS residue sets and intra-chain SS pairs
    chain_ss = defaultdict(set)
    intrachain_ss = defaultdict(list)  # chain_id -> [(resseq1, resseq2)]
    for ch1, res1, ch2, res2 in ss_bonds:
        chain_ss[ch1].add(res1)
        chain_ss[ch2].add(res2)
        if ch1 == ch2:
            intrachain_ss[ch1].append((res1, res2))

    if ss_bonds and args.verbose:
        print(f"Disulfide bonds: {len(ss_bonds)}")
        for ch1, res1, ch2, res2 in ss_bonds:
            print(f"  {ch1}:{res1} - {ch2}:{res2}")

    # Apply protonation overrides (--protonate and --his)
    # --protonate without args: all ASP->ASPP, GLU->GLUP, HIS->HSP
    # --protonate with args: specific CHAIN:NUM[:STATE] overrides
    # Default protonated forms per FF
    if args.ff == 'charmm':
        _PROT_DEFAULTS = {'ASP': 'ASPP', 'GLU': 'GLUP', 'HIS': 'HSP',
                          'HIE': 'HSP', 'HID': 'HSP', 'HSE': 'HSP', 'HSD': 'HSP'}
    else:
        _PROT_DEFAULTS = {'ASP': 'ASH', 'GLU': 'GLH', 'HIS': 'HIP',
                          'HIE': 'HIP', 'HID': 'HIP', 'HSE': 'HIP', 'HSD': 'HIP'}

    protonate_all = args.protonate == 'all'
    prot_overrides = {}
    if args.protonate and args.protonate != 'all':
        # Validate: must contain ':' (CHAIN:NUM format), not a filename
        if ':' not in args.protonate:
            print(f"ERROR: Invalid --protonate value '{args.protonate}'. "
                  f"Use --protonate (all) or --protonate CHAIN:NUM[:STATE],... "
                  f"Note: place --protonate after the input file.",
                  file=sys.stderr)
            sys.exit(1)
        for spec in args.protonate.split(','):
            parts = spec.split(':')
            if len(parts) == 3:
                prot_overrides[(parts[0], int(parts[1]))] = parts[2]
            elif len(parts) == 2:
                # CHAIN:NUM without STATE — use default protonated form
                prot_overrides[(parts[0], int(parts[1]))] = None  # resolve later

    his_overrides = {}
    for his_spec in args.his:
        parts = his_spec.split(':')
        if len(parts) == 3:
            his_overrides[(parts[0], int(parts[1]))] = parts[2]

    # Map each protonation STATE to the residue family it's valid for.
    # A user passing --protonate X:N:HIP must point at a HIS-family residue,
    # not a VAL — catch this before OpenMM does and emit a clear error.
    _PROT_PARENT = {
        # ASH family — must target ASP
        'ASH': 'ASP', 'ASPP': 'ASP', 'ASPH': 'ASP',
        # GLH family — must target GLU
        'GLH': 'GLU', 'GLUP': 'GLU', 'GLUH': 'GLU',
        # HIS variants — must target HIS
        'HIP': 'HIS', 'HSP': 'HIS', 'HISH': 'HIS',
        'HIE': 'HIS', 'HSE': 'HIS', 'HISE': 'HIS',
        'HID': 'HIS', 'HSD': 'HIS', 'HISD': 'HIS',
        # CYS variants — must target CYS
        'CYX': 'CYS', 'CYM': 'CYS',
        # LYS variant — must target LYS
        'LYN': 'LYS', 'LSN': 'LYS',
        # Vanilla standard names are also legal "targets" for themselves
        'ASP': 'ASP', 'GLU': 'GLU', 'HIS': 'HIS', 'CYS': 'CYS', 'LYS': 'LYS',
    }

    # Pre-validate every --protonate target: residue name vs requested state.
    # Build a (chain_id, resseq) -> resname lookup from the parsed PDB.
    seen_residues = {(r.chain_id, r.resseq): r.resname
                     for chain in protein_chains for r in chain.residues}
    bad_targets = []  # (cid, rseq, requested_state, actual_resname_or_missing)
    for (cid, rseq), state in prot_overrides.items():
        actual = seen_residues.get((cid, rseq))
        if actual is None:
            bad_targets.append((cid, rseq, state, None))
            continue
        if state is None:
            # No explicit STATE — only valid if residue is in _PROT_DEFAULTS keys
            if actual not in _PROT_DEFAULTS:
                bad_targets.append((cid, rseq, '(default)', actual))
            continue
        expected_parent = _PROT_PARENT.get(state.upper())
        if expected_parent is None:
            bad_targets.append((cid, rseq, state, actual))
            continue
        actual_parent = _PROT_PARENT.get(actual.upper(), actual.upper())
        if actual_parent != expected_parent:
            bad_targets.append((cid, rseq, state, actual))
    if bad_targets:
        # Pre-index residues by (chain, parent_family) for nearby-suggestion
        # output: when a target is wrong, show the closest residues in the same
        # chain that ARE valid for the requested state.
        chain_family_residues = defaultdict(list)  # (cid, parent) -> [(rseq, resname)]
        for chain in protein_chains:
            for r in chain.residues:
                parent = _PROT_PARENT.get(r.resname.upper(), r.resname.upper())
                chain_family_residues[(r.chain_id, parent)].append(
                    (r.resseq, r.resname))

        def _nearest(cid, rseq, parent, n=5):
            candidates = chain_family_residues.get((cid, parent), [])
            if not candidates:
                return []
            return sorted(candidates, key=lambda rr: abs(rr[0] - rseq))[:n]

        print("ERROR: --protonate targets that don't match the actual residue:",
              file=sys.stderr)
        for cid, rseq, state, actual in bad_targets:
            if actual is None:
                # No residue at all — show neighbouring resseqs that DO exist
                existing = sorted({r.resseq for chain in protein_chains
                                   for r in chain.residues
                                   if chain.chain_id == cid})
                if existing:
                    print(f"  {cid}:{rseq}:{state}  →  no residue at that "
                          f"chain/resnum. Chain {cid} has resseq "
                          f"{existing[0]}..{existing[-1]} "
                          f"({len(existing)} residues).",
                          file=sys.stderr)
                else:
                    print(f"  {cid}:{rseq}:{state}  →  chain {cid} not found "
                          f"in the input.",
                          file=sys.stderr)
            elif state not in _PROT_PARENT and state != '(default)':
                print(f"  {cid}:{rseq}:{state}  →  unknown protonation state "
                      f"(valid: ASH/ASPP, GLH/GLUP, HIE/HID/HIP/HSE/HSD/HSP, "
                      f"CYX/CYM, LYN/LSN). Residue at this position is {actual}.",
                      file=sys.stderr)
            else:
                expected = _PROT_PARENT.get(state.upper(), '?')
                nearby = _nearest(cid, rseq, expected)
                hint = ''
                if nearby:
                    nearby_str = ', '.join(
                        f"{cid}:{rs}({rn})" for rs, rn in nearby)
                    hint = f" Nearest {expected} in chain {cid}: {nearby_str}."
                else:
                    hint = f" Chain {cid} has no {expected} residues at all."
                print(f"  {cid}:{rseq}:{state}  →  residue at that position is "
                      f"{actual}, but {state} is only valid for {expected}.{hint}",
                      file=sys.stderr)
        print("Check that the chain IDs and residue numbers match your input "
              "PDB. Use `grep '^ATOM' input.pdb | awk '{print $5,$6,$4}' | "
              "sort -u` to list (chain, resnum, resname) triples.",
              file=sys.stderr)
        sys.exit(1)

    # Apply overrides or auto-detect from PDB atoms
    for chain in protein_chains:
        for res in chain.residues:
            key = (res.chain_id, res.resseq)
            if key in prot_overrides:
                state = prot_overrides[key]
                if state is None:
                    state = _PROT_DEFAULTS.get(res.resname)
                if state:
                    old_name = res.resname
                    res.resname = state
                    if args.verbose:
                        print(f"  {old_name} {res.chain_id}:{res.resseq} -> "
                              f"{res.resname} (--protonate)")
            elif protonate_all and res.resname in _PROT_DEFAULTS:
                old_name = res.resname
                res.resname = _PROT_DEFAULTS[old_name]
                if args.verbose:
                    print(f"  {old_name} {res.chain_id}:{res.resseq} -> "
                          f"{res.resname} (--protonate all)")
            elif key in his_overrides:
                res.resname = his_overrides[key]
            elif res.resname == 'HIS':
                # Auto-detect protonation from H atoms present
                atom_names = {a[0] for a in res.atoms}
                has_hd1 = 'HD1' in atom_names
                has_he2 = 'HE2' in atom_names
                if has_hd1 and has_he2:
                    res.resname = 'HIP'  # doubly protonated
                elif has_hd1:
                    res.resname = 'HID'  # delta protonated
                else:
                    res.resname = 'HIE'  # epsilon protonated (default)
                if args.verbose and res.resname != 'HIE':
                    print(f"  HIS {res.chain_id}:{res.resseq} -> {res.resname} "
                          f"(auto-detected from H atoms)")

    # Add missing protonation H atoms with proper geometry via PDBFixer
    if args.protonate or args.his:
        _add_protonation_hydrogens(protein_chains, args.input, args.ff,
                                   verbose=args.verbose)

    # Ignore hydrogens if requested
    if args.ignh:
        for chain in protein_chains:
            for res in chain.residues:
                res.atoms = [(n, x, y, z) for n, x, y, z in res.atoms
                             if not n.startswith('H') and n not in ('1H', '2H', '3H')]

    # Build topologies (builder already created above for chain filtering)
    chain_tops = []

    for chain in protein_chains:
        if args.verbose:
            print(f"\nBuilding topology for chain {chain.chain_id} "
                  f"({len(chain.residues)} residues)")

        ct = builder.build_chain(chain, chain_ss.get(chain.chain_id, set()),
                                 intrachain_ss.get(chain.chain_id, []))
        if ct is not None:
            chain_tops.append(ct)
            n_bonds = len(ct.bonds)
            n_angles = len(ct.angles)
            n_dihedrals = len(ct.dihedrals)
            print(f"Chain {chain.chain_id}: {len(ct.atoms)} atoms, {n_bonds} bonds, "
                  f"{n_angles} angles, {n_dihedrals} dihedrals")
        else:
            print(f"WARNING: Failed to build topology for chain {chain.chain_id}",
                  file=sys.stderr)

    # Detect glycan/glycolipid links early to know which residues are handled
    glycan_links = detect_glycan_links(chains, pdb_path=input_path)
    glycolipid_ceramide_resseqs = set()  # (chain_id, resseq) of ceramides in glycolipids
    glycolipid_sugar_resseqs = set()     # sugars in glycolipid trees

    if glycan_links:
        trees = build_glycan_trees(glycan_links, chains)
        for tree, link_atoms, prot_links, cer_links in trees:
            if cer_links:
                for cl in cer_links:
                    glycolipid_ceramide_resseqs.add((cl[0], cl[1]))
                for ch, rs in tree:
                    glycolipid_sugar_resseqs.add((ch, rs))

    charmm_sugar_names = set(PDB_TO_CARB.values())

    # Build non-protein chains (no terminal patches, no SS bonds)
    # Skip chains that are entirely sugars or glycolipids — handled later
    for chain in other_chains:
        all_sugar_or_lipid = all(
            r.resname in sugar_names or r.resname in charmm_sugar_names
            or (chain.chain_id, r.resseq) in glycolipid_ceramide_resseqs
            for r in chain.residues
        )
        if all_sugar_or_lipid:
            if args.verbose:
                print(f"\nSkipping chain {chain.chain_id} (sugar/glycolipid residues, "
                      f"handled by glycan/glycolipid detection)")
            continue

        if args.verbose:
            print(f"\nBuilding topology for non-protein chain {chain.chain_id} "
                  f"({len(chain.residues)} residues)")

        ct = builder.build_chain(chain)
        if ct is not None:
            # Rename from default Protein_ to Other_
            ct.name = ct.name.replace('Protein_', 'Other_')
            chain_tops.append(ct)
            print(f"Chain {chain.chain_id}: {len(ct.atoms)} atoms, "
                  f"{len(ct.bonds)} bonds (non-protein)")
        else:
            print(f"WARNING: Failed to build topology for chain {chain.chain_id}",
                  file=sys.stderr)

    if not chain_tops and not glycan_links and not small_mol_counts:
        print("Error: No topologies built", file=sys.stderr)
        sys.exit(1)

    # Detect and build glycan/glycolipid chains
    glycan_tops = []
    protein_sugar_links = []  # for intermolecular interactions
    if glycan_links:
        # Build residue lookup for glycolipid ceramide
        res_lookup = {}
        for chain in chains:
            for res in chain.residues:
                res_lookup[(chain.chain_id, res.resseq)] = res

        for tree, link_atoms, prot_links, cer_links in trees:
            if not tree:
                continue
            # Filter glycan_links relevant to this tree
            tree_set = set(tree)
            tree_links = [
                gl for gl in glycan_links
                if (gl[0], gl[1]) in tree_set or (gl[3], gl[4]) in tree_set
            ]
            # Only sugar-sugar links for building the glycan molecule
            sugar_sugar_links = [
                gl for gl in tree_links
                if (gl[0], gl[1]) in tree_set and (gl[3], gl[4]) in tree_set
            ]

            if cer_links:
                # Glycolipid: build ceramide + sugar tree as one moleculetype
                cl = cer_links[0]  # use first ceramide link
                cer_res = res_lookup.get((cl[0], cl[1]))
                if cer_res is not None:
                    ct = builder.build_glycolipid_chain(
                        cer_res, tree, link_atoms, chains,
                        sugar_sugar_links, cl)
                    if ct is not None:
                        glycan_tops.append(ct)
                        print(f"Glycolipid {ct.name}: {len(ct.atoms)} atoms, "
                              f"{len(ct.bonds)} bonds")
                    else:
                        print("WARNING: Failed to build glycolipid topology",
                              file=sys.stderr)
            else:
                # Pure glycan tree
                ct = builder.build_glycan_chain(tree, link_atoms, chains,
                                                sugar_sugar_links)
                if ct is not None:
                    glycan_tops.append(ct)
                    print(f"Glycan {ct.name}: {len(ct.atoms)} atoms, "
                          f"{len(ct.bonds)} bonds")

            # Collect protein-sugar links relevant to this tree
            for pl in prot_links:
                if (pl[3], pl[4]) in tree_set and pl not in protein_sugar_links:
                    protein_sugar_links.append(pl)

    chain_tops.extend(glycan_tops)

    # Build small molecule topologies (single-residue CGenFF molecules)
    small_mol_tops = []  # list of (ChainTopology, count)
    if small_mol_counts:
        for resname, count in small_mol_counts.items():
            # Build a single-residue chain for this molecule type
            dummy_chain = PDBChain(chain_id='X')
            # Find a representative residue from the original parsed chains
            rep_res = None
            for ch in read_pdb_chains(input_path):
                for r in ch.residues:
                    if r.resname == resname:
                        rep_res = r
                        break
                if rep_res is not None:
                    break
            if rep_res is None:
                print(f"WARNING: Cannot find {resname} residue for topology",
                      file=sys.stderr)
                continue
            rep_res.resseq = 1  # reset to 1 (GRO may have global numbering)
            dummy_chain.residues = [rep_res]
            ct = builder.build_chain(dummy_chain)
            if ct is not None:
                ct.name = resname
                small_mol_tops.append((ct, count))
                print(f"Small molecule {resname}: {len(ct.atoms)} atoms, "
                      f"{count} copies")

    if not chain_tops and not small_mol_tops:
        print("Error: No topologies built", file=sys.stderr)
        sys.exit(1)

    # Merge chains if requested
    if args.merge and len(chain_tops) > 1:
        merged = _merge_chains(chain_tops)
        chain_tops = [merged]
        print(f"Merged into single moleculetype: {merged.name}")

    # Output paths
    if args.output:
        top_path = Path(args.output)
    else:
        top_path = orig_input_path.parent / 'topol.top'

    out_dir = top_path.parent

    # Collect inter-chain SS bonds (both chains different)
    interchain_ss = []
    for ch1, res1, ch2, res2 in ss_bonds:
        if ch1 != ch2:
            interchain_ss.append((ch1, res1, ch2, res2))

    # Build per-chain bonded_types list for write_top
    bonded_types_list = []
    for ct in chain_tops:
        bt = (builder.carb_bonded_types
              if (ct.name.startswith('Glycan_') or ct.name.startswith('Glycolipid_'))
              and builder.carb_bonded_types
              else builder.bonded_types)
        bonded_types_list.append(bt)

    # Write position restraint files (still separate, used with #ifdef POSRES)
    for ct in chain_tops:
        posre_path = out_dir / f"posre_{ct.name}.itp"
        write_posre(ct, posre_path)
        print(f"Wrote {posre_path}")

    # Position restraints for small molecules too
    for ct, count in small_mol_tops:
        posre_path = out_dir / f"posre_{ct.name}.itp"
        write_posre(ct, posre_path)
        print(f"Wrote {posre_path}")

    # Write small molecule .itp files
    for ct, count in small_mol_tops:
        itp_path = out_dir / f"{ct.name}.itp"
        bt = builder.bonded_types
        if builder.carb_bonded_types:
            bt = builder.carb_bonded_types
        with open(itp_path, 'w') as f:
            f.write(f"; Moleculetype: {ct.name}\n")
            f.write("; Generated by dvbfixer top\n\n")
            _write_moleculetype(f, ct, bt)
        print(f"Wrote {itp_path}")

    # Write inter-chain bond file if needed (SS bonds + protein-glycan bonds)
    has_interchain_bonds = interchain_ss or protein_sugar_links
    if has_interchain_bonds:
        ss_path = out_dir / "interchain_ss.itp"
        _write_interchain_ss(interchain_ss, chain_tops, protein_chains, ss_path,
                             args.ff, protein_sugar_links)
        n_ss = len(interchain_ss)
        n_pg = len(protein_sugar_links)
        parts = []
        if n_ss:
            parts.append(f"{n_ss} inter-chain SS bond(s)")
        if n_pg:
            parts.append(f"{n_pg} protein-glycan bond(s)")
        print(f"Wrote {ss_path} ({', '.join(parts)})")
        print("WARNING: interchain_ss.itp must stay at the end of topol.top, after [ molecules ].")
        print("         After gmx solvate/genion, move the #include line below SOL/ion entries.")

    # Detect ions/BUF, small molecules, and water in PDB (preserving PDB order)
    ions_path = ff_dir / 'ions.itp'
    ion_names = set()
    if ions_path.exists():
        ion_names = _parse_ion_names(ions_path)
    small_mol_names_set = {ct.name for ct, _ in small_mol_tops}
    # Count all extra molecules preserving PDB order
    countable_names = ion_names | small_mol_names_set | _WATER_RESNAMES
    extra_molecules = _count_molecules(input_path, countable_names)
    # Fix water count: _count_molecules uses (chain, resseq) dedup which breaks
    # for large systems where PDB wraps resseq at 9999. Use atom count / 3 instead.
    water_count = _count_water(input_path)
    extra_molecules = [
        ('SOL', water_count) if name in _WATER_RESNAMES else (name, count)
        for name, count in extra_molecules
    ]
    if extra_molecules:
        for mol_name, mol_count in extra_molecules:
            print(f"Found {mol_count} {mol_name} molecule(s) in PDB")

    # Write TOP file with modular .itp includes
    small_mol_names = [ct.name for ct, _ in small_mol_tops]
    system_name = orig_input_path.stem
    # ion_set is None for CHARMM (use bundled ions.itp) and the resolved set name
    # for AMBER (emit water-matched ion atom types + moleculetypes).
    write_top_ion_set = None if args.ff == 'charmm' else args.ion_set
    write_top(chain_tops, top_path, ff_dir, args.ff, bonded_types_list,
              args.water, system_name,
              has_interchain_ss=has_interchain_bonds,
              extra_molecules=extra_molecules,
              small_mol_itps=small_mol_names,
              ion_set=write_top_ion_set)
    print(f"Wrote {out_dir / 'ffparams.itp'}")
    for ct in chain_tops:
        print(f"Wrote {out_dir / ct.name}.itp")
    print(f"Wrote {out_dir / 'water.itp'}")
    print(f"Wrote {out_dir / 'ions.itp'}")
    print(f"Wrote {top_path}")

    # Write output PDB with topology-matched atom names
    if args.pdb:
        pdb_path = Path(args.pdb)
    else:
        pdb_path = out_dir / 'conf.pdb'
    # Extract CRYST1 (box vectors) from input PDB
    cryst1_line = None
    with open(input_path) as f:
        for line in f:
            if line.startswith('CRYST1'):
                cryst1_line = line
                break

    # Collect ion/BUF/water PDB lines for output
    extra_pdb_lines = []
    if extra_molecules:
        extra_mol_names = ion_names | _WATER_RESNAMES | {ct.name for ct, _ in small_mol_tops}
        extra_pdb_lines = _extract_molecule_lines(input_path, extra_mol_names)
    write_pdb(chain_tops, pdb_path, extra_pdb_lines=extra_pdb_lines,
              cryst1=cryst1_line)
    print(f"Wrote {pdb_path}")

    # Clean up temp PDB from GRO conversion
    if tmp_pdb is not None:
        tmp_pdb.unlink(missing_ok=True)


def _write_interchain_ss(ss_bonds, chain_tops, protein_chains, path, ff_type,
                         protein_sugar_links=None):
    """Write inter-chain bond topology (SS bonds + protein-glycan bonds).

    Uses [ intermolecular_interactions ] which must be #included in .top
    after [ molecules ].
    """
    # Compute global atom offsets per chain topology (in [ molecules ] order)
    chain_offsets = {}  # ct.name -> offset
    offset = 0
    for ct in chain_tops:
        chain_offsets[ct.name] = offset
        offset += len(ct.atoms)

    def find_atom(chain_id, resseq, atomname):
        """Find global atom index by chain_id, resseq, atomname."""
        for ct in chain_tops:
            ct_offset = chain_offsets[ct.name]
            for atom in ct.atoms:
                if (atom.chain_id == chain_id and
                    atom.orig_resseq == resseq and
                    atom.atomname == atomname):
                    return ct_offset + atom.index
        return None

    with open(path, 'w') as f:
        f.write("; Inter-chain bonds (SS + protein-glycan)\n")
        f.write("; Generated by dvbfixer top\n\n")
        f.write("[ intermolecular_interactions ]\n\n")
        f.write("[ bonds ]\n")
        f.write(";  ai    aj  funct   r0 (nm)   k (kJ/mol/nm^2)\n")

        # SS bonds
        for ch1, res1, ch2, res2 in ss_bonds:
            idx1 = find_atom(ch1, res1, 'SG')
            idx2 = find_atom(ch2, res2, 'SG')
            if idx1 is not None and idx2 is not None:
                f.write(f"{idx1:5d} {idx2:5d}     6    0.204   250000\n")
                f.write(f"; {ch1}:{res1}:SG - {ch2}:{res2}:SG\n")

        # Protein-glycan bonds (ASN ND2 - NAG C1)
        if protein_sugar_links:
            for prot_ch, prot_rs, prot_atom, sug_ch, sug_rs in protein_sugar_links:
                idx1 = find_atom(prot_ch, prot_rs, prot_atom)
                idx2 = find_atom(sug_ch, sug_rs, 'C1')
                if idx1 is not None and idx2 is not None:
                    f.write(f"{idx1:5d} {idx2:5d}     6    0.1430   250000\n")
                    f.write(f"; {prot_ch}:{prot_rs}:{prot_atom} - "
                            f"{sug_ch}:{sug_rs}:C1\n")


def _merge_chains(chain_tops):
    """Merge multiple chain topologies into one."""
    merged = ChainTopology(
        name="Protein",
        nrexcl=chain_tops[0].nrexcl,
    )

    offset = 0
    resnr_offset = 0
    for ct in chain_tops:
        for atom in ct.atoms:
            new_atom = AtomEntry(
                index=atom.index + offset,
                atom_type=atom.atom_type,
                resnr=atom.resnr + resnr_offset,
                resname=atom.resname,
                atomname=atom.atomname,
                cgnr=atom.cgnr,
                charge=atom.charge,
                mass=atom.mass,
                x=atom.x,
                y=atom.y,
                z=atom.z,
                chain_id=atom.chain_id,
                orig_resseq=atom.orig_resseq,
                orig_resname=atom.orig_resname,
            )
            merged.atoms.append(new_atom)

        for i, j in ct.bonds:
            merged.bonds.append((i + offset, j + offset))
        for i, j in ct.pairs:
            merged.pairs.append((i + offset, j + offset))
        for i, j, k in ct.angles:
            merged.angles.append((i + offset, j + offset, k + offset))
        for dih in ct.dihedrals:
            if len(dih) == 5:
                i, j, k, l, t = dih
                merged.dihedrals.append((i + offset, j + offset, k + offset, l + offset, t))
            else:
                i, j, k, l = dih
                merged.dihedrals.append((i + offset, j + offset, k + offset, l + offset))
        for i, j, k, l in ct.impropers:
            merged.impropers.append((i + offset, j + offset, k + offset, l + offset))
        for cm in ct.cmap:
            merged.cmap.append(tuple(x + offset for x in cm))

        max_idx = max(a.index for a in ct.atoms) if ct.atoms else 0
        offset += max_idx
        max_resnr = max((a.resnr for a in ct.atoms), default=0)
        resnr_offset += max_resnr

    return merged


if __name__ == '__main__':
    main()
