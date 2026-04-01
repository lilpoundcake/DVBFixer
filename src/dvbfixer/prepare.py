"""Fix missing atoms and residues in a PDB structure using PDBFixer.

Adds missing residues, missing heavy atoms, and hydrogens. Writes a .dat file
recording which atoms were added, for use by 'dvbfixer minimize' to apply
selective restraints.
"""

import argparse
import json
import sys
from pathlib import Path

from openmm.app import Modeller, PDBFile
from pdbfixer import PDBFixer


DEFAULT_PH = 7.0

# Residues that form glycosidic bonds through ND2 (N-linked glycosylation)
GLYCOSYLATED_RESIDUES = {'ASN'}
SUGAR_RESNAMES = {'NAG', 'NDG', 'BMA', 'MAN', 'FUC', 'FUL', 'GAL', 'BGC', 'GLC', 'SIA'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer prepare",
        description="Fix missing atoms and residues in a PDB structure using PDBFixer. "
        "Writes a .dat file recording added atoms for selective restraints "
        "during minimization.",
    )
    p.add_argument("input", help="Input PDB file")
    p.add_argument("-o", "--output", help="Output PDB file (default: <input>_prepared.pdb)")
    p.add_argument("--dat", help="Restraint data file path (default: <output>.dat)")
    p.add_argument("--ph", type=float, default=DEFAULT_PH,
                   help=f"pH for adding hydrogens (default: {DEFAULT_PH})")
    p.add_argument("--keep-water", action="store_true",
                   help="Keep crystallographic waters")
    p.add_argument("--keep-heterogens", action="store_true",
                   help="Keep all heterogens (ligands, ions, etc.)")
    p.add_argument("--mutate", action="append", default=[],
                   metavar="CHAIN:RESNUM:NEW_AA",
                   help="Mutate a residue (e.g. A:39:ALA). Can be used multiple times.")
    p.add_argument("--rename", action="store_true",
                   help="Rename non-canonical residues (AMBER/CHARMM) to standard names before processing")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print detailed progress")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# .dat file: records what PDBFixer added
# ---------------------------------------------------------------------------

def build_dat(fixer, new_atom_indices):
    """Build a dict describing what PDBFixer added, suitable for JSON export."""
    atoms_list = list(fixer.topology.atoms())

    added_atoms = []
    added_residues_summary = {}
    for idx in sorted(new_atom_indices):
        atom = atoms_list[idx]
        res = atom.residue
        entry = {
            "chain": res.chain.id,
            "resid": res.id,
            "icode": res.insertionCode,
            "resname": res.name,
            "atom": atom.name,
            "element": atom.element.symbol,
        }
        added_atoms.append(entry)

        rkey = f"{res.chain.id}/{res.name}{res.id}"
        if rkey not in added_residues_summary:
            added_residues_summary[rkey] = {"heavy": 0, "hydrogen": 0}
        if atom.element.symbol == 'H':
            added_residues_summary[rkey]["hydrogen"] += 1
        else:
            added_residues_summary[rkey]["heavy"] += 1

    return {
        "description": "PDBFixer restraint data. Added atoms get weak/no restraints during minimization. "
                       "Edit 'added_atoms' list to change which atoms are treated as 'new'.",
        "total_added": len(added_atoms),
        "residue_summary": added_residues_summary,
        "added_atoms": added_atoms,
    }


def write_dat(dat, path):
    with open(path, 'w') as f:
        json.dump(dat, f, indent=2)
    print(f"Saved restraint data: {path}")


# ---------------------------------------------------------------------------
# Glycosylation-aware hydrogen fix
# ---------------------------------------------------------------------------

def find_glycosylated_atoms(input_path):
    """Parse CONECT records to find protein atoms bonded to sugars.

    Returns set of (chain_id, resid, atom_name) for protein atoms with
    glycosidic bonds (e.g. ASN ND2 bonded to NAG C1).
    """
    with open(input_path) as f:
        lines = f.readlines()

    serials = {}
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            serial = int(line[6:11])
            chain = line[21]
            resname = line[17:20].strip()
            resid = line[22:26].strip()
            atomname = line[12:16].strip()
            serials[serial] = (chain, resname, resid, atomname)

    glycosylated = set()
    for line in lines:
        if not line.startswith('CONECT'):
            continue
        parts = []
        s = line[6:]
        while len(s) >= 5:
            chunk = s[:5].strip()
            if chunk:
                parts.append(int(chunk))
            s = s[5:]
        if len(parts) < 2:
            continue
        src = parts[0]
        for dst in parts[1:]:
            if src not in serials or dst not in serials:
                continue
            s_info, d_info = serials[src], serials[dst]
            # Protein atom bonded to sugar
            if s_info[1] in GLYCOSYLATED_RESIDUES and d_info[1] in SUGAR_RESNAMES:
                glycosylated.add((s_info[0], s_info[2], s_info[3]))
            elif d_info[1] in GLYCOSYLATED_RESIDUES and s_info[1] in SUGAR_RESNAMES:
                glycosylated.add((d_info[0], d_info[2], d_info[3]))

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


# ---------------------------------------------------------------------------
# PDBFixer
# ---------------------------------------------------------------------------

# Map AMBER/CHARMM protonation variant names to standard PDB names
# (PDBFixer only accepts standard 3-letter codes for mutations)
_VARIANT_TO_STANDARD = {
    'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS',
    'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS',
    'ASH': 'ASP', 'ASPP': 'ASP',
    'GLH': 'GLU', 'GLUP': 'GLU',
    'CYX': 'CYS', 'CYM': 'CYS',
    'LYN': 'LYS',
}


def parse_mutations(mutate_args):
    """Parse --mutate arguments into PDBFixer mutation format.

    Input format: ['A:39:ALA', 'B:100:GLY', 'A:83:HIP']
    Handles AMBER protonation variants (HIP, ASH, GLH, etc.).
    Returns: (mutations_by_chain, variant_overrides)
      mutations_by_chain: dict chain_id -> [(resnum, standard_aa)]
      variant_overrides: dict (chain_id, resnum) -> variant_name
    """
    from collections import defaultdict
    mutations_by_chain = defaultdict(list)
    variant_overrides = {}  # (chain, resnum) -> user-requested variant name
    for spec in mutate_args:
        parts = spec.split(':')
        if len(parts) != 3:
            print(f"Error: invalid --mutate format '{spec}' (expected CHAIN:RESNUM:NEW_AA)",
                  file=sys.stderr)
            sys.exit(1)
        chain, resnum, new_aa = parts
        new_aa = new_aa.upper()

        # If it's a protonation variant, record it and use standard name for PDBFixer
        standard_aa = _VARIANT_TO_STANDARD.get(new_aa, new_aa)
        if new_aa != standard_aa:
            variant_overrides[(chain, resnum)] = new_aa

        mutations_by_chain[chain].append((resnum, standard_aa))
    return mutations_by_chain, variant_overrides


def apply_mutations(fixer, mutations_by_chain, verbose=False):
    """Apply point mutations using PDBFixer."""
    res_lookup = {}
    for res in fixer.topology.residues():
        res_lookup[(res.chain.id, res.id)] = res.name

    for chain_id, muts in mutations_by_chain.items():
        pdbfixer_muts = []
        for resnum, new_aa in muts:
            old_aa = res_lookup.get((chain_id, resnum))
            if old_aa is None:
                print(f"Warning: residue {chain_id}:{resnum} not found, skipping mutation")
                continue
            mut_str = f"{old_aa}-{resnum}-{new_aa}"
            pdbfixer_muts.append(mut_str)
            if verbose:
                print(f"  Mutation: {chain_id}:{old_aa}{resnum} -> {new_aa}")
        if pdbfixer_muts:
            fixer.applyMutations(pdbfixer_muts, chain_id)


def run_pdbfixer(input_path, ph, keep_water, keep_heterogens, verbose, mutations=None):
    """Run PDBFixer to add missing atoms/residues. Returns (fixer, new_atom_indices)."""
    # Detect glycosylated atoms before PDBFixer modifies anything
    glycosylated_atoms = find_glycosylated_atoms(input_path)
    if glycosylated_atoms and verbose:
        print(f"Detected {len(glycosylated_atoms)} glycosylated atom(s)")

    with open(input_path) as f:
        fixer = PDBFixer(pdbfile=f)

    # Apply mutations if requested
    variant_overrides = {}  # (chain, resnum) -> variant name (HIP, ASH, etc.)
    if mutations:
        mutations_by_chain, variant_overrides = parse_mutations(mutations)
        if mutations_by_chain:
            print(f"Applying {sum(len(v) for v in mutations_by_chain.values())} mutation(s)...")
            apply_mutations(fixer, mutations_by_chain, verbose)
            if variant_overrides and verbose:
                for (ch, rn), var in variant_overrides.items():
                    print(f"  Protonation override: {ch}:{rn} → {var}")

    # Capture AMBER/CHARMM protonation variant names from input PDB before
    # PDBFixer normalizes them. These are re-applied after replaceNonstandardResidues.
    for res in fixer.topology.residues():
        key = (res.chain.id, res.id)
        if key not in variant_overrides and res.name in _VARIANT_TO_STANDARD:
            variant_overrides[key] = res.name

    # Strip hydrogens and reload so findMissingAtoms gets clean template matching.
    # Wrong H (e.g. from mutations or external tools) confuse PDBFixer.
    # Hydrogens are re-added at the end by addMissingHydrogens anyway.
    modeller = Modeller(fixer.topology, fixer.positions)
    h_atoms = [a for a in modeller.topology.atoms() if a.element.symbol == 'H']
    if h_atoms:
        modeller.delete(h_atoms)
        if verbose:
            print(f"Stripped {len(h_atoms)} H for clean template matching (will re-add)")
        import tempfile as _tf
        with _tf.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as _tmp:
            PDBFile.writeFile(modeller.topology, modeller.positions, _tmp, keepIds=True)
            _tmp_path = _tmp.name
        try:
            with open(_tmp_path) as _f:
                fixer = PDBFixer(pdbfile=_f)
        finally:
            Path(_tmp_path).unlink()

    # Record original atoms before any modifications
    original_atoms = set()
    for atom in fixer.topology.atoms():
        original_atoms.add((atom.residue.chain.id, atom.residue.id,
                            atom.residue.insertionCode, atom.name))

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()

    if verbose:
        if fixer.missingResidues:
            print("Missing residues:")
            for (chain_idx, res_idx), resnames in sorted(fixer.missingResidues.items()):
                chain_id = list(fixer.topology.chains())[chain_idx].id
                print(f"  Chain {chain_id} position {res_idx}: {' '.join(resnames)}")
        if fixer.missingAtoms:
            print("Missing atoms:")
            for res, atoms in fixer.missingAtoms.items():
                print(f"  {res.chain.id}/{res.name}{res.id}: {', '.join(a.name for a in atoms)}")
        if fixer.missingTerminals:
            print("Missing terminals:")
            for res, atoms in fixer.missingTerminals.items():
                print(f"  {res.chain.id}/{res.name}{res.id}: {', '.join(a if isinstance(a, str) else a.name for a in atoms)}")
        if fixer.nonstandardResidues:
            print("Nonstandard residues:")
            for res, std in fixer.nonstandardResidues:
                print(f"  {res.chain.id}/{res.name}{res.id} -> {std}")

    n_missing_res = sum(len(v) for v in fixer.missingResidues.values())
    n_missing_atoms = sum(len(v) for v in fixer.missingAtoms.items())
    print(f"PDBFixer: {n_missing_res} missing residues, {n_missing_atoms} residues with missing atoms")

    if not keep_heterogens:
        fixer.removeHeterogens(keepWater=keep_water)
    fixer.replaceNonstandardResidues()
    fixer.addMissingAtoms()

    # Restore protonation variant names that PDBFixer's replaceNonstandardResidues
    # reverted to standard names (HIP→HIS, ASH→ASP, GLH→GLU, CYX→CYS, etc.).
    # Also rename any remaining HIS to explicit variants (HIE default) to prevent
    # OpenMM's unreliable auto-detection from crashing.
    for res in fixer.topology.residues():
        key = (res.chain.id, res.id)
        if key in variant_overrides:
            old = res.name
            res.name = variant_overrides[key]
            if verbose:
                print(f"  {old} {res.chain.id}:{res.id} → {res.name} (variant override)")
        elif res.name == 'HIS':
            # No explicit override — detect from H atoms or default to HIE
            atom_names = {a.name for a in res.atoms()}
            if 'HD1' in atom_names and 'HE2' in atom_names:
                res.name = 'HIP'
            elif 'HD1' in atom_names:
                res.name = 'HID'
            elif 'HE2' in atom_names:
                res.name = 'HIE'
            else:
                res.name = 'HIE'  # default
            if verbose:
                print(f"  HIS {res.chain.id}:{res.id} → {res.name}")

    fixer.addMissingHydrogens(ph)

    # Remove extra hydrogens from glycosylated residues
    if glycosylated_atoms:
        remove_extra_glycan_hydrogens(fixer, glycosylated_atoms, verbose)

    # Identify which atoms in the new topology are newly added
    new_atom_indices = set()
    added_residue_ids = set()
    for atom in fixer.topology.atoms():
        key = (atom.residue.chain.id, atom.residue.id,
               atom.residue.insertionCode, atom.name)
        if key not in original_atoms:
            new_atom_indices.add(atom.index)
            added_residue_ids.add((atom.residue.chain.id, atom.residue.id))

    if verbose:
        n_new_heavy = sum(1 for idx in new_atom_indices
                          for atom in [list(fixer.topology.atoms())[idx]]
                          if atom.element.symbol != 'H')
        print(f"Added {len(new_atom_indices)} atoms ({n_new_heavy} heavy), "
              f"across {len(added_residue_ids)} residues")

    return fixer, new_atom_indices


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_prepared")
    dat_path = Path(args.dat) if args.dat else output_path.with_suffix(".dat")

    if args.rename:
        from dvbfixer.rename import canonicalize_pdb
        import tempfile as _tf
        _tmp = Path(_tf.mktemp(suffix='.pdb'))
        n = canonicalize_pdb(input_path, _tmp, args.verbose)
        if n > 0:
            print(f"Canonicalized {n} non-canonical residue(s)")
            input_path = _tmp
        elif _tmp.exists():
            _tmp.unlink()

    print(f"=== PDBFixer: {input_path} ===")
    fixer, new_atom_indices = run_pdbfixer(
        input_path, args.ph, args.keep_water, args.keep_heterogens, args.verbose,
        mutations=args.mutate
    )

    with open(output_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    print(f"Saved prepared structure: {output_path}")

    dat = build_dat(fixer, new_atom_indices)

    # Merge with upstream .dat (e.g. from model step) if it exists
    upstream_dat = input_path.with_suffix('.dat')
    if upstream_dat.exists():
        with open(upstream_dat) as f:
            prev = json.load(f)
        # Carry forward upstream added_atoms that still exist in output
        # (keyed by chain+resid+icode+atom to avoid duplicates)
        existing_keys = {(a["chain"], a["resid"], a["icode"], a["atom"])
                         for a in dat["added_atoms"]}
        carried = 0
        for a in prev["added_atoms"]:
            key = (a["chain"], a["resid"], a["icode"], a["atom"])
            if key not in existing_keys:
                dat["added_atoms"].append(a)
                existing_keys.add(key)
                rkey = f"{a['chain']}/{a['resname']}{a['resid']}"
                if rkey not in dat["residue_summary"]:
                    dat["residue_summary"][rkey] = {"heavy": 0, "hydrogen": 0}
                if a.get("element") == 'H':
                    dat["residue_summary"][rkey]["hydrogen"] += 1
                else:
                    dat["residue_summary"][rkey]["heavy"] += 1
                carried += 1
        dat["total_added"] = len(dat["added_atoms"])
        if carried > 0:
            print(f"Merged {carried} atoms from upstream {upstream_dat.name}")

    write_dat(dat, dat_path)
