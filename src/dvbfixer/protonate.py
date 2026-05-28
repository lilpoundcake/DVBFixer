"""Set residue protonation states in a PDB file based on PROPKA3 pKa predictions.

Runs PROPKA3 to predict per-residue pKa values, then renames titratable residues
to their correct protonation state at the target pH. Designed for AMBER force
field naming conventions (HID/HIE/HIP, ASH, GLH, CYM, LYN).

Protonation logic at a given pH:
  - HIS: pKa > pH -> HIP (doubly protonated); pKa < pH -> HIE (default neutral)
          or HID (if Nd1 is the donor based on local H-bond network)
  - ASP: pKa > pH -> ASH (protonated); otherwise ASP (charged)
  - GLU: pKa > pH -> GLH (protonated); otherwise GLU (charged)
  - CYS: pKa < pH -> CYM (deprotonated thiolate); otherwise CYS
         CYS in disulfide bonds (pKa=99.99 from PROPKA) -> CYX
  - LYS: pKa < pH -> LYN (neutral); otherwise LYS (charged)
  - TYR: pKa < pH -> rename not needed (standard TYR handles both)
"""

import argparse
import io
import sys
from pathlib import Path

WATER_RESNAMES = {'HOH', 'WAT', 'TIP3', 'TIP', 'SOL', 'T3P', 'T4P', 'T5P'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer protonate",
        description="Predict pKa values with PROPKA3 and set correct protonation "
        "state names in a PDB file for a given pH. Uses AMBER residue naming "
        "(HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN)."
    )
    p.add_argument("input", help="Input PDB file")
    p.add_argument("-o", "--output", help="Output PDB file (default: <input>_prot.pdb)")
    p.add_argument(
        "--ph", type=float, default=7.0,
        help="Target pH for protonation assignment (default: 7.0)"
    )
    p.add_argument(
        "--his-default", choices=["HIE", "HID"], default="HIE",
        help="Default neutral HIS tautomer when pKa < pH (default: HIE = Ne2 protonated)"
    )
    p.add_argument(
        "--cys-disulfide-pka", type=float, default=90.0,
        help="PROPKA pKa threshold above which CYS is assumed to be in a disulfide "
             "bond and renamed to CYX (default: 90.0)"
    )
    p.add_argument(
        "--summary", action="store_true",
        help="Print pKa summary table for all titratable residues"
    )
    p.add_argument(
        "--no-hydrogens", action="store_true",
        help="Only rename residues, do not add/fix hydrogen atoms"
    )
    p.add_argument(
        "--ff", nargs='+', default=['amber19/protein.ff19SB.xml', 'amber19/tip3p.xml'],
        help="Force field XML files for hydrogen addition (default: amber19/protein.ff19SB.xml amber19/tip3p.xml)"
    )
    p.add_argument(
        "--keep-water", action="store_true",
        help="Keep water molecules (HOH, WAT, TIP3, SOL) in output (default: remove)"
    )
    p.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print only residues that get non-standard protonation"
    )
    return p.parse_args(argv)


# GLYCAM glycoprotein residues that must be renamed to their standard parents
# before PROPKA3 sees them. PROPKA3 only recognizes the 20 canonical amino acids.
_PROPKA_RENAME = {'NLN': 'ASN', 'OLS': 'SER', 'OLT': 'THR'}


def _sanitize_for_propka(input_path):
    """Build a temp-file PDB with GLYCAM names renamed to standard parents
    and all heterogens (sugars, ligands) stripped. PROPKA3 only knows
    standard amino acids and will silently drop or error on NLN/UYB/etc.

    Returns the temp file path (caller must delete).
    """
    import tempfile
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    keep_resnames = PROTEIN_RESIDUES | SOLVENT_IONS

    fd, tmp_path = tempfile.mkstemp(suffix='.pdb', prefix='propka_')
    with open(input_path) as inp, open(tmp_path, 'w') as out:
        for line in inp:
            if line.startswith(('ATOM', 'HETATM', 'ANISOU')):
                resname = line[17:20].strip()
                # Rename GLYCAM glycoprotein residues so PROPKA processes them
                # as their standard parent for pKa calculation.
                if resname in _PROPKA_RENAME:
                    parent = _PROPKA_RENAME[resname]
                    line = line[:17] + f"{parent:<3s}" + line[20:]
                    resname = parent
                # Drop heterogens (sugars, ligands) — PROPKA doesn't need them.
                if resname not in keep_resnames:
                    continue
            elif line.startswith('TER'):
                # Keep TER but rewrite resname if it's a GLYCAM protein.
                if len(line) > 20:
                    resname = line[17:20].strip()
                    if resname in _PROPKA_RENAME:
                        line = line[:17] + f"{_PROPKA_RENAME[resname]:<3s}" + line[20:]
            elif line.startswith('CONECT'):
                # CONECT records reference HETATMs we stripped; skip.
                continue
            out.write(line)
    import os
    os.close(fd)
    return tmp_path


def run_propka(input_path):
    """Run PROPKA3 on the input PDB and return the MolecularContainer.

    Pre-sanitizes input by renaming GLYCAM glycoprotein residues
    (NLN→ASN, OLS→SER, OLT→THR) and stripping heterogens so PROPKA3
    can process glycoprotein inputs.
    """
    from propka.run import single

    sanitized = _sanitize_for_propka(input_path)
    try:
        # Suppress PROPKA warnings to stderr
        old_stderr = sys.stderr
        sys.stderr = io.StringIO()
        try:
            mc = single(sanitized, write_pka=False)
        finally:
            sys.stderr = old_stderr
    finally:
        Path(sanitized).unlink(missing_ok=True)

    return mc


def get_pka_results(mc):
    """Extract pKa predictions from PROPKA MolecularContainer.

    Returns list of dicts with keys: restype, chain, resnum, icode, pka, model_pka
    """
    results = []
    conformation = mc.conformations.get('AVR')
    if conformation is None:
        # Fall back to first available conformation
        conformation = next(iter(mc.conformations.values()))

    for grp in conformation.groups:
        if grp.pka_value is None:
            continue
        results.append({
            "restype": grp.residue_type,
            "chain": grp.atom.chain_id,
            "resnum": grp.atom.res_num,
            "icode": grp.atom.icode.strip(),
            "pka": grp.pka_value,
            "model_pka": grp.model_pka,
        })

    return results


def decide_protonation(pka_results, ph, his_default, cys_ss_pka):
    """Decide protonation state for each titratable residue.

    Returns dict: (chain, resnum, icode) -> new_resname
    Only includes residues that need renaming (non-standard protonation).
    """
    renames = {}

    for r in pka_results:
        key = (r["chain"], r["resnum"], r["icode"])
        rt = r["restype"]
        pka = r["pka"]

        if rt == "HIS":
            if pka > ph:
                renames[key] = "HIP"  # doubly protonated (charged)
            else:
                renames[key] = his_default  # neutral tautomer
        elif rt == "ASP":
            if pka > ph:
                renames[key] = "ASH"  # protonated (neutral)
            # else: standard ASP (deprotonated, charged) — no rename needed
        elif rt == "GLU":
            if pka > ph:
                renames[key] = "GLH"  # protonated (neutral)
        elif rt == "CYS":
            if pka >= cys_ss_pka:
                renames[key] = "CYX"  # disulfide bridge
            elif pka < ph:
                renames[key] = "CYM"  # deprotonated thiolate
        elif rt == "LYS":
            if pka < ph:
                renames[key] = "LYN"  # neutral lysine

    return renames


def rename_residues(lines, renames):
    """Rename residues in PDB lines according to the renames dict.

    Matches on (chain_id, resnum, icode) and replaces resname at cols 17-19.
    Applies to ATOM, HETATM, TER, ANISOU lines.
    """
    output = []
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM", "ANISOU") and len(line) > 26:
            chain = line[21]
            seq_str = line[22:26].strip()
            if seq_str and seq_str.lstrip('-').isdigit():
                resnum = int(seq_str)
                icode = line[26].strip() if len(line) > 26 else ""
                key = (chain, resnum, icode)
                if key in renames:
                    new_name = renames[key]
                    line = line[:17] + f"{new_name:<3s}" + line[20:]
        elif rec == "TER" and len(line) > 26:
            chain = line[21]
            seq_str = line[22:26].strip()
            if seq_str and seq_str.lstrip('-').isdigit():
                resnum = int(seq_str)
                icode = line[26].strip() if len(line) > 26 else ""
                key = (chain, resnum, icode)
                if key in renames:
                    new_name = renames[key]
                    line = line[:17] + f"{new_name:<3s}" + line[20:]
        output.append(line)
    return output


def _strip_hydrogens(topology, positions):
    """Remove all hydrogen atoms from topology/positions."""
    from openmm.app import Modeller
    h_atoms = [a for a in topology.atoms() if a.element.symbol == 'H']
    if not h_atoms:
        return topology, positions
    modeller = Modeller(topology, positions)
    modeller.delete(h_atoms)
    return modeller.topology, modeller.positions


_VARIANT_TO_PARENT_FOR_ADDHYDROGENS = {
    'LYN': 'LYS', 'CYX': 'CYS', 'CYM': 'CYS',
    'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS',
    'ASH': 'ASP', 'GLH': 'GLU',
}


def _rename_variants_to_parent(top):
    """Rename variant residues (LYN/HID/HIE/CYX/...) to their standard parent
    names so OpenMM's addHydrogens can find them in hydrogens.xml (keyed by
    parent name only). Returns dict {(chain_id, res_id): original_name}.
    Residue references in the OLD topology are stale after addHydrogens
    rebuilds the topology — use this dict for the post-pass lookup.
    """
    saved = {}
    for res in top.residues():
        if res.name in _VARIANT_TO_PARENT_FOR_ADDHYDROGENS:
            saved[(res.chain.id, res.id)] = res.name
            res.name = _VARIANT_TO_PARENT_FOR_ADDHYDROGENS[res.name]
    return saved


def _restore_variants_post_addhydrogens(top, saved):
    """Walk the NEW topology after addHydrogens and restore variant names."""
    for res in top.residues():
        key = (res.chain.id, res.id)
        if key in saved:
            res.name = saved[key]


def _fix_lyn_hz_naming(top, saved, renames):
    """Rename HZ1 → HZ3 on every LYN residue (chemically equivalent — HZ2
    and HZ3 in the AMBER LYN template share the same charge and bond
    topology). hydrogens.xml gates HZ3 by variant="LYS" so addHydrogens with
    variant=LYN produces HZ1+HZ2 — the OPPOSITE of what the AMBER LYN
    template expects (HZ2+HZ3).

    LYN identification merges two sources: residues whose pre-addHydrogens
    name was LYN (from `saved`), and residues whose PROPKA-assigned variant
    is LYN (from `renames`).
    """
    lyn_keys = {k for k, v in saved.items() if v == 'LYN'}
    for (ch, rs, ic), variant in renames.items():
        if variant == 'LYN':
            lyn_keys.add((ch, str(rs)))
    for res in top.residues():
        if (res.chain.id, res.id) not in lyn_keys:
            continue
        for atom in res.atoms():
            if atom.name == 'HZ1':
                atom.name = 'HZ3'
                break


def _add_hydrogens_to_output(input_path, output_path, args, renames):
    """Load original PDB, strip H, add H with correct protonation variants, write output."""
    from openmm.app import ForceField, Modeller, PDBFile
    from dvbfixer.ffutils import (PROTEIN_RESIDUES, SOLVENT_IONS,
                                   detect_glycam_input,
                                   create_forcefield_with_openff,
                                   fix_atom_hetatm_records)

    pdb = PDBFile(str(input_path))

    # Detect GLYCAM residues. If present, we keep ALL residues in the system
    # and use AMBER14+GLYCAM (which has NLN/OLS/OLT/sugar templates) instead
    # of ff19SB (which has none).
    info = detect_glycam_input(pdb.topology)
    glycam_present = bool(info['glycam_proteins'] or info['glycam_sugars'])

    # Strip existing hydrogens
    topology, positions = _strip_hydrogens(pdb.topology, pdb.positions)

    known = PROTEIN_RESIDUES | SOLVENT_IONS
    # GLYCAM mode: keep heterogens in the topology (GLYCAM FF parametrizes
    # NLN/OLS/OLT + sugar residues directly). Otherwise strip heterogens and
    # re-append them verbatim at write time.
    if glycam_present:
        to_delete = []
    else:
        to_delete = [res for res in topology.residues() if res.name not in known]
    has_hetatm = len(to_delete) > 0

    if has_hetatm:
        modeller = Modeller(topology, positions)
        modeller.delete(to_delete)
        stripped_top, stripped_pos = modeller.topology, modeller.positions
    else:
        stripped_top, stripped_pos = topology, positions

    # Build variants list: one entry per residue in stripped topology
    # Maps PROPKA renames to OpenMM variant names
    variants = []
    for res in stripped_top.residues():
        chain = res.chain.id
        try:
            resnum = int(res.id)
        except ValueError:
            variants.append(None)
            continue
        icode = res.insertionCode.strip() if hasattr(res, 'insertionCode') else ''
        key = (chain, resnum, icode)
        if key in renames:
            variants.append(renames[key])
        else:
            variants.append(None)

    print("Adding hydrogens for assigned protonation states...")
    modeller = Modeller(stripped_top, stripped_pos)

    if glycam_present:
        # Build GLYCAM-aware FF. create_forcefield_with_openff loads the
        # provided XMLs (caller upgraded args.ff to amber14+GLYCAM) and
        # registers SMIRNOFF for any unknown ligands.
        forcefield = create_forcefield_with_openff(
            args.ff, modeller.topology, verbose=args.verbose,
        )
        # Load GLYCAM-specific H definitions so addHydrogens knows where to
        # place H atoms on UYB/4YB/VMB/NLN/OLS/OLT/etc.
        try:
            Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')
        except Exception:
            pass
        # Add intra-residue bonds for GLYCAM residues — OpenMM's PDBFile
        # doesn't infer them for non-standard residues, but addHydrogens
        # needs them to place H correctly.
        try:
            from dvbfixer.acpype_export import add_glycam_bonds
            add_glycam_bonds(modeller.topology, forcefield, args.verbose,
                              positions=modeller.positions)
        except Exception as e:
            if args.verbose:
                print(f"  add_glycam_bonds skipped: {e}")
    else:
        forcefield = ForceField(*args.ff)

    # Use PDBFixer to add any missing heavy atoms first
    from pdbfixer import PDBFixer
    import tempfile as _tempfile
    with _tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as _tmp:
        PDBFile.writeFile(modeller.topology, modeller.positions, _tmp, keepIds=True)
        _tmp_path = _tmp.name
    try:
        with open(_tmp_path) as _f:
            fixer = PDBFixer(pdbfile=_f)
        fixer.findMissingResidues()
        fixer.missingResidues = {}  # Don't add missing residues, only atoms
        # In GLYCAM mode, prevent PDBFixer from "fixing" NLN/OLS/OLT (it
        # doesn't recognize them as standard residues and may strip/replace).
        if glycam_present:
            from dvbfixer.ffutils import GLYCAM_PROTEIN_RESIDUES
            # PDBFixer's substitutions[chain_idx] is a list of (res_idx, new_name).
            # Just clear it — we don't want any substitutions.
            try:
                fixer.findNonstandardResidues()
                fixer.nonstandardResidues = []
            except Exception:
                pass
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        modeller = Modeller(fixer.topology, fixer.positions)
        # Rebuild variants list for potentially reordered topology
        variants = []
        for res in modeller.topology.residues():
            chain = res.chain.id
            try:
                resnum = int(res.id)
            except ValueError:
                variants.append(None)
                continue
            icode = res.insertionCode.strip() if hasattr(res, 'insertionCode') else ''
            key = (chain, resnum, icode)
            if key in renames:
                variants.append(renames[key])
            else:
                variants.append(None)
    finally:
        Path(_tmp_path).unlink()

    # Rename variant residues to parent so addHydrogens can find them in
    # hydrogens.xml (keyed by parent name). After addHydrogens, restore the
    # variant names AND patch HZ1→HZ3 on LYN residues (AMBER LYN expects
    # HZ2+HZ3; addHydrogens with variant=LYN produces HZ1+HZ2 because of an
    # OpenMM hydrogens.xml vs AMBER template inconsistency).
    _saved = _rename_variants_to_parent(modeller.topology)
    try:
        if glycam_present:
            # ignoreExternalBonds=True: the protein-glycan ND2-C1 bond doesn't
            # match any single template, so addHydrogens must tolerate it.
            modeller.addHydrogens(forcefield, pH=args.ph, variants=variants)
        else:
            modeller.addHydrogens(forcefield, pH=args.ph, variants=variants)
    finally:
        _fix_lyn_hz_naming(modeller.topology, _saved, renames)
        _restore_variants_post_addhydrogens(modeller.topology, _saved)

    # Rename residues in the final topology to match AMBER names
    for res in modeller.topology.residues():
        chain = res.chain.id
        try:
            resnum = int(res.id)
        except ValueError:
            continue
        icode = res.insertionCode.strip() if hasattr(res, 'insertionCode') else ''
        key = (chain, resnum, icode)
        if key in renames:
            res.name = renames[key]

    if has_hetatm:
        with open(output_path, 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)

        # Read back and insert original HETATM/CONECT before END
        with open(output_path) as f:
            protein_lines = f.readlines()

        with open(str(input_path)) as f:
            orig_lines = f.readlines()
        hetatm_lines = [l for l in orig_lines if l.startswith("HETATM")]
        conect_lines = [l for l in orig_lines if l.startswith("CONECT")]

        with open(output_path, 'w') as f:
            for line in protein_lines:
                if line.startswith("END"):
                    f.writelines(hetatm_lines)
                    f.writelines(conect_lines)
                f.write(line)
    else:
        with open(output_path, 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)

    # Rewrite HETATM→ATOM for AMBER protonation variants (HID/HIE/HIP/ASH/
    # GLH/CYX/CYM/LYN) and GLYCAM glycoprotein residues (NLN/OLS/OLT).
    # PDBFile.writeFile defaults non-standard names to HETATM.
    fix_atom_hetatm_records(output_path)

    n_h = sum(1 for a in modeller.topology.atoms() if a.element.symbol == 'H')
    print(f"Added {n_h} hydrogen atoms")


def _scan_glycam_residues(input_path):
    """Quick text scan of input PDB. Returns (glycam_positions, has_glycam)
    where glycam_positions is a set of (chain, resnum, icode) tuples for
    NLN/OLS/OLT (and any GLYCAM sugar) in the file.
    """
    from dvbfixer.ffutils import is_glycam_residue
    positions = set()
    has_glycam = False
    with open(input_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            resname = line[17:20].strip()
            if is_glycam_residue(resname):
                has_glycam = True
                chain = line[21]
                seq_str = line[22:26].strip()
                if seq_str and seq_str.lstrip('-').isdigit():
                    resnum = int(seq_str)
                    icode = line[26].strip() if len(line) > 26 else ''
                    positions.add((chain, resnum, icode))
    return positions, has_glycam


# Default FF args (in argparse) so we can detect whether the user changed them
_DEFAULT_FF = ['amber19/protein.ff19SB.xml', 'amber19/tip3p.xml']


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_prot")

    # GLYCAM detection: if input has NLN/OLS/OLT or GLYCAM-named sugars,
    # switch FF to AMBER14+GLYCAM (ff19SB has no GLYCAM templates → crash).
    glycam_positions, has_glycam = _scan_glycam_residues(input_path)
    if has_glycam and args.ff == _DEFAULT_FF:
        args.ff = ['amber14-all.xml', 'amber14/GLYCAM_06j-1.xml',
                   'amber14/tip3pfb.xml']
        if args.verbose:
            print(f"GLYCAM residues detected → using FF: {' '.join(args.ff)}")

    print(f"Running PROPKA3 on {input_path} at pH {args.ph}...")
    mc = run_propka(input_path)
    pka_results = get_pka_results(mc)

    if not pka_results:
        print("No titratable residues found.", file=sys.stderr)
        sys.exit(1)

    if args.summary:
        print(f"\n{'Type':<4s} {'Chain':<5s} {'ResNum':<8s} {'pKa':>6s}  {'Model':>6s}  State at pH {args.ph}")
        print("-" * 55)
        for r in sorted(pka_results, key=lambda x: (x["chain"], x["resnum"])):
            icode_str = r["icode"] if r["icode"] else " "
            resid_str = f"{r['resnum']}{icode_str.strip()}"
            # Determine state label
            state = "standard"
            key = (r["chain"], r["resnum"], r["icode"])
            renames = decide_protonation(pka_results, args.ph, args.his_default, args.cys_disulfide_pka)
            if key in renames:
                state = renames[key]
            print(f"{r['restype']:<4s} {r['chain']:<5s} {resid_str:<8s} {r['pka']:6.2f}  {r['model_pka']:6.2f}  {state}")
        print()

    renames = decide_protonation(pka_results, args.ph, args.his_default, args.cys_disulfide_pka)

    # PROPKA saw the GLYCAM glycoprotein residues as ASN/SER/THR (after our
    # temp-PDB rename). Don't apply those protonation renames to the actual
    # NLN/OLS/OLT positions — they have different chemistry (the sidechain N
    # or O is bonded to a sugar, not protonated).
    if glycam_positions:
        renames = {k: v for k, v in renames.items() if k not in glycam_positions}

    # Preserve AMBER protonation variants already present in the input. PROPKA
    # ran on the sanitized PDB (which had its residues renamed to canonical
    # parents), so it might predict "standard" protonation for a residue that
    # the user had explicitly labeled HID/HIE/HIP/CYX/etc. in the input.
    # If `renames` has no entry for that position, carry the input variant
    # name forward so the output preserves it.
    _AMBER_VARIANT_NAMES = {'HID', 'HIE', 'HIP', 'ASH', 'GLH',
                             'CYX', 'CYM', 'LYN'}
    input_variants = {}
    with open(input_path) as _f:
        for _ln in _f:
            if not _ln.startswith(('ATOM', 'HETATM')):
                continue
            _rn = _ln[17:20].strip()
            if _rn not in _AMBER_VARIANT_NAMES:
                continue
            _ch = _ln[21]
            _ss = _ln[22:26].strip()
            if _ss and _ss.lstrip('-').isdigit():
                _rs = int(_ss)
                _ic = _ln[26].strip() if len(_ln) > 26 else ''
                input_variants[(_ch, _rs, _ic)] = _rn
    for key, name in input_variants.items():
        if key not in renames and key not in glycam_positions:
            renames[key] = name

    if args.verbose or args.summary:
        if renames:
            print(f"Non-standard protonation at pH {args.ph}:")
            for (chain, resnum, icode), new_name in sorted(renames.items()):
                ic_str = icode if icode else ""
                # Find original restype
                orig = next((r["restype"] for r in pka_results
                             if r["chain"] == chain and r["resnum"] == resnum
                             and r["icode"] == icode), "???")
                pka = next((r["pka"] for r in pka_results
                            if r["chain"] == chain and r["resnum"] == resnum
                            and r["icode"] == icode), 0.0)
                print(f"  {chain}/{orig} {resnum}{ic_str} -> {new_name}  (pKa={pka:.2f})")
        else:
            print(f"All residues have standard protonation at pH {args.ph}")

    # Read and rename
    with open(input_path) as f:
        lines = f.readlines()

    output_lines = rename_residues(lines, renames)

    if not args.keep_water:
        filtered = []
        for line in output_lines:
            if line.startswith(("ATOM  ", "HETATM", "ANISOU")):
                resname = line[17:20].strip()
                if resname in WATER_RESNAMES:
                    continue
            elif line.startswith("TER") and len(line) > 20:
                resname = line[17:20].strip()
                if resname in WATER_RESNAMES:
                    continue
            filtered.append(line)
        output_lines = filtered

    if not args.no_hydrogens:
        _add_hydrogens_to_output(input_path, output_path, args, renames)
    else:
        with open(output_path, 'w') as f:
            f.writelines(output_lines)

    n_his = sum(1 for v in renames.values() if v in ("HIP", "HIE", "HID"))
    n_other = len(renames) - n_his
    print(f"Wrote {output_path} ({n_his} HIS renamed, {n_other} other residues renamed)")
