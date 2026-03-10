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


def run_propka(input_path):
    """Run PROPKA3 on the input PDB and return the MolecularContainer."""
    from propka.run import single

    # Suppress PROPKA warnings to stderr
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        mc = single(str(input_path), write_pka=False)
    finally:
        sys.stderr = old_stderr

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


def _add_hydrogens_to_output(input_path, output_path, args, renames):
    """Load original PDB, strip H, add H with correct protonation variants, write output."""
    from openmm.app import ForceField, Modeller, PDBFile

    pdb = PDBFile(str(input_path))

    # Strip existing hydrogens
    topology, positions = _strip_hydrogens(pdb.topology, pdb.positions)

    # Strip non-protein residues that the force field can't handle
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    known = PROTEIN_RESIDUES | SOLVENT_IONS
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
    forcefield = ForceField(*args.ff)
    modeller = Modeller(stripped_top, stripped_pos)

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

    modeller.addHydrogens(forcefield, pH=args.ph, variants=variants)

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

    n_h = sum(1 for a in modeller.topology.atoms() if a.element.symbol == 'H')
    print(f"Added {n_h} hydrogen atoms")


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_prot")

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
