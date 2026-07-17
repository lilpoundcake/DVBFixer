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

from dvbfixer.ffutils.variants import (
    fix_lyn_hz_naming as _fix_lyn_hz_naming_shared,
)
from dvbfixer.ffutils.variants import (
    rename_variants_to_parent_in_topology as _rename_variants_to_parent,
)
from dvbfixer.ffutils.variants import (
    restore_variants_post_addhydrogens as _restore_variants_post_addhydrogens_shared,
)
from dvbfixer.ffutils.variants import (
    text_rename_variants_to_parent as _text_rename_variants_to_parent_shared,
)

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
        "--ff", nargs='+', default=['auto'],
        help="Force field for hydrogen addition. Accepts a short name "
             "(auto, amber, amber+glycam, charmm, ...) or an explicit list "
             "of OpenMM XML paths. Default: 'auto' — detect from residue "
             "names in the input. See docs/force-fields.md."
    )
    p.add_argument(
        "--keep-water", action="store_true",
        help="Keep water molecules (HOH, WAT, TIP3, SOL) in output (default: remove)"
    )
    p.add_argument(
        "--protassign", action=argparse.BooleanOptionalAction, default=True,
        help="Run MolProbity Reduce to optimise HIS tautomers (HID/HIE/HIP) "
             "and detect ASN/GLN side-chain flips based on local H-bond "
             "network. **Default ON** — gives every protonate run the same "
             "higher-quality H-network without remembering a flag. Pass "
             "--no-protassign to disable (PROPKA-only pH-driven decisions). "
             "Requires the `reduce` binary (bundled with AmberTools in the "
             "dvbfixer env)."
    )
    p.add_argument(
        "--protassign-binary", dest="protassign_binary", default=None,
        help="Override the `reduce` binary path (default: search PATH, then "
             "the dvbfixer env's bin dir)."
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


# AMBER-style variant → CHARMM-style equivalent. Used when the resolved
# FF is CHARMM: PROPKA + decide_protonation always produce AMBER names
# (because that's what OpenMM's hydrogens.xml recognises for addHydrogens);
# after addHydrogens completes, the OUTPUT PDB gets rewritten to
# CHARMM-native names so downstream CHARMM tools recognise them.
#
# Only the variants that OpenMM's shipped `charmm36.xml` actually contains
# as templates are listed here. ASH/GLH/LYN have NO OpenMM-CHARMM36
# equivalent (CHARMM-GUI uses ASPP/GLUP/LSN as 4-char names but those
# aren't in the OpenMM XML), so on `--ff charmm` those residues fall back
# to their standard parent (ASP/GLU/LYS). Warned in the log.
_AMBER_TO_CHARMM_VARIANT = {
    'HID': 'HSD',
    'HIE': 'HSE',
    'HIP': 'HSP',
    'CYX': 'CYS',   # CHARMM tracks disulfide via SSBOND/DISU patch, not resname
    'CYM': 'CYM',   # CHARMM36 also uses CYM for the deprotonated thiolate
}

# AMBER variants that CHARMM36 has no template for — mapped to their
# standard parent so the pipeline still succeeds (charge state may not
# be what PROPKA asked for; user warned).
_AMBER_UNSUPPORTED_IN_CHARMM = {
    'ASH': 'ASP',
    'GLH': 'GLU',
    'LYN': 'LYS',
}


def _is_charmm_ff(ff_xmls):
    """True when the resolved --ff list is CHARMM (any charmm*.xml present)."""
    return any('charmm' in str(x).lower() for x in (ff_xmls or ()))


def decide_protonation(pka_results, ph, his_default, cys_ss_pka):
    """Decide protonation state for each titratable residue.

    Returns dict: (chain, resnum, icode) -> new_resname
    Only includes residues that need renaming (non-standard protonation).

    Names always use AMBER-style variants (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN)
    because those are the only names OpenMM's `hydrogens.xml` recognises
    in its `variants=` argument. For CHARMM outputs, the caller remaps
    these to CHARMM-native names (HSD/HSE/HSP/ASPP/GLUP/CYS/CYM/LSN) via
    `_AMBER_TO_CHARMM_VARIANT` AFTER addHydrogens completes.

    Note: PROPKA also reports pKas for TYR and ARG (and terminal N+/C-),
    but neither AMBER14/19 nor CHARMM36 ship a `TYD`/`TYN`/`ARN` variant
    template; there is no way to add a "deprotonated tyrosinate" or
    "neutral arginine" hydrogen set via the stock addHydrogens flow.
    Callers should log those pKas but not rename.
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
        # TYR / ARG: no FF-supported deprotonated/neutral variant. The
        # caller logs the pKa in --verbose mode but skips the rename.

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


def _text_rename_variants_to_parent(pdb_path, verbose=False):
    """Delegate to ``ffutils.variants.text_rename_variants_to_parent``.

    Kept as a thin wrapper so the module-level call sites and any external
    tests importing this name remain valid.
    """
    return _text_rename_variants_to_parent_shared(pdb_path, verbose=verbose)


def _restore_variants_post_addhydrogens(top, saved):
    """Delegate to ``ffutils.variants.restore_variants_post_addhydrogens``."""
    _restore_variants_post_addhydrogens_shared(top, saved)


def _fix_lyn_hz_naming(top, saved, renames):
    """Delegate to ``ffutils.variants.fix_lyn_hz_naming``.

    ``renames`` here is protonate.py's PROPKA-produced
    ``{(chain, resseq, icode): variant}`` map — pass it through as the
    ``extra_lyn_keys`` argument.
    """
    _fix_lyn_hz_naming_shared(top, saved, extra_lyn_keys=renames)


def _add_hydrogens_to_output(input_path, output_path, args, renames):
    """Load original PDB, strip H, add H with correct protonation variants, write output."""
    from openmm.app import ForceField, Modeller, PDBFile

    from dvbfixer.ffutils import (
        PROTEIN_RESIDUES,
        SOLVENT_IONS,
        create_forcefield_with_openff,
        detect_glycam_input,
        fix_atom_hetatm_records,
    )

    # Rename AMBER + CHARMM variant residues (LYN/HID/HIE/HIP/ASH/GLH/
    # CYX/CYM/HSD/HSE/HSP/ASPP/GLUP/LSN) to their standard parent in the
    # PDB TEXT before OpenMM's PDBFile parses it. OpenMM's parser
    # recognises only standard AA names for peptide-bond inference: a
    # residue named "LYN" mid-chain is treated as unknown, its peptide
    # bond to the previous residue is NOT inferred, and downstream
    # addHydrogens fails with the confusing "residue X (ASN) missing 1 C
    # atom externally bonded" (the missing bond is on the LYN-neighbour
    # ASN, not on ASN itself). Renaming to LYS/HIS/… in the raw text
    # lets OpenMM parse the chain correctly; `_pretext_saved` tracks
    # the original names so we can restore them on the topology after
    # addHydrogens completes.
    input_path, _pretext_saved = _text_rename_variants_to_parent(
        input_path, verbose=args.verbose)

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
        # suppresses GLYCAM sugar templates that would fuzzy-match PDB
        # sugar residues to the wrong entry.
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
    import tempfile as _tempfile

    from pdbfixer import PDBFixer
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
    # Two sources of variant-rename tracking:
    #   1. `_pretext_saved` — populated by the raw-text pre-rename above.
    #      Keys are (chain, str(resid), icode.strip()) as they appear in
    #      the PDB text; values are original variant names (LYN, HSE, ...).
    #   2. `_rename_variants_to_parent(topology)` — legacy topology-level
    #      rename. If the pre-rename already handled everything, this
    #      finds nothing new. Kept for defence in depth (e.g. GLYCAM
    #      topology-only paths that might introduce variants after load).
    _saved = _rename_variants_to_parent(modeller.topology)
    # Merge the pre-rename map into _saved so restore covers both.
    # Topology-level keys use res.chain.id + res.id (both strings) —
    # normalise the pre-rename keys accordingly.
    for (ch, rs, ic), orig in _pretext_saved.items():
        _saved.setdefault((ch, rs), orig)
    try:
        try:
            if glycam_present:
                # ignoreExternalBonds=True: the protein-glycan ND2-C1 bond
                # doesn't match any single template, so addHydrogens must
                # tolerate it.
                modeller.addHydrogens(forcefield, pH=args.ph,
                                       variants=variants)
            else:
                modeller.addHydrogens(forcefield, pH=args.ph,
                                       variants=variants)
        except Exception as e:
            # OpenMM's raw error uses topology INDEX ("residue 117 (ASN)"),
            # not the PDB resseq the user knows. `explain_template_error`
            # translates + shows chain/resseq/icode/atom-set diffs and
            # the neighbouring residues so the user can pinpoint the
            # actual problematic residue in the input PDB.
            from dvbfixer.ffutils import explain_template_error
            diag = explain_template_error(e, modeller.topology, forcefield)
            print(f"\nERROR: protonate addHydrogens failed:\n  {e}",
                  file=sys.stderr)
            if diag:
                for line in diag.split('\n'):
                    print(f"  {line}", file=sys.stderr)
            print(
                "\n  Common causes:\n"
                "    - a protein residue is a HETATM (peptide bond not\n"
                "      inferred → 'missing external C atom'). protonate\n"
                "      already runs `sanitize_protein_hetatm` on the\n"
                "      input; if the message persists, an upstream tool\n"
                "      is producing an unusual atom layout.\n"
                "    - a chain break between two consecutive residues\n"
                "      (prev-C to this-N distance too large for OpenMM).\n"
                "    - an FF-external atom left over from a previous step.\n"
                "  Re-run the upstream step with -v to inspect.\n",
                file=sys.stderr)
            raise
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

    # FF-aware output naming: if the resolved FF is CHARMM, remap the
    # AMBER-style variant names we wrote (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN)
    # to their CHARMM-native equivalents (HSD/HSE/HSP/ASPP/GLUP/CYS/CYM/
    # LSN) so downstream CHARMM tools recognise them. addHydrogens itself
    # ran with the AMBER names because that's the only vocabulary
    # `hydrogens.xml` speaks.
    if _is_charmm_ff(getattr(args, 'ff', None)):
        _remap_amber_variants_to_charmm_in_pdb(output_path,
                                                verbose=args.verbose)

    n_h = sum(1 for a in modeller.topology.atoms() if a.element.symbol == 'H')
    print(f"Added {n_h} hydrogen atoms")


def _remap_amber_variants_to_charmm_in_pdb(pdb_path, verbose=False):
    """Rewrite AMBER-style variant resnames in an already-written PDB
    file to their CHARMM36 equivalents:
      HID→HSD, HIE→HSE, HIP→HSP, CYX→CYS, CYM→CYM (idempotent).
    ASH/GLH/LYN have NO OpenMM-CHARMM36 template and get folded to their
    STANDARD parent (ASP/GLU/LYS) — the protonation state PROPKA
    requested cannot be expressed with the current CHARMM36 XML. Emits
    an INFO line noting the dropped protonation whenever this happens.

    All resnames written stay in the 3-char slot (cols 17-19), so column
    alignment is safe. In-place. Idempotent.
    """
    with open(pdb_path) as f:
        lines = f.readlines()
    n_rewrite = 0
    n_dropped = 0
    dropped_names = set()
    out = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')) and len(line) >= 20:
            rn = line[17:20].strip()
            new_rn = _AMBER_TO_CHARMM_VARIANT.get(rn)
            if new_rn is not None and new_rn != rn:
                # All CHARMM equivalents are 3-char, safe in-place.
                line = line[:17] + f"{new_rn:<3s}" + line[20:]
                n_rewrite += 1
            elif rn in _AMBER_UNSUPPORTED_IN_CHARMM:
                # Drop to the standard parent — CHARMM36 has no template
                # for this protonation state.
                parent = _AMBER_UNSUPPORTED_IN_CHARMM[rn]
                line = line[:17] + f"{parent:<3s}" + line[20:]
                n_dropped += 1
                dropped_names.add(rn)
        out.append(line)
    if n_rewrite or n_dropped:
        with open(pdb_path, 'w') as f:
            f.writelines(out)
        if verbose:
            if n_rewrite:
                print(f"  [protonate] remapped {n_rewrite} AMBER variant "
                      f"lines to CHARMM36 equivalents in {pdb_path}")
            if n_dropped:
                names = ', '.join(sorted(dropped_names))
                print(f"  WARNING: {n_dropped} residue(s) with AMBER-only "
                      f"variants ({names}) folded back to standard parent "
                      f"— OpenMM's charmm36.xml has no template for "
                      f"protonated ASP/GLU or neutral LYS. Charge state "
                      f"differs from PROPKA's prediction.")


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


# Kept for backwards-compat imports; FF selection now flows through
# ffutils.resolve_ff (see main()).
_DEFAULT_FF = ['amber19/protein.ff19SB.xml', 'amber19/tip3p.xml']


# ---------------------------------------------------------------------------
# ProtAssign-style optimisation via MolProbity Reduce
#
# `reduce -build -flip` runs the Word-Lovell-Richardson 1999 algorithm:
# local H-bond network + van-der-Waals clash optimisation over HIS tautomers
# (HID/HIE/HIP) and ASN/GLN side-chain 180° flips.
#
# We wrap reduce, then parse its output to extract:
#   - HIS tautomer per residue (from which H atoms reduce placed)
#   - ASN/GLN flip decisions (from heavy-atom coordinate diff vs input)
# The decisions are merged with PROPKA's pH-driven renames (PROPKA wins
# HIP — that's pKa-driven, more reliable than reduce's local count).
# ---------------------------------------------------------------------------

def _find_reduce_binary(override=None):
    """Locate the `reduce` binary. Returns Path or None.

    Search order: explicit override → PATH → directory of the running
    Python interpreter (env's bin dir, where conda installs binaries).
    """
    import os
    import shutil
    from pathlib import Path
    if override:
        p = Path(override)
        return p if p.is_file() and os.access(p, os.X_OK) else None
    hit = shutil.which('reduce')
    if hit:
        return Path(hit)
    env_bin = Path(sys.executable).parent / 'reduce'
    if env_bin.is_file() and os.access(env_bin, os.X_OK):
        return env_bin
    return None


def _run_reduce(pdb_path, binary):
    """Run `reduce -build -flip -quiet -noheterogens` on pdb_path.

    Returns the path to a temp output PDB. Caller deletes it.
    """
    import subprocess
    import tempfile
    fd, out_path = tempfile.mkstemp(suffix='.reduce.pdb', prefix='dvbfixer_')
    cmd = [str(binary), '-build', '-flip', '-quiet', '-noheterogens',
           str(pdb_path)]
    with open(out_path, 'w') as outf:
        result = subprocess.run(cmd, stdout=outf,
                                stderr=subprocess.PIPE, text=True)
    # Reduce returns 0 on success but uses non-zero codes for some warnings;
    # check that the output file is non-empty instead.
    import os as _os
    if _os.path.getsize(out_path) == 0:
        _os.unlink(out_path)
        raise RuntimeError(
            f"reduce produced empty output. stderr:\n{result.stderr}")
    return out_path


def _read_atom_positions(pdb_path):
    """Parse ATOM/HETATM records → {(chain, resseq, icode, atomname): (x, y, z), resname}.

    Returns two dicts keyed by (chain, resseq, icode, atomname):
      positions: (x, y, z) tuple in Å
      resnames: residue name string
    """
    positions = {}
    resnames = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            atom = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            ss = line[22:26].strip()
            if not ss or not ss.lstrip('-').isdigit():
                continue
            resseq = int(ss)
            icode = line[26].strip() if len(line) > 26 else ''
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            key = (chain, resseq, icode, atom)
            positions[key] = (x, y, z)
            resnames[key] = resname
    return positions, resnames


def _parse_reduce_decisions(reduce_pdb, input_pdb):
    """Extract HIS tautomer + ASN/GLN flip decisions from reduce's output.

    HIS tautomer is inferred from which H atoms reduce placed on each HIS:
      HD1 + HE2 → HIP
      HD1 only  → HID
      HE2 only  → HIE
      neither   → None (caller falls back to --his-default)

    ASN flip: in the input, OD1 and ND2 are at specific coordinates. In
    reduce's output, the atom NAMES are the same but COORDINATES may have
    been swapped (reduce flips by reassigning labels). We detect by
    checking whether |OD1_reduce - ND2_input| < 0.1 Å.

    Same for GLN: OE1, NE2.

    Returns (his_picks dict, asn_flip set, gln_flip set), all keyed by
    (chain, resseq, icode).
    """
    in_pos, in_res = _read_atom_positions(input_pdb)
    out_pos, out_res = _read_atom_positions(reduce_pdb)

    # Build a HIS residue index: which H atoms reduce placed?
    his_picks = {}
    his_residues = set()
    for (ch, rs, ic, atom), resname in out_res.items():
        if resname == 'HIS':
            his_residues.add((ch, rs, ic))
    for key in his_residues:
        ch, rs, ic = key
        has_hd1 = (ch, rs, ic, 'HD1') in out_pos
        has_he2 = (ch, rs, ic, 'HE2') in out_pos
        if has_hd1 and has_he2:
            his_picks[key] = 'HIP'
        elif has_hd1:
            his_picks[key] = 'HID'
        elif has_he2:
            his_picks[key] = 'HIE'
        # else: skip — let downstream decide

    def _swapped(in_a, in_b, out_a, tol=0.1):
        """True if out_a is at the position where in_b was (within tol Å)."""
        if in_a is None or in_b is None or out_a is None:
            return False
        d = sum((out_a[i] - in_b[i]) ** 2 for i in range(3)) ** 0.5
        return d < tol

    asn_flip = set()
    gln_flip = set()
    asn_residues = {(ch, rs, ic) for (ch, rs, ic, _), rn in in_res.items()
                    if rn == 'ASN'}
    gln_residues = {(ch, rs, ic) for (ch, rs, ic, _), rn in in_res.items()
                    if rn == 'GLN'}

    for key in asn_residues:
        ch, rs, ic = key
        in_od1 = in_pos.get((ch, rs, ic, 'OD1'))
        in_nd2 = in_pos.get((ch, rs, ic, 'ND2'))
        out_od1 = out_pos.get((ch, rs, ic, 'OD1'))
        if _swapped(in_od1, in_nd2, out_od1):
            asn_flip.add(key)

    for key in gln_residues:
        ch, rs, ic = key
        in_oe1 = in_pos.get((ch, rs, ic, 'OE1'))
        in_ne2 = in_pos.get((ch, rs, ic, 'NE2'))
        out_oe1 = out_pos.get((ch, rs, ic, 'OE1'))
        if _swapped(in_oe1, in_ne2, out_oe1):
            gln_flip.add(key)

    return his_picks, asn_flip, gln_flip


def _apply_flips_to_pdb_text(lines, asn_flips, gln_flips):
    """Swap heavy-atom coordinates for flipped ASN/GLN residues in PDB text.

    Atom names are kept the same — only the coordinate triples in cols
    30-54 are exchanged between OD1/ND2 (ASN) or OE1/NE2 (GLN). Hydrogens
    HD21/HD22 (resp HE21/HE22) are also swapped if present, for
    completeness — but they get stripped before addHydrogens anyway.
    """
    flip_pairs = {
        ('ASN', 'OD1'): 'ND2', ('ASN', 'ND2'): 'OD1',
        ('ASN', 'HD21'): 'HD22', ('ASN', 'HD22'): 'HD21',
        ('GLN', 'OE1'): 'NE2', ('GLN', 'NE2'): 'OE1',
        ('GLN', 'HE21'): 'HE22', ('GLN', 'HE22'): 'HE21',
    }
    # Build a coord lookup: (chain, resseq, icode, atomname) -> coord str
    coords = {}
    for line in lines:
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        resname = line[17:20].strip()
        if resname not in ('ASN', 'GLN'):
            continue
        chain = line[21]
        ss = line[22:26].strip()
        if not ss or not ss.lstrip('-').isdigit():
            continue
        resseq = int(ss)
        icode = line[26].strip() if len(line) > 26 else ''
        key = (chain, resseq, icode)
        flip_set = asn_flips if resname == 'ASN' else gln_flips
        if key not in flip_set:
            continue
        atom = line[12:16].strip()
        coords[(chain, resseq, icode, atom)] = line[30:54]

    # Rewrite lines, swapping the coord field for paired atoms.
    out = []
    for line in lines:
        if not line.startswith(('ATOM', 'HETATM')):
            out.append(line)
            continue
        resname = line[17:20].strip()
        if resname not in ('ASN', 'GLN'):
            out.append(line)
            continue
        chain = line[21]
        ss = line[22:26].strip()
        if not ss or not ss.lstrip('-').isdigit():
            out.append(line)
            continue
        resseq = int(ss)
        icode = line[26].strip() if len(line) > 26 else ''
        key = (chain, resseq, icode)
        flip_set = asn_flips if resname == 'ASN' else gln_flips
        if key not in flip_set:
            out.append(line)
            continue
        atom = line[12:16].strip()
        partner_atom = flip_pairs.get((resname, atom))
        if partner_atom is None:
            out.append(line)
            continue
        partner_coords = coords.get((chain, resseq, icode, partner_atom))
        if partner_coords is None:
            out.append(line)  # partner atom not present
            continue
        new_line = line[:30] + partner_coords + line[54:]
        out.append(new_line)
    return out


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_prot")

    # Sanitize protein-residue HETATM records + drop spurious mid-chain
    # TER records BEFORE any downstream OpenMM load. If a protein residue
    # comes in as HETATM (some upstream tools produce this), OpenMM's
    # PDBFile parser won't infer its peptide bond to neighbours and
    # `addHydrogens` fails with a confusing "missing 1 C atom externally
    # bonded" template error. No-op on clean input.
    from dvbfixer.ffutils import sanitize_protein_hetatm
    sanitized_path = sanitize_protein_hetatm(input_path, verbose=args.verbose)
    if sanitized_path != str(input_path):
        # Rewritten to a temp file — use the sanitized version for the
        # rest of protonate.
        input_path = Path(sanitized_path)

    # Position scan is still needed downstream (to filter GLYCAM residues
    # out of PROPKA-driven renames). FF selection now goes through the
    # shared resolver so short-names like --ff amber / charmm work.
    glycam_positions, has_glycam = _scan_glycam_residues(input_path)
    from dvbfixer.ffutils import print_ff_selection, resolve_ff
    args.ff, _ff_alias, _ff_reason = resolve_ff(
        args.ff, input_path, verbose=args.verbose)
    print_ff_selection(_ff_alias, _ff_reason, args.ff)

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

    # ProtAssign-style optimisation (opt-in via --protassign):
    # Wrap MolProbity Reduce to pick HIS tautomers and detect ASN/GLN flips
    # from local H-bond network + clash analysis. Merge into PROPKA renames.
    protassign_asn_flips = set()
    protassign_gln_flips = set()
    if args.protassign:
        reduce_binary = _find_reduce_binary(args.protassign_binary)
        if reduce_binary is None:
            print("ERROR: --protassign (default ON) requires the `reduce` "
                  "binary. Install AmberTools "
                  "(`conda install -c conda-forge ambertools`), pass "
                  "--protassign-binary PATH, or pass --no-protassign to "
                  "skip the optimisation (PROPKA-only mode).",
                  file=sys.stderr)
            sys.exit(1)
        print(f"Running MolProbity Reduce ({reduce_binary}) for HIS tautomers "
              f"+ ASN/GLN flip optimisation...")
        import os as _os
        try:
            reduce_out = _run_reduce(input_path, reduce_binary)
            his_picks, asn_flips, gln_flips = _parse_reduce_decisions(
                reduce_out, input_path)
        finally:
            try:
                _os.unlink(reduce_out)
            except (OSError, NameError, UnboundLocalError):
                pass
        protassign_asn_flips = asn_flips
        protassign_gln_flips = gln_flips
        # GLYCAM glycoprotein residues at NLN/OLS/OLT positions: skip — their
        # sidechain N/O is sugar-bonded, not a normal amide/hydroxyl.
        his_picks = {k: v for k, v in his_picks.items()
                     if k not in glycam_positions}
        asn_flips = {k for k in asn_flips if k not in glycam_positions}
        gln_flips = {k for k in gln_flips if k not in glycam_positions}
        # Overlay Reduce's HIS picks onto PROPKA's renames. PROPKA wins HIP
        # (pKa-driven, more reliable for charged-state decisions); Reduce
        # wins HID vs HIE (local H-bond geometry).
        applied_his = 0
        for key, variant in his_picks.items():
            current = renames.get(key)
            if current == 'HIP':
                continue  # PROPKA's HIP wins
            if current != variant:
                renames[key] = variant
                applied_his += 1
        if args.verbose:
            print(f"  --protassign: {applied_his} HIS tautomer override(s) from Reduce")
            print(f"  --protassign: {len(asn_flips)} ASN flip(s), "
                  f"{len(gln_flips)} GLN flip(s)")
            for k in sorted(asn_flips):
                print(f"    ASN flip at {k[0]}:{k[1]}{k[2]}")
            for k in sorted(gln_flips):
                print(f"    GLN flip at {k[0]}:{k[1]}{k[2]}")

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

        # PROPKA also computes pKas for TYR (deprotonated tyrosinate) and
        # ARG (neutral arginine). Neither AMBER14/19 nor CHARMM36 ships a
        # `TYD`/`TYN`/`ARN` variant template — there's no way to add a
        # "deprotonated tyrosinate" or "neutral arginine" hydrogen set via
        # the stock addHydrogens flow. Log any hits so the user knows to
        # handle manually (custom template, CpHMD, etc.).
        _unsupported = []
        for r in pka_results:
            rt = r["restype"]
            pka = r["pka"]
            if rt == "TYR" and pka < args.ph:
                _unsupported.append((r["chain"], r["resnum"],
                                     r["icode"], rt, pka,
                                     "deprotonated tyrosinate — no "
                                     "TYD/TYN variant in AMBER14/19 or "
                                     "CHARMM36; skipping rename"))
            elif rt == "ARG" and pka < args.ph:
                _unsupported.append((r["chain"], r["resnum"],
                                     r["icode"], rt, pka,
                                     "neutral arginine — no ARN variant "
                                     "in AMBER14/19 or CHARMM36; "
                                     "skipping rename"))
        if _unsupported:
            print(f"Unsupported PROPKA-titratable groups at pH {args.ph} "
                  f"(no FF variant available):")
            for ch, rs, ic, rt, pk, msg in _unsupported:
                ic_str = ic if ic else ""
                print(f"  {ch}/{rt} {rs}{ic_str}  (pKa={pk:.2f}) — {msg}")

    # Read and rename
    with open(input_path) as f:
        lines = f.readlines()

    # Apply ASN/GLN coord swaps from --protassign BEFORE the residue rename
    # (rename doesn't touch coordinates, so order doesn't really matter, but
    # this keeps the data flow easy to follow).
    if protassign_asn_flips or protassign_gln_flips:
        lines = _apply_flips_to_pdb_text(
            lines, protassign_asn_flips, protassign_gln_flips)

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
        # If we applied ASN/GLN flips, _add_hydrogens_to_output must see the
        # flipped coordinates. Write the flipped+renamed lines to a temp file
        # and use it as the OpenMM input.
        h_input_path = input_path
        h_tmp_path = None
        if protassign_asn_flips or protassign_gln_flips:
            import tempfile as _tempfile
            _fd, h_tmp_path = _tempfile.mkstemp(
                suffix='.pdb', prefix='dvbfixer_protassign_')
            with open(h_tmp_path, 'w') as _tf:
                _tf.writelines(output_lines)
            h_input_path = Path(h_tmp_path)
        try:
            _add_hydrogens_to_output(h_input_path, output_path, args, renames)
        finally:
            if h_tmp_path is not None:
                import os as _os
                try:
                    _os.unlink(h_tmp_path)
                except OSError:
                    pass
    else:
        with open(output_path, 'w') as f:
            f.writelines(output_lines)

    n_his = sum(1 for v in renames.values() if v in ("HIP", "HIE", "HID"))
    n_other = len(renames) - n_his
    print(f"Wrote {output_path} ({n_his} HIS renamed, {n_other} other residues renamed)")
