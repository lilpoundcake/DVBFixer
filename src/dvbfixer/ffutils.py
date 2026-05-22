"""Shared force field utilities for dvbfixer.

Provides OpenFF-based parametrization for non-standard residues (glycans,
ligands) that lack templates in standard AMBER force fields. Uses
SMIRNOFFTemplateGenerator from openmmforcefields to auto-parametrize
unknown residues on the fly.
"""

from openmm.app import ForceField

# SMILES for common glycan residues (from PDB Chemical Component Dictionary).
# These are used to create OpenFF Molecule objects for parametrization.
KNOWN_GLYCAN_SMILES = {
    'NAG': 'CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1O',   # N-acetyl-D-glucosamine
    'NDG': 'CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1O',   # N-acetyl-D-glucosamine (alt)
    'BMA': 'OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O',              # beta-D-mannose
    'MAN': 'OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O',              # alpha-D-mannose
    'FUC': 'C[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O',               # L-fucose
    'FUL': 'C[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O',               # L-fucose (alt)
    'GAL': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',               # D-galactose
    'BGC': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',               # beta-D-glucose
    'GLC': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',               # D-glucose
    'SIA': 'OC[C@@H](O)[C@@H]1OC(=O)[C@H](O)[C@@H](O)[C@@H]1NC(C)=O',  # sialic acid (simplified)
}

# Standard protein/water/ion residues that don't need OpenFF parametrization
PROTEIN_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX', 'CYM', 'ASH', 'GLH',
    'LYN', 'MSE', 'ACE', 'NME', 'NHE',
    # GLYCAM protein residues
    'NLN', 'OLS', 'OLT',
}
SOLVENT_IONS = {
    'HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN',
    'NA+', 'CL-', 'K+', 'MG2+', 'CA2+',
}

# GLYCAM force field naming detection. GLYCAM uses 3-char codes:
#   [linkage][sugar][anomer] (e.g. UYB, 4YB, VMB, 0YA)
# Plus glycoprotein residues (NLN/OLS/OLT) and reducing-end caps.
_GLYCAM_LINKAGE_CHARS = set('0123456789VWUZXYTSRQPvwuzxytsr')
_GLYCAM_ANOMER_CHARS = {'A', 'B'}
GLYCAM_PROTEIN_RESIDUES = {'NLN', 'OLS', 'OLT'}
GLYCAM_CAPS = {'ROH', 'OME', 'TBT', 'CMET'}

# Residue names that should be written as ATOM (not HETATM) in PDB output.
# OpenMM's PDBFile.writeFile defaults non-standard names to HETATM; this set
# is used by fix_atom_hetatm_records() to rewrite them after writing.
FORCE_ATOM_RESIDUES = frozenset(PROTEIN_RESIDUES)


def is_glycam_sugar(name):
    """True if `name` is a GLYCAM sugar code (3-char linkage+sugar+anomer or cap)."""
    if name in GLYCAM_CAPS:
        return True
    return (len(name) == 3
            and name[0] in _GLYCAM_LINKAGE_CHARS
            and name[2] in _GLYCAM_ANOMER_CHARS)


def is_glycam_residue(name):
    """True if `name` is any GLYCAM-named residue (sugar OR glycoprotein)."""
    return name in GLYCAM_PROTEIN_RESIDUES or is_glycam_sugar(name)


def detect_glycam_input(topology):
    """Scan topology for GLYCAM and PDB sugar residues.

    Returns dict with keys:
      - glycam_proteins: set of (chain_id, res_id) for NLN/OLS/OLT
      - glycam_sugars:   set of (chain_id, res_id) for GLYCAM-named sugars
      - pdb_sugars:      set of (chain_id, res_id) for known PDB sugars
                          (NAG, BMA, MAN, FUC, ...) — present in KNOWN_GLYCAN_SMILES
      - unknown_hets:    set of (chain_id, res_id) for anything non-protein
                          non-solvent that's not in the above
    """
    known_prot_solv = PROTEIN_RESIDUES | SOLVENT_IONS
    info = {
        'glycam_proteins': set(),
        'glycam_sugars': set(),
        'pdb_sugars': set(),
        'unknown_hets': set(),
    }
    for res in topology.residues():
        key = (res.chain.id, res.id)
        name = res.name
        if name in GLYCAM_PROTEIN_RESIDUES:
            info['glycam_proteins'].add(key)
        elif is_glycam_sugar(name):
            info['glycam_sugars'].add(key)
        elif name in KNOWN_GLYCAN_SMILES:
            info['pdb_sugars'].add(key)
        elif name not in known_prot_solv:
            info['unknown_hets'].add(key)
    return info


def fix_atom_hetatm_records(pdb_path):
    """Rewrite HETATM→ATOM for protein residues that OpenMM's PDBFile.writeFile
    incorrectly emitted as HETATM (AMBER protonation variants HID/HIE/HIP/
    ASH/GLH/CYX/CYM/LYN and GLYCAM glycoprotein residues NLN/OLS/OLT).

    Reads pdb_path, rewrites in place. Idempotent.
    """
    try:
        with open(pdb_path) as f:
            lines = f.readlines()
    except FileNotFoundError:
        return
    changed = False
    out = []
    for line in lines:
        if line.startswith('HETATM') and len(line) >= 20:
            resname = line[17:20].strip()
            if resname in FORCE_ATOM_RESIDUES:
                line = 'ATOM  ' + line[6:]
                changed = True
        out.append(line)
    if changed:
        with open(pdb_path, 'w') as f:
            f.writelines(out)


def _find_unknown_residue_names(topology):
    """Find residue names in topology that aren't protein/water/ions."""
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    unknown = set()
    for res in topology.residues():
        if res.name not in known:
            unknown.add(res.name)
    return unknown


def explain_template_error(exc, topology, forcefield=None):
    """Turn an opaque OpenMM template-match error into a useful diagnostic.

    OpenMM's `_matchAllResiduesToTemplates` reports failures by TOPOLOGY
    INDEX (0-based position in `topology.residues()` iteration order),
    not by the PDB resseq the user wrote. Error messages like

        "No template found for residue 181 (PHE). The set of heavy atoms
         matches PHE, but the residue is missing 2 H atoms."

    are misleading — there might be zero PHE residues with resseq 181 in
    the input PDB; "181" is just the position of THIS residue in OpenMM's
    iteration over the topology.

    This helper extracts the residue index from the error message, looks
    up the actual residue in the topology, and returns a multi-line
    string identifying the (chain, resseq, icode, resname) plus — when a
    forcefield is provided — the specific atom-set mismatch (missing /
    extra atoms) against the named template.

    Returns None if the error message doesn't match the OpenMM template-
    error format (caller should fall back to the original `str(exc)`).
    """
    import re

    msg = str(exc)
    m = re.search(r'residue\s+(\d+)\s+\(([A-Za-z0-9_]+)\)', msg)
    if not m:
        return None
    try:
        res_idx = int(m.group(1))
    except ValueError:
        return None
    expected_resname = m.group(2)

    # Look up the actual residue by topology iteration order.
    residues = list(topology.residues())
    if res_idx < 0 or res_idx >= len(residues):
        return None
    res = residues[res_idx]

    chain_id = res.chain.id if res.chain else '?'
    res_id = res.id
    icode = ''
    if hasattr(res, 'insertionCode') and res.insertionCode:
        icode = res.insertionCode.strip()
    res_name = res.name

    lines = [
        f"Failed residue (topology index {res_idx}, NOT PDB resseq):",
        f"  chain = {chain_id}    resseq = {res_id}{icode}    resname = {res_name}",
    ]
    if res_name != expected_resname:
        lines.append(
            f"  (OpenMM error said '{expected_resname}' — this is the template "
            f"name it tried to fit, not the input resname)"
        )

    # List the atom set we currently have for this residue.
    atom_names = [a.name for a in res.atoms()]
    lines.append(f"  atoms in topology ({len(atom_names)}): {' '.join(atom_names)}")

    # If a forcefield was provided, try to compute the actual atom-set
    # mismatch against the matching template.
    if forcefield is not None:
        tpl_name = None
        if res_name in forcefield._templates:
            tpl_name = res_name
        elif expected_resname in forcefield._templates:
            tpl_name = expected_resname
        if tpl_name is not None:
            template = forcefield._templates[tpl_name]
            tpl_atoms = {a.name for a in template.atoms}
            cur_atoms = set(atom_names)
            missing = sorted(tpl_atoms - cur_atoms)
            extra = sorted(cur_atoms - tpl_atoms)
            lines.append(f"  template '{tpl_name}' expects {len(tpl_atoms)} atoms")
            if missing:
                lines.append(f"  MISSING from input vs template: {' '.join(missing)}")
            if extra:
                lines.append(f"  EXTRA in input not in template: {' '.join(extra)}")
            if not missing and not extra:
                lines.append(
                    "  (atom names match — failure is likely from external-bond "
                    "expectations, not atom set)"
                )

    # Neighbour residues — often the real source of the problem (e.g. an
    # NLN whose adjacent ASN is missing its peptide bond, or a sugar tree
    # missing a glycosidic bond to a sibling).
    if res_idx > 0:
        prev = residues[res_idx - 1]
        lines.append(
            f"  prev residue (idx {res_idx-1}): {prev.chain.id}:{prev.name}{prev.id}"
        )
    if res_idx + 1 < len(residues):
        nxt = residues[res_idx + 1]
        lines.append(
            f"  next residue (idx {res_idx+1}): {nxt.chain.id}:{nxt.name}{nxt.id}"
        )

    return '\n'.join(lines)


def create_forcefield_with_openff(ff_xmls, topology, small_mol_ff='openff-2.2.0',
                                  extra_molecules=None, verbose=False):
    """Create OpenMM ForceField with automatic OpenFF parametrization for unknown residues.

    When PDB-named sugars (NAG, FUC, GAL, MAN, etc.) are present in the
    topology and GLYCAM_06j-1.xml was loaded, this function deletes GLYCAM's
    sugar/nucleic acid templates (keeping only NLN/OLS/OLT/ROH/etc.) so
    SMIRNOFF can handle the sugars cleanly. Without this, the 1400+ GLYCAM
    templates fuzzy-match PDB sugars to wrong templates (NAG → UVA, etc.).

    Args:
        ff_xmls: List of force field XML files (e.g., ['amber19/protein.ff19SB.xml', ...])
        topology: OpenMM Topology to scan for unknown residues
        small_mol_ff: OpenFF force field name (default: openff-2.2.0 Sage)
        extra_molecules: Optional list of additional openff.toolkit.Molecule objects
        verbose: Print info about parametrized residues

    Returns:
        ForceField object with SMIRNOFFTemplateGenerator registered for unknown residues
    """
    ff = ForceField(*ff_xmls)

    # If GLYCAM xml was loaded AND PDB-named sugars are in topology, suppress
    # GLYCAM sugar/NA templates (keep glycoprotein/cap templates).
    pdb_sugars = {r.name for r in topology.residues()
                  if r.name in KNOWN_GLYCAN_SMILES}
    glycam_loaded = any('GLYCAM' in str(x) for x in ff_xmls)
    if glycam_loaded and pdb_sugars:
        amber_only_xmls = [x for x in ff_xmls if 'GLYCAM' not in str(x)]
        if amber_only_xmls:
            amber_only = ForceField(*amber_only_xmls)
            glycam_extra = set(ff._templates) - set(amber_only._templates)
            _KEEP = {'NLN', 'OLS', 'OLT', 'ROH', 'OME', 'TBT', 'CMET'}
            removed = [n for n in glycam_extra if n not in _KEEP]
            for n in removed:
                del ff._templates[n]
            if verbose and removed:
                print(f"Suppressed {len(removed)} GLYCAM sugar/NA templates "
                      f"(PDB sugars detected → SMIRNOFF will handle them)")

    unknown = _find_unknown_residue_names(topology)
    if not unknown and not extra_molecules:
        return ff

    # Build Molecule objects for known glycans
    from openff.toolkit import Molecule
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator

    molecules = []
    matched = set()
    for resname in unknown:
        if resname in KNOWN_GLYCAN_SMILES:
            mol = Molecule.from_smiles(KNOWN_GLYCAN_SMILES[resname],
                                       allow_undefined_stereo=True)
            molecules.append(mol)
            matched.add(resname)

    if extra_molecules:
        molecules.extend(extra_molecules)

    unmatched = unknown - matched
    if unmatched and verbose:
        print(f"Warning: no SMILES for residues: {', '.join(sorted(unmatched))}")
        print("  These residues may cause errors. Use --sdf to provide molecule definitions.")

    if molecules:
        smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield=small_mol_ff)
        ff.registerTemplateGenerator(smirnoff.generator)
        if verbose:
            print(f"OpenFF parametrization registered for: {', '.join(sorted(matched))}")

        # Pre-generate templates for each PDB sugar residue in topology so they
        # get registered by exact name. Without this, OpenMM's fuzzy template
        # matcher tries to fit nucleotide templates (CN, A, etc.) to sugars
        # before falling back to the SMIRNOFF generator.
        seen_resnames = set()
        for res in topology.residues():
            if res.name not in matched or res.name in seen_resnames:
                continue
            seen_resnames.add(res.name)
            try:
                smirnoff.generator(ff, res)
            except Exception as e:
                if verbose:
                    print(f"  Pre-register SMIRNOFF template {res.name}: {e}")

    return ff
