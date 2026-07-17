"""Shared force field utilities for dvbfixer.

FF selection (short-name aliases + auto-detection) via `resolve_ff`, shared
`fix_atom_hetatm_records`, GLYCAM residue detection, and a small helper
`create_forcefield_with_glycam_suppression` that prunes GLYCAM's 1400+
sugar / nucleic-acid templates when the input carries PDB-standard sugar
names (which otherwise fuzzy-match to the wrong GLYCAM template).

For arbitrary unknown ligands (cofactors, drug-like molecules) that lack a
template in the resolved FF, use `--parametrize-ligands` on `minimize` /
`zbs` — that runs `parametrize`'s GAFF2 + AM1-BCC pipeline and registers a
real OpenMM template per ligand. See `docs/force-fields.md`.
"""

from openmm.app import ForceField

# PDB-standard sugar residue names — the ambiguous set from the auto-detect
# path (a hit alone does not identify an FF; user is warned to convert first).
# Kept as a module-level constant because detect_ff_from_pdb consults it.
_PDB_SUGAR_NAMES = {
    'NAG', 'NDG', 'BMA', 'MAN', 'FUC', 'FUL', 'GAL', 'BGC', 'GLC', 'SIA',
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


# ---------------------------------------------------------------------------
# Shared FF selection: short-name aliases, auto-detection, and resolver.
# Used by prepare / minimize / protonate / pull / zbs. `top.py` uses a
# different (GROMACS FF-dir) namespace; see docs/force-fields.md.
# ---------------------------------------------------------------------------

FF_ALIASES = {
    'amber':          ['amber19/protein.ff19SB.xml', 'amber19/tip3p.xml'],
    'amber19':        ['amber19/protein.ff19SB.xml', 'amber19/tip3p.xml'],
    'amber14':        ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml'],
    'amber+glycam':   ['amber14-all.xml', 'amber14/GLYCAM_06j-1.xml',
                       'amber14/tip3pfb.xml'],
    'amber14+glycam': ['amber14-all.xml', 'amber14/GLYCAM_06j-1.xml',
                       'amber14/tip3pfb.xml'],
    'amber+lipid':    ['amber14-all.xml', 'amber14/lipid17.xml',
                       'amber14/tip3p.xml'],
    'amber+nucleic':  ['amber14-all.xml', 'amber14/DNA.OL15.xml',
                       'amber14/RNA.OL3.xml', 'amber14/tip3p.xml'],
    'charmm':         ['charmm36.xml', 'charmm36/water.xml'],
    'charmm36':       ['charmm36.xml', 'charmm36/water.xml'],
    'charmm2024':     ['charmm36_2024.xml', 'charmm36/water.xml'],
}

# Marker residue names that unambiguously identify a force field.
# CHARMM markers: CHARMM-specific protonation names + CHARMM-GUI 4-char sugar
# names + ceramide names. See top.CERAMIDE_RTP / glycam._CHARMM_4CHAR_RESNAMES.
_CHARMM_PROTONATION_MARKERS = {'HSD', 'HSE', 'HSP', 'ASPP', 'GLUP', 'LSN'}
_CHARMM_SUGAR_MARKERS = {
    'BGLC', 'AGLC', 'BMAN', 'AMAN', 'BGAL', 'AGAL', 'BFUC', 'AFUC',
    'BGLCNA', 'AGLCNA', 'BGALNA', 'AGALNA',
    'ANE5', 'BNE5', 'ANE5AC', 'BNE5AC',
    'AIDO', 'BIDO',
}
_CHARMM_CERAMIDE_MARKERS = {
    'CER1', 'CER160', 'CER180', 'CER181',
    'CER2', 'CER200', 'CER220', 'CER240', 'CER241', 'CER3E',
}
_CHARMM_MARKERS = (_CHARMM_PROTONATION_MARKERS
                   | _CHARMM_SUGAR_MARKERS
                   | _CHARMM_CERAMIDE_MARKERS)

# Ambiguous names — standard PDB Chemical Component Dictionary sugar names.
# Present in raw crystal PDBs / GLYCAM-pre-rename / CHARMM-pre-rename inputs.
# Neither amber14+GLYCAM (needs 3-char linkage codes like UYB, VMB) nor
# charmm36 (needs BGLCNA, BMAN, ...) has templates for these bare names.
# Detected only to emit a "convert first" warning, NOT to auto-pick an FF.
_AMBIGUOUS_SUGAR_MARKERS = frozenset(_PDB_SUGAR_NAMES)


def _scan_resnames(pdb_path):
    """Return set of resnames present in a PDB file (ATOM/HETATM lines).

    Reads both 3-char (cols 18-20) and CHARMM 4-char (cols 18-21) forms
    so CHARMM-GUI-style names like ASPP/BGLCNA are detected.
    """
    names = set()
    try:
        with open(pdb_path) as f:
            for line in f:
                if not (line.startswith('ATOM') or line.startswith('HETATM')):
                    continue
                # Try 4-char first (CHARMM-GUI); if the 4-char slice hits a
                # known CHARMM 4-char name we take it, otherwise fall back
                # to the standard 3-char slice.
                cand4 = line[17:21].strip() if len(line) >= 21 else ''
                if cand4 in _CHARMM_SUGAR_MARKERS or cand4 in _CHARMM_CERAMIDE_MARKERS:
                    names.add(cand4)
                    continue
                cand3 = line[17:20].strip()
                if cand3:
                    names.add(cand3)
    except (FileNotFoundError, OSError):
        pass
    return names


def detect_ff_from_pdb(pdb_path):
    """Scan a PDB file and return (alias, reason).

    - `('charmm', reason)` if any CHARMM-specific marker is found.
    - `('amber+glycam', reason)` if any GLYCAM-specific marker is found.
    - `('amber', reason)` otherwise (with a warning-worthy reason if only
      ambiguous PDB sugar names are present).

    CHARMM wins over GLYCAM when both markers appear (CHARMM protonation
    names are unambiguous FF-prep signals).
    """
    names = _scan_resnames(pdb_path)

    charmm_hits = names & _CHARMM_MARKERS
    if charmm_hits:
        sample = ', '.join(sorted(charmm_hits)[:3])
        return 'charmm', f'CHARMM residue(s) detected ({sample})'

    glycam_protein = names & GLYCAM_PROTEIN_RESIDUES
    glycam_caps = names & GLYCAM_CAPS
    glycam_sugars = {n for n in names if is_glycam_sugar(n)}
    glycam_hits = glycam_protein | glycam_caps | glycam_sugars
    if glycam_hits:
        sample = ', '.join(sorted(glycam_hits)[:3])
        return 'amber+glycam', f'GLYCAM residue(s) detected ({sample})'

    ambig_hits = names & _AMBIGUOUS_SUGAR_MARKERS
    if ambig_hits:
        sample = ', '.join(sorted(ambig_hits)[:3])
        return 'amber', (
            f'PDB-standard sugar name(s) detected ({sample}) with no '
            f'FF-specific markers — cannot auto-select. Run '
            f'`dvbfixer convert --to-amber` (→ GLYCAM UYB/VMB/...) or '
            f'`dvbfixer convert --to-charmm` (→ BGLCNA/BMAN/...) before '
            f'this step, or pass --ff explicitly'
        )

    return 'amber', 'no non-standard residues detected'


def _looks_like_xml_path(item):
    """True if `item` looks like an OpenMM FF XML path, not a short-name alias."""
    return item.endswith('.xml') or '/' in item or '\\' in item


def resolve_ff(user_ff, pdb_path, *, verbose=False):
    """Resolve a user's --ff argument into (xml_list, alias_name, reason).

    Accepts:
      - None or 'auto' or ['auto']: run detect_ff_from_pdb, expand alias.
      - A single short name in FF_ALIASES: expand; also run auto-detect to
        emit an *upgrade* if the input clearly needs a different FF (user
        said 'amber' but input has GLYCAM markers → upgrade to
        'amber+glycam' with a log line explaining why).
      - A list containing any '.xml' path: pass through unchanged
        (backward compat with the old `--ff a.xml b.xml` UX).

    Returns (xml_list, alias_name, upgrade_reason_or_None). The caller
    typically prints one line at startup:

        FF: {alias_name}  ({upgrade_reason} if upgrade_reason else "")
          → {" ".join(xml_list)}
    """
    # Normalise to a list
    if user_ff is None:
        items = ['auto']
    elif isinstance(user_ff, str):
        items = [user_ff]
    else:
        items = list(user_ff)

    # Explicit XML list? Pass through.
    if items and any(_looks_like_xml_path(x) for x in items):
        return items, 'custom', None

    if len(items) != 1:
        # Multiple short names — take the first, ignore the rest with a warning.
        if verbose:
            print(f"WARN: --ff got {len(items)} short-names; using first "
                  f"({items[0]}), ignoring {items[1:]}")
        items = items[:1]

    tok = items[0].lower().strip() if items else 'auto'

    if tok in ('auto', ''):
        alias, reason = detect_ff_from_pdb(pdb_path)
        if alias not in FF_ALIASES:
            alias = 'amber'
        return FF_ALIASES[alias], alias, reason

    if tok not in FF_ALIASES:
        raise ValueError(
            f"Unknown --ff short-name '{items[0]}'. "
            f"Valid: {', '.join(sorted(FF_ALIASES))} or pass explicit "
            f"XML paths (see docs/force-fields.md)."
        )

    # User asked for a specific short name; check if input clearly warrants
    # an upgrade (e.g. user said 'amber' but input has GLYCAM markers).
    detected, detect_reason = detect_ff_from_pdb(pdb_path)
    if detected != tok and detected in ('charmm', 'amber+glycam'):
        # Only auto-upgrade when detected FF is more specific than the plain
        # 'amber' family. Don't downgrade user's explicit CHARMM/GLYCAM choice.
        if tok in ('amber', 'amber19', 'amber14'):
            reason = (f"upgraded from '{tok}' → '{detected}' because "
                      f"{detect_reason}")
            return FF_ALIASES[detected], detected, reason

    return FF_ALIASES[tok], tok, None


def print_ff_selection(alias, reason, xml_list, prefix=""):
    """Emit the standard two-line FF selection banner used by every tool."""
    tag = f"{prefix}FF: {alias}"
    if reason:
        tag += f"  ({reason})"
    print(tag)
    print(f"{prefix}  → {' '.join(xml_list)}")


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
      - pdb_sugars:      set of (chain_id, res_id) for standard PDB sugar
                          names (NAG, BMA, MAN, FUC, ... — see _PDB_SUGAR_NAMES)
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
        elif name in _PDB_SUGAR_NAMES:
            info['pdb_sugars'].add(key)
        elif name not in known_prot_solv:
            info['unknown_hets'].add(key)
    return info


def sanitize_protein_hetatm(pdb_path, verbose=False):
    """Rewrite `pdb_path` so protein/GLYCAM-glycoprotein residues are
    guaranteed to be `ATOM` and any spurious mid-chain `TER` records are
    dropped. Both issues break OpenMM's peptide-bond inference.

    Two fixes applied (both no-ops on clean inputs):

    1. **HETATM → ATOM** for any residue name in `FORCE_ATOM_RESIDUES`
       (standard AAs + AMBER protonation variants HID/HIE/HIP/ASH/GLH/
       CYX/CYM/LYN + GLYCAM glycoprotein residues NLN/OLS/OLT). OpenMM's
       PDBFile parser only infers peptide bonds between `ATOM` records;
       a lone `HETATM ASN` sits in isolation with no bond to its
       neighbours, and downstream `addHydrogens` fails with the confusing
       "missing 1 C atom externally bonded" template error.

    2. **Spurious `TER` records** between two protein residues on the
       same chain. A TER forces OpenMM to start a new chain, breaking
       the polymer.

    Returns a temp-file path if any rewrite happened, or the original
    `pdb_path` unchanged otherwise.
    """
    import tempfile as _tf

    with open(pdb_path) as f:
        lines = f.readlines()

    res_at_pos = {}
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')) and len(line) >= 27:
            ch = line[21]
            try:
                rs = int(line[22:26])
            except ValueError:
                continue
            ic = line[26] if len(line) > 26 else ' '
            rn = line[17:20].strip()
            res_at_pos.setdefault((ch, rs, ic), rn)

    n_hetatm_fix = 0
    n_ter_drop = 0
    out_lines = []
    last_res_key = None
    pending_ter = None

    def _flush_pending(flush_list):
        nonlocal pending_ter
        if pending_ter is not None:
            flush_list.append(pending_ter)
            pending_ter = None

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')) and len(line) >= 27:
            ch = line[21]
            try:
                rs = int(line[22:26])
                ic = line[26] if len(line) > 26 else ' '
            except ValueError:
                _flush_pending(out_lines)
                out_lines.append(line)
                continue
            rn = line[17:20].strip()

            if line.startswith('HETATM') and rn in FORCE_ATOM_RESIDUES:
                line = 'ATOM  ' + line[6:]
                n_hetatm_fix += 1

            if pending_ter is not None:
                prev_rn = res_at_pos.get(last_res_key) if last_res_key else None
                same_chain = last_res_key and last_res_key[0] == ch
                both_protein = (
                    prev_rn in FORCE_ATOM_RESIDUES
                    and rn in FORCE_ATOM_RESIDUES
                )
                if same_chain and both_protein:
                    n_ter_drop += 1
                    pending_ter = None
                else:
                    out_lines.append(pending_ter)
                    pending_ter = None

            out_lines.append(line)
            last_res_key = (ch, rs, ic)

        elif line.startswith('TER'):
            _flush_pending(out_lines)
            pending_ter = line

        else:
            _flush_pending(out_lines)
            out_lines.append(line)

    _flush_pending(out_lines)

    if n_hetatm_fix == 0 and n_ter_drop == 0:
        return str(pdb_path)

    tmp_path = _tf.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False).name
    with open(tmp_path, 'w') as f:
        f.writelines(out_lines)
    if verbose:
        if n_hetatm_fix:
            print(f"  [sanitize] rewrote {n_hetatm_fix} HETATM → ATOM lines "
                  f"for protein/GLYCAM glycoprotein residues")
        if n_ter_drop:
            print(f"  [sanitize] dropped {n_ter_drop} spurious TER record(s) "
                  f"between same-chain protein residues")
    return tmp_path


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


def create_forcefield_with_openff(ff_xmls, topology,
                                  extra_generators=None, verbose=False,
                                  **_legacy_kwargs):
    """Build an OpenMM ForceField with GLYCAM template suppression.

    When `GLYCAM_06j-1.xml` is loaded AND the input topology carries
    PDB-standard sugar names (NAG / BMA / MAN / …), the ~1400 GLYCAM sugar
    and nucleic-acid templates fuzzy-match the wrong residues (NAG → UVA,
    …). This helper drops those GLYCAM templates so the PDB-named sugars
    fall through to whatever handler you want (a `--parametrize-ligands`
    generator, an explicit XML, etc.). GLYCAM's glycoprotein and cap
    templates (NLN/OLS/OLT/ROH/OME/TBT/CMET) are kept.

    Note: dvbfixer used to auto-register a SMIRNOFF generator here for
    unknown ligands. That path was removed — SMIRNOFF doesn't handle
    cross-residue bonds (glycosidic / protein-glycan linkages have no
    parameters, geometry blows up on minimize). Use `--parametrize-ligands`
    for real per-ligand GAFF2 + AM1-BCC via antechamber. See
    docs/force-fields.md.

    Args:
        ff_xmls: OpenMM FF XML paths (from ffutils.resolve_ff).
        topology: OpenMM Topology (scanned for PDB-sugar names).
        extra_generators: Optional iterable of already-built template
            generators to register on the returned ForceField (used by
            `--parametrize-ligands` to plug in GAFF2 templates).
        verbose: Extra logging.
        _legacy_kwargs: Silently swallowed (was `small_mol_ff`,
            `extra_molecules`); kept for callers that hadn't been updated.
    """
    ff = ForceField(*ff_xmls)

    # Suppress GLYCAM sugar/NA templates when PDB-named sugars are present.
    pdb_sugars = {r.name for r in topology.residues()
                  if r.name in _PDB_SUGAR_NAMES}
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
                      f"(PDB sugars detected)")

    for gen in (extra_generators or ()):
        # openmmforcefields' *TemplateGenerator objects expose a bound
        # `.generator` method; OpenMM's registerTemplateGenerator wants
        # that method, not the object itself. Callables are passed through
        # unchanged so callers can also hand in raw functions.
        hook = getattr(gen, 'generator', gen)
        ff.registerTemplateGenerator(hook)

    return ff
