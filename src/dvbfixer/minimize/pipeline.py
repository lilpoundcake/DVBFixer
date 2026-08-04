"""Energy-minimize a PDB structure with OpenMM using selective restraints.

Reads a .dat file (from 'dvbfixer prepare') to determine which atoms were
added by PDBFixer. Original atoms get strong positional restraints while
newly added atoms are free to relax into physically reasonable positions.

Three-tier restraint system:
  - Original heavy atoms: strong restraint (default 100 kcal/mol/A^2)
  - New backbone (N, CA, C, O, CB): weak restraint (default 5 kcal/mol/A^2)
  - New sidechain + all hydrogens: free (no restraint)

Minimization runs in two phases:
  1. Full restraints (default 1000 iterations)
  2. Restraints reduced 10x (default 1000 iterations)
"""

import sys
from pathlib import Path

from openmm import CustomExternalForce, LangevinMiddleIntegrator
from openmm.app import PME, ForceField, Modeller, PDBFile, Simulation
from openmm.unit import kelvin, nanometer, picosecond

from dvbfixer.ffutils.geometry import collect_ss_pairs as _collect_ss_pairs
from dvbfixer.ffutils.geometry import (
    drop_spurious_inter_aa_bonds as _drop_spurious_inter_aa_bonds,
)
from dvbfixer.ffutils.variants import AMBER_VARIANTS, scan_variant_names
from dvbfixer.minimize.cli import BACKBONE_NAMES, parse_args
from dvbfixer.minimize.refine import (
    _rigid_track_glycan_trees,
    refine_with_obminimize,
    refine_with_xtb,
)

# ---------------------------------------------------------------------------
# .dat file loading
# ---------------------------------------------------------------------------

def load_dat(path):
    """Load a ``.dat`` file, return set of (chain, resid, icode, atom_name).

    Delegates to :func:`dvbfixer.ffutils.dat.load_added_keys` — schema
    lives there.
    """
    from dvbfixer.ffutils.dat import load_added_keys
    return load_added_keys(path)


def resolve_new_atom_indices(topology, added_keys, verbose=False):
    """Match (chain, resid, icode, atom_name) keys from .dat against a topology.
    Returns set of atom indices that are 'new' (added by PDBFixer).
    """
    normalized_keys = {(c, r, ic.strip(), a) for c, r, ic, a in added_keys}

    indices = set()
    matched_keys = set()
    for atom in topology.atoms():
        res = atom.residue
        key = (res.chain.id, res.id, res.insertionCode.strip(), atom.name)
        if key in normalized_keys:
            indices.add(atom.index)
            matched_keys.add(key)

    unmatched = normalized_keys - matched_keys
    if unmatched and verbose:
        print(f"  {len(unmatched)} atoms from .dat not in topology "
              f"(hydrogens are re-added during minimization)")

    return indices


# ---------------------------------------------------------------------------
# Restraints & minimization
# ---------------------------------------------------------------------------

def build_restraint_force(topology, positions, new_atom_indices, strong_k, weak_k):
    """Build a CustomExternalForce with positional restraints.

    - Protein original heavy atoms: strong restraint to initial position
    - Newly added backbone atoms: weak restraint (keeps loop shape reasonable)
    - Newly added sidechain atoms: no restraint (free to move)
    - Heterogen heavy atoms (sugars, ligands): weak restraint (0.7.9+) — see
      below for why this is no longer "free"
    - All hydrogens: no restraint (free to move)

    Heterogens used to get ZERO restraint ("BioLuminate-style: protein is
    fixed-ish, ligands relax"), which is fine for a small, torsion-poor
    ligand but was found to let a multi-residue glycan tree (many free
    glycosidic torsions, no restoring force, worse under `--no-solvent`'s
    unscreened electrostatics) drift a mean of 4+ Å and up to 10+ Å from its
    covalent anchor during minimize, collapsing into unrelated, clashing
    parts of the protein surface. This matches established MD/structure-
    prep practice (CHARMM-GUI equilibration protocols, published
    glycoprotein MD setups): protein AND glycan heavy atoms are restrained
    together during initial minimization, commonly at ~1-10 kcal/mol/Å² —
    squarely in the range of this restraint scheme's existing `weak_k`
    tier, reused here rather than inventing a new one. Multiple studies
    found fully unrestrained minimization "almost always led to worse
    final structures" with "significant off-pathway sampling."
    """
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    known = PROTEIN_RESIDUES | SOLVENT_IONS

    # 1 kcal/mol/A^2 = 418.4 kJ/mol/nm^2
    conv = 418.4

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    n_strong = 0
    n_weak = 0
    n_free = 0

    for atom in topology.atoms():
        pos = positions[atom.index]
        x0 = pos[0].value_in_unit(nanometer)
        y0 = pos[1].value_in_unit(nanometer)
        z0 = pos[2].value_in_unit(nanometer)

        is_hydrogen = atom.element.symbol == 'H'
        is_new = atom.index in new_atom_indices
        is_heterogen = atom.residue.name not in known

        if is_hydrogen:
            # All H atoms are free — including newly-placed heterogen H,
            # which needs to relax into proper geometry.
            n_free += 1
            continue
        elif is_heterogen:
            k = weak_k * conv
            force.addParticle(atom.index, [k, x0, y0, z0])
            n_weak += 1
        elif not is_new:
            k = strong_k * conv
            force.addParticle(atom.index, [k, x0, y0, z0])
            n_strong += 1
        elif atom.name in BACKBONE_NAMES:
            k = weak_k * conv
            force.addParticle(atom.index, [k, x0, y0, z0])
            n_weak += 1
        else:
            n_free += 1

    return force, n_strong, n_weak, n_free


def _has_hydrogens(topology):
    """Check if topology already contains hydrogen atoms."""
    return any(a.element.symbol == 'H' for a in topology.atoms())


# Disulfide-bond snapshot/cleanup (_collect_ss_pairs / _drop_spurious_inter_aa_bonds)
# now lives in ffutils.geometry (imported above) — shared with protonate.py
# so the "snapshot true SS pairing before a rebuild/load, restore after"
# sequence can't drift out of sync across call sites again. It already had
# once: the `keep_h` branch below was missing it while the `rebuild_h`
# branch had it, and `protonate.py` had neither.


def _read_amber_renames(pdb_path):
    """Read PDB text to find AMBER protonation names before OpenMM normalizes them.

    Wraps the shared ``ffutils.variants.scan_variant_names`` scanner,
    keeping its ``(chain, resseq, icode)`` keys as-is — this module's
    callers used to project them down to a legacy ``(chain, resseq)``
    shape, which silently conflated two residues sharing a resSeq via
    insertion code (e.g. Kabat CDR-loop ``H:82``/``H:82A``). Filters out
    CHARMM-only variants (HSD/HSE/HSP/ASPP/GLUP/LSN) — minimize's
    downstream ``addHydrogens`` path only knows AMBER variant names.
    """
    saved = scan_variant_names(pdb_path)
    renames = {}
    for (chain, resseq, icode), name in saved.items():
        if name in AMBER_VARIANTS:
            renames[(chain, resseq, icode)] = name
    return renames


def _build_variants(topology, amber_renames):
    """Build variants list from AMBER renames dict.

    Returns list of variant names (or None) for each residue, suitable for
    Modeller.addHydrogens(variants=...).

    ``amber_renames`` keys are ``(chain, resid, icode)``. Matches the
    residue's real icode first; only falls back to an icode-less entry
    when the residue itself has no icode — an icode-bearing residue must
    never inherit a same-resSeq sibling's variant.
    """
    if not amber_renames:
        return None
    variants = []
    has_any = False
    for res in topology.residues():
        ic = res.insertionCode.strip() if hasattr(res, 'insertionCode') else ''
        key = (res.chain.id, res.id, ic)
        if key not in amber_renames and not ic:
            key = (res.chain.id, res.id, '')
        if key in amber_renames:
            variants.append(amber_renames[key])
            has_any = True
        else:
            variants.append(None)
    return variants if has_any else None


def _get_known_residues():
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    return PROTEIN_RESIDUES | SOLVENT_IONS


def _strip_hetatm(topology, positions):
    """Remove non-protein/solvent/ion residues. Returns (new_top, new_pos).

    Also renames GLYCAM glycoprotein residues (NLN/OLS/OLT) back to their
    standard parents (ASN/SER/THR) — they have no template in the standard
    AMBER19 protein FF. The protein-glycan bond would be missing afterward
    but the residue stays mid-chain; addHydrogens/PDBFixer will add the
    missing HD22/HG/HG1 to make the standard template match.
    """
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    _GLYCAM_BACK = {'NLN': 'ASN', 'OLS': 'SER', 'OLT': 'THR'}
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    to_delete = [res for res in topology.residues() if res.name not in known]

    # Rename NLN/OLS/OLT in topology so the standard AMBER template matches.
    for res in topology.residues():
        if res.name in _GLYCAM_BACK:
            res.name = _GLYCAM_BACK[res.name]

    if not to_delete:
        return topology, positions

    modeller = Modeller(topology, positions)
    modeller.delete(to_delete)
    return modeller.topology, modeller.positions


def minimize(topology, positions, new_atom_indices, args, amber_renames=None):
    """Set up system with solvent, apply selective restraints, minimize.

    Two modes:
    - keep_heterogens (default): minimize the whole system (protein + sugars +
      ligands) using the resolved --ff (AMBER + GLYCAM for glycoproteins,
      CHARMM36 for CHARMM inputs, plus any per-ligand GAFF2 templates
      registered via --parametrize-ligands). Heterogens not in .dat get no
      restraint (free to relax).
    - strip-heterogens (legacy): strip non-protein/non-solvent residues
      before parametrization and restore with original coordinates
      afterward.
    """
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    if args.keep_heterogens:
        stripped_top = topology
        stripped_pos = positions
        n_stripped = 0
        stripped_new_indices = set(new_atom_indices)
    else:
        # Legacy: strip non-protein HETATM
        stripped_top, stripped_pos = _strip_hetatm(topology, positions)
        n_stripped = sum(1 for _ in topology.residues()) - sum(1 for _ in stripped_top.residues())
        if n_stripped > 0:
            print(f"Stripped {n_stripped} non-protein residues for minimization")
            # The stripped heterogens will be spliced back at their ORIGINAL
            # coords after the protein-only minimize. Warn about the
            # inevitable interface strain and point at the higher-fidelity
            # path (--parametrize-ligands) so users don't take the artifact
            # for a real signal.
            print(
                "WARNING: heterogens will be restored at their INPUT "
                "coordinates while the protein has moved during minimization. "
                "Protein-ligand interface geometry (H-bonds, contacts) may "
                "be strained. For proper interface geometry re-run with "
                "`dvbfixer minimize --parametrize-ligands` instead (real "
                "GAFF2 templates for the ligand; whole system optimised "
                "together)."
            )
        # Remap new_atom_indices to stripped topology
        known = PROTEIN_RESIDUES | SOLVENT_IONS
        old_to_new = {}
        new_idx = 0
        for atom in topology.atoms():
            if atom.residue.name in known:
                old_to_new[atom.index] = new_idx
                new_idx += 1
        stripped_new_indices = {old_to_new[i] for i in new_atom_indices if i in old_to_new}

    # Capture the identity of every "new" atom as a (chain, resid, icode,
    # atom_name) key — NOT just its integer index — before any further
    # rebuild happens. `Modeller.addHydrogens()`/PDBFixer's
    # `addMissingAtoms()` (both called further down, possibly more than
    # once) each build an entirely new Topology, shifting every atom's
    # index whenever one is inserted. Re-resolving these keys against the
    # FINAL topology right before `build_restraint_force` (via the same
    # `resolve_new_atom_indices` helper the `.dat`-driven path already
    # uses) keeps the restraint tier assignment correct instead of
    # silently applying strong/weak restraints to the wrong atoms.
    new_atom_keys = {
        (a.residue.chain.id, a.residue.id, a.residue.insertionCode.strip(), a.name)
        for a in stripped_top.atoms() if a.index in stripped_new_indices
    }

    print("Loading force field...")
    # Check if there are heterogens — if so and keep_heterogens, use glycan-aware FF
    has_heterogens = args.keep_heterogens and any(
        res.name not in (PROTEIN_RESIDUES | SOLVENT_IONS)
        for res in stripped_top.residues()
    )

    # Build GAFF2 templates for any unknown-to-FF ligand and register them
    # on the ForceField below via extra_generators. `build_ligand_generator`
    # is a cheap no-op (returns None immediately) when every heterogen
    # already has a template, so this always runs rather than only under
    # --parametrize-ligands: without it, an unparametrizable ligand (e.g.
    # DAN when the user didn't think to pass the flag) silently fails
    # whole-system createSystem, triggers the legacy strip-and-splice
    # fallback below, and comes back at its ORIGINAL pre-minimize
    # coordinates with no tracking mechanism (unlike covalently-bonded
    # glycans, which _rigid_track_glycan_trees can re-anchor) — producing
    # exactly the clash reported against 1VCU's DANA. Explicit
    # --parametrize-ligands keeps strict=True (any failure raises
    # LigandParamError — the user asked for it, so don't silently
    # degrade); the auto path uses strict=False so environments without
    # AmberTools/openff still fall through to today's strip-and-splice
    # behaviour unchanged, just without a proactive attempt.
    if has_heterogens:
        from openmm.app import ForceField as _FF
        base_templates = set(_FF(*args.ff)._templates)
        from dvbfixer.lig_params import build_ligand_generator
        _explicit_parametrize = getattr(args, 'parametrize_ligands', False)
        _lig_gen = build_ligand_generator(
            stripped_top, stripped_pos, base_templates,
            strict=_explicit_parametrize, verbose=args.verbose,
        )
        args._ligand_generators = [_lig_gen] if _lig_gen else None

    if has_heterogens:
        from dvbfixer.ffutils import create_forcefield_with_openff
        # Loads the resolved --ff XMLs and prunes GLYCAM sugar/NA templates
        # that would fuzzy-match PDB-named sugars to the wrong entry. Any
        # GAFF2 ligand templates from --parametrize-ligands are appended
        # via extra_generators (see below).
        forcefield = create_forcefield_with_openff(
            args.ff, stripped_top,
            extra_generators=getattr(args, '_ligand_generators', None),
            verbose=args.verbose,
        )
        # Load glycam-hydrogens.xml so Modeller.addHydrogens (called below
        # to re-place H on variant residues) can also (re)place H on the
        # GLYCAM sugar + NLN/OLS/OLT residues. Without this, addHydrogens'
        # internal createSystem fails on every sugar with "missing N H
        # atoms" because the default hydrogens.xml has no GLYCAM entries.
        # The build_glycam_system wrapper does this for prepare/protonate;
        # minimize builds the FF directly so must load it here. Safe to
        # call repeatedly (no-op if already loaded).
        try:
            Modeller.loadHydrogenDefinitions("glycam-hydrogens.xml")
        except Exception:
            pass
    else:
        forcefield = ForceField(*args.ff)

    modeller = Modeller(stripped_top, stripped_pos)

    # Filter spurious inter-residue bonds BEFORE any addHydrogens call.
    # OpenMM's addHydrogens internally invokes createSystem, so we must
    # drop the "1 C-O bond too many"-style edges from CONECT inference
    # here to avoid template mismatches inside addHydrogens itself.
    _dropped_early = _drop_spurious_inter_aa_bonds(modeller.topology,
                                                    verbose=args.verbose)
    if _dropped_early:
        print(f"  Dropped {_dropped_early} spurious inter-residue bond(s) "
              f"before H placement")

    # Detect SS bonds from CONECT in input PDB and force CYX template for
    # those residues. Without this, GLYCAM FF triggers CYS/CYX/CYM ambiguity
    # since CYS in SS bonds has no HG (looks like CYM).
    ss_cys = set()
    if has_heterogens:
        try:
            from dvbfixer.acpype_export import detect_ss_bonds
            # We need the input PDB path for CONECT inspection.
            ss_cys = detect_ss_bonds(args.input)
            if ss_cys and args.verbose:
                print(f"  Detected {len(ss_cys)} CYS in SS bonds → CYX")
        except Exception as e:
            if args.verbose:
                print(f"  SS detection skipped: {e}")
    for res in modeller.topology.residues():
        if res.name == 'CYS' and (res.chain.id, int(res.id)) in ss_cys:
            res.name = 'CYX'

    # Rename all HIS to explicit variants before any OpenMM operation.
    # OpenMM's template matching and auto-detection crash on ambiguous HIS.
    for res in modeller.topology.residues():
        if res.name == 'HIS':
            atom_names = {a.name for a in res.atoms()}
            if 'HD1' in atom_names and 'HE2' in atom_names:
                res.name = 'HIP'
            elif 'HD1' in atom_names:
                res.name = 'HID'
            elif 'HE2' in atom_names:
                res.name = 'HIE'
            elif 'ND1' in atom_names and 'NE2' not in atom_names:
                res.name = 'HID'
            elif 'NE2' in atom_names and 'ND1' not in atom_names:
                res.name = 'HIE'
            else:
                res.name = 'HIE'  # default

    # OpenMM's Modeller.addHydrogens looks up `res.name` in its internal
    # `_residueHydrogens` dict — keyed by the standard residue name ("LYS",
    # "HIS", "CYS"), not by AMBER variant names ("LYN", "HID", "HIE", ...).
    # A topology residue with `res.name = "LYN"` is silently SKIPPED by
    # addHydrogens — no H atoms get added — and the next createSystem then
    # fails with "No template found for residue (LYN). ... missing N H atoms".
    # Rename variant residues to their parent names before addHydrogens and
    # restore the original names afterwards.
    _VARIANT_TO_PARENT = {
        'LYN': 'LYS', 'CYX': 'CYS', 'CYM': 'CYS',
        'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS',
        'ASH': 'ASP', 'GLH': 'GLU',
    }

    def _rename_variants_to_parent(top, amber_renames=None):
        """Rename variant residues to parent so addHydrogens recognises them.
        Returns dict {(chain_id, res_id, icode): original_name} for the
        post-pass — residue refs in the OLD topology are stale after
        addHydrogens (which rebuilds the topology entirely).

        Ground truth for the variant name is ``amber_renames`` (populated
        by ``_read_amber_renames`` from raw PDB text before any OpenMM
        parsing). Falling back to the current topology name misses
        HID/HIE/HIP/ASH/GLH/CYX/CYM because OpenMM's intermediate
        ``PDBFile.writeFile`` normalizes them to canonical HIS/ASP/GLU/CYS
        — only LYN survives that normalization (not in OpenMM's standard
        residue map). That asymmetry is why zbs pre-0.6.3 preserved LYN
        but canonicalized every other variant.

        ``amber_renames`` keys are ``(chain, resid, icode)`` — matched on
        the residue's real icode first, falling back to an icode-less
        entry only when the residue itself has none, so two residues
        sharing a resSeq via insertion code never share a variant.
        """
        saved = {}
        for res in top.residues():
            ic = res.insertionCode.strip() if hasattr(res, 'insertionCode') else ''
            key = (res.chain.id, res.id, ic)
            # Prefer the raw-text-derived variant name; fall back to the
            # topology name (only reliable for LYN).
            variant = None
            lookup_key = key
            if amber_renames and lookup_key not in amber_renames and not ic:
                lookup_key = (res.chain.id, res.id, '')
            if amber_renames and lookup_key in amber_renames:
                candidate = amber_renames[lookup_key]
                if candidate in _VARIANT_TO_PARENT:
                    variant = candidate
            elif res.name in _VARIANT_TO_PARENT:
                variant = res.name
            if variant:
                saved[key] = variant
                res.name = _VARIANT_TO_PARENT[variant]
        return saved

    def _restore_variants_in_topology(top, saved):
        for res in top.residues():
            ic = res.insertionCode.strip() if hasattr(res, 'insertionCode') else ''
            key = (res.chain.id, res.id, ic)
            if key in saved:
                res.name = saved[key]

    def _fix_lyn_hz_naming(top, saved):
        """OpenMM's hydrogens.xml gates HZ3 by variant="LYS", so addHydrogens
        with variant=LYN produces HZ1+HZ2 — but the AMBER LYN template expects
        HZ2+HZ3. Rename HZ1 → HZ3 on every LYN residue (chemically equivalent).
        Uses `saved` to identify which residues started life as LYN, since
        `res.name` may still be "LYS" at call time.
        """
        for res in top.residues():
            key = (res.chain.id, res.id)
            if saved.get(key) != 'LYN':
                continue
            for atom in res.atoms():
                if atom.name == 'HZ1':
                    atom.name = 'HZ3'
                    break

    # Default: keep existing hydrogens. --rebuild-h strips and re-adds via OpenMM.
    #
    # When the input has AMBER protonation variants (LYN / ASH / GLH / CYX /
    # CYM / HID / HIE / HIP), the topology loaded through PDBFile has
    # already been text-renamed to parent names above so OpenMM inferred
    # proper intra-residue bonds. But the ATOM SET may not match either
    # the parent or variant template (e.g. LYN input is missing HZ1,
    # ASH input has an extra HD2 that isn't in the ASP template).
    # `Modeller.addHydrogens` with ``variants=[...]`` rebuilds each H set
    # per variant template — but only if we strip existing H first, since
    # it never REMOVES extra H atoms. Force the strip-and-rebuild path
    # whenever variants are present.
    keep_h = not args.rebuild_h and not amber_renames
    if amber_renames and not args.rebuild_h:
        print(f"  {len(amber_renames)} AMBER variant residue(s) present — "
              f"stripping H and re-adding via variant templates")
    if keep_h:
        # Always run PDBFixer to detect and fix missing atoms (heavy + H).
        # Even in keep_h mode, mutated residues may have incomplete sidechains
        # (e.g. LYS from VAL mutation missing CE + H atoms).
        import tempfile as _tf

        from pdbfixer import PDBFixer as _PDBFixer
        with _tf.NamedTemporaryFile(mode='w', suffix='.pdb',
                                    delete=False) as _tmp:
            PDBFile.writeFile(modeller.topology, modeller.positions,
                              _tmp, keepIds=True)
            _tmp_path = _tmp.name
        try:
            with open(_tmp_path) as _f:
                _fixer = _PDBFixer(pdbfile=_f)
        finally:
            Path(_tmp_path).unlink()
        _fixer.findMissingResidues()
        _fixer.findMissingAtoms()
        n_missing = sum(len(v) for v in _fixer.missingAtoms.values())
        n_terminals = sum(len(v) for v in _fixer.missingTerminals.values())
        if n_missing or n_terminals:
            print(f"Fixing {n_missing} missing atom(s), {n_terminals} terminal(s)...")
            # Snapshot the true disulfide pairing BEFORE the rebuild — PDBFixer's
            # addMissingAtoms() re-derives SG-SG bonds by a pure distance scan
            # (Topology.createDisulfideBonds(), no 1:1 matching) and can pick up
            # spurious extra pairs on tightly-packed CYS clusters, crashing
            # createSystem with "N S atom(s) too many". The sibling `rebuild_h`
            # branch below already does this; this branch previously didn't.
            _true_ss_pairs = _collect_ss_pairs(_fixer.topology)
            from dvbfixer.ffutils.geometry import rebuild_missing_atoms_with_retry as _rebuild
            _rebuild(_fixer, verbose=args.verbose, log_prefix="[minimize] ")
            from dvbfixer.ffutils.geometry import fix_ca_chirality as _fix_chir
            _fix_chir(_fixer.topology, _fixer.positions, verbose=args.verbose)
            modeller = Modeller(_fixer.topology, _fixer.positions)
            # PDBFixer rebuilds bonds — re-drop spurious bonds before addHydrogens.
            _drop_spurious_inter_aa_bonds(modeller.topology, verbose=args.verbose,
                                           valid_ss_pairs=_true_ss_pairs)
            # PDBFixer's addMissingHydrogens ignores AMBER variant names
            # (its _describeVariant only recognises standard PDB names) and
            # canonicalises HIE/HID/HIP/ASH/GLH/CYX/CYM back to HIS/ASP/
            # GLU/CYS before placing H — violates the CLAUDE.md hard rule.
            # Follow the same variant-aware pattern as the branches below:
            # rename → Modeller.addHydrogens(variants=...) → restore.
            variants = _build_variants(modeller.topology, amber_renames)
            _saved = _rename_variants_to_parent(modeller.topology, amber_renames)
            try:
                if variants:
                    modeller.addHydrogens(forcefield, pH=min(args.ph, 9.99), variants=variants)
                else:
                    modeller.addHydrogens(forcefield, pH=min(args.ph, 9.99))
                from dvbfixer.ffutils.geometry import repair_misplaced_hydrogens
                repair_misplaced_hydrogens(
                    modeller.topology, modeller.positions, verbose=args.verbose,
                )
            finally:
                _fix_lyn_hz_naming(modeller.topology, _saved)
                _restore_variants_in_topology(modeller.topology, _saved)
        elif not _has_hydrogens(modeller.topology):
            print("No hydrogens found, adding them...")
            variants = _build_variants(modeller.topology, amber_renames)
            _saved = _rename_variants_to_parent(modeller.topology, amber_renames)
            try:
                if variants:
                    modeller.addHydrogens(forcefield, pH=min(args.ph, 9.99), variants=variants)
                else:
                    modeller.addHydrogens(forcefield, pH=min(args.ph, 9.99))
                from dvbfixer.ffutils.geometry import repair_misplaced_hydrogens
                repair_misplaced_hydrogens(
                    modeller.topology, modeller.positions, verbose=args.verbose,
                )
            finally:
                _fix_lyn_hz_naming(modeller.topology, _saved)
                _restore_variants_in_topology(modeller.topology, _saved)
        else:
            # Per-residue H check: PDBFixer's findMissingAtoms only counts
            # heavy atoms, so a protein residue with all its heavy atoms but
            # zero H (e.g. a user-edited LYN with no hydrogens; common after
            # text-level renaming tools like `dvbfixer glycam` or after a
            # mutation that didn't run addHydrogens) won't trigger the
            # branch above. createSystem then fails with a "template X does
            # not match residue Y" error. Detect this case and run
            # addHydrogens with variants — OpenMM only adds missing H, so
            # existing H on other residues are untouched.
            from dvbfixer.ffutils import PROTEIN_RESIDUES as _PR
            no_h_residues = []
            for res in modeller.topology.residues():
                if res.name not in _PR:
                    continue
                if res.name == 'PRO':
                    # Proline has no backbone H — skip
                    expected_min = 1   # at least HA / HB
                else:
                    expected_min = 1
                h_count = sum(1 for a in res.atoms() if a.element.symbol == 'H')
                heavy_count = sum(1 for a in res.atoms() if a.element.symbol != 'H')
                # A protein residue with heavy atoms but ZERO hydrogens is
                # clearly de-protonated by accident; addHydrogens will fix
                # only those (preserves H on already-protonated residues).
                if heavy_count > 1 and h_count < expected_min:
                    no_h_residues.append(res)
            if no_h_residues:
                print(f"  {len(no_h_residues)} protein residue(s) missing all "
                      f"H atoms (e.g. {no_h_residues[0].chain.id}:"
                      f"{no_h_residues[0].name}{no_h_residues[0].id}) — "
                      f"running addHydrogens to fill them in")
                variants = _build_variants(modeller.topology, amber_renames)
                _saved = _rename_variants_to_parent(modeller.topology, amber_renames)
                try:
                    if variants:
                        modeller.addHydrogens(forcefield, pH=min(args.ph, 9.99),
                                              variants=variants)
                    else:
                        modeller.addHydrogens(forcefield, pH=min(args.ph, 9.99))
                    from dvbfixer.ffutils.geometry import repair_misplaced_hydrogens
                    repair_misplaced_hydrogens(
                        modeller.topology, modeller.positions, verbose=args.verbose,
                    )
                finally:
                    _fix_lyn_hz_naming(modeller.topology, _saved)
                    _restore_variants_in_topology(modeller.topology, _saved)
            else:
                print("Keeping existing hydrogens from input")
    else:
        # Strip H, fix any missing heavy atoms via PDBFixer, then re-add correct H.
        # addHydrogens adds correct H for every PROTEIN residue, so stripping
        # those is safe — but scope the strip to protein residues only. This
        # whole branch exists to fix AMBER protonation VARIANTS (HIE/HID/
        # ASH/GLH/...), which is a protein-only concept; addHydrogens has no
        # hydrogens.xml entry for arbitrary heterogens (drug-like ligands
        # with no GLYCAM/GAFF2-aware placement rule) and stripping their H
        # here is unrecoverable — it silently produces a heavy-atom-only
        # residue that then fails createSystem/GAFF2 template matching
        # downstream (confirmed: DAN in 1VCU lost all 10 H this way even
        # after a GAFF2 template was built for it, since the strip ran
        # first). GLYCAM sugars are excluded too — their H is already
        # correct from `prepare`'s dedicated GLYCAM-aware H-addition pass;
        # re-running addHydrogens on them here would be unnecessary churn
        # for an unrelated (protein-side) protonation fix.
        from dvbfixer.ffutils import PROTEIN_RESIDUES as _PR_STRIP
        h_to_delete = [a for a in modeller.topology.atoms()
                       if a.element.symbol == 'H'
                       and a.residue.name in _PR_STRIP]
        if h_to_delete:
            modeller.delete(h_to_delete)

        # Use PDBFixer to detect and add any missing heavy atoms
        import tempfile as _tf

        from pdbfixer import PDBFixer as _PDBFixer
        with _tf.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as _tmp:
            PDBFile.writeFile(modeller.topology, modeller.positions, _tmp, keepIds=True)
            _tmp_path = _tmp.name
        try:
            with open(_tmp_path) as _f:
                _fixer = _PDBFixer(pdbfile=_f)
        finally:
            Path(_tmp_path).unlink()
        _fixer.findMissingResidues()
        _fixer.findMissingAtoms()
        n_missing = sum(len(v) for v in _fixer.missingAtoms.values())
        n_terminals = sum(len(v) for v in _fixer.missingTerminals.values())
        if n_missing or n_terminals:
            print(f"Fixing {n_missing} missing atom(s), {n_terminals} terminal(s)")
        # PDBFixer.addMissingAtoms() rebuilds the topology from scratch and
        # re-derives disulfides via OpenMM's Topology.createDisulfideBonds()
        # — a pure distance-cutoff scan with no 1:1 matching. On inputs
        # where multiple CYS SG atoms sit close together (e.g. two chain
        # copies whose N-termini pack near each other) this adds spurious
        # extra SG-SG bonds on top of the real ones, and createSystem then
        # fails with "N S atom(s) too many" for the CYX template. Capture
        # the correct pairing (already established upstream by prepare's
        # PROPKA+Reduce disulfide detection) before the rebuild so it can
        # be restored after.
        _true_ss_pairs = _collect_ss_pairs(modeller.topology)
        from dvbfixer.ffutils.geometry import rebuild_missing_atoms_with_retry as _rebuild
        _rebuild(_fixer, verbose=args.verbose, log_prefix="[minimize] ")
        from dvbfixer.ffutils.geometry import fix_ca_chirality as _fix_chir
        _fix_chir(_fixer.topology, _fixer.positions, verbose=args.verbose)
        modeller = Modeller(_fixer.topology, _fixer.positions)
        # PDBFixer rebuilds bonds — re-drop spurious bonds before addHydrogens.
        _drop_spurious_inter_aa_bonds(modeller.topology, verbose=args.verbose,
                                       valid_ss_pairs=_true_ss_pairs)

        # When keeping heterogens with GLYCAM, OpenMM's PDBFile doesn't infer
        # intra-residue bonds for GLYCAM residues (NLN/OLS/OLT + sugars).
        # add_glycam_bonds populates them from FF templates. MUST be called
        # BEFORE addHydrogens so the topology is correct when hydrogens are
        # placed. Also must be called BEFORE renaming GLYCAM residues, since
        # add_glycam_bonds expects the original GLYCAM names.
        if has_heterogens:
            try:
                from dvbfixer.acpype_export import add_glycam_bonds
                add_glycam_bonds(modeller.topology, forcefield, args.verbose,
                                  positions=modeller.positions)
            except Exception as e:
                if args.verbose:
                    print(f"  add_glycam_bonds skipped: {e}")

        print(f"Adding hydrogens (pH {args.ph})...")
        variants = _build_variants(modeller.topology, amber_renames)
        if variants:
            n_var = sum(1 for v in variants if v is not None)
            print(f"  {n_var} AMBER protonation variants detected")
        _saved = _rename_variants_to_parent(modeller.topology, amber_renames)
        try:
            try:
                if variants:
                    modeller.addHydrogens(forcefield, pH=min(args.ph, 9.99), variants=variants)
                else:
                    modeller.addHydrogens(forcefield, pH=min(args.ph, 9.99))
            except (ValueError, KeyError) as e:
                # Whole-topology addHydrogens can fail on a heterogen whose
                # forcefield template requires an external bond (e.g. a
                # glycosylated NLN/OLS/OLT) that geometry doesn't actually
                # support in this input (no sugar close enough to link to)
                # — a ValueError. It can also raise KeyError deep inside
                # openmmforcefields' GAFFTemplateGenerator: addHydrogens
                # builds its OWN temporary topology/residueTemplates
                # internally and invokes the registered GAFF2 generator
                # callback against that, which can hit an unregistered
                # GAFF atom type (confirmed: `KeyError: 'c'` for DAN,
                # even though DAN already has a working GAFF2 template for
                # the OUTER createSystem call this file makes directly —
                # addHydrogens's internal matching is a different, more
                # fragile code path in openmmforcefields/OpenMM). Mirrors
                # prepare's identical fallback (prepare/pipeline.py's
                # "falling back to protein-only"): retry without a
                # forcefield — addHydrogens then does plain geometric H
                # placement and skips createSystem/template matching
                # entirely, so a template mismatch elsewhere can't abort
                # the whole run.
                print(f"  WARNING: whole-topology addHydrogens failed "
                      f"({type(e).__name__}: {e}); falling back to "
                      f"protein-only.")
                if variants:
                    modeller.addHydrogens(pH=min(args.ph, 9.99), variants=variants)
                else:
                    modeller.addHydrogens(pH=min(args.ph, 9.99))
            from dvbfixer.ffutils.geometry import repair_misplaced_hydrogens
            repair_misplaced_hydrogens(
                modeller.topology, modeller.positions, verbose=args.verbose,
            )
        finally:
            _fix_lyn_hz_naming(modeller.topology, _saved)
            _restore_variants_in_topology(modeller.topology, _saved)

    # Unconditional Cα chirality safety net. Every path above may leave
    # residues in D configuration:
    # - keep_h Branch 1: rebuilt sidechains via PDBFixer templates
    # - keep_h Branch 2/3: no chirality fix at all
    # - rebuild_h: PDBFixer rebuild before H placement
    # - Input arriving from model/prepare/pull may already carry D
    # Idempotent: returns 0 when nothing is D.
    from dvbfixer.ffutils.geometry import (
        assert_all_l,
    )
    from dvbfixer.ffutils.geometry import (
        fix_ca_chirality as _fix_chir_final,
    )
    _chir_fixed = _fix_chir_final(modeller.topology, modeller.positions,
                                   verbose=args.verbose)
    if _chir_fixed and not args.verbose:
        print(f"  Repaired {_chir_fixed} Cα chirality inversion(s)")
    # Reflection is deterministic: any non-degenerate D-CA becomes L
    # after fix_ca_chirality, so assert_all_l here only fires when the
    # geometry is genuinely degenerate (N-CA-C near-colinear) or when
    # a new inversion was introduced by an earlier bug. Fail fast
    # rather than build restraints against unfixable coords.
    assert_all_l(modeller.topology, modeller.positions)

    # Drop any spurious inter-residue bond between two standard AAs that
    # isn't the canonical peptide C-N. CONECT inference or an upstream
    # over-eager bond guesser can add e.g. a HID sidechain C bonded to a
    # neighbouring residue's O — then residueTemplates rejects the
    # residue because "the set of atoms matches HID, but the residue has
    # 1 C-O bond too many". Same filter as pdbutils.inference applies
    # to CONECT records; we replay it at topology level for defence.
    _dropped = _drop_spurious_inter_aa_bonds(modeller.topology, verbose=args.verbose)
    if _dropped:
        print(f"  Dropped {_dropped} spurious inter-residue bond(s) "
              f"before createSystem")

    # Build residueTemplates for protein variants only (CYX, HIE/HID/HIP, LYN).
    # Avoids GLYCAM FF template ambiguity (e.g. CYS without HG matches both
    # CYM and CYX). Sugars are left to auto-match — forcing them is too risky
    # since PDB sugar names like BGL may exist in GLYCAM with different atoms.
    # Built before pre-solvent check so the check uses the same disambiguation
    # as the real createSystem (otherwise pre-check fails on CYX/CYM and
    # triggers an unnecessary fallback to legacy strip-and-splice).
    res_templates = {}
    if has_heterogens:
        n_term_keys = set()
        c_term_keys = set()
        for chain in modeller.topology.chains():
            rl = list(chain.residues())
            if rl:
                n_term_keys.add((chain.id, int(rl[0].id)))
                c_term_keys.add((chain.id, int(rl[-1].id)))
        for res in modeller.topology.residues():
            if res.name not in ('CYX', 'HIE', 'HID', 'HIP', 'LYN'):
                continue
            try:
                key = (res.chain.id, int(res.id))
            except (ValueError, TypeError):
                continue
            is_terminal = key in n_term_keys or key in c_term_keys
            if is_terminal:
                term_name = ('N' if key in n_term_keys else 'C') + res.name
                if term_name in forcefield._templates:
                    res_templates[res] = term_name
                    continue
            if res.name in forcefield._templates:
                res_templates[res] = res.name

    # If heterogens are present, pre-validate they can be parametrized
    # before addSolvent (which internally calls createSystem and would
    # propagate an uncaught NAG-template error). On failure, fall through
    # to the legacy strip-and-splice path before any solvent is placed.
    if has_heterogens and not args.no_solvent:
        try:
            from openmm.app import NoCutoff
            forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                                    ignoreExternalBonds=True,
                                    residueTemplates=res_templates)
        except (ValueError, KeyError) as e:
            # ValueError: OpenMM's normal "no template found"/"missing N
            # atoms" report. KeyError: its temporary re-matching can raise
            # this instead when a GAFF2 generator (registered above via
            # `_ligand_generators`) is involved — see CLAUDE.md's hard rule
            # on this exact `createSystem`/GAFF2-generator interaction.
            # Anything else is a real bug, not an expected template
            # mismatch — let it propagate instead of masquerading as one.
            from dvbfixer.ffutils import explain_template_error
            diag = explain_template_error(e, modeller.topology, forcefield)
            print(f"\nWARNING: pre-solvent parametrization failed:\n  {e}\n")
            if diag:
                for line in diag.split('\n'):
                    print(f"  {line}")
                print()
            print("  → falling back to --strip-heterogens flow.")
            print("    NOTE: heterogens will be restored at their INPUT "
                  "coords while the protein moves during minimize, so the "
                  "interface geometry may be strained. For proper interface "
                  "geometry re-run with `--parametrize-ligands` (real GAFF2 "
                  "templates for the ligand; whole system optimised "
                  "together).\n")
            import argparse as _ap
            legacy_args = _ap.Namespace(**vars(args))
            legacy_args.keep_heterogens = False
            legacy_args.rebuild_h = True
            out_top, out_pos = minimize(topology, positions, new_atom_indices,
                                         legacy_args, amber_renames=amber_renames)
            # Rigid-tracking inside the inner minimize already aligned the
            # glycan tree to the post-min anchor and snapped the linkage atom
            # to its ideal trigonal-planar position. A prep-coordinate
            # restore pass over the glycosidic amide here would overwrite
            # that post-min alignment with stale coordinates — don't add one.
            return out_top, out_pos

    if not args.no_solvent:
        print("Adding solvent...")
        modeller.addSolvent(forcefield, model='tip3p',
                            padding=args.padding * nanometer)

    print("Creating system...")
    from openmm.app import NoCutoff
    nbm = NoCutoff if args.no_solvent else PME
    try:
        if has_heterogens:
            system = forcefield.createSystem(modeller.topology, nonbondedMethod=nbm,
                                             ignoreExternalBonds=True,
                                             residueTemplates=res_templates)
        else:
            system = forcefield.createSystem(modeller.topology, nonbondedMethod=nbm)
    except (ValueError, KeyError) as e:
        # ValueError: OpenMM's normal "no template found"/"missing N
        # atoms" report. KeyError: its temporary re-matching can raise
        # this instead when a GAFF2 generator is involved — see
        # CLAUDE.md's hard rule on this exact createSystem/GAFF2-generator
        # interaction. `(ValueError, Exception)` used to be caught here,
        # which is exactly `except Exception` (ValueError already IS an
        # Exception) — any real, unrelated bug (a typo raising
        # AttributeError, say) silently masqueraded as an expected
        # template mismatch and got the fallback/diagnostic treatment
        # below instead of surfacing as a crash.
        if not has_heterogens:
            # Protein-only path — no auto-fallback available. Surface
            # `explain_template_error` so the user sees the specific
            # chain/resseq/atom-set (OpenMM's raw message uses topology
            # index which is useless for debugging), plus a hint about
            # the two most common causes on this path.
            from dvbfixer.ffutils import explain_template_error
            diag = explain_template_error(e, modeller.topology, forcefield)
            print(f"\nERROR: minimize createSystem failed:\n  {e}\n",
                  file=sys.stderr)
            if diag:
                for line in diag.split('\n'):
                    print(f"  {line}", file=sys.stderr)
            print(
                "\n  Common causes on the protein-only path:\n"
                "    - CONECT record (in the input PDB or added by\n"
                "      dvbfixer's auto-inference) linked an unrelated\n"
                "      atom into the residue's external bonds.\n"
                "    - An upstream tool changed the topology so a\n"
                "      standard AA now has an unexpected neighbour.\n"
                "  Re-run with --no-infer-conect to bypass the\n"
                "  CONECT-inference layer and confirm which side is\n"
                "  at fault.\n",
                file=sys.stderr)
            raise
        # Heterogen parametrization failed (e.g. PDB-named sugars without H).
        # Auto-fall back to legacy strip-and-restore flow with --rebuild-h
        # forced so addHydrogens fills in any missing H on renamed NLN→ASN etc.
        from dvbfixer.ffutils import explain_template_error
        diag = explain_template_error(e, modeller.topology, forcefield)
        print(f"\nWARNING: whole-system parametrization failed:\n  {e}\n")
        if diag:
            for line in diag.split('\n'):
                print(f"  {line}")
            print()
        print("  → falling back to --strip-heterogens flow "
              "(heterogens kept in output with original coords).")
        print("    NOTE: heterogens will be restored at their INPUT coords "
              "while the protein moves during minimize, so the interface "
              "geometry may be strained. For proper interface geometry "
              "re-run with `--parametrize-ligands` (real GAFF2 templates "
              "for the ligand; whole system optimised together).\n")
        import argparse as _ap
        legacy_args = _ap.Namespace(**vars(args))
        legacy_args.keep_heterogens = False
        legacy_args.rebuild_h = True
        out_top, out_pos = minimize(topology, positions, new_atom_indices, legacy_args,
                                     amber_renames=amber_renames)
        # Rigid-tracking inside the inner minimize already aligned the
        # glycan tree to the post-min anchor and snapped the linkage atom
        # to its ideal trigonal-planar position. A prep-coordinate restore
        # pass over the glycosidic amide (CG/OD1/ND2/HD21) here would
        # overwrite that with stale coordinates, breaking the bond to the
        # rigid-tracked glycan — don't add one.
        return out_top, out_pos

    # Re-resolve by identity against the FINAL topology — `stripped_new_indices`
    # (raw ints captured at function entry) no longer correspond to the same
    # atoms if any addHydrogens/PDBFixer rebuild happened in between.
    final_new_indices = resolve_new_atom_indices(
        modeller.topology, new_atom_keys, verbose=args.verbose,
    )
    restraint, n_strong, n_weak, n_free = build_restraint_force(
        modeller.topology, modeller.positions,
        final_new_indices, args.restraint_k, args.weak_k
    )
    system.addForce(restraint)
    print(f"Restraints: {n_strong} strong (original), {n_weak} weak (new backbone), {n_free} free")

    # Chirality is enforced *positionally* — reflect any D residue via
    # ``fix_ca_chirality`` (already ran above at the safety net and
    # asserted L via ``assert_all_l``) and re-verify in the iterative
    # phase-2 loop below. NO CustomTorsionForce is added: field
    # consensus (VMD/NAMD chirality plugin, pdb4amber, CHARMM-GUI) is
    # reflect-then-minimize-then-verify; a stiff runtime bias would
    # produce NaN forces on broken inputs and fight physics on the
    # rare cases where reflection genuinely needs a second pass.

    integrator = LangevinMiddleIntegrator(300 * kelvin, 1.0 / picosecond, 0.002 * picosecond)

    if args.platform:
        from openmm import Platform
        platform = Platform.getPlatformByName(args.platform)
        simulation = Simulation(modeller.topology, system, integrator, platform)
    else:
        simulation = Simulation(modeller.topology, system, integrator)

    simulation.context.setPositions(modeller.positions)

    state = simulation.context.getState(getEnergy=True)
    _e_pre = state.getPotentialEnergy()
    print(f"Energy before minimization: {_e_pre}")
    # Astronomical starting energy = severe atom clashes upstream
    # (e.g. PDBFixer / addHydrogens put an atom at another's position).
    # Warn early: no minimize configuration can recover from this and
    # any additional force term risks NaN via double-precision overflow.
    if _e_pre.value_in_unit(_e_pre.unit) > 1e10:
        print(f"WARNING: pre-minimize energy is astronomically high "
              f"({_e_pre}). The input has severe atom clashes "
              f"(likely coincident atoms from the prepare step). Run "
              f"`dvbfixer diagnose {args.input}` to identify the "
              f"offending residue(s). minimize will attempt to proceed "
              f"but may crash with an inf/NaN force at the first step.",
              file=sys.stderr)

    # Phase 1: full restraints
    print(f"Minimizing (phase 1: {args.max_iter} iterations, original atoms restrained)...")
    simulation.minimizeEnergy(maxIterations=args.max_iter)

    state = simulation.context.getState(getEnergy=True)
    print(f"Energy after phase 1: {state.getPotentialEnergy()}")

    # Phase 2: 10x reduced restraints
    print("Minimizing (phase 2: reduced restraints)...")
    for i in range(restraint.getNumParticles()):
        idx, params = restraint.getParticleParameters(i)
        restraint.setParticleParameters(i, idx, [params[0] / 10.0, params[1], params[2], params[3]])
    restraint.updateParametersInContext(simulation.context)

    # Phase 2 minimize with weakened restraints. Do NOT wrap this in an
    # iterative reflect-and-re-minimize loop: reflecting a sidechain
    # after minimize has equilibrated it swings CB across the CA-N-C
    # plane by ~2 Å, which the weakened phase-2 restraints cannot
    # recover from — resulting in CA-CB bonds stretched to 3 Å and
    # neighbour-residue steric clashes. Reflection is safe only *before*
    # OpenMM sees the geometry (prepare, and the pre-createSystem safety
    # net at line 655), where the sidechain hasn't been equilibrated
    # against its neighbours yet.
    simulation.minimizeEnergy(maxIterations=args.max_iter)

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"Energy after phase 2: {state.getPotentialEnergy()}")

    # Post-minimize chirality check + reflect-and-re-minimize loop.
    # After phase-2 minimize, some residues may have drifted into D-Cα
    # geometry (real FF energy minimum for that local packing). We
    # reflect any D-Cα back to L via fix_ca_chirality and run a short
    # follow-up minimize under the (already-weakened) phase-2 restraints
    # so the reflected CB relaxes into a compatible position without
    # swinging to 2 Å bond lengths. Bounded loop: at most a few
    # iterations to keep runtime under control.
    from dvbfixer.ffutils.geometry import find_d_residues, fix_ca_chirality
    min_topology = simulation.topology
    min_positions = state.getPositions()
    for _reflect_iter in range(3):
        offenders = find_d_residues(min_topology, min_positions)
        if not offenders:
            break
        print(f"  {len(offenders)} residue(s) drifted into D-Cα during "
              f"minimize (iter {_reflect_iter + 1}); reflecting and "
              f"re-minimizing.", file=sys.stderr)
        # Reflect in-place on a mutable list of Vec3.
        pos_list = list(min_positions)
        n_fixed = fix_ca_chirality(min_topology, pos_list,
                                    verbose=args.verbose)
        if not n_fixed:
            # Reflect returned 0 — nothing repairable (all near-degenerate).
            # Break to avoid an infinite loop.
            break
        simulation.context.setPositions(pos_list)
        simulation.minimizeEnergy(maxIterations=max(200, args.max_iter // 4))
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        min_positions = state.getPositions()
        print(f"  Energy after reflect+re-minimize iter {_reflect_iter + 1}: "
              f"{state.getPotentialEnergy()}", file=sys.stderr)
    else:
        # Loop exhausted without breaking → the FF's local minimum for
        # these residues genuinely sits on the D side (rare, exotic
        # packing). Chirality invariant is non-negotiable: reflect
        # unconditionally and accept the small packing strain. Any
        # downstream tool that cares can re-relax; we WILL NOT emit a D
        # residue.
        final_offenders = find_d_residues(min_topology, min_positions)
        if final_offenders:
            pos_list = list(min_positions)
            n_forced = fix_ca_chirality(min_topology, pos_list,
                                         verbose=args.verbose)
            print(f"WARNING: {len(final_offenders)} residue(s) preferred "
                  f"D-Cα geometry after 3 reflect+re-minimize iterations. "
                  f"Forcing L via unconditional reflection (may leave "
                  f"minor local strain — manual inspection recommended):",
                  file=sys.stderr)
            for ch, rid, name, tri in final_offenders[:20]:
                print(f"  {ch}/{name}{rid}: triple={tri:+.5f} nm³",
                      file=sys.stderr)
            if len(final_offenders) > 20:
                print(f"  ...and {len(final_offenders) - 20} more",
                      file=sys.stderr)
            if n_forced:
                # A plain follow-up minimize here would pull CB back
                # toward the FF's preferred D minimum, undoing the fix
                # — UNLESS the restraint anchors are moved first (see
                # below). Sidechain internal geometry is preserved by
                # fix_ca_chirality (reflects the whole sidechain, not
                # just CB), so bond lengths within the sidechain are
                # unchanged; only CB's position relative to backbone
                # neighbours changes.
                # Round-trip through the simulation context so
                # min_positions is a Quantity(list, unit) — the raw
                # pos_list is list[Quantity] which trips PDBFile's
                # np.isnan check downstream.
                simulation.context.setPositions(pos_list)
                min_positions = simulation.context.getState(
                    getPositions=True,
                ).getPositions()
                # Verify.
                still_d = find_d_residues(min_topology, min_positions)
                if still_d:
                    print(f"  {len(still_d)} residue(s) STILL D after "
                          f"forced reflection — geometry near-degenerate "
                          f"(fix_ca_chirality skipped them):",
                          file=sys.stderr)
                    for ch, rid, name, tri in still_d[:5]:
                        print(f"    {ch}/{name}{rid}: triple={tri:+.5f} nm³",
                              file=sys.stderr)

                # Rigid reflection swings the whole sidechain to its
                # mirror image, which can land a hydrogen inside a
                # neighbouring residue (hydrogens carry no restraint at
                # all in `restraint`, so they were never free to dodge
                # during the reflect itself). Resolve this with one
                # bounded LOCAL minimize — but first re-anchor every
                # restrained atom of the just-reflected residues to
                # their NEW (already-corrected) position. Backbone
                # atoms are untouched by fix_ca_chirality, so
                # re-anchoring them is a no-op; only the reflected
                # sidechain's anchor actually moves. Because the
                # anchor now points at the corrected geometry (not the
                # old D-favoring one), this minimize has no energetic
                # or restraint-driven path back to D — it can only let
                # the genuinely unrestrained neighbours (hydrogens)
                # relax out of the way.
                target_keys = {(ch, rid, name)
                               for ch, rid, name, _tri in final_offenders}
                target_atom_indices = {
                    atom.index
                    for res in min_topology.residues()
                    if (res.chain.id, str(res.id), res.name) in target_keys
                    for atom in res.atoms()
                }
                n_reanchored = 0
                for i in range(restraint.getNumParticles()):
                    atom_idx, params = restraint.getParticleParameters(i)
                    if atom_idx in target_atom_indices:
                        p = min_positions[atom_idx]
                        x0 = p[0].value_in_unit(nanometer)
                        y0 = p[1].value_in_unit(nanometer)
                        z0 = p[2].value_in_unit(nanometer)
                        restraint.setParticleParameters(
                            i, atom_idx, [params[0], x0, y0, z0])
                        n_reanchored += 1
                if n_reanchored:
                    restraint.updateParametersInContext(simulation.context)
                    simulation.minimizeEnergy(
                        maxIterations=max(200, args.max_iter // 4))
                    state = simulation.context.getState(
                        getEnergy=True, getPositions=True)
                    min_positions = state.getPositions()
                    print(f"  Local relax after forced reflection "
                          f"(anchors moved to corrected geometry): "
                          f"{state.getPotentialEnergy()}", file=sys.stderr)
                    # Sanity check only — the re-anchored restraint
                    # should make this impossible, but the chirality
                    # invariant is non-negotiable, so re-reflect rather
                    # than ever emit a D residue.
                    post_relax_d = find_d_residues(min_topology, min_positions)
                    if post_relax_d:
                        print(f"  {len(post_relax_d)} residue(s) reverted to "
                              f"D-Cα during post-reflect local relax; "
                              f"reflecting again unconditionally.",
                              file=sys.stderr)
                        pos_list = list(min_positions)
                        fix_ca_chirality(min_topology, pos_list,
                                          verbose=args.verbose)
                        simulation.context.setPositions(pos_list)
                        min_positions = simulation.context.getState(
                            getPositions=True,
                        ).getPositions()

    # Restore non-protein residues with original coordinates (legacy path only).
    # When keep_heterogens is on, n_stripped==0 and we write the minimized topology directly.
    if n_stripped > 0:
        # Strip solvent from minimized result first
        if not args.no_solvent:
            min_topology, min_positions = strip_solvent(min_topology, min_positions)

        # Merge: use minimized positions for protein atoms, original for HETATM
        # Match by (chain, resid, parent_resname, atomname) since
        # addHydrogens may have changed indices. `parent_resname` (resname
        # normalised through `_VARIANT_TO_PARENT`, so e.g. HIS/HIE/HID/HIP
        # all compare equal) is included — not just chain+resid — as
        # defense in depth: `min_topology` only ever contains protein/
        # solvent-ion residues (stripped by `_strip_hetatm` above), so this
        # never excludes a legitimate match (variant renames included) —
        # but a HETATM residue that ends up sharing a (chain, resid) with a
        # real protein residue (confirmed on a real structure: a
        # resSeq-collision bug in `model/renumber.py`'s HETATM placement,
        # now fixed at its source) must never have its own atoms
        # overwritten by that protein residue's minimized coordinates just
        # because an atom name like `CA`/`CB` happens to coincide.
        import numpy as np
        from openmm.unit import nanometer as nm_unit

        def _parent_resname(name):
            return _VARIANT_TO_PARENT.get(name, name)

        # Build position lookup from minimized protein
        min_pos_map = {}
        for atom in min_topology.atoms():
            res = atom.residue
            key = (res.chain.id, res.id, _parent_resname(res.name), atom.name)
            min_pos_map[key] = min_positions[atom.index].value_in_unit(nm_unit)

        # Build result: original positions, overwritten by minimized where available
        n_atoms = len(positions)
        result = np.zeros((n_atoms, 3))
        for i in range(n_atoms):
            result[i] = positions[i].value_in_unit(nm_unit)

        n_updated = 0
        for atom in topology.atoms():
            res = atom.residue
            key = (res.chain.id, res.id, _parent_resname(res.name), atom.name)
            if key in min_pos_map:
                result[atom.index] = min_pos_map[key]
                n_updated += 1

        if args.verbose:
            print(f"Restored: {n_updated} protein atoms updated, "
                  f"{n_atoms - n_updated} HETATM atoms kept original")

        # Rigid-body transform each glycan tree to follow its protein anchor.
        # Without this, the protein ASN/SER/THR moves during minimization but
        # the glycan stays at prep coordinates → ND2-C1 bond stretches and
        # clashes appear between glycan and moved protein. We compute the
        # Kabsch transform from prep→post-min anchor residue coords and apply
        # it to every atom in the bonded glycan tree.
        result = _rigid_track_glycan_trees(
            topology, positions, result, min_pos_map, verbose=args.verbose,
        )

        return topology, result * nm_unit

    return min_topology, min_positions


def strip_solvent(topology, positions):
    """Remove water and ions from the final structure."""
    modeller = Modeller(topology, positions)
    modeller.deleteWater()

    ions_to_delete = []
    for res in modeller.topology.residues():
        if res.name.upper() in ('NA', 'CL', 'NA+', 'CL-', 'K', 'K+', 'MG', 'MG2+', 'CA2+'):
            ions_to_delete.append(res)
    if ions_to_delete:
        modeller.delete(ions_to_delete)

    return modeller.topology, modeller.positions


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_minimized")

    if not args.no_infer_conect:
        from dvbfixer.pdbutils import _materialise_inferred_pdb
        input_path = Path(_materialise_inferred_pdb(
            input_path, verbose=args.verbose))

    # Resolve --ff short-name / auto-detect from input residue names.
    from dvbfixer.ffutils import (
        FF_ALIASES,
        has_pdb_standard_sugars,
        print_ff_selection,
        resolve_ff,
    )
    args.ff, _ff_alias, _ff_reason = resolve_ff(
        args.ff, args.input, verbose=args.verbose)

    # CHARMM + PDB-standard sugars: process under amber+glycam
    # (charmm36.xml has no sugar templates); rewrite sugar residue
    # names back to CHARMM convention post-run.
    _charmm_output_requested = False
    if _ff_alias == 'charmm' and has_pdb_standard_sugars(input_path):
        print("  [ff] charmm requested but input has PDB-standard "
              "sugars charmm36.xml can't parametrise; processing "
              "with amber+glycam and rewriting output residue names "
              "to CHARMM convention.")
        _charmm_output_requested = True
        _ff_alias = 'amber+glycam'
        args.ff = FF_ALIASES['amber+glycam']
        _ff_reason = (
            'auto-switched from charmm (PDB sugars need GLYCAM); '
            'output residue names will be rewritten back to CHARMM'
        )

    # amber+glycam + PDB-standard sugars: auto-convert to GLYCAM
    # canonical naming so templates match.
    if _ff_alias == 'amber+glycam' and has_pdb_standard_sugars(input_path):
        import atexit

        from dvbfixer.glycam import convert_to_glycam
        from dvbfixer.tempfiles import make_temp_path
        _conv_out = make_temp_path(suffix='.pdb', prefix='dvbfixer_glycam_')
        # Registered unconditionally (not only on success) — used as
        # `input_path` for the rest of this run, so it can't be deleted
        # right away; atexit is the same pattern
        # `pdbutils._materialise_inferred_pdb` already uses for this exact
        # "temp file becomes the working input" situation.
        atexit.register(lambda p=_conv_out: p.unlink(missing_ok=True))
        try:
            convert_to_glycam(str(input_path), str(_conv_out),
                              add_roh=True, verbose=args.verbose)
            print("  [ff] auto-converted PDB-standard sugar names → "
                  "GLYCAM canonical for amber+glycam FF matching.")
            input_path = Path(_conv_out)
        except Exception as e:
            print(f"  [ff] convert_to_glycam failed ({e}); continuing "
                  f"with original names — createSystem may fail on "
                  f"unrecognised sugar residues.")

    print_ff_selection(_ff_alias, _ff_reason, args.ff)

    # .dat is optional: if --dat given, use it; otherwise auto-detect next
    # to the USER'S input (args.input) — not the mutated `input_path`,
    # which by now points at a /tmp/dvbfixer_conect_*.pdb temp file that
    # obviously has no .dat sibling.
    if args.dat:
        dat_path = Path(args.dat)
    else:
        dat_path = Path(args.input).with_suffix(".dat")

    if args.rename:
        import atexit

        from dvbfixer.rename import canonicalize_pdb
        from dvbfixer.tempfiles import make_temp_path
        _tmp = make_temp_path(suffix='.pdb')
        n = canonicalize_pdb(input_path, _tmp, args.verbose)
        if n > 0:
            print(f"Canonicalized {n} non-canonical residue(s)")
            input_path = _tmp
            # Used as `input_path` for the rest of this run — can't
            # delete it right away; defer to process exit (previously
            # never cleaned up at all in this branch).
            atexit.register(lambda p=_tmp: p.unlink(missing_ok=True))
        elif _tmp.exists():
            _tmp.unlink()

    print(f"=== Minimization: {input_path} ===")

    # Read AMBER protonation names before OpenMM normalizes them
    amber_renames = _read_amber_renames(str(input_path))
    if amber_renames:
        print(f"Detected {len(amber_renames)} AMBER protonation variants in input")

    # Snapshot original GLYCAM residue names from the input PDB text. If the
    # legacy strip-and-splice fallback fires inside minimize(), it renames
    # NLN→ASN/OLS→SER/OLT→THR in place on the topology. Restore these on
    # the final output so the GLYCAM input names survive end-to-end.
    _glycam_orig_names = {}
    try:
        with open(input_path) as _gf:
            for _line in _gf:
                if not _line.startswith(('ATOM', 'HETATM')):
                    continue
                _rn = _line[17:20].strip()
                if _rn in ('NLN', 'OLS', 'OLT'):
                    _ch = _line[21]
                    _seq = _line[22:26].strip()
                    if _seq and _seq.lstrip('-').isdigit():
                        _ic = _line[26].strip() if len(_line) > 26 else ''
                        # Topology Residue.id is a string of the resSeq
                        _glycam_orig_names[(_ch, _seq, _ic)] = _rn
    except Exception:
        pass

    # Pre-rename AMBER/CHARMM variant residue names to their standard
    # parents in a temp copy BEFORE PDBFile parses the file. Rationale:
    # OpenMM's `PDBFile._standardResidues` set covers only the 20
    # canonical AAs, so variants (LYN / ASH / GLH / CYX / CYM / HID /
    # HIE / HIP) get loaded with zero intra-residue bonds — every
    # downstream template match then fails ("residue has no bonds
    # between its atoms"). Renaming to parent lets OpenMM infer proper
    # bonds; ``amber_renames`` (already captured above from the raw
    # text) drives the eventual restoration to variant names in the
    # topology and output PDB.
    from dvbfixer.ffutils.variants import text_rename_variants_to_parent
    _pdb_load_path, _pre_saved = text_rename_variants_to_parent(
        str(input_path), verbose=args.verbose,
    )
    try:
        pdb = PDBFile(_pdb_load_path)
    finally:
        if _pdb_load_path != str(input_path):
            Path(_pdb_load_path).unlink(missing_ok=True)
    topology = pdb.topology
    positions = pdb.positions

    # OpenMM's PDBFile.__init__ unconditionally calls
    # Topology.createDisulfideBonds() on every load (pdbfile.py), a pure
    # distance-cutoff scan with no 1:1 matching — it ignores CONECT
    # entirely for this, so even a perfectly clean input file's true
    # disulfide pairing gets discarded and re-derived, potentially with
    # spurious extra SG-SG bonds when several CYS cluster close together
    # (e.g. two chain copies whose N-termini sit near each other).
    # Resolve immediately after load via nearest-distance 1:1 matching so
    # the rest of the pipeline never sees the contaminated version.
    _drop_spurious_inter_aa_bonds(topology, verbose=args.verbose,
                                   positions=positions)

    # Load restraint data from .dat if available
    if args.dat and not dat_path.exists():
        print(f"Error: .dat file not found: {dat_path}", file=sys.stderr)
        sys.exit(1)

    if dat_path.exists():
        added_keys = load_dat(dat_path)
        new_atom_indices = resolve_new_atom_indices(topology, added_keys, args.verbose)

        # Also load variant overrides from .dat (saved by prepare) — use the
        # shared DatRecord to avoid re-parsing.
        from dvbfixer.ffutils.dat import DatRecord

        _dat = DatRecord.load(dat_path)
        for key_str, var_name in (_dat.variant_overrides or {}).items():
            # New format: "chain:resid:icode" (3 fields, icode may be
            # empty). Older .dat files on disk from before this fix used
            # "chain:resid" (2 fields, no icode) — still accepted.
            _parts = key_str.split(':', 2)
            ch, rn = _parts[0], _parts[1]
            ic = _parts[2] if len(_parts) > 2 else ''
            key = (ch, rn, ic)
            if key not in amber_renames:
                amber_renames[key] = var_name
        if amber_renames:
            print(f"  Total protonation variants (PDB + .dat): {len(amber_renames)}")
    else:
        print("No .dat file — all atoms get uniform strong restraints")
        new_atom_indices = set()

    final_topology, final_positions = minimize(
        topology, positions, new_atom_indices, args,
        amber_renames=amber_renames,
    )

    # Strip solvent from output. When keep_heterogens is on, the splice-back
    # inside minimize was skipped, so solvent is still in the topology and
    # must be stripped here. In legacy strip mode, solvent stripping happened
    # in the splice path only when n_stripped > 0; we redo it here as a
    # safety net for purely protein inputs.
    if not args.no_solvent:
        final_topology, final_positions = strip_solvent(final_topology, final_positions)

    # Optional refinement passes — auto-parametrize any heterogen via xtb
    # GFN-FF or OpenBabel MMFF94. These don't need user-provided FF params.
    if args.xtb_refine:
        final_positions = refine_with_xtb(
            final_topology, final_positions,
            cycles=args.xtb_cycles,
            heterogens_only=args.refine_heterogens_only,
            verbose=args.verbose,
        )
    if args.obminimize_refine:
        final_positions = refine_with_obminimize(
            final_topology, final_positions,
            ff=args.obminimize_ff,
            steps=args.obminimize_steps,
            heterogens_only=args.refine_heterogens_only,
            verbose=args.verbose,
        )

    # Restore GLYCAM glycoprotein residue names that the legacy strip path
    # may have renamed in-place (NLN→ASN, OLS→SER, OLT→THR). Keyed by the
    # (chain_id, res_id_str, icode) tuple captured from the raw input PDB.
    if _glycam_orig_names:
        _restored = 0
        for _res in final_topology.residues():
            _key = (_res.chain.id, str(_res.id),
                    _res.insertionCode.strip()
                    if hasattr(_res, 'insertionCode') else '')
            _orig = _glycam_orig_names.get(_key)
            if _orig and _res.name != _orig:
                _res.name = _orig
                _restored += 1
        if _restored and args.verbose:
            print(f"  Restored {_restored} GLYCAM residue names "
                  f"(NLN/OLS/OLT) renamed during legacy fallback")

    with open(output_path, 'w') as f:
        PDBFile.writeFile(final_topology, final_positions, f, keepIds=True)
    # Rewrite HETATM→ATOM for AMBER protonation variants (HID/HIE/HIP/ASH/
    # GLH/CYX/CYM/LYN) and GLYCAM glycoprotein residues (NLN/OLS/OLT) that
    # PDBFile.writeFile emits as HETATM.
    from dvbfixer.ffutils import fix_atom_hetatm_records
    fix_atom_hetatm_records(output_path)
    # OpenMM's PDBFile canonicalises HIE/HID/HIP → HIS, ASH → ASP,
    # GLH → GLU, CYX → CYS on write (empirically verified on 8.5.1).
    # amber_renames holds the ground truth from raw input text; apply
    # it as a text-level rewrite so the final output carries the
    # variant names the FF actually needs. Also folds in the
    # GROMACS amber99sb-ildn LYN HZ3→HZ1 shift.
    _target_ff = 'charmm' if any('charmm' in x for x in args.ff) else 'amber'
    from dvbfixer.ffutils.ff_names import apply_variants_to_pdb_text
    _n_var = apply_variants_to_pdb_text(
        output_path, amber_renames or {},
        target_ff=_target_ff, verbose=args.verbose,
        include_gromacs_shifts=(getattr(args, 'atom_naming', 'gromacs')
                                 == 'gromacs'),
    )
    if _n_var and args.verbose:
        print(f"  Applied {_n_var} variant name/atom rewrites to output")

    # CHARMM output was requested but the pipeline ran under
    # amber+glycam (PDB-standard sugars needed GLYCAM templates).
    # Rewrite sugar residue names GLYCAM → CHARMM in the output PDB
    # so the user gets what they asked for at the CLI. Coordinates
    # and topology unchanged; only sugar residue names change.
    if _charmm_output_requested:
        from dvbfixer.glycam import convert_to_charmm
        from dvbfixer.tempfiles import make_temp_path
        _tmp = make_temp_path(suffix='.pdb', prefix='dvbfixer_charmm_out_')
        try:
            convert_to_charmm(str(output_path), str(_tmp), verbose=args.verbose)
            import shutil as _sh
            _sh.move(str(_tmp), str(output_path))
            print("  [ff] rewrote sugar residue names GLYCAM → CHARMM "
                  "for output (user asked for --ff charmm).")
        except Exception as e:
            print(f"  [ff] convert_to_charmm on output failed ({e}); "
                  f"output has GLYCAM sugar names — run "
                  f"`dvbfixer convert --to-charmm` manually.")
        finally:
            # On success `shutil.move` already relocated this path (so
            # there's nothing left to remove); on failure a partial or
            # nonexistent file may remain — clean it up either way.
            Path(_tmp).unlink(missing_ok=True)

    print(f"\nSaved minimized structure: {output_path}")
