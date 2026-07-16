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

import argparse
import json
import sys
from pathlib import Path

from openmm import CustomExternalForce, LangevinMiddleIntegrator
from openmm.app import ForceField, Modeller, PDBFile, PME, Simulation
from openmm.unit import kelvin, nanometer, picosecond


DEFAULT_FF = 'auto'
DEFAULT_PH = 7.0
DEFAULT_PADDING = 1.0  # nm
DEFAULT_RESTRAINT_K = 100.0  # kcal/mol/A^2 for original atoms
DEFAULT_WEAK_K = 5.0  # kcal/mol/A^2 for added atoms (backbone only)
DEFAULT_MAX_ITER = 1000

BACKBONE_NAMES = {'N', 'CA', 'C', 'O', 'CB'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer minimize",
        description="Energy-minimize a PDB structure with OpenMM. Uses selective "
        "restraints: original atoms are restrained, newly added atoms (from "
        "PDBFixer .dat file) are free to relax.",
    )
    p.add_argument("input", help="Input PDB file")
    p.add_argument("-o", "--output", help="Output minimized PDB (default: <input>_minimized.pdb)")
    p.add_argument("--dat", help="Restraint data file from 'dvbfixer prepare' (default: <input>.dat)")
    p.add_argument("--ph", type=float, default=DEFAULT_PH,
                   help=f"pH for hydrogen addition if needed (default: {DEFAULT_PH})")
    p.add_argument("--ff", nargs='+', default=[DEFAULT_FF],
                   help="Force field selection. Accepts a short name "
                        "(auto, amber, amber+glycam, charmm, ...) or an "
                        "explicit list of OpenMM XML paths. Default: 'auto' — "
                        "detect from residue names in the input. "
                        "See docs/force-fields.md.")
    p.add_argument("--padding", type=float, default=DEFAULT_PADDING,
                   help=f"Solvent padding in nm (default: {DEFAULT_PADDING})")
    p.add_argument("--restraint-k", type=float, default=DEFAULT_RESTRAINT_K,
                   help=f"Restraint force constant for original atoms in kcal/mol/A^2 (default: {DEFAULT_RESTRAINT_K})")
    p.add_argument("--weak-k", type=float, default=DEFAULT_WEAK_K,
                   help=f"Restraint force constant for added backbone atoms in kcal/mol/A^2 (default: {DEFAULT_WEAK_K})")
    p.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER,
                   help=f"Max minimization iterations per phase (default: {DEFAULT_MAX_ITER})")
    p.add_argument("--rebuild-h", action="store_true",
                   help="Strip and re-add hydrogens via OpenMM (default: keep existing)")
    p.add_argument("--strip-heterogens", dest="keep_heterogens",
                   action="store_false", default=True,
                   help="Strip heterogens before minimization, restore coords after "
                        "(protein-only mode). Default: minimize the whole system.")
    p.add_argument("--no-solvent", action="store_true",
                   help="Minimize in vacuum (no solvent box)")
    p.add_argument("--xtb-refine", action="store_true",
                   help="After OpenMM minimization, run xtb GFN-FF universal "
                        "force field as a refinement pass. Auto-parametrizes any "
                        "organic molecule (sugars, ligands) without templates. "
                        "Requires `xtb` binary in PATH. Slower but higher quality.")
    p.add_argument("--xtb-cycles", type=int, default=200,
                   help="Max xtb optimization cycles (default: 200)")
    p.add_argument("--obminimize-refine", action="store_true",
                   help="After OpenMM minimization, run OpenBabel obminimize "
                        "(MMFF94) as a refinement pass. Auto-typing for any "
                        "organic molecule. Faster than xtb, lower quality.")
    p.add_argument("--obminimize-ff", default="UFF",
                   choices=["MMFF94", "MMFF94s", "UFF", "GAFF", "Ghemical"],
                   help="OpenBabel force field for --obminimize-refine (default: UFF — "
                        "handles N-glycosidic linkages correctly; MMFF94s mistypes the "
                        "anomeric C as sp2 giving 120° angles instead of 109°)")
    p.add_argument("--obminimize-steps", type=int, default=500,
                   help="OpenBabel minimization steps (default: 500)")
    p.add_argument("--refine-heterogens-only", action="store_true",
                   help="Restrict xtb/obminimize refinement to heterogen residues "
                        "(protein heavy atoms frozen). Default: refine whole system.")
    p.add_argument("--platform", choices=["CPU", "CUDA", "OpenCL", "Reference"],
                   help="OpenMM platform (default: auto-select fastest)")
    p.add_argument("--rename", action="store_true",
                   help="Rename non-canonical residues (AMBER/CHARMM) to standard names before processing")
    p.add_argument("--no-infer-conect", dest="no_infer_conect",
                   action="store_true",
                   help="Skip automatic CONECT inference (default: infer missing "
                        "SS / glycosidic / glycosylation bonds before minimize).")
    p.add_argument("--parametrize-ligands", action="store_true",
                   help="For each heterogen residue that lacks a template in "
                        "the resolved --ff, run GAFF2 + AM1-BCC parametrisation "
                        "via antechamber/parmchk2 and register the resulting "
                        "GAFF template with OpenMM before createSystem. "
                        "Requires openmmforcefields + openff-toolkit + "
                        "AmberTools (antechamber, parmchk2). Cached to "
                        "~/.cache/dvbfixer/lig_params/ (override with "
                        "$DVBFIXER_LIG_CACHE). See docs/force-fields.md.")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print detailed progress")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# .dat file loading
# ---------------------------------------------------------------------------

def load_dat(path):
    """Load .dat file, return set of (chain, resid, icode, atom_name) for added atoms."""
    with open(path) as f:
        dat = json.load(f)

    added_keys = set()
    for entry in dat["added_atoms"]:
        added_keys.add((entry["chain"], entry["resid"],
                        entry["icode"], entry["atom"]))

    print(f"Loaded restraint data: {path} ({len(added_keys)} added atoms)")
    return added_keys


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
    - Newly added sidechain atoms & all hydrogens: no restraint (free to move)
    - Heterogen heavy atoms (sugars, ligands): NO restraint — free to refine
      (BioLuminate-style: protein is fixed-ish, ligands relax)
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

        if is_hydrogen or is_heterogen:
            # All H atoms and all heterogen heavy atoms are free
            n_free += 1
            continue
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


# AMBER protonation variant names that need explicit variants in addHydrogens
AMBER_VARIANTS = {'HIE', 'HID', 'HIP', 'ASH', 'GLH', 'CYM', 'CYX', 'LYN'}


def _read_amber_renames(pdb_path):
    """Read PDB text to find AMBER protonation names before OpenMM normalizes them.

    OpenMM's PDBFile converts GLH->GLU, HIE->HIS, CYX->CYS on read.
    This reads the raw PDB to capture the original names.

    Returns dict: (chain_id, resid_str) -> amber_name
    """
    renames = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip()
            if resname in AMBER_VARIANTS:
                chain = line[21]
                resid = line[22:26].strip()
                renames[(chain, resid)] = resname
    return renames


def _build_variants(topology, amber_renames):
    """Build variants list from AMBER renames dict.

    Returns list of variant names (or None) for each residue, suitable for
    Modeller.addHydrogens(variants=...).
    """
    if not amber_renames:
        return None
    variants = []
    has_any = False
    for res in topology.residues():
        key = (res.chain.id, res.id)
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
        # Remap new_atom_indices to stripped topology
        known = PROTEIN_RESIDUES | SOLVENT_IONS
        old_to_new = {}
        new_idx = 0
        for atom in topology.atoms():
            if atom.residue.name in known:
                old_to_new[atom.index] = new_idx
                new_idx += 1
        stripped_new_indices = {old_to_new[i] for i in new_atom_indices if i in old_to_new}

    print("Loading force field...")
    # Check if there are heterogens — if so and keep_heterogens, use glycan-aware FF
    has_heterogens = args.keep_heterogens and any(
        res.name not in (PROTEIN_RESIDUES | SOLVENT_IONS)
        for res in stripped_top.residues()
    )

    # --parametrize-ligands: build GAFF2 templates for unknown ligands and
    # register them on the ForceField below via extra_generators. Strict mode:
    # any failure raises LigandParamError (propagated to main); we do NOT
    # silently fall back to strip-heterogens when the user explicitly asked
    # for parametrisation.
    if has_heterogens and getattr(args, 'parametrize_ligands', False):
        from openmm.app import ForceField as _FF
        base_templates = set(_FF(*args.ff)._templates)
        from dvbfixer.lig_params import build_ligand_generator
        _lig_gen = build_ligand_generator(
            stripped_top, stripped_pos, base_templates,
            strict=True, verbose=args.verbose,
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
    else:
        forcefield = ForceField(*args.ff)

    modeller = Modeller(stripped_top, stripped_pos)

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

    def _rename_variants_to_parent(top):
        """Rename variant residues to parent so addHydrogens recognises them.
        Returns dict {(chain_id, res_id): original_name} for the post-pass —
        residue refs in the OLD topology are stale after addHydrogens (which
        rebuilds the topology entirely).
        """
        saved = {}
        for res in top.residues():
            if res.name in _VARIANT_TO_PARENT:
                saved[(res.chain.id, res.id)] = res.name
                res.name = _VARIANT_TO_PARENT[res.name]
        return saved

    def _restore_variants_in_topology(top, saved):
        for res in top.residues():
            key = (res.chain.id, res.id)
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
    keep_h = not args.rebuild_h
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
            _fixer.addMissingAtoms()
            _fixer.addMissingHydrogens(args.ph)
            modeller = Modeller(_fixer.topology, _fixer.positions)
        elif not _has_hydrogens(modeller.topology):
            print("No hydrogens found, adding them...")
            variants = _build_variants(modeller.topology, amber_renames)
            _saved = _rename_variants_to_parent(modeller.topology)
            try:
                if variants:
                    modeller.addHydrogens(forcefield, pH=args.ph, variants=variants)
                else:
                    modeller.addHydrogens(forcefield, pH=args.ph)
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
                _saved = _rename_variants_to_parent(modeller.topology)
                try:
                    if variants:
                        modeller.addHydrogens(forcefield, pH=args.ph,
                                              variants=variants)
                    else:
                        modeller.addHydrogens(forcefield, pH=args.ph)
                finally:
                    _fix_lyn_hz_naming(modeller.topology, _saved)
                    _restore_variants_in_topology(modeller.topology, _saved)
            else:
                print("Keeping existing hydrogens from input")
    else:
        # Strip H, fix any missing heavy atoms via PDBFixer, then re-add correct H.
        # addHydrogens adds correct H for every residue, so stripping is safe.
        h_to_delete = [a for a in modeller.topology.atoms() if a.element.symbol == 'H']
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
        _fixer.addMissingAtoms()
        modeller = Modeller(_fixer.topology, _fixer.positions)

        print(f"Adding hydrogens (pH {args.ph})...")
        variants = _build_variants(modeller.topology, amber_renames)
        if variants:
            n_var = sum(1 for v in variants if v is not None)
            print(f"  {n_var} AMBER protonation variants detected")
        _saved = _rename_variants_to_parent(modeller.topology)
        try:
            if variants:
                modeller.addHydrogens(forcefield, pH=args.ph, variants=variants)
            else:
                modeller.addHydrogens(forcefield, pH=args.ph)
        finally:
            _fix_lyn_hz_naming(modeller.topology, _saved)
            _restore_variants_in_topology(modeller.topology, _saved)

    # When keeping heterogens with GLYCAM, OpenMM's PDBFile doesn't infer
    # intra-residue bonds for GLYCAM residues (NLN/OLS/OLT + sugars).
    # add_glycam_bonds populates them from FF templates.
    if has_heterogens:
        try:
            from dvbfixer.acpype_export import add_glycam_bonds
            add_glycam_bonds(modeller.topology, forcefield, args.verbose,
                              positions=modeller.positions)
        except Exception as e:
            if args.verbose:
                print(f"  add_glycam_bonds skipped: {e}")

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
        except Exception as e:
            from dvbfixer.ffutils import explain_template_error
            diag = explain_template_error(e, modeller.topology, forcefield)
            print(f"\nWARNING: pre-solvent parametrization failed:\n  {e}\n")
            if diag:
                for line in diag.split('\n'):
                    print(f"  {line}")
                print()
            print(f"  → falling back to --strip-heterogens flow.\n")
            import argparse as _ap
            legacy_args = _ap.Namespace(**vars(args))
            legacy_args.keep_heterogens = False
            legacy_args.rebuild_h = True
            out_top, out_pos = minimize(topology, positions, new_atom_indices,
                                         legacy_args, amber_renames=amber_renames)
            # Rigid-tracking inside the inner minimize already aligned the
            # glycan tree to the post-min anchor and snapped the linkage atom
            # to its ideal trigonal-planar position. Calling
            # _restore_glycosylated_h here would overwrite the post-min amide
            # with stale prep coords and break that alignment.
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
    except (ValueError, Exception) as e:
        if not has_heterogens:
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
        print(f"  → falling back to --strip-heterogens flow "
              f"(heterogens kept in output with original coords).\n")
        import argparse as _ap
        legacy_args = _ap.Namespace(**vars(args))
        legacy_args.keep_heterogens = False
        legacy_args.rebuild_h = True
        out_top, out_pos = minimize(topology, positions, new_atom_indices, legacy_args,
                                     amber_renames=amber_renames)
        # Rigid-tracking inside the inner minimize already aligned the
        # glycan tree to the post-min anchor and snapped the linkage atom
        # to its ideal trigonal-planar position. Calling
        # _restore_glycosylated_h here would overwrite the post-min amide
        # (CG/OD1/ND2/HD21) with stale prep coords, breaking the bond
        # to the rigid-tracked glycan. Skip it.
        return out_top, out_pos

    restraint, n_strong, n_weak, n_free = build_restraint_force(
        modeller.topology, modeller.positions,
        stripped_new_indices, args.restraint_k, args.weak_k
    )
    system.addForce(restraint)
    print(f"Restraints: {n_strong} strong (original), {n_weak} weak (new backbone), {n_free} free")

    integrator = LangevinMiddleIntegrator(300 * kelvin, 1.0 / picosecond, 0.002 * picosecond)

    if args.platform:
        from openmm import Platform
        platform = Platform.getPlatformByName(args.platform)
        simulation = Simulation(modeller.topology, system, integrator, platform)
    else:
        simulation = Simulation(modeller.topology, system, integrator)

    simulation.context.setPositions(modeller.positions)

    state = simulation.context.getState(getEnergy=True)
    print(f"Energy before minimization: {state.getPotentialEnergy()}")

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

    simulation.minimizeEnergy(maxIterations=args.max_iter)

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"Energy after phase 2: {state.getPotentialEnergy()}")

    min_topology = simulation.topology
    min_positions = state.getPositions()

    # Restore non-protein residues with original coordinates (legacy path only).
    # When keep_heterogens is on, n_stripped==0 and we write the minimized topology directly.
    if n_stripped > 0:
        # Strip solvent from minimized result first
        if not args.no_solvent:
            min_topology, min_positions = strip_solvent(min_topology, min_positions)

        # Merge: use minimized positions for protein atoms, original for HETATM
        # Match by (chain, resid, atomname) since addHydrogens may have changed indices
        from openmm.unit import nanometer as nm_unit
        import numpy as np

        # Build position lookup from minimized protein
        min_pos_map = {}
        for atom in min_topology.atoms():
            res = atom.residue
            key = (res.chain.id, res.id, atom.name)
            min_pos_map[key] = min_positions[atom.index].value_in_unit(nm_unit)

        # Build result: original positions, overwritten by minimized where available
        n_atoms = len(positions)
        result = np.zeros((n_atoms, 3))
        for i in range(n_atoms):
            result[i] = positions[i].value_in_unit(nm_unit)

        n_updated = 0
        for atom in topology.atoms():
            res = atom.residue
            key = (res.chain.id, res.id, atom.name)
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


def _find_binary(name):
    """Locate a binary, checking PATH and the current Python env's bin dir."""
    import shutil as _sh
    import os as _os
    import sys as _sys
    found = _sh.which(name)
    if found:
        return found
    # Fall back to Python env's bin directory (handles direct-executable env)
    py_bin = _os.path.dirname(_sys.executable)
    candidate = _os.path.join(py_bin, name)
    if _os.access(candidate, _os.X_OK):
        return candidate
    return None


def _write_xyz(path, topology, positions, comment=""):
    """Write topology + positions to XYZ format (used by xtb)."""
    from openmm.unit import angstrom, nanometer
    atoms = list(topology.atoms())
    with open(path, 'w') as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for atom in atoms:
            p = positions[atom.index].value_in_unit(nanometer)
            x = float(p[0]) * 10.0  # nm → Å
            y = float(p[1]) * 10.0
            z = float(p[2]) * 10.0
            sym = atom.element.symbol
            f.write(f"{sym:<3s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def _read_xyz_coords(path):
    """Read coordinates from XYZ file. Returns list of Vec3 in nm."""
    from openmm import Vec3
    coords = []
    with open(path) as f:
        n = int(f.readline().strip())
        f.readline()  # comment line
        for _ in range(n):
            parts = f.readline().split()
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            # Å → nm
            coords.append(Vec3(x / 10.0, y / 10.0, z / 10.0))
    return coords


def _build_frozen_atom_list(topology, heterogen_only):
    """Return list of atom indices (1-based) to freeze during xtb opt.

    If heterogen_only=True, freezes protein atoms; otherwise freezes nothing.
    """
    if not heterogen_only:
        return []
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    frozen = []
    for atom in topology.atoms():
        if atom.residue.name in known:
            frozen.append(atom.index + 1)  # 1-based for xtb
    return frozen


def refine_with_xtb(topology, positions, cycles=200, heterogens_only=False,
                    verbose=False):
    """Refine geometry with xtb GFN-FF universal force field.

    GFN-FF auto-parametrizes any organic molecule from connectivity rules.
    For whole-protein systems use heterogens_only=True (extract sub-system),
    otherwise xtb tries to handle the full system which is slow.
    """
    from openmm import Vec3
    from openmm.unit import nanometer, Quantity

    xtb_bin = _find_binary('xtb')
    if xtb_bin is None:
        print("WARNING: xtb binary not found in PATH — skipping xtb refinement")
        return positions

    if heterogens_only:
        sub_top, sub_pos, idx_map, anchor_indices = _extract_heterogen_subsystem(
            topology, positions
        )
        n_sub = sum(1 for _ in sub_top.atoms())
        if n_sub == 0:
            return positions
        print(f"\n=== xtb GFN-FF refinement ({cycles} cycles, "
              f"heterogens-only — {n_sub} atoms, "
              f"{len(anchor_indices)} protein anchors frozen) ===")
        new_sub_pos = _run_xtb(sub_top, sub_pos, xtb_bin, cycles, verbose,
                                frozen_indices=anchor_indices)
        if new_sub_pos is None:
            return positions
        coords = []
        for p in positions:
            v = p.value_in_unit(nanometer)
            coords.append(Vec3(float(v[0]), float(v[1]), float(v[2])))
        for full_idx, sub_idx in idx_map.items():
            if sub_idx in anchor_indices:
                continue  # anchors didn't move; skip
            p = new_sub_pos[sub_idx].value_in_unit(nanometer)
            coords[full_idx] = Vec3(float(p[0]), float(p[1]), float(p[2]))
        return Quantity(coords, nanometer)
    else:
        n = sum(1 for _ in topology.atoms())
        if n > 5000:
            print(f"\nINFO: full-system xtb on {n} atoms takes hours — "
                  f"auto-switching to --refine-heterogens-only (use that flag "
                  f"explicitly to silence this notice)")
            return refine_with_xtb(
                topology, positions, cycles=cycles,
                heterogens_only=True, verbose=verbose,
            )
        print(f"\n=== xtb GFN-FF refinement ({cycles} cycles, "
              f"whole system — {n} atoms) ===")
        new_pos = _run_xtb(topology, positions, xtb_bin, cycles, verbose)
        return new_pos if new_pos is not None else positions


def _run_xtb(topology, positions, xtb_bin, cycles, verbose, frozen_indices=None):
    """Internal: write XYZ, run xtb --opt --gfnff, return refined Quantity or None.

    frozen_indices: set of sub-topology atom indices (0-based) to freeze.
    Written to xcontrol with $fix atoms (1-based).
    """
    import shutil as _sh
    import subprocess as _sp
    import tempfile as _tf
    import os as _os
    from openmm.unit import nanometer, Quantity

    workdir = _tf.mkdtemp(prefix='dvbfixer_xtb_')
    old_cwd = _os.getcwd()
    try:
        input_xyz = _os.path.join(workdir, 'input.xyz')
        _write_xyz(input_xyz, topology, positions, "dvbfixer xtb input")

        cmd = [xtb_bin, 'input.xyz', '--opt', '--gfnff',
               '--cycles', str(cycles), '--norestart']

        # Freeze protein anchor atoms via xcontrol $fix block (1-based indices)
        if frozen_indices:
            xc_path = _os.path.join(workdir, 'xcontrol')
            sorted_idx = sorted(i + 1 for i in frozen_indices)
            # Group consecutive indices into ranges for compactness
            ranges = []
            start = prev = sorted_idx[0]
            for i in sorted_idx[1:]:
                if i == prev + 1:
                    prev = i
                else:
                    ranges.append((start, prev))
                    start = prev = i
            ranges.append((start, prev))
            range_str = ",".join(
                f"{s}" if s == e else f"{s}-{e}" for s, e in ranges
            )
            with open(xc_path, 'w') as f:
                f.write("$fix\n")
                f.write(f"  atoms: {range_str}\n")
                f.write("$end\n")
            cmd.extend(['--input', 'xcontrol'])
        _os.chdir(workdir)
        if verbose:
            print(f"  Running: {' '.join(cmd)}")
        try:
            result = _sp.run(cmd, capture_output=True, text=True, timeout=3600)
        except _sp.TimeoutExpired:
            print("WARNING: xtb timeout")
            return None
        if result.returncode != 0:
            print(f"WARNING: xtb failed (returncode {result.returncode})")
            if verbose:
                print("--- stderr ---")
                print((result.stderr or "")[-3000:])
                print("--- stdout (last 30 lines) ---")
                print('\n'.join(result.stdout.splitlines()[-30:]))
            return None

        opt_xyz = _os.path.join(workdir, 'xtbopt.xyz')
        if not _os.path.exists(opt_xyz):
            print("WARNING: xtb produced no xtbopt.xyz")
            return None

        new_coords = _read_xyz_coords(opt_xyz)
        for line in result.stdout.splitlines():
            if 'TOTAL ENERGY' in line:
                print(f"  xtb {line.strip()}")
                break
        print(f"  xtb refined {len(new_coords)} atoms")
        return Quantity(new_coords, nanometer)
    finally:
        _os.chdir(old_cwd)
        _sh.rmtree(workdir, ignore_errors=True)


def _rigid_track_glycan_trees(in_top, in_pos, result_array, min_pos_map,
                               verbose=False):
    """Rigid-body transform each glycan tree to follow its protein anchor.

    The legacy strip-and-splice path minimizes protein-only, then restores
    HETATM coords verbatim. But the protein anchor (ASN/SER/THR side chain)
    moves during minimization, leaving the glycan in the wrong place — the
    ND2-C1 bond stretches and the glycan can clash with the moved protein.

    Fix: for each protein-heterogen bond, compute the Kabsch (rotation +
    translation) transform from prep anchor heavy atoms → post-min anchor
    heavy atoms, then apply it to every atom in the BFS-connected glycan
    tree. Result: glycan keeps its relative orientation to the anchor amide;
    bond length and stereochemistry are preserved.

    in_top/in_pos: original (prep) topology + positions.
    result_array: numpy (n_atoms, 3) — protein already overwritten with
                  minimized positions, HETATM still at prep positions.
    min_pos_map: dict (chain, resid, atomname) → minimized position (nm).
    """
    import numpy as np
    from openmm.unit import nanometer as nm_unit
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    known = PROTEIN_RESIDUES | SOLVENT_IONS
    atoms = list(in_top.atoms())

    # Find glycosylated protein residues + which heterogen atom they connect to.
    het_atom_set = {a.index for a in atoms if a.residue.name not in known}
    anchor_to_het_atoms = {}  # (chain, rid) → set of bonded heterogen atom indices
    for b in in_top.bonds():
        ai, bi = b[0].index, b[1].index
        if ai in het_atom_set and bi not in het_atom_set:
            r = b[1].residue
            anchor_to_het_atoms.setdefault((r.chain.id, r.id), set()).add(ai)
        elif bi in het_atom_set and ai not in het_atom_set:
            r = b[0].residue
            anchor_to_het_atoms.setdefault((r.chain.id, r.id), set()).add(bi)

    if not anchor_to_het_atoms:
        return result_array

    # Build heterogen adjacency for BFS through the glycan tree.
    het_adj = {idx: set() for idx in het_atom_set}
    for b in in_top.bonds():
        ai, bi = b[0].index, b[1].index
        if ai in het_atom_set and bi in het_atom_set:
            het_adj[ai].add(bi)
            het_adj[bi].add(ai)

    # Geometry recipe per (residue, anchor_atom). For each recipe, we
    # compute the IDEAL linkage atom (e.g. NAG C1) position from the post-min
    # protein side-chain geometry — bypassing the AMBER-placed HD21/HD22 which
    # may have been put on the wrong side.
    #
    # N-linked (ASN/NLN ND2): trans amide. C1 lies in the amide plane (defined
    # by CG, OD1, ND2), at 120° from the CG-ND2 axis, on the OPPOSITE side of
    # the plane from OD1. This is the canonical E,Z amide configuration.
    #
    # O-linked (SER/THR): sp3 tetrahedral. C1 substitutes for the hydroxyl H,
    # at the position the H currently occupies (anti-periplanar to CB).
    _IDEAL_GEOMETRY = {
        ('ASN', 'ND2'): {'bond_nm': 0.145, 'type': 'amide_trans',
                          'axis': 'CG', 'plane_neighbor': 'OD1', 'angle_deg': 120.0},
        ('NLN', 'ND2'): {'bond_nm': 0.145, 'type': 'amide_trans',
                          'axis': 'CG', 'plane_neighbor': 'OD1', 'angle_deg': 120.0},
        ('SER', 'OG'):  {'bond_nm': 0.143, 'type': 'replace_H', 'h_atom': 'HG'},
        ('OLS', 'OG'):  {'bond_nm': 0.143, 'type': 'replace_H', 'h_atom': 'HG'},
        ('THR', 'OG1'): {'bond_nm': 0.143, 'type': 'replace_H', 'h_atom': 'HG1'},
        ('OLT', 'OG1'): {'bond_nm': 0.143, 'type': 'replace_H', 'h_atom': 'HG1'},
    }

    def _ideal_linkage_pos(anchor_atom_name, anchor_res_name, anchor_ch, anchor_rid):
        """Compute the IDEAL linkage atom position from the post-min anchor
        side-chain geometry. Returns ideal_pos (numpy, nm) or None on missing
        reference atoms.
        """
        geom = _IDEAL_GEOMETRY.get((anchor_res_name, anchor_atom_name))
        if geom is None:
            return None
        akey = (anchor_ch, anchor_rid, anchor_atom_name)
        if akey not in min_pos_map:
            return None
        a_pos = np.array(min_pos_map[akey])

        if geom['type'] == 'replace_H':
            hkey = (anchor_ch, anchor_rid, geom['h_atom'])
            if hkey not in min_pos_map:
                return None
            h_pos = np.array(min_pos_map[hkey])
            v = h_pos - a_pos
            n = np.linalg.norm(v)
            if n < 1e-6:
                return None
            return a_pos + geom['bond_nm'] * (v / n)

        if geom['type'] == 'amide_trans':
            ax_key = (anchor_ch, anchor_rid, geom['axis'])
            pn_key = (anchor_ch, anchor_rid, geom['plane_neighbor'])
            if ax_key not in min_pos_map or pn_key not in min_pos_map:
                return None
            ax_pos = np.array(min_pos_map[ax_key])         # e.g. CG
            pn_pos = np.array(min_pos_map[pn_key])         # e.g. OD1
            # u_axis = unit vector from anchor (ND2) toward axis atom (CG)
            u_axis = ax_pos - a_pos
            n = np.linalg.norm(u_axis)
            if n < 1e-6:
                return None
            u_axis = u_axis / n
            # OD1 relative to anchor (ND2)
            v_pn = pn_pos - a_pos
            # In-plane perpendicular: component of OD1-ND2 perpendicular to
            # CG-ND2 axis. This points TOWARD OD1 in the amide plane.
            v_perp = v_pn - np.dot(v_pn, u_axis) * u_axis
            n_perp = np.linalg.norm(v_perp)
            if n_perp < 1e-6:
                return None
            u_perp_toward_pn = v_perp / n_perp
            # Trans-amide C1 direction:
            #   angle between u_axis and u_C1 = 120° (CG-ND2-C1)
            #   u_C1 in amide plane, on OPPOSITE side from OD1 across CG-ND2.
            # In-plane decomposition:
            #   u_C1 = cos(angle) * (-u_axis_from_anchor_to_axis) (so angle from u_axis = 120°)
            # equivalently: u_C1 = cos(180° - angle) * u_axis + sin(angle) * (-u_perp_toward_pn)
            # but simpler:
            theta = np.radians(geom['angle_deg'])
            # u_C1 is at angle theta from u_axis (measured at ND2), in the amide
            # plane, on the side AWAY from OD1.
            #   u_C1 = cos(theta) * u_axis + sin(theta) * (-u_perp_toward_pn)
            u_c1 = np.cos(theta) * u_axis + np.sin(theta) * (-u_perp_toward_pn)
            return a_pos + geom['bond_nm'] * u_c1

        return None

    n_trees = 0
    n_atoms_moved = 0
    for (anchor_ch, anchor_rid), root_het_atoms in anchor_to_het_atoms.items():
        # Anchor residue in the input topology
        anchor_res = None
        for r in in_top.residues():
            if r.chain.id == anchor_ch and r.id == anchor_rid:
                anchor_res = r
                break
        if anchor_res is None:
            continue

        # Collect heavy atoms of the anchor side chain + linkage neighbor for
        # the Kabsch fit. We exclude H atoms (their positions are dominated by
        # AMBER addHydrogens, not by minimization itself).
        prep_pts = []
        post_pts = []
        for a in anchor_res.atoms():
            if a.element.symbol == 'H':
                continue
            key = (a.residue.chain.id, a.residue.id, a.name)
            if key not in min_pos_map:
                continue
            pp = in_pos[a.index].value_in_unit(nm_unit)
            mp = min_pos_map[key]
            prep_pts.append([float(pp[0]), float(pp[1]), float(pp[2])])
            post_pts.append([float(mp[0]), float(mp[1]), float(mp[2])])

        if len(prep_pts) < 3:
            continue

        # Inject the IDEAL linkage atom position as an extra Kabsch reference.
        # Without this, the bad prep-side glycosidic geometry (e.g. CG-ND2-C1
        # ≈ 167° instead of 122° after Modeller loop modeling) is faithfully
        # reproduced at the post-min anchor — the Kabsch is dominated by the
        # anchor-side coordinates and preserves the relative orientation.
        # Injecting the ideal C1 (from post-min amide trigonal-planar geometry)
        # pulls the fit toward correct stereochemistry.
        anchor_atoms_by_name = {a.name: a for a in anchor_res.atoms()}
        # Use the first bonded heterogen atom as the linkage atom (e.g. NAG C1)
        link_het_idx = next(iter(root_het_atoms))
        link_het_atom = atoms[link_het_idx]
        # Find the protein anchor atom name (e.g. ND2 for ASN)
        anchor_link_atom_name = None
        for b in in_top.bonds():
            ai, bi = b[0].index, b[1].index
            if ai == link_het_idx and bi in {a.index for a in anchor_res.atoms()}:
                anchor_link_atom_name = atoms[bi].name
                break
            if bi == link_het_idx and ai in {a.index for a in anchor_res.atoms()}:
                anchor_link_atom_name = atoms[ai].name
                break
        # Resolve the anchor residue name. The legacy path renames NLN→ASN
        # in topology BEFORE AMBER addHydrogens, so anchor_res.name may now
        # be 'ASN' even though the prep input had 'NLN'. Both keys point to
        # the same trigonal-planar geometry recipe (CG+HD21), so either name
        # works — we just need a registered entry.
        post_anchor_name = anchor_res.name
        ideal_pos = None
        if anchor_link_atom_name:
            for try_name in (anchor_res.name,
                             {'ASN': 'NLN', 'SER': 'OLS', 'THR': 'OLT',
                              'NLN': 'ASN', 'OLS': 'SER', 'OLT': 'THR'}.get(anchor_res.name)):
                if try_name is None:
                    continue
                if (try_name, anchor_link_atom_name) in _IDEAL_GEOMETRY:
                    post_anchor_name = try_name
                    break
            ideal_pos = _ideal_linkage_pos(anchor_link_atom_name,
                                            post_anchor_name,
                                            anchor_ch, anchor_rid)
        if ideal_pos is not None:
            prep_link = in_pos[link_het_idx].value_in_unit(nm_unit)
            prep_pts.append([float(prep_link[0]), float(prep_link[1]),
                              float(prep_link[2])])
            post_pts.append(ideal_pos.tolist())

        prep_arr = np.array(prep_pts)
        post_arr = np.array(post_pts)
        # Kabsch
        prep_c = prep_arr.mean(axis=0)
        post_c = post_arr.mean(axis=0)
        P = prep_arr - prep_c
        Q = post_arr - post_c
        H = P.T @ Q
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        t = post_c - R @ prep_c

        # BFS to collect glycan tree atoms (heterogen atoms reachable from
        # the anchor's bonded heterogen atoms).
        seen = set()
        queue = list(root_het_atoms)
        while queue:
            cur = queue.pop()
            if cur in seen:
                continue
            seen.add(cur)
            for nb in het_adj.get(cur, ()):
                if nb not in seen:
                    queue.append(nb)

        # Apply transform to each glycan atom
        for idx in seen:
            p = in_pos[idx].value_in_unit(nm_unit)
            v = np.array([float(p[0]), float(p[1]), float(p[2])])
            v_new = R @ v + t
            result_array[idx] = v_new
            n_atoms_moved += 1

        # Post-Kabsch correction: snap the linkage atom exactly to the ideal
        # post-min trigonal-planar position by translating the whole glycan
        # tree by (ideal - tracked). The Kabsch fit only approximately satisfies
        # the ideal C1 target (it gets averaged with the 5 ASN reference
        # points), but for correct stereochemistry the C1 needs to land on
        # the exact amide plane direction. Translation preserves intra-glycan
        # geometry; the next obminimize pass relaxes any residual strain.
        if anchor_link_atom_name and 'ideal_pos' in locals() and ideal_pos is not None:
            tracked_link = result_array[link_het_idx]
            correction = ideal_pos - tracked_link
            for idx in seen:
                result_array[idx] = result_array[idx] + correction

            # Also snap the amide H (HD21 for ASN/NLN) to its canonical
            # trigonal-planar position — cis to OD1, in the amide plane.
            # AMBER places HD21 at some valid sp2 position but not always at
            # the canonical cis-OD1 spot; this enforces CG-ND2-HD21 ≈ 120°
            # and the H in the amide plane.
            if (post_anchor_name, anchor_link_atom_name) == ('ASN', 'ND2') or \
               (post_anchor_name, anchor_link_atom_name) == ('NLN', 'ND2'):
                # Find HD21 atom in in_top (prep NLN topology has HD21)
                hd21_in_idx = None
                for a in anchor_res.atoms():
                    if a.name == 'HD21':
                        hd21_in_idx = a.index
                        break
                if hd21_in_idx is not None:
                    # Compute ideal HD21 position (cis to OD1, in amide plane).
                    # Same formula as for C1 but with opposite perp sign.
                    nd2_key = (anchor_ch, anchor_rid, 'ND2')
                    cg_key = (anchor_ch, anchor_rid, 'CG')
                    od1_key = (anchor_ch, anchor_rid, 'OD1')
                    if all(k in min_pos_map for k in (nd2_key, cg_key, od1_key)):
                        nd2_p = np.array(min_pos_map[nd2_key])
                        cg_p = np.array(min_pos_map[cg_key])
                        od1_p = np.array(min_pos_map[od1_key])
                        u_ax = cg_p - nd2_p
                        u_ax = u_ax / np.linalg.norm(u_ax)
                        v_pn = od1_p - nd2_p
                        v_perp = v_pn - np.dot(v_pn, u_ax) * u_ax
                        n_perp = np.linalg.norm(v_perp)
                        if n_perp > 1e-6:
                            u_perp_to_od1 = v_perp / n_perp
                            theta = np.radians(120.0)
                            # HD21: trigonal-planar, cis to OD1 → same perp sign as OD1
                            u_h = np.cos(theta) * u_ax + np.sin(theta) * u_perp_to_od1
                            result_array[hd21_in_idx] = nd2_p + 0.101 * u_h  # N-H 1.01 Å
        n_trees += 1

    if verbose and n_trees:
        print(f"  Rigid-tracked {n_trees} glycan tree(s) "
              f"to follow protein anchor ({n_atoms_moved} atoms)")
    return result_array


def _restore_glycosylated_h(out_top, out_pos, in_top, in_pos):
    """Restore side-chain positions of glycosylated ASN/SER/THR (NLN/OLS/OLT)
    to their input values to preserve the protein-glycan stereochemistry.

    After legacy strip-and-splice minimize, the AMBER ASN template treats the
    glycosylated residue as a normal amide. The minimization can:
    - Flip HD21/HD22 to the wrong side of ND2 (collision with linked C1)
    - Rotate the side chain so the amide plane no longer aligns with the sugar

    Restoring the entire side chain (CB, HB*, CG, OD1, ND2, HD21 for ASN;
    CB, OG, HG for SER; CB, CG2, OG1, HG1 for THR) snaps the protein-glycan
    interface back to the prepared geometry. Backbone atoms remain refined.
    """
    from openmm import Vec3
    from openmm.unit import nanometer, Quantity
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    known = PROTEIN_RESIDUES | SOLVENT_IONS

    glycosylated_residues = set()
    het_atom_set = {a.index for a in in_top.atoms() if a.residue.name not in known}
    for b in in_top.bonds():
        ai, bi = b[0].index, b[1].index
        if ai in het_atom_set and bi not in het_atom_set:
            r = b[1].residue
            glycosylated_residues.add((r.chain.id, r.id))
        elif bi in het_atom_set and ai not in het_atom_set:
            r = b[0].residue
            glycosylated_residues.add((r.chain.id, r.id))
    for r in in_top.residues():
        if r.name in ('NLN', 'OLS', 'OLT'):
            glycosylated_residues.add((r.chain.id, r.id))

    if not glycosylated_residues:
        return out_top, out_pos

    # Atoms to restore: ONLY the amide/hydroxyl group at the glycosidic linkage
    # (NOT HB2/HB3 or other tetrahedral H — those need AMBER-quality minimization
    # for proper sp3 angles). Restoring CG-OD1-ND2-HD21 preserves amide planarity
    # and HD21 orientation w.r.t. the linked sugar. CB and HB* are NOT touched —
    # they get the minimized sp3 geometry.
    _SIDECHAIN = {
        'ASN': {'CG', 'OD1', 'ND2', 'HD21', 'HD22'},
        'NLN': {'CG', 'OD1', 'ND2', 'HD21'},
        'SER': {'OG', 'HG', 'HG1'},
        'OLS': {'OG'},
        'THR': {'OG1', 'HG1'},
        'OLT': {'OG1'},
    }

    input_positions = {}
    for atom in in_top.atoms():
        key = (atom.residue.chain.id, atom.residue.id, atom.name)
        p = in_pos[atom.index].value_in_unit(nanometer)
        input_positions[key] = (float(p[0]), float(p[1]), float(p[2]))

    n_restored = 0
    out_pos_list = []
    for atom in out_top.atoms():
        p = out_pos[atom.index].value_in_unit(nanometer)
        v = Vec3(float(p[0]), float(p[1]), float(p[2]))
        res_key = (atom.residue.chain.id, atom.residue.id)
        if res_key in glycosylated_residues:
            # Determine which side-chain atom set to use. Output may have
            # ASN (renamed from NLN); input may have NLN — try both.
            sc_atoms = (_SIDECHAIN.get(atom.residue.name, set())
                        | _SIDECHAIN.get('NLN', set())
                        | _SIDECHAIN.get('OLS', set())
                        | _SIDECHAIN.get('OLT', set()))
            if atom.name in sc_atoms:
                ipos = input_positions.get(
                    (atom.residue.chain.id, atom.residue.id, atom.name)
                )
                if ipos is not None:
                    v = Vec3(*ipos)
                    n_restored += 1
        out_pos_list.append(v)
    if n_restored:
        print(f"  Restored {n_restored} side-chain atoms on glycosylated residues")
    return out_top, Quantity(out_pos_list, nanometer)


def _extract_heterogen_subsystem(topology, positions, padding_residues=False):
    """Build a sub-topology of heterogen residues PLUS protein anchor atoms.

    For each protein-heterogen bond (e.g. ASN ND2 → NAG C1), the protein
    atom is included as a single "anchor" atom (without its full residue).
    This preserves the glycosidic bond constraint during refinement; downstream
    the anchor atoms are frozen so the protein doesn't move.

    Returns (sub_topology, sub_positions, atom_index_map, anchor_sub_indices).
    atom_index_map[full_atom_index] = sub_atom_index for all included atoms.
    anchor_sub_indices = set of sub_topology atom indices that are protein
                        anchors (should be frozen during refinement).
    """
    from openmm.app import Topology, element
    from openmm import Vec3
    from openmm.unit import nanometer, Quantity
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    known = PROTEIN_RESIDUES | SOLVENT_IONS

    # Find protein atoms bonded to heterogens, then include their ENTIRE
    # residues as anchors. Including just the single linkage atom (e.g. ASN
    # ND2 alone) leaves xtb confused about bond chemistry — it perceives a
    # lone N and pulls the C-N bond to ~1 Å. Including the full ASN gives
    # xtb the proper sp3 N context.
    het_atom_indices = set()
    for res in topology.residues():
        if res.name not in known:
            for atom in res.atoms():
                het_atom_indices.add(atom.index)
    anchor_residues = set()  # residues whose atoms anchor a heterogen
    linkage_het_atoms = set()  # heterogen atoms at the cross-residue bond
    for b in topology.bonds():
        ai, bi = b[0].index, b[1].index
        if ai in het_atom_indices and bi not in het_atom_indices:
            anchor_residues.add(b[1].residue.index)
            linkage_het_atoms.add(ai)
        elif bi in het_atom_indices and ai not in het_atom_indices:
            anchor_residues.add(b[0].residue.index)
            linkage_het_atoms.add(bi)
    anchor_protein_atoms = set()
    for res in topology.residues():
        if res.index in anchor_residues:
            for atom in res.atoms():
                anchor_protein_atoms.add(atom.index)
    # ALSO freeze the heterogen linkage atom itself (e.g. NAG C1 bonded to
    # ASN ND2). Reason: UFF/MMFF don't have explicit amide planarity terms,
    # so during refinement the linkage C can rotate ~15° out of the amide
    # plane. Freezing it preserves the OpenMM-AMBER geometry of the
    # protein-glycan interface; the rest of the sugar refines freely.

    sub_top = Topology()
    sub_pos_list = []
    atom_index_map = {}
    anchor_sub_indices = set()

    # First pass: add heterogen residues
    for chain in topology.chains():
        new_chain = None
        for res in chain.residues():
            if res.name in known:
                continue
            if new_chain is None:
                new_chain = sub_top.addChain(chain.id)
            new_res = sub_top.addResidue(res.name, new_chain, res.id,
                                          res.insertionCode)
            for atom in res.atoms():
                new_atom = sub_top.addAtom(atom.name, atom.element, new_res)
                atom_index_map[atom.index] = new_atom.index
                # Linkage heterogen atom (e.g. NAG C1) is INTENTIONALLY LEFT
                # FREE so UFF can minimize the protein-glycan bond length and
                # angle. UFF's sp2 N type (N_2) has equilibrium angle 120° and
                # an improper torsion that keeps the amide near-planar, so the
                # geometry settles naturally instead of being snapped rigid.
                p = positions[atom.index].value_in_unit(nanometer)
                sub_pos_list.append(Vec3(float(p[0]), float(p[1]), float(p[2])))

    # Second pass: add full anchor residues with their original names/atoms.
    # Use a separate "Z" chain so they don't merge with heterogen residues.
    if anchor_protein_atoms:
        anchor_chain = sub_top.addChain("Z")
        from collections import defaultdict
        anchors_by_res = defaultdict(list)
        full_atoms = list(topology.atoms())
        for ai in anchor_protein_atoms:
            anchors_by_res[full_atoms[ai].residue].append(ai)
        for orig_res, atom_indices in anchors_by_res.items():
            anc_res = sub_top.addResidue(orig_res.name, anchor_chain,
                                         str(orig_res.id),
                                         orig_res.insertionCode)
            # Freeze only the BACKBONE + first sidechain carbon of the anchor
            # residue. The amide group (CG, OD1, ND2, HD21, HD22 for ASN; OG/HG
            # for SER; OG1/HG1 for THR) and the linkage atom are LEFT FREE so
            # UFF/MMFF can actually minimize the protein-glycan bond length,
            # angle, and amide planarity. This is the answer to "why doesn't
            # the tool minimize the ASN-glycan bond" — previously the whole
            # anchor was frozen and the bond was only rigid-tracked, not
            # minimized by any FF.
            _FREEZE_ANCHOR_ATOMS = {
                'ASN': {'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O'},
                'NLN': {'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O'},
                'SER': {'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O'},
                'OLS': {'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O'},
                'THR': {'N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG2',
                        'HG21', 'HG22', 'HG23', 'C', 'O'},
                'OLT': {'N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG2',
                        'HG21', 'HG22', 'HG23', 'C', 'O'},
            }
            freeze_names = _FREEZE_ANCHOR_ATOMS.get(
                orig_res.name,
                # fallback: freeze everything for unfamiliar anchors
                {a.name for a in [full_atoms[ai] for ai in atom_indices]}
            )
            for ai in atom_indices:
                a = full_atoms[ai]
                new_atom = sub_top.addAtom(a.name, a.element, anc_res)
                atom_index_map[ai] = new_atom.index
                if a.name in freeze_names:
                    anchor_sub_indices.add(new_atom.index)
                p = positions[ai].value_in_unit(nanometer)
                sub_pos_list.append(Vec3(float(p[0]), float(p[1]), float(p[2])))

    # Carry all bonds where both atoms made it into the sub-topology
    sub_atoms = list(sub_top.atoms())
    for b in topology.bonds():
        a1 = atom_index_map.get(b[0].index)
        a2 = atom_index_map.get(b[1].index)
        if a1 is not None and a2 is not None:
            sub_top.addBond(sub_atoms[a1], sub_atoms[a2])

    sub_pos = Quantity(sub_pos_list, nanometer)
    return sub_top, sub_pos, atom_index_map, anchor_sub_indices


def refine_with_obminimize(topology, positions, ff='MMFF94s', steps=500,
                            heterogens_only=False, verbose=False):
    """Refine geometry with OpenBabel obminimize (MMFF94/UFF/GAFF).

    OpenBabel auto-types any organic molecule via SMARTS rules.
    Returns refined positions (Quantity).
    """
    import shutil as _sh
    import subprocess as _sp
    import tempfile as _tf
    import os as _os
    from openmm.app import PDBFile
    from openmm.unit import nanometer, Quantity
    from openmm import Vec3

    if _find_binary('obminimize') is None:
        print("WARNING: obminimize binary not found in PATH — skipping refinement")
        return positions

    if heterogens_only:
        sub_top, sub_pos, idx_map, anchor_indices = _extract_heterogen_subsystem(
            topology, positions
        )
        n_sub = sum(1 for _ in sub_top.atoms())
        if n_sub == 0:
            return positions
        print(f"\n=== OpenBabel obminimize refinement ({ff}, {steps} steps, "
              f"heterogens-only — {n_sub} atoms, "
              f"{len(anchor_indices)} protein anchors frozen) ===")
        new_sub_pos = _run_obminimize(sub_top, sub_pos, ff, steps, verbose,
                                       frozen_indices=anchor_indices)
        if new_sub_pos is None:
            return positions
        coords = []
        for p in positions:
            v = p.value_in_unit(nanometer)
            coords.append(Vec3(float(v[0]), float(v[1]), float(v[2])))
        for full_idx, sub_idx in idx_map.items():
            if sub_idx in anchor_indices:
                continue  # anchors didn't move
            p = new_sub_pos[sub_idx].value_in_unit(nanometer)
            coords[full_idx] = Vec3(float(p[0]), float(p[1]), float(p[2]))
        return Quantity(coords, nanometer)
    else:
        n_atoms = sum(1 for _ in topology.atoms())
        if n_atoms > 5000:
            print(f"\nINFO: full-system obminimize on {n_atoms} atoms is "
                  f"memory-intensive — auto-switching to --refine-heterogens-only "
                  f"(use that flag explicitly to silence this notice)")
            return refine_with_obminimize(
                topology, positions, ff=ff, steps=steps,
                heterogens_only=True, verbose=verbose,
            )
        print(f"\n=== OpenBabel obminimize refinement ({ff}, {steps} steps, "
              f"whole system — {n_atoms} atoms) ===")
        new_pos = _run_obminimize(topology, positions, ff, steps, verbose)
        return new_pos if new_pos is not None else positions


def _run_obminimize_pybel(topology, positions, ff, steps, frozen_indices, verbose):
    """Run OpenBabel via Python API with atom freezing.

    OBFFConstraints.AddAtomConstraint(idx) freezes atoms by 1-based index.
    Used when frozen_indices is non-empty (CLI obminimize has no freeze flag).
    """
    import tempfile as _tf
    import os as _os
    from openmm.app import PDBFile
    from openmm.unit import nanometer, Quantity
    from openmm import Vec3

    try:
        from openbabel import openbabel as ob, pybel
    except ImportError:
        print("WARNING: openbabel Python bindings missing — skipping refinement")
        return None

    workdir = _tf.mkdtemp(prefix='dvbfixer_obmin_')
    try:
        in_pdb = _os.path.join(workdir, 'in.pdb')
        with open(in_pdb, 'w') as f:
            PDBFile.writeFile(topology, positions, f, keepIds=True)

        # Try requested FF, fall back to UFF on setup failure
        ffs_to_try = [ff] if ff == 'UFF' else [ff, 'UFF']
        for try_ff in ffs_to_try:
            mol = next(pybel.readfile('pdb', in_pdb))
            obff = ob.OBForceField.FindForceField(try_ff)
            if obff is None:
                continue
            constraints = ob.OBFFConstraints()
            for idx in sorted(frozen_indices):
                constraints.AddAtomConstraint(idx + 1)  # 1-based
            if not obff.Setup(mol.OBMol, constraints):
                if verbose:
                    print(f"  {try_ff} Setup failed; trying next FF")
                continue
            obff.SteepestDescent(steps)
            obff.GetCoordinates(mol.OBMol)
            used_ff = try_ff
            if try_ff != ff:
                print(f"  ({ff} setup failed; refined with UFF instead)")
            print(f"  obminimize refined {mol.OBMol.NumAtoms()} atoms "
                  f"({used_ff}, {len(frozen_indices)} frozen)")
            # Extract coords
            n_atoms = sum(1 for _ in topology.atoms())
            if mol.OBMol.NumAtoms() != n_atoms:
                print("WARNING: atom count mismatch after obminimize")
                return None
            coords = []
            for ai in range(1, mol.OBMol.NumAtoms() + 1):
                ob_a = mol.OBMol.GetAtom(ai)
                # OpenBabel stores Å; convert to nm
                coords.append(Vec3(ob_a.GetX() / 10.0,
                                   ob_a.GetY() / 10.0,
                                   ob_a.GetZ() / 10.0))
            return Quantity(coords, nanometer)
        print(f"WARNING: obminimize Python API: no FF could be set up")
        return None
    finally:
        import shutil as _sh
        _sh.rmtree(workdir, ignore_errors=True)


def _run_obminimize(topology, positions, ff, steps, verbose, frozen_indices=None):
    """Internal: run obminimize via CLI (no freeze) or pybel API (with freeze).

    When the requested FF (MMFF94/MMFF94s) lacks parameters for some atoms,
    automatically retries with UFF (universal FF, covers any element).
    """
    if frozen_indices:
        return _run_obminimize_pybel(
            topology, positions, ff, steps, frozen_indices, verbose
        )
    import shutil as _sh
    import subprocess as _sp
    import tempfile as _tf
    import os as _os
    from openmm.app import PDBFile
    from openmm.unit import nanometer, Quantity
    from openmm import Vec3

    workdir = _tf.mkdtemp(prefix='dvbfixer_obmin_')
    try:
        in_pdb = _os.path.join(workdir, 'in.pdb')
        out_pdb = _os.path.join(workdir, 'out.pdb')
        with open(in_pdb, 'w') as f:
            PDBFile.writeFile(topology, positions, f, keepIds=True)
        obmin_bin = _find_binary('obminimize')

        # Try requested FF, fall back to UFF if param lookup fails
        ffs_to_try = [ff]
        if ff != 'UFF':
            ffs_to_try.append('UFF')

        result = None
        used_ff = None
        for try_ff in ffs_to_try:
            cmd = [obmin_bin, '-ff', try_ff, '-n', str(steps), in_pdb]
            try:
                result = _sp.run(cmd, capture_output=True, text=True, timeout=3600)
            except _sp.TimeoutExpired:
                print("WARNING: obminimize timeout")
                return None
            stderr = result.stderr or ""
            stdout = result.stdout or ""
            param_missing = (
                'could not setup force field' in stderr.lower()
                or 'could not find van der waals' in stdout.lower()
                or 'could not find van der waals' in stderr.lower()
            )
            if result.returncode == 0 and stdout.strip():
                used_ff = try_ff
                if try_ff != ff:
                    print(f"  ({ff} lacked parameters; refined with UFF instead)")
                break
            if param_missing and try_ff != ffs_to_try[-1]:
                print(f"  {try_ff} parameters incomplete — retrying with UFF...")
                continue
            # Hard failure or unexpected output — give up
            print(f"WARNING: obminimize failed with {try_ff} "
                  f"(rc={result.returncode}). Output left unrefined.")
            if verbose:
                print((stderr or stdout)[-1500:])
            return None

        if used_ff is None:
            return None
        with open(out_pdb, 'w') as f:
            f.write(result.stdout)
        new_pdb = PDBFile(out_pdb)
        if sum(1 for _ in new_pdb.topology.atoms()) != sum(1 for _ in topology.atoms()):
            print("WARNING: obminimize output atom count mismatch")
            return None
        coords = []
        for p in new_pdb.positions:
            v = p.value_in_unit(nanometer)
            coords.append(Vec3(float(v[0]), float(v[1]), float(v[2])))
        print(f"  obminimize refined {len(coords)} atoms ({used_ff})")
        return Quantity(coords, nanometer)
    finally:
        _sh.rmtree(workdir, ignore_errors=True)


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
    from dvbfixer.ffutils import resolve_ff, print_ff_selection
    args.ff, _ff_alias, _ff_reason = resolve_ff(
        args.ff, args.input, verbose=args.verbose)
    print_ff_selection(_ff_alias, _ff_reason, args.ff)

    # .dat is optional: if --dat given, use it; otherwise auto-detect; otherwise skip
    if args.dat:
        dat_path = Path(args.dat)
    else:
        dat_path = input_path.with_suffix(".dat")

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

    pdb = PDBFile(str(input_path))
    topology = pdb.topology
    positions = pdb.positions

    # Load restraint data from .dat if available
    if args.dat and not dat_path.exists():
        print(f"Error: .dat file not found: {dat_path}", file=sys.stderr)
        sys.exit(1)

    if dat_path.exists():
        added_keys = load_dat(dat_path)
        new_atom_indices = resolve_new_atom_indices(topology, added_keys, args.verbose)

        # Also load variant overrides from .dat (saved by prepare)
        import json
        with open(dat_path) as _df:
            _dat = json.load(_df)
        for key_str, var_name in _dat.get('variant_overrides', {}).items():
            ch, rn = key_str.split(':', 1)
            if (ch, rn) not in amber_renames:
                amber_renames[(ch, rn)] = var_name
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
    print(f"\nSaved minimized structure: {output_path}")
