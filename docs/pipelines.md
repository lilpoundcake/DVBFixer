# Pipelines

[← README](../README.md)

## Quick: Full pipeline with `zbs`

```bash
# Run everything: renumber -> model -> prepare -> minimize
# (PROPKA + MolProbity Reduce run inside prepare, 0.7.7+ — no separate protonate step)
dvbfixer zbs 6DDV.pdb -v
# -> 6DDV_zbs.pdb

# Skip terminal modeling for structures with missing N/C termini
dvbfixer zbs 6DDV.pdb --no-terminal -v

# Skip steps as needed
dvbfixer zbs input.pdb --skip-model --skip-minimize -v
```

See [`zbs`](commands/zbs.md) for all flags.

## Manual: Step by step

For a GROMACS MD output PDB file with no chain IDs:

```bash
# 1. Split chains
dvbfixer split md_output.pdb -v

# 2. Renumber (requires SEQRES)
dvbfixer renumber structure_with_seqres.pdb -v

# 3. Rebuild missing loops (writes .dat for restraints)
dvbfixer model renumbered.pdb -v

# 4. Fix missing atoms/residues (merges model .dat); PROPKA + MolProbity Reduce
#    run inside this step by default, picking AMBER protonation variants and
#    placing final H accordingly
dvbfixer prepare modeled.pdb -v

# 5. Minimize the whole system; preserves AMBER variant names on write
dvbfixer minimize prepared.pdb -v
```

Standalone `dvbfixer protonate` is a post-hoc re-protonation tool, not a
required pipeline stage — use it if you need to re-run PROPKA/Reduce on an
already-prepared PDB (e.g. to switch pH after the fact):

```bash
dvbfixer protonate minimized.pdb -o reprotonated.pdb --ph 6.5 -v
```

An older 7-step pattern (`minimize` → `protonate --no-hydrogens` (rename
only) → `minimize` → `protonate --no-hydrogens` (restore names)) was
dropped in 0.7.7: `--no-hydrogens` leaves existing H atoms in the wrong
positions relative to the newly-assigned variant names, and it's no longer
needed now that `prepare` picks variants up front and `minimize` preserves
AMBER names on its own.

For a PDB database file with Kabat antibody numbering and missing loops:

```bash
# All at once (recommended)
dvbfixer zbs 1HZH.pdb -v

# Or step by step with custom options
dvbfixer renumber 1HZH.pdb -v
dvbfixer model 1HZH_renum.pdb --num-loops 4 -v
dvbfixer prepare 1HZH_renum_model.pdb -v
dvbfixer minimize 1HZH_renum_model_prepared.pdb -v
```

## Glycoprotein preparation with GLYCAM

**End-to-end pipeline** (recommended for glycoproteins with PDB-named sugars):

```bash
# 1. Convert PDB glycan names + atom names to GLYCAM convention.
#    Detects glycosidic bonds from CONECT, renames sugars (NAG→UYB/0YB,
#    BMA→VMB, MAN→2MA, etc.) and glycoprotein residues (ASN→NLN,
#    SER→OLS, THR→OLT). Renames atoms (C7→C2N, O7→O2N, HO3→H3O, etc.).
dvbfixer convert crystal.pdb -o glycam.pdb -v

# 2. Add hydrogens with AMBER14+GLYCAM_06j-1 templates. Also runs PROPKA +
#    MolProbity Reduce by default, picking AMBER protonation variants
#    (ASH/GLH/HIP/LYN/CYM/CYX, HIS tautomer, ASN/GLN flips) up front.
#    Auto-detects GLYCAM residues; runs add_glycam_bonds(positions=...)
#    to populate intra-residue + peptide + sugar-sugar glycosidic bonds
#    before addHydrogens. RDKit/OpenBabel H polish is skipped (GLYCAM
#    already provides correct H placement). Output: NLN/UYB/4YB/VMB
#    preserved with H atoms.
dvbfixer prepare glycam.pdb -o prep.pdb -v

# 3. Energy-minimize the whole system (protein + glycans) with
#    AMBER14+GLYCAM. Glycan heavy atoms get the same weak restraint
#    tier as newly-modeled backbone atoms (0.7.10 — prevents a
#    multi-residue glycan tree drifting off its covalent anchor);
#    protein heavy atoms strongly restrained. NLN/OLS/OLT names
#    preserved.
dvbfixer minimize prep.pdb -o min.pdb --no-solvent -v

# 4. (Optional) Re-run PROPKA3 standalone if you want to re-derive
#    protonation state after minimize's geometry changes — prepare
#    (step 2) already assigned variants once, so this is a refinement,
#    not a required step. Auto-renames NLN→ASN internally for the pKa
#    calc, then maps results back. FF auto-switches to amber14+GLYCAM
#    for the H addition.
dvbfixer protonate min.pdb -o prot.pdb -v

# 5. Generate GROMACS topology (AMBER14+GLYCAM, mixed 1-4 scaling
#    via [ pairs_nb ]).
dvbfixer top prot.pdb --acpype -o gmx/ -v
# -> gmx/topol.top, gmx/prot.gro, gmx/posre_prot.itp
```

GLYCAM names (NLN/UYB/4YB/VMB/0YB/0fA/0LA/2MA/...) and protonation
variants (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN) survive end-to-end through
all five steps. After the pipeline, run a normal GROMACS workflow:
`gmx editconf` → `gmx solvate` → `gmx genion` → `gmx grompp` → `gmx mdrun`.

### Switching the same structure to CHARMM36

Use [`dvbfixer convert --to-charmm`](commands/convert.md) to reverse the
GLYCAM naming back to standard PDB / CHARMM-compatible names, then use
the RTP CHARMM path in [`top`](commands/top.md):

```bash
# After step 1 (or any later step with GLYCAM names), reverse to CHARMM
dvbfixer convert glycam.pdb --to-charmm -o charmm.pdb -v
# -> charmm.pdb has NAG/NDG/BMA/MAN/GAL/FUL/SIA (no GLYCAM codes),
#    ASN (no NLN), and standard PDB atom names

# CHARMM36 GROMACS topology (modular .itp files)
dvbfixer top charmm.pdb --ff charmm -o gmx_charmm/topol.top -v
```

The same source structure can drive both AMBER14+GLYCAM (via `--acpype`)
and CHARMM36 (via `--ff charmm`) — pick the FF by which renaming pass
you use.

### Alternative: GLYCAM-Web workflow (for adding glycans from scratch)

When the input has no glycans and you need to add them via the GLYCAM-Web
glycan builder:

```bash
# 1. Prepare the deglycosylated protein
dvbfixer zbs antibody.pdb -v

# 2. Extract glycosylation site residues into donor.pdb (PyMOL select)
# 3. Submit donor.pdb to GLYCAM-Web -> download glycam_output.pdb

# 4. Transplant GLYCAM output back, with relaxation
dvbfixer transplant antibody_zbs.pdb --donor donor.pdb \
    --graft glycam_output.pdb --relax -v

# 5. Continue through prepare -> minimize -> top (prepare already runs
#    PROPKA + Reduce internally; standalone protonate below is an
#    optional re-derivation after minimize's geometry changes, not a
#    required step)
dvbfixer prepare   antibody_zbs_transplant_relaxed.pdb -o prep.pdb -v
dvbfixer minimize  prep.pdb -o min.pdb --no-solvent -v
dvbfixer protonate min.pdb -o prot.pdb -v
dvbfixer top       prot.pdb --acpype -o gmx/ -v
```

## GROMACS topology generation

```bash
# RTP-based (fast, modular .itp output — no FF dir needed)
dvbfixer top input.pdb --ff amber                   # AMBER + TIP3P + JC-TIP3P ions (default)
dvbfixer top input.pdb --ff amber --water opc       # AMBER + OPC + Li-Merz HFE-OPC ions
dvbfixer top input.pdb --ff amber --water tip4pew   # AMBER + TIP4P-Ew + JC-TIP4P-Ew ions
dvbfixer top input.pdb --ff charmm                  # CHARMM + CHARMM-TIP3P + CHARMM ions

# Override ion parameter set (AMBER only)
dvbfixer top input.pdb --water opc --ion-set lm-iod-opc      # IOD-fit instead of HFE
dvbfixer top input.pdb --water tip3p --ion-set dang-legacy   # bundled pre-2008 values

# Glycolipid from CHARMM-GUI (auto-detected, CHARMM36 only)
dvbfixer top glycolipid_charmm.pdb --ff charmm -o gmx_top/

# ACPYPE-based (proteins + GLYCAM glycans, handles mixed 1-4 scaling)
# Note: --acpype hardcodes TIP3P; --water/--ion-set are ignored.
dvbfixer top input.pdb --acpype
```

Ion atom-type set is auto-selected from `--water` and covers Na⁺/K⁺/Cl⁻/Ca²⁺/Mg²⁺/Zn²⁺.
See [`top`](commands/top.md#water-matched-ions-amber) for the full mapping.

## Glycan conformational analysis

```bash
# Cluster glycan conformations from MD trajectory
dvbfixer cluster topol.tpr md.xtc --plot -v
# -> md_representatives.pdb, md_summary.json, interactive HTML plots
```

## Homology modeling (antibody engineering)

```bash
# Combine Fv from Fab template + constant domains from IgG template
dvbfixer homology target.fasta --template fab.pdb --template igg.pdb --antibody -v

# With full pipeline
dvbfixer homology target.fasta --template fab.pdb --template igg.pdb --minimize -v
```

## Small molecule parametrization

```bash
# Parametrize a buffer component for GROMACS
dvbfixer parametrize acetate.pdb -n ACET --net-charge -1
# -> ACET.itp, ACET.gro, posre_ACET.itp
```

## See also

- [Command index](commands/index.md)
- [Known issues](known-issues.md)
- [BEST_PRACTICES.md](../BEST_PRACTICES.md) — additional recipes and gotchas
