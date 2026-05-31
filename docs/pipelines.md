# Pipelines

[← README](../README.md)

## Quick: Full pipeline with `zbs`

```bash
# Run everything: renumber -> model -> prepare -> minimize -> protonate -> minimize
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

# 4. Fix missing atoms/residues (merges model .dat)
dvbfixer prepare modeled.pdb -v

# 5. Minimize (pass 1 — standard protonation)
dvbfixer minimize prepared.pdb -v

# 6. Set protonation state names (rename only)
dvbfixer protonate minimized.pdb --no-hydrogens -v

# 7. Minimize (pass 2 — detects AMBER names, adds correct H)
dvbfixer minimize protonated.pdb --no-solvent -v

# 8. Re-apply AMBER names (OpenMM reverts them)
dvbfixer protonate minimized2.pdb --no-hydrogens -v
```

For a PDB database file with Kabat antibody numbering and missing loops:

```bash
# All at once (recommended)
dvbfixer zbs 1HZH.pdb -v

# Or step by step with custom options
dvbfixer renumber 1HZH.pdb -v
dvbfixer model 1HZH_renum.pdb --num-loops 4 -v
dvbfixer prepare 1HZH_renum_model.pdb -v
dvbfixer minimize 1HZH_renum_model_prepared.pdb -v
dvbfixer protonate 1HZH_renum_model_prepared_minimized.pdb --no-hydrogens -v
dvbfixer minimize 1HZH_renum_model_prepared_minimized_prot.pdb -v
dvbfixer protonate 1HZH_renum_model_prepared_minimized_prot_minimized.pdb --no-hydrogens -v
```

## Glycoprotein preparation with GLYCAM

**End-to-end pipeline** (recommended for glycoproteins with PDB-named sugars):

```bash
# 1. Convert PDB glycan names + atom names to GLYCAM convention.
#    Detects glycosidic bonds from CONECT, renames sugars (NAG→UYB/0YB,
#    BMA→VMB, MAN→2MA, etc.) and glycoprotein residues (ASN→NLN,
#    SER→OLS, THR→OLT). Renames atoms (C7→C2N, O7→O2N, HO3→H3O, etc.).
dvbfixer convert crystal.pdb -o glycam.pdb -v

# 2. Add hydrogens with AMBER14+GLYCAM_06j-1 templates.
#    Auto-detects GLYCAM residues; runs add_glycam_bonds(positions=...)
#    to populate intra-residue + peptide + sugar-sugar glycosidic bonds
#    before addHydrogens. RDKit/OpenBabel H polish is skipped (GLYCAM
#    already provides correct H placement). Output: NLN/UYB/4YB/VMB
#    preserved with H atoms.
dvbfixer prepare glycam.pdb -o prep.pdb -v

# 3. Energy-minimize the whole system (protein + glycans) with
#    AMBER14+GLYCAM. Glycans relax freely (no .dat restraint); protein
#    heavy atoms strongly restrained. NLN/OLS/OLT names preserved.
dvbfixer minimize prep.pdb -o min.pdb --no-solvent -v

# 4. Assign protonation states with PROPKA3. Auto-renames NLN→ASN
#    internally for the pKa calc, then maps results back. FF auto-
#    switches to amber14+GLYCAM for the H addition.
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

# 5. Continue through the same prepare → minimize → protonate → top pipeline
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
