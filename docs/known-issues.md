# Known issues

[← README](../README.md)

- **N-terminal ASH/GLH in ACPYPE mode**: AMBER14 has no N/C-terminal protonated ASP/GLU templates (NASH/NGLH — never parameterized via RESP in any AMBER version). When `--acpype` encounters ASH or GLH at chain termini, it strips the protonation hydrogen (HD2/HE2) and uses the standard deprotonated template (NASP/NGLU). A `UserWarning` is emitted. Internal (non-terminal) ASH/GLH residues are preserved correctly.

- **Chain ID mismatch in .dat workflow**: The `.dat` file stores chain IDs from PDBFixer. If the prepared PDB is saved through a tool that reassigns chain IDs (PyMOL, VMD), the `.dat` entries won't match the new chain letters. Workaround: ensure chain IDs remain consistent between prepare and minimize steps, or manually edit the `.dat` file.

- **Hydrogen handling in minimize**: By default, existing hydrogens are kept. Use `--rebuild-h` to strip and re-add via OpenMM (needed when protonation state changes). When AMBER protonation names (GLH, HIE, CYX, etc.) are detected in the input PDB, they are passed as `variants` to `addHydrogens` to ensure correct protonation hydrogens.

- **OpenMM normalizes AMBER names**: `PDBFile` reader converts GLH→GLU, HIE→HIS, CYX→CYS. The minimize tool reads raw PDB text first to capture original names. `PDBFile.writeFile` also writes standard names, so a final protonate text-based rename is needed to restore AMBER names.

- **Pull valence checking**: The `pull` tool validates bonds before and after pulling. Pre-pull: checks valence (bond count vs element max), warns about unusual element pairs. Post-pull: checks convergence (distance vs target), bond length range, and steric clashes within the pulling residues. All checks are warnings only — they do not prevent the operation.

- **Glycoprotein minimization in `minimize`**: Default minimizes the WHOLE system (protein + glycans + ligands) with AMBER14 + GLYCAM_06j-1 + SMIRNOFF. Strip-and-splice mode (`--strip-heterogens`) is an opt-in protein-only flow with HETATM coords spliced back from the input. When the GLYCAM full-system path can't parametrize PDB-named sugars (NAG/BMA/MAN without GLYCAM templates), the tool auto-falls back to strip-and-splice and `_rigid_track_glycan_trees` does Kabsch tracking + canonical trans-amide C1/HD21 placement to preserve the glycan geometry relative to the moved protein.

- **Mixed 1-4 scaling (AMBER+GLYCAM)**: AMBER uses fudgeLJ=0.5/fudgeQQ=0.8333, GLYCAM uses 1.0/1.0. GROMACS only supports one global value. The `--acpype` flag on `top` and `--gromacs` on `transplant` solve this via ACPYPE's `[ pairs_nb ]` directive with per-pair LJ/Coulomb parameters.

- **AMBER14 has no terminal protonated ASP/GLU**: AMBER14 lacks NASH/NGLH/CASH/CGLH templates (no RESP charges were ever computed for terminal protonated ASP/GLU — a 15+ year gap). Affects both `dvbfixer top --acpype` and `dvbfixer top --ff amber --protonate`. When ASH/GLH is requested at a terminus, the protonation H is dropped, the residue is converted to standard ASP/GLU (using the existing NASP/CASP/NGLU/CGLU templates), and a `UserWarning` is emitted. HIS variants (HID/HIE/HIP) are unaffected — terminal templates exist (NHIE/CHIE etc.). CHARMM is unaffected — it uses TDB patches that combine cleanly with ASPP/GLUP.

- **Modeller terminal alignment**: `align2d` can misplace terminal gaps (e.g. matching last template residue to last target residue). This is auto-corrected by `_fix_terminal_alignment` which forces gaps to the actual N/C termini.

- **FASTA chain IDs required**: `dvbfixer model --fasta` matches sequences to PDB chains by chain ID embedded in the FASTA header. Accepted forms: `>chain_X`, `>PDBID_X` (e.g. `>1abc_A`), or bare `>X`. Sequences are NOT matched by file order. Headers without a parseable chain ID produce a clear error.

- **HIS tautomer selection**: PROPKA only predicts the overall pKa, not which nitrogen is protonated. The `--his-default` flag sets a global default (HIE or HID). For accurate per-residue tautomer assignment, use tools like MolProbity's Reduce or Schrodinger's ProtAssign.

- **Water + ion mismatch causes LINCS failure**: Prior to the `--ion-set` flag, `dvbfixer top --water` only changed the water moleculetype while keeping bundled Aqvist Na⁺/Dang Cl⁻ ions regardless of water choice. Combining OPC water with Dang Cl⁻ caused Cl⁻ to over-attract to protein cations; in a real user case a 4× trastuzumab + OPC system saw atomic pressure crash to −9000 bar in 10 ps of NPT and LINCS died at step 8027. Now `--ion-set auto` (default) picks the matched set: TIP3P→JC-TIP3P, SPC/E→JC-SPCE, TIP4P-Ew→JC-TIP4P-Ew, OPC→Li-Merz HFE-OPC. Pass `--ion-set dang-legacy` only when reproducing pre-flag runs.

- **CHARMM water restriction**: CHARMM ions (SOD/CLA/POT/CAL/MGA) are fitted to CHARMM-TIP3P. `dvbfixer top --ff charmm` only accepts `--water tip3p|spc|spce`; `--water opc|tip4p|tip4pew` is rejected at the CLI level. To use OPC water with this protein, switch to `--ff amber`. `--ion-set` is a no-op with `--ff charmm`.

- **`--acpype` mode is TIP3P-locked**: The `--acpype` pipeline (OpenMM → ParmEd → ACPYPE) hardcodes TIP3P water + AMBER14+GLYCAM ions and ignores `--water`/`--ion-set`. A future enhancement could add OPC support there; today, use the RTP-based `dvbfixer top --water opc` path if you need OPC.

## See also

- [Command index](commands/index.md)
- [Pipelines](pipelines.md)
- [BEST_PRACTICES.md](../BEST_PRACTICES.md) — additional gotchas and recipes
