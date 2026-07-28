# Changelog

All notable changes to dvbfixer are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Backfilled from git history — commits before v0.3.0 are grouped by
feature area rather than by strict release. Older entries are
best-effort summaries; consult `git log` for exact provenance.

## [Unreleased]

### Documented

- **propka 3.5.1 × Python 3.14 incompatibility.** propka reads
  `self.__annotations__` (instance-level) inside its `Parameters`
  dataclass; Python 3.14's PEP 649/749 lazy-annotations change makes
  instance `__annotations__` raise `AttributeError` instead of falling
  through to the class dict, crashing `dvbfixer protonate` (and the
  PROPKA step inside `prepare` / `zbs`) at
  `propka.parameters.parse_line` with
  `AttributeError: 'Parameters' object has no attribute '__annotations__'`.
  `environment.yml` already pins `python >=3.11,<3.14` to keep the env
  on 3.12/3.13 where propka works — do not loosen it. Documented in
  [CLAUDE.md](../CLAUDE.md) and [docs/known-issues.md](known-issues.md).

- **micromamba env creation fails on macOS Docker host bind mounts.**
  When `MAMBA_ROOT_PREFIX` lives under a host bind mount such as
  `/home/agent` (a `fakeowner` / VirtioFS bind of macOS `/Users`),
  `micromamba create` aborts during `Linking 'ncurses'` with
  `filesystem error: cannot copy symlink: Invalid argument` on the
  case-variant terminfo pair `share/terminfo/32/2621A` vs `2621a`.
  The mount is case-insensitive and rejects libmamba's `copy_symlink`
  for these case-colliding entries; `always_copy` / `--copy` do not
  help (copy mode still recreates in-package symlinks rather than
  dereferencing them), and no micromamba flag dereferences or skips
  them. Workaround: create the env + package cache on the container's
  native overlay filesystem (`export MAMBA_ROOT_PREFIX=/opt/mamba`).
  Documented in [docs/known-issues.md](known-issues.md) and
  [docs/installation.md](installation.md).

## [0.7.7] — 2026-07

### Added

- **PROPKA + MolProbity Reduce now run INSIDE legacy prepare.**
  Prior to 0.7.7 the default pipeline (`dvbfixer prepare` /
  `dvbfixer zbs`) built its `Modeller.addHydrogens(variants=[...])`
  list only from user `--mutate` overrides + HIS HD1/HE2 presence.
  Neither PROPKA nor Reduce was invoked, so pKa-driven ASH / GLH /
  HIP / LYN / CYM variants and per-residue HIS tautomer picks
  never appeared in output. Only the opt-in `--backend tleap-reduce`
  ran them. Fixed via new helper
  `dvbfixer.prepare.pipeline._run_propka_reduce_variants` which
  mirrors the tleap-reduce backend's PROPKA-decide + Reduce-tautomer
  overlay logic (fixed in 0.7.4) and is called by
  `dvbfixer.prepare` before addHydrogens. PROPKA drives ASP/GLU/
  LYS/CYS/HIS→HIP; Reduce fills in HID vs HIE for neutral HIS; user
  `--mutate` overrides win on collision; CONECT-detected SS pairs
  force CYX regardless of PROPKA's CYS pKa.

- **New CLI flags on `prepare` and `zbs`**:
  - `--propka / --no-propka` (default: on).
  - `--protassign / --no-protassign` (default: on).
  - `--his-default {HIE,HID}` (default: HIE) — fallback tautomer
    when both PROPKA and Reduce are ambiguous.
  - `--cys-ss-pka` (default: 8.0) — PROPKA pKa cutoff for CYX.

  zbs propagates the flags to prepare.

- New tests `tests/test_prepare_propka_integration.py` (6 tests).

### Changed

- **`zbs` docstring** updated: "renumber → model → prepare
  (with PROPKA + Reduce) → minimize". The stale "protonate step"
  language from pre-0.7.0 finally removed.

- **`--skip-protonate` on zbs is deprecated**. Kept for backward
  compat: now maps to `--no-propka --no-protassign` on prepare
  (with a stderr warning). Standalone `dvbfixer protonate` command
  is unchanged and can still be used as a post-hoc re-protonation
  tool (e.g. to switch pH on an already-prepared PDB).

## [0.7.6] — 2026-07

### Changed

- **`detect_ff_from_pdb` no longer surrenders on PDB-standard sugars.**
  Previously an input with IUPAC/PDB sugar codes (BGL / BMA / AMA /
  NAG / …) fell back to `amber` with a WARNING telling the user to
  run `dvbfixer convert` manually first. Now returns `amber+glycam`
  with a reason string; downstream prepare + minimize invoke
  `convert_to_glycam` under the hood so the tool Just Works on
  downloaded RCSB glycoproteins.

- **`zbs` no longer force-overrides minimize `--ff` to
  `amber14-all.xml`.** Previously when the user passed `--ff auto`,
  zbs auto-detected the FF for prepare but forced minimize back to
  plain `amber14`, breaking glycoprotein pipelines. Now the same
  `--ff` value flows to both steps; each runs its own auto-detection
  independently.

### Added

- **CHARMM + PDB-sugars hybrid path.** When `--ff charmm` is
  requested (or auto-selected) AND the input has PDB-standard sugar
  names `charmm36.xml` can't parametrise, prepare + minimize
  process the run under `amber+glycam` (which has sugar templates)
  and rewrite output residue names to CHARMM convention (`BGLCNA`,
  `BMAN`, …) via `convert_to_charmm` after the pipeline completes.
  Coordinates and topology unchanged; only sugar residue names
  differ between amber and charmm output.

- **`--atom-naming {gromacs,standard}` CLI flag** on `prepare`,
  `minimize`, `protonate`, and `zbs`. Default `gromacs` preserves
  the current behaviour (GROMACS amber99sb-ildn shifts: HB3→HB1
  keeping HB2, HZ3→HZ1 on LYN, O/OXT→OC2/OC1, H→HN for CHARMM).
  Pass `standard` to keep IUPAC/AMBER-native names (HB2/HB3,
  HZ1/HZ2/HZ3, O/OXT, plain H) — matches ff14SB / VMD / most PDB
  downloaders.

- New helper `dvbfixer.ffutils.has_pdb_standard_sugars(pdb_path)`
  — text-level scan for the auto-convert trigger.

- New tests `tests/test_ff_auto_detect_sugar.py` covering
  auto-detect on glycoproteins, amber→amber+glycam upgrade, and
  the has_pdb_standard_sugars helper.

## [0.7.5] — 2026-07

### Changed

- **Default prep backend flipped from `tleap-reduce` back to `legacy`
  (Modeller+PDBFixer).** The tleap-reduce backend hard-fails on
  glycoproteins, ligands, PTMs, and any covalent-HETATM input
  (tleap has no template for those). Legacy prep handles them via
  the existing `build_glycam_system` + `--parametrize-ligands`
  paths. The chirality invariant that motivated the switch
  originally is now enforced downstream in
  `minimize/pipeline.py`'s post-phase-2 unconditional force-reflect
  fallback (added in 0.7.4), which is prep-backend-agnostic — so
  legacy prep's PDBFixer.addMissingAtoms D-Cα risk is neutralised
  where it matters. Applies to `prepare`, `protonate`, and
  transitively to `zbs`. `tleap-reduce` remains fully functional
  via explicit `--backend tleap-reduce` for pure-protein inputs.

### Added

- Comprehensive integration test suite covering every subcommand
  (except `homology`) against real PDB inputs across the input
  classes dvbfixer supports: pure protein, antibody, glycoprotein,
  protein+ligand, multi-MODEL, and the curated `test/shit/`
  regression set.

## [0.7.4] — 2026-07

### Fixed

- **PROPKA-driven variant assignment on the tleap+reduce backend was
  broken end-to-end.** The 0.7.3 code stuffed
  `propka_dict[(chain, resnum)] = (restype, pka)` and iterated the PDB.
  PROPKA emits multiple `Group` records per residue (side-chain
  acid/base + terminal N+/C-), so last-write-wins silently clobbered
  side-chain pKas with terminal pKas — ASH/GLH/LYN were correct only
  when PROPKA's iteration order happened to land the right entry
  last. HIS→HIP was NEVER assigned because HIP was inferred from
  HD1+HE2 atom presence, but `reduce -build` picks a single tautomer.
  Rewrote to route through `dvbfixer.protonate.decide_protonation`
  (group-type-filtered, `(chain, resnum, icode)`-keyed) then overlay
  Reduce's HID/HIE tautomer for neutral HIS. `_patch_variant_hydrogens`
  gained an HIP branch that places the missing imidazole H. Deleted
  dead `_STD_PKA` and `_classify_variant`.

- **AMBER variant residue names silently broke minimize.** OpenMM's
  `PDBFile._standardResidues` set covers only the 20 canonical AAs;
  LYN/ASH/GLH/CYX/CYM/HID/HIE/HIP load with zero intra-residue bonds
  → `createSystem` fails "no bonds between its atoms". Additionally
  `Modeller.addHydrogens` copies input bonds but does NOT rebuild
  missing ones from templates, and only ADDS missing H, never
  REMOVES extras. Fix in `minimize/pipeline.py`:
  1. `text_rename_variants_to_parent` runs BEFORE `PDBFile.load` so
     OpenMM infers proper intra-residue bonds from parent templates.
  2. When `amber_renames` is non-empty, force strip-H + variant-aware
     addHydrogens rebuild.
  3. Cap the pH passed to `addHydrogens` at 9.99 to bypass the
     `hydrogens.xml` `maxph="10.0"` HZ3 gate that broke terminal LYS
     at high pH.
  4. `pdbutils/inference.py::_apply_filter` keeps intra-residue
     bonds for variant-named residues in emitted CONECT.

- **Disulfide bonds dropped during minimize.**
  `_drop_spurious_inter_aa_bonds` treated every inter-residue bond
  between two protein residues that wasn't a C-N peptide bond as
  spurious. After `text_rename_variants_to_parent` folded CYX → CYS,
  SG-SG disulfides were between two "CYS"-named residues → dropped.
  Added an SG-SG exception for the CYS family (CYS/CYX/CYM).

- **D-Cα residues occasionally survived minimize.** The prior
  WARN-only design left D residues in the output when the FF's local
  minimum for a residue genuinely sat on the D side (rare, exotic
  packing). Replaced with two-tier enforcement:
  1. Bounded reflect-and-re-minimize loop (max 3 iterations).
  2. If any residue is still D after the loop, do a final
     unconditional `fix_ca_chirality` reflection and skip any
     follow-up minimize. Since `fix_ca_chirality` mirrors the whole
     sidechain through the CA-N-C plane, all internal bond lengths
     and angles are preserved (CA-CB ≈ 0.154 nm, CB-HB ≈ 0.109 nm);
     only CB's position relative to backbone neighbours changes. The
     chirality invariant is now non-negotiable: zero D-Cα in output.

### Added

- `tests/test_prep_backend_variants.py` — 16 focused tests for the
  PROPKA + variant-H paths, covering ASH, GLH, LYN, HIP, HID, HIE,
  CYX-from-SS, terminal skip.

## [0.7.0] — 2026-07

### Fixed
- **Full GROMACS-canonical atom naming on every user-visible PDB
  output.** Prior 0.6.5 fix only handled LYN NZ hydrogens. Full
  sub-agent survey (Jul 2026) identified pervasive naming
  differences OpenMM ff14SB (IUPAC) vs GROMACS amber99sb-ildn
  (older AMBER numbering) that GROMACS's own `aminoacids.arn`
  does NOT rewrite at `pdb2gmx` read time. User reported
  `pdb2gmx -ff amber99sb-ildn` errors on dvbfixer output for
  many residues. Now handled by `apply_variants_to_pdb_text`:
  - **Methylene H shift** on all β/γ/δ/ε methylenes: `HB3 → HB1`
    (HB2 stays; result `{HB1, HB2}` matches GROMACS). Same for
    HG, HD, HE. Applied to every AA that has the corresponding
    methylene (LYS/ARG/PRO get all four; SER/CYS/etc get just
    HB; GLN/GLU/MET get β+γ; etc.).
  - **GLY α-methylene**: `HA3 → HA1`
  - **ILE γ-methylene**: `HG13 → HG11`
  - **LYN NZ H**: `HZ3 → HZ1` (already in 0.6.5, kept)
  - **ACE/NME cap H**: `HH31/HH32/HH33 → H1/H2/H3` (amber14→amber19/GROMACS)
  - **N-terminal residues**: `H → H1` on first residue of each
    protein chain (only if a bare H exists; skip if Modeller
    already placed H1/H2/H3 for a charged terminus).
  - **N-terminal proline (NPRO)**: `H3 → H1` on the ring N.
  - **C-terminal residues**: `O → OC2`, `OXT → OC1` on last
    residue of each protein chain. Applied atomically.
  - **CHARMM output**: universal backbone `H → HN` on every
    protein residue when `target_ff='charmm'`.
  - **DNA / RNA**: `H2' / H2'' → H2'1 / H2'2`, `H5' / H5'' →
    H5'1 / H5'2`, `HO'2 → HO2'` (RNA). Map exposed but
    end-to-end testing deferred pending column-perfect
    fixture handling.

  All renames are single-source (no shift-pair) so they're
  naturally idempotent — running `apply_variants_to_pdb_text`
  twice produces the same output.

### Renamed
- `apply_variants_to_pdb_text` parameter `include_gromacs_lyn` →
  `include_gromacs_shifts`. Old name kept as a deprecated alias
  for one release.

## [0.6.6] — 2026-07

### Added
- **`dvbfixer model --strip-heterogens`** — remove HETATM records
  (ligands, sugars, ions, cofactors) before Modeller runs. Off by
  default (Modeller usually benefits from heterogen context for
  loop refinement) but useful when heterogen geometry causes loop
  artifacts. Waters are preserved when `--keep-water` is also
  passed. Orphan CONECT records referencing dropped serials are
  filtered too.
- **zbs propagation**: when `dvbfixer zbs --strip-heterogens` is
  invoked, the flag is now forwarded to the model step (previously
  it was only forwarded to prepare and both minimize passes, so
  Modeller saw all heterogens even when the user had explicitly
  asked to strip them for the pipeline).

## [0.6.5] — 2026-07

### Fixed
- **AMBER protonation variants (GLH/ASH/HID/HIE/HIP/CYX/CYM)
  finally survive `dvbfixer zbs`.** User verified failure across
  0.6.3 and 0.6.4 with multiple fixtures (2VLP, 1TM7). Empirical
  test on OpenMM 8.5.1 showed `PDBFile` canonicalises variant
  residue names on both read AND write — the topology-level
  `_rename_variants_to_parent → addHydrogens → _restore_variants_in_topology`
  dance was a no-op because the "restore" walked a topology that
  only ever had canonical names.
  Fix: text-level rewrite of every user-visible output PDB using
  `amber_renames` (populated from the RAW input text via
  `scan_variant_names`, which OpenMM never sees). Verified on
  2VLP at pH 5: GLH=48, ASH=13, HIP=126 all in final output.
- **`dvbfixer zbs --rename` no longer destroys protonate's work.**
  The flag was being propagated to minimize step 2 (which runs
  AFTER protonate); canonicalising there threw away every AMBER
  variant PROPKA assigned. Removed the propagation — `--rename`
  now only applies to prepare + minimize step 1 (both run BEFORE
  protonate), where its "strip non-canonical names" behaviour is
  actually useful.

### Added
- **New shared FF-name module `dvbfixer.ffutils.ff_names`** —
  central place for AMBER↔CHARMM protonation-variant residue
  and atom-name maps + a text-level PDB rewrite primitive
  (`apply_variants_to_pdb_text`). Reused by `prepare`,
  `minimize`, `protonate`. `convert` (glycam.py) is unchanged;
  a future PR will migrate it to use the shared module too.

## [0.6.4] — 2026-07

Three bugs surfaced by a `dvbfixer zbs test/stereo_bug/1TM1.pdb
--no-solvent --strip-heterogens` reproduction.

### Fixed
- **GLH/ASH/HIE/HID/HIP/CYX STILL canonicalized after 0.6.3.**
  The 0.6.3 fix wired `amber_renames` into Branches 2 and 3 of
  minimize's `keep_h=True` path but missed **Branch 1** (the
  `if n_missing or n_terminals:` branch, which fires on virtually
  every real input). Branch 1 called `PDBFixer.addMissingHydrogens`
  — a hard-rule violation from CLAUDE.md, because PDBFixer's
  `_describeVariant` only recognises standard PDB names and
  rewrites variant residues to canonical HIS/ASP/GLU/CYS before
  placing H. Fixed by making Branch 1 follow the same
  `_rename → Modeller.addHydrogens(variants=...) → _restore`
  pattern as Branches 2/3.
- **`fix_ca_chirality` didn't run on inputs with all atoms
  present.** SER224 in 1TM1 came out as D-chirality after zbs.
  Root cause: `fix_ca_chirality` was gated by
  `if n_missing or n_terminals:` — if PDBFixer detected no missing
  atoms, chirality was never checked. Also: reflecting only CB
  (leaving OG/CG/etc on the wrong side) broke sidechain geometry.
  Fixed: `fix_ca_chirality` now reflects EVERY sidechain atom
  through the CA-N-C plane (not just CB), and runs unconditionally
  at the end of minimize before `createSystem`.
- **"1 C-O bond too many" createSystem failure.** Over-eager
  CONECT inference wired backbone amide H atoms to CA/HA/O
  simultaneously, exceeding hydrogen valence and making residue
  templates fail to match. Added
  `minimize.pipeline._drop_spurious_inter_aa_bonds` which:
  (1) drops non-peptide inter-residue bonds between two standard
  AAs, and (2) enforces H valence = 1 by dropping every second
  bond involving each H atom. Called after every `Modeller(...)`
  reconstruction and before every `createSystem`.

### Known limitation
`dvbfixer zbs --parametrize-ligands` on inputs with badly
corrupted CONECT records (like 1TM1) still fails downstream
because OpenMM's `Modeller.addHydrogens` hits a zero-division on
degenerate geometry after the H-valence filter reshapes the bond
graph. Investigation continues; use `--strip-heterogens` in the
meantime.

## [0.6.3] — 2026-07

### Fixed
- **zbs canonicalized all AMBER variants except LYN.** User reported
  that after `dvbfixer zbs`, output PDB had GLU (was GLH), ASP (was
  ASH), HIS (was HIE/HID/HIP), CYS (was CYX) — only LYN survived.
  Root cause: `minimize/pipeline.py:_rename_variants_to_parent`
  walked the topology, but OpenMM's intermediate
  `PDBFile.writeFile` normalizes HID/HIE/HIP/ASH/GLH/CYX/CYM to
  canonical names (LYN is the only variant OpenMM doesn't
  normalize). The `_saved` dict therefore captured only LYN, and
  the restore step only put LYN back. Fixed by making
  `_rename_variants_to_parent` consult `amber_renames` (the
  raw-text-derived variant dict populated by `_read_amber_renames`
  before any OpenMM parsing) as the ground truth.

### Changed
- **LYN output PDBs now use GROMACS-compatible `HZ1` + `HZ2`.** User
  reported that feeding dvbfixer's PDB into `pdb2gmx -ff amber99sb-ildn`
  failed on LYN residues because ff14SB's LYN template names its two
  NZ hydrogens `HZ2 + HZ3`, but GROMACS amber99sb-ildn's
  `aminoacids.hdb` H-add rule expects `HZ1 + HZ2`. Sub-agent survey
  (Jul 2026) confirmed LYN is the ONLY protein atom-name difference
  between ff14SB and GROMACS amber99sb-ildn for standard residues +
  HID/HIE/HIP/CYX/CYM/ASH/GLH/LYN. Added
  `ffutils.variants.rename_lyn_hz_for_gromacs` (topology-level) and
  `rename_lyn_hz_for_gromacs_in_pdb_text` (file-level) helpers.
  Wired into the final user-visible PDB write step of `prepare`,
  `minimize`, and `protonate` (both H-adding and text-rename paths).
  Intermediate temp writes remain ff14SB so OpenMM's `createSystem`
  still matches the LYN template.

## [0.6.2] — 2026-07

### Fixed
- **Cα chirality inversion on PDBFixer sidechain rebuild.** User
  reported `dvbfixer zbs` producing `B/VAL98:CA` as a D-amino acid
  (triple product (N×C)·CB = -0.0017 nm³). Root cause:
  `PDBFixer.addMissingAtoms()` rebuilds missing sidechain heavy
  atoms from ideal AMBER templates; on branched-Cβ residues
  (VAL/ILE/THR) the template alignment can pick the D face.
  Nothing in the pipeline detected this. Added
  `dvbfixer.ffutils.geometry.fix_ca_chirality` which reflects any
  D-configured CB through the CA-N-C plane. Wired into `prepare`
  and `minimize` immediately after `addMissingAtoms`. Diagnose
  now reports 0 chirality findings after zbs.

### Docs
- `docs/known-issues.md` — new sections documenting the PDBFixer
  chirality behaviour and the LYN hydrogen naming FF-vs-hydrogens.xml
  quirk (the latter is not a bug — ff14SB template ground truth is
  `HZ2 + HZ3`, dvbfixer output already matches).
- `ffutils.variants.fix_lyn_hz_naming` docstring polish making it
  explicit that ff14SB is authoritative.

## [0.6.1] — 2026-07

### Changed
- **`diagnose`: default clash overlap thresholds now match ChimeraX**
  (WARN 0.6 Å / ERROR 0.9 Å), not MolProbity (0.4 / 0.5 Å). User
  reported that BioLuminate and ChimeraX report no clashes on inputs
  where diagnose reported many — MolProbity's clashscore floor is a
  strict-validation setting, not a "does this look OK" setting.

### Added
- **`--clash-mode {chimerax,molprobity,bioluminate}`** preset flag.
  - `chimerax` — 0.6 / 0.9 Å (new default; matches ChimeraX `clashes`)
  - `molprobity` — 0.4 / 0.5 Å (restores pre-0.6.1 strict floor)
  - `bioluminate` — 0.75 / 1.0 Å (matches BioLuminate "Bad" / "Ugly")
- **`--clash-cutoff WARN,ERROR`** escape hatch for explicit tuning
  (overrides `--clash-mode`). E.g. `--clash-cutoff 0.35,0.45` for
  extra-strict validation.

## [0.6.0] — 2026-07

Follow-up on a cross-validated assessment (two independent agents: one
code-review, one functional test on real PDB fixtures). Addresses the
noisiest failure modes plus adds test coverage for previously
untested code paths.

### Added
- **Multi-MODEL detection.** Multi-model PDBs used to silently produce
  meaningless reports (chain breaks across MODEL boundaries, spurious
  valence violations from CONECT superposition). Now diagnose extracts
  MODEL 1 to a temp file, emits a WARNING banner, and reports on
  MODEL 1 only.
- **SEQRES-vs-ATOM check.** Terminal truncations (e.g. Fab His-tag in
  SEQRES but not resolved in ATOM) are now flagged as WARNING. PDBFixer's
  `findMissingResidues` only catches internal gaps.
- **Disulfide-geometry check.** SS bond length (2.05 ± 0.10 Å),
  Cα-Cα distance (5.5 – 7.0 Å), CB-SG-SG-CB dihedral (χ_ss = 60 – 120°)
  reported as `disulfide_geometry` findings.
- **`--format json` output.** Machine-readable findings + summary
  for CI gating and scripting.
- **`-v/--verbose` implementation.** Prints per-check-family timing
  and finding counts to stderr.

### Fixed
- **Waters no longer generate massive chain-break noise.** Every
  crystallographic water used to trigger a "chain break after HOH..."
  WARNING. Waters (HOH/WAT/SOL/TIP3/TIP4/TIP5/SPC/SPCE/DOD/H2O) are
  now excluded from chain-break AND steric checks by default. Add
  `--include-water` to keep the old behaviour.
- **Heavy-atom H-bonds no longer flagged as clashes.** Structures
  without explicit hydrogens (unprotonated inputs) had every salt
  bridge and backbone amide interaction reported as an ERROR clash
  ("LYS:NZ clashes with GLU:OE1 — overlap 0.4 Å"). Any pair of
  H-bond-capable heavy atoms (N/O/S/F ↔ N/O/S/F) at 2.5 – 3.4 Å is
  now skipped as a heavy-atom H-bond, matching CCP4's envelope.
- **Intra-residue vdW pairs skipped from clash detection.** Anything
  within the same residue is FF-template geometry — 1-4 exclusion
  already covered most cases, but glycan / non-standard-AA close
  contacts within a single residue no longer surface as clashes.

### Test coverage
Added 16 new tests covering previously-untested paths:
- Dihedral computation (trans-peptide → 180°, cis → 0°)
- Cα chirality (L-Ala pass, D-Ala flagged)
- Valence (5-bond C flagged, 4-bond OK)
- Bond length (canonical, sub-WARN stretched, above-ERROR stretched)
- Disulfide geometry (canonical, stretched, non-bonded pair)
- Multi-model input handling
- Water suppression (both directions of `--include-water`)
- JSON output format

## [0.5.2] — 2026-07

### Fixed
- **`diagnose`: hydrogen bonds no longer reported as clashes.** Every
  backbone N-H..O=C amide interaction, sidechain O-H..O, and
  ARG/LYS..carboxylate salt bridge in a well-folded structure showed
  up as an ERROR ("O clashes with H — overlap 0.7 Å"). MolProbity's
  ``probe`` handles these by classifying them as ``hbond`` contacts;
  our Python engine now applies the same filter — a pair is
  suppressed when one atom is a polar H (covalently bonded to N/O/S),
  the other is an acceptor (N/O/S/F), and their distance falls in
  the 1.4 – 2.6 Å H-bond envelope. Fixes the "α-helix reported as
  15 ERRORs" complaint on 0.5.1.

## [0.5.1] — 2026-07

### Fixed
- **`diagnose`: bond exclusion now covers HETATMs and non-standard
  residues.** `steric._build_bond_exclusion` previously used only
  `topology.bonds()`, which is empty for every glycan, ligand,
  AMBER protonation variant (HIE / CYX / LYN / ASH / GLH), and
  CHARMM variant (HSD / HSE / HSP) that OpenMM's PDBFile can't
  match to a standard template. Their real covalent bonds were
  reported as ERROR clashes (bonded C-C at 1.53 Å read as 1.87 Å
  overlap). Fixed by augmenting the topology bond set with a
  distance-inferred graph (heavy-heavy ≤ 1.9 Å, X-H ≤ 1.3 Å,
  S-S ≤ 2.25 Å) built from the same cKDTree used for clash search.

### Changed
- **`diagnose`: thresholds recalibrated for lower false-positive
  rate.** Matches BioLuminate's Protein Report and MolProbity's
  official clashscore floor.
  - Bond-length WARNING at 20 % deviation (was 10 %); ERROR at 50 %
    (was 30 %). The SER HG-on-OXT bug case sits at 75 %, so genuine
    breakage is still caught.
  - Steric clash WARNING at 0.4 Å overlap (was 0.2 Å) — matches
    MolProbity's official clashscore threshold. ERROR unchanged at
    0.5 Å.

## [0.5.0] — 2026-07

### Added
- **`dvbfixer diagnose`** — new standalone subcommand. Inspects a PDB
  file and emits a plain-text per-residue findings report; **never
  mutates the input**. Inspired by BioLuminate's Protein Report
  widget and MolProbity's all-atom validation. Three check families:
  - **Structural integrity** — missing atoms / residues / terminals
    (via PDBFixer), coincident atoms + misplaced hydrogens (reusing
    the 0.4.1 / 0.4.2 detection extracted into shared
    `ffutils.geometry.detect_coincident_atoms` +
    `detect_misplaced_hydrogens`), altLoc conflicts, chain breaks,
    insertion codes.
  - **Chemistry / bond geometry** — valence violations (reusing
    `pull.MAX_BONDS`), bond-length outliers (atom-name-aware canonical
    table distinguishes backbone C=O 1.23 Å from sidechain C-O
    1.43 Å, backbone C-N 1.33 Å from sidechain C-N 1.47 Å), cis
    peptides (INFO for cis-PRO, WARNING elsewhere), non-planar
    amides, Cα chirality via triple-product sign check.
  - **Steric analysis** — all-atom clash detection. Python engine
    (scipy cKDTree + MolProbity-standard vdW radii; excludes 1-2,
    1-3, and 1-4 bonded pairs) by default; shells to MolProbity's
    `probe` binary if it's on PATH. WARNING for ≥ 0.2 Å overlap,
    ERROR for ≥ 0.5 Å.
  - 3-tier ERROR / WARNING / INFO classification; exit-code convention
    (0 clean, 1 has ERROR, 2 I/O error) makes it usable in shell
    gates: `dvbfixer diagnose input.pdb --severity ERROR && dvbfixer prepare …`.
  - `--only {all,structural,chemistry,steric}` restricts categories;
    `--severity` filters minimum severity; `-o` writes report to file.
  - Documentation at `docs/commands/diagnose.md`; auto-generated
    reference at `docs/reference/diagnose.md`.

### Refactor
- **`ffutils.geometry`** — extracted `detect_misplaced_hydrogens` and
  `detect_coincident_atoms` as pure detection helpers (return findings
  as new `MisplacedHydrogen` / `CoincidentAtoms` dataclasses).
  `repair_misplaced_hydrogens` now composes `detect + place`, and
  `prepare.pipeline`'s pre-strip and the new `diagnose.structural`
  share the same detection code.

## [0.4.2] — 2026-07

### Fixed
- **`prepare` / `protonate` / `minimize`** — post-`addHydrogens` geometry
  sanity check. The 0.4.1 pre-strip only caught inputs where the (H,
  heavy-atom) pair was already coincident. The real failure mode — HG
  MISSING and OXT PRESENT on a C-terminal SER — bypassed it because the
  input had no coincident pair to detect; OpenMM's `Modeller.addHydrogens`
  then placed the newly-added HG right on top of the existing OXT via
  its CSER template path. dvbfixer now runs
  `dvbfixer.ffutils.geometry.repair_misplaced_hydrogens` immediately after
  every `addHydrogens` invocation across `prepare`, `protonate`, and
  `minimize`. The helper walks every hydrogen, verifies distance to its
  bonded heavy-atom parent, and repairs any H that landed > 1.5 Å from
  its parent OR within 0.5 Å of another atom in the same residue by
  re-placing it in linear-anti direction at the canonical O-H / N-H /
  C-H bond length. Verified on the reported reproducer: `HG` moves from
  0.001 Å apart from `OXT` (broken) to 0.97 Å from `OG` (canonical).

### Added
- `src/dvbfixer/ffutils/geometry.py` — new module hosting
  `repair_misplaced_hydrogens(topology, positions, verbose=False)`.

## [0.4.1] — 2026-07

### Fixed
- **`prepare`** — coincident-atom detection before hydrogen strip. When
  an input file placed a hydrogen at exactly another atom's position
  (reported against `test/broken_SER/SER.pdb`: SER 126's `HG` sitting
  0.001 Å from `OXT`), `Modeller.addHydrogens` was re-placing `HG` at
  the same coincident position, leaving it 1.7 Å from its own `OG`
  after prepare — a broken sp3 hydroxyl geometry that then corrupted
  the downstream `zbs` protonate step. Prepare now scans every protein
  residue for `(H, heavy_atom)` pairs within 0.5 Å and, when found,
  strips BOTH atoms so PDBFixer's `addMissingAtoms` re-adds the heavy
  atom in its correct position (typically the OXT terminal-atom slot)
  and `addHydrogens` places `H` without the interfering coincident
  atom. Verified by a new regression test at
  `tests/test_prepare_broken_geom.py` that runs on the reported fixture.

## [0.4.0] — 2026-07

### Refactor
- **Guardrails (Phase 0)** — added `tests/` with pytest scaffolding,
  22 seed regression tests for the lightweight subcommands, `[tool.ruff]`
  / `[tool.mypy]` / `[tool.pytest.ini_options]` in `pyproject.toml`,
  `.github/workflows/ci.yml` (fast + full lanes), and a mechanical
  `ruff --fix` pass across `src/dvbfixer/`.
- **Dedup shared helpers (Phase 1)**:
  - `ffutils/variants.py` — consolidates AMBER/CHARMM variant capture,
    text-level rename, and post-`addHydrogens` restore that
    prepare/minimize/protonate/transplant each re-implemented.
  - `ffutils/dat.py` — `DatRecord` dataclass codifying the pipeline's
    `.dat` schema (added_atoms, variant_overrides, removed_residues,
    residue_summary, templates, target_chains). `total_added` is
    derived on save. Six modules were hand-rolling `json.load` /
    `json.dump`; all now go through this.
  - `ffutils.build_glycam_system` — one-call wrapper for the
    `create_forcefield_with_openff` + `loadHydrogenDefinitions` +
    `add_glycam_bonds` ritual.
  - `pdbutils/` split into `inference.py` (OpenBabel-based CONECT
    perception + domain overrides) and `io.py` (line-level serial
    remap / append-before-end helpers).
- **God-module splits (Phase 2)** — `minimize`, `model`, `prepare`,
  `top` converted from flat modules to packages. The follow-up wave
  after the initial argparse extraction cleaned the rest along the
  boundaries the plan called out:
  - `minimize/` → `cli.py` / `pipeline.py` / `refine.py`.
  - `model/` → `cli.py` / `pipeline.py` / `renumber.py` (the
    three-tier residue-numbering strategy) / `modeller_run.py`
    (`_PinnedLoopModel` + `run_modeller` + PIR helpers +
    `_fix_terminal_alignment` + `_explain_modeller_error`).
  - `prepare/` → `cli.py` / `pipeline.py` / `mutations.py`
    (`parse_mutations` + `apply_deletions_to_pdb_text` + SSBOND
    repair + glycan-walk BFS) / `glycan.py` (glycosylation detection,
    NLN/OLS/OLT rename, RDKit/OpenBabel heterogen-H).
  - `top/` → `cli.py` / `pipeline.py` / `ff_data.py` (data-only
    constants: PDB_TO_*, ION_PARAMS, `_GLYCAN_LINKAGE_PARAMS`) /
    `writers.py` (moleculetype/PDB/posre + FF-content splicing).
    Remaining `TopologyBuilder` + `--acpype` extractions
    documented as follow-ups in the package `__init__.py`.
- **Types + hardening (Phase 3)** — annotated the API-surface modules
  (`cli.py`, `ffutils/`, `pdbutils/`, `align.py`); flipped mypy CI
  gate from advisory to required for those paths.
- **Documentation (Phase 4)**:
  - `scripts/gen_cli_reference.py` generates `docs/reference/{cmd}.md`
    from each subcommand's `parse_args()`. CI runs
    `gen_cli_reference.py --check` on the full lane so committed
    reference pages never drift from the actual `--help`.
  - `CLAUDE.md` shrunk from ~85 KB to ~3.6 KB (pure index of the
    docs tree + hard rules). Historical design notes moved verbatim
    to `docs/DESIGN_NOTES.md`, then split into their targeted homes:
    the 4 architecture-adjacent Notes sections went to
    `ARCHITECTURE.md`; the 17 per-subcommand "algorithm" sections
    went to `docs/commands/{cmd}.md` under a `## How it works`
    heading. `DESIGN_NOTES.md` now a ~20-line stub pointer.
  - `CHANGELOG.md` (this file) created and backfilled.

## [0.3.0] — 2026-07

### Added
- **`protonate`** — FF-aware output naming: renames AMBER protonation
  variants to CHARMM equivalents (HID → HSD, HIE → HSE, HIP → HSP,
  CYX → CYS via SSBOND, CYM → CYM) when `--ff charmm` is resolved.
  ASH/GLH/LYN fall back to standard names with a warning (no CHARMM XML
  templates).
- **`protonate`** — `--protassign` (default ON) wraps MolProbity Reduce
  for HIS tautomer picking and ASN/GLN side-chain flip detection.
- **`model`** — fast Python-only no-gap pre-check: string-compares
  ATOM-derived sequence against SEQRES per chain; on match, skips
  Modeller entirely (~2 min → ~0.25 s on gap-free TCR/MHC inputs).
- **`model`** — `--pin-input` (default ON) freezes flankers during
  Modeller's loop-refinement MD via `_PinnedLoopModel` subclass. The
  initial automodel CG still runs on all atoms.
- **`model`** — `--num-output N` saves the top-N ranked candidates by
  molpdf; filenames get `_N` suffix.
- **`minimize`** — `--parametrize-ligands`: on-the-fly GAFF2 + AM1-BCC
  template generation for unknown ligands via antechamber /
  `GAFFTemplateGenerator`. Cached to `~/.cache/dvbfixer/lig_params/`
  (override with `$DVBFIXER_LIG_CACHE`).
- **`zbs`** — `--align-to-input` (default ON): line-level Kabsch
  rewrite after each pipeline step keeps interim outputs in the
  original input's Cartesian frame. Preserves SEQRES / CONECT / all
  non-ATOM records.

### Changed
- **`glycam`** subcommand renamed to **`convert`**, with a
  deprecation-warning alias. Now supports `--to-charmm` in addition
  to the default `--to-amber`; also renames AMBER↔CHARMM protonation
  variants and stale protonation H atoms.
- **`ffutils`** — new `--ff` short-name aliases (`auto`, `amber`,
  `amber+glycam`, `charmm`, ...) accepted by every OpenMM-using tool.
  Auto-detect scans for CHARMM protonation names, GLYCAM sugar codes,
  and CHARMM-GUI 4-char sugars; ambiguous PDB sugar names emit a
  "convert first" warning.
- **`prepare`** — heterogen hydrogen addition switched from the broken
  SMIRNOFF path to real AMBER14+GLYCAM_06j-1 templates via
  `ffutils.create_forcefield_with_openff`. SMIRNOFF failed on
  cross-residue glycosidic bonds.
- **`minimize`** — default is now full-system minimization (protein +
  glycans + ligands). Legacy strip-and-splice behaviour available via
  `--strip-heterogens`; auto-fallback when GLYCAM parametrization
  fails on PDB-named sugars.
- **`top`** — `--water` now drives ion LJ parameter selection via
  new `--ion-set` flag. TIP3P → Joung-Cheatham; SPC-E → JC-SPCE;
  TIP4P-Ew → JC-TIP4P-Ew; OPC → Li-Merz HFE. Prevents the pre-flag
  bug where OPC water + Dang Cl⁻ ions crashed LINCS in ~8 ps.
- **`zbs`** — reduced from 7 to 6 steps; the two `protonate
  --no-hydrogens` calls were replaced by one full `protonate` between
  the minimize passes.

### Fixed
- **`protonate`** — text-level pre-rename of AMBER/CHARMM variants to
  standard parents before OpenMM parses. Fixes the LYN-blocks-
  peptide-bond crash where a mid-chain LYN would cause `addHydrogens`
  to fail on the PREVIOUS residue with "missing 1 externally bonded
  C atom".
- **`model`** — three-tier deterministic residue-number strategy
  (K-finder → Needleman-Wunsch → align2d mask). Fixes FcgRI chain A
  regression where 5 modeled residues at 219-223 ended up numbered
  283-287. HETATM resseqs are preserved verbatim.
- **`model`** — N-terminal gaps extend BACKWARD from the first
  template residue instead of starting at 1.
- **`model`** — `_reorder_chains_for_modeller` groups each chain ID
  into one contiguous block, fixing silent NAG drops (3ry6 chain C)
  and BLK alignment crashes (FcgRI).

## [0.2.0] and earlier

Historical, best-effort roll-up of the pre-0.3.0 milestones. Consult
git log for exact commit provenance.

### Added
- **`transplant`** — Kabsch-aligned graft of GLYCAM-Web output into
  an acceptor PDB; `--relax` runs AMBER14+GLYCAM minimization;
  `--gromacs` exports via ACPYPE with `[ pairs_nb ]` for mixed 1-4
  scaling.
- **`top`** — RTP-based GROMACS topology generation (AMBER99SB-ILDN
  + CHARMM36) with modular `.itp` output; `--acpype` alternative
  mode via OpenMM → ParmEd → ACPYPE for the mixed 1-4 case.
- **`top`** — glycan and glycolipid support (CHARMM36) with distance-
  based glycosidic-bond detection, charge redistribution, and
  ceramide-sugar linkage parameters.
- **`parametrize`** — GAFF2 small-molecule pipeline (antechamber →
  parmchk2 → tleap → ParmEd) with AM1-BCC default. RESP available
  via `--qm-engine gaussian` (existing two-step flow),
  `--qm-engine psi4` (free, subprocess to a separate conda env), or
  `--qm-engine pyscf` (recommended, pure-Python in the main env).
- **`homology`** — multi-template Modeller with automatic chain
  mapping; `--antibody` mode uses ANARCI for Kabat / IMGT numbering
  and CDR detection.
- **`cluster`** — glycan conformational clustering from MD
  trajectories via glycosidic-torsion RMSD (GFDB method); global vs
  per-linkage modes; interactive Plotly output.
- **`renumber`** — antibody scheme support (Kabat / Chothia / IMGT /
  Martin / Aho / EU) via ANARCI for V-domains + bundled human IgG1 /
  Cκ / Cλ EU references for C-domains.
- **`prepare`** — `--mutate` with substitution (`A:39:CYX`) and
  deletion (`A:39:del`) support, including glycan-walk BFS and
  disulfide-partner repair.
- **`conect`** — standalone subcommand + shared inference engine
  auto-called by prepare / top / minimize / transplant / convert.
- **`puppet`** — strip PDB to backbone-only polyglycine (with
  `--keep` for residues to preserve intact).
- **`split`** — multi-MODEL detection: preserves per-MODEL chain-ID
  consistency instead of cascading letters (A B C / D E F / …).

### Environment
- Migrated from separate CLI scripts to a single `dvbfixer` command
  with subcommand dispatch (`dvbfixer.cli`).
- Bundled AMBER99SB-ILDN and CHARMM36 force-field directories in
  `FF/` for RTP parsing; no external GROMACS FF dir needed at
  runtime.
- Added optional runtime deps: `plotly` (cluster), `hmmer` + `anarci`
  (antibody modes), `xtb` + `openbabel` (minimize refinement), and
  `pyscf` (parametrize RESP).
