# Changelog

All notable changes to dvbfixer are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Backfilled from git history — commits before v0.3.0 are grouped by
feature area rather than by strict release. Older entries are
best-effort summaries; consult `git log` for exact provenance.

## [Unreleased]

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
