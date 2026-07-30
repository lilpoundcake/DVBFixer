# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## What this project is

**dvbfixer** — a Python package providing CLI tools for preparing PDB
(Protein Data Bank) structural biology files. Installed as a single
`dvbfixer` command with 17 subcommands.

## Ongoing work: `feat/tleap-reduce-backend` branch

An in-progress replacement of the H+heavy-atom repair pipeline lives on
branch `feat/tleap-reduce-backend`. It swaps `PDBFixer.addMissingAtoms`
+ `Modeller.addHydrogens` for **AmberTools `tleap` + MolProbity
`reduce`** (subprocess-based, deterministic, L-only by construction).
Fixes the D-Cα (openmm/pdbfixer#145) and coincident-H
(Modeller.addHydrogens bug) issues that surface on gap-filled model
outputs. Verified against `test/shit/{1EMV,1FR2,2VLN,2VLQ}_original.pdb`
— all four yield zero D-Cα, zero coincident atoms.

On main the legacy Modeller+PDBFixer path is still the only backend;
its known-good boundaries are the hard rules below. If you need a
D-Cα-free deterministic pipeline for standard proteins, use the branch.

## Where to look

| For… | Read… |
|---|---|
| Command inventory + quick start | [`README.md`](README.md) |
| Per-command "how to use" prose | [`docs/commands/`](docs/commands/) |
| Exact `--help` output (auto-generated, matches code) | [`docs/reference/`](docs/reference/) |
| Module map, data flow, `.dat` schema | [`ARCHITECTURE.md`](ARCHITECTURE.md) |
| End-to-end recipes | [`docs/pipelines.md`](docs/pipelines.md), [`BEST_PRACTICES.md`](BEST_PRACTICES.md) |
| Known issues + active workarounds | [`docs/known-issues.md`](docs/known-issues.md) |
| Historical design notes, gotchas | [`docs/DESIGN_NOTES.md`](docs/DESIGN_NOTES.md) |
| Installation + Modeller license | [`docs/installation.md`](docs/installation.md) |
| Force-field selection matrix | [`docs/force-fields.md`](docs/force-fields.md) |

## Environment

```bash
micromamba create -f environment.yml
micromamba run -n dvbfixer pip install -e ".[dev]"
```

`.[dev]` pulls in pytest / ruff / mypy. Without dev extras the tools
still install and run.

Modeller needs a free academic license from
<https://salilab.org/modeller/registration.html>. Set the key in
`<env>/lib/modeller-10.8/modlib/modeller/config.py`.

**Python is pinned `>=3.11,<3.14`** in `environment.yml`. Do not loosen
this. propka 3.5.1 reads `self.__annotations__` (instance-level) inside
its `Parameters` dataclass; Python 3.14's PEP 649/749 annotation change
makes instance `__annotations__` raise `AttributeError`, which crashes
`dvbfixer protonate` / `prepare` at the PROPKA step. Repro:
`python -c "from propka.parameters import Parameters; Parameters().__annotations__"`
returns a dict on <3.14, raises on 3.14.

**On macOS Docker / VirtioFS, create the env on the container's native
overlay filesystem, not on a host bind mount.** `/home/agent` here is a
`fakeowner` (macOS `/Users`) bind mount that is case-insensitive and
rejects libmamba's `copy_symlink` on ncurses's case-variant terminfo
pairs (`share/terminfo/32/2621A` vs `2621a`), failing the install with
`filesystem error: cannot copy symlink: Invalid argument`. `always_copy`
/ `--copy` do NOT help — copy mode still recreates symlinks rather than
dereferencing them, and there is no micromamba flag to skip/deref them.
Fix: point the root at an overlay path, e.g.
`sudo mkdir -p /opt/mamba && sudo chown -R agent:agent /opt/mamba`,
then `export MAMBA_ROOT_PREFIX=/opt/mamba` before `micromamba create`,
and keep `/opt/mamba/envs/dvbfixer/bin` on `PATH` so subprocesses that
call bare `dvbfixer` resolve it. See
[docs/known-issues.md](docs/known-issues.md) for the full write-up.

## Working in this codebase — hard rules

- **Default prep backend is `legacy`** (Modeller+PDBFixer). Handles
  every input class dvbfixer supports (glycans, ligands, PTMs,
  covalent HETATM links). The chirality invariant is guaranteed by
  minimize's post-phase-2 unconditional force-reflect (0.7.4+), so
  PDBFixer's D-Cα risk is neutralised downstream — legacy is safe
  by construction now. **`tleap-reduce`** is retained as opt-in
  (`--backend tleap-reduce`) for pure-protein inputs where the user
  specifically wants deterministic L-only heavy atoms from tleap
  itself; it rejects non-canonical residues. When the tleap-reduce
  backend IS selected, all its rules below still apply.
- **Do not remove `keepIds=True`** from any `PDBFile.writeFile` call.
  Losing chain IDs mid-pipeline breaks the `.dat` handoff that
  `minimize` uses for tiered restraints.
- **PROPKA + MolProbity Reduce run INSIDE the legacy prepare backend**
  (since 0.7.7). The helper
  `dvbfixer.prepare.pipeline._run_propka_reduce_variants` produces
  the variant map for `Modeller.addHydrogens(variants=[...])`;
  disable per-flag via `--no-propka` / `--no-protassign`. Standalone
  `dvbfixer protonate` still exists as a post-hoc re-protonation
  tool.
- **Do not call `PDBFixer.addMissingHydrogens(pH)`.** On the new
  backend, H placement is Reduce's job (subprocess). On the legacy
  backend, use `modeller.addHydrogens(forcefield, pH=..., variants=[...])`.
  PDBFixer's `addMissingHydrogens` uses its own `_describeVariant`
  that only recognises standard PDB names (HIS / ASP / GLU / CYS /
  LYS) — it ignores AMBER variant labels (HIE / HID / HIP / ASH /
  GLH / CYX / CYM / LYN). On the legacy path, follow every
  `Modeller.addHydrogens` call with
  `dvbfixer.ffutils.geometry.repair_misplaced_hydrogens(topology,
  positions)` (widened 2026-07-24 to also break sibling-H clashes
  via proper sp3 tetrahedral placement).
- **Chirality invariant**: `dvbfixer.ffutils.geometry.assert_all_l`
  must succeed after every heavy-atom repair. `find_d_residues`
  and `fix_ca_chirality` are the detector and reflector primitives.
  tleap is L-only by construction so `assert_all_l` should never
  trip on the new backend — it fires only if a downstream bug
  regresses. Post-minimize enforcement is bounded reflect+re-minimize
  (3 iters) then UNCONDITIONAL force-reflect: the output MUST have
  zero D-Cα, even if it means accepting minor local packing strain
  on the rare residue whose FF minimum genuinely lies on the D side.
  Do NOT re-introduce a WARN-only path — the chirality invariant is
  non-negotiable.
- **Do not read a `.dat` file with hand-rolled `json.load`.** Use
  `dvbfixer.ffutils.dat.DatRecord` — the schema (added_atoms,
  variant_overrides, removed_residues, residue_summary,
  templates, target_chains) lives there, `total_added` is derived
  on save, and merge semantics ("downstream wins on collision") are
  codified.
- **Do not re-implement AMBER/CHARMM variant handling.** The four
  functions in `dvbfixer.ffutils.variants`
  (`scan_variant_names`, `text_rename_variants_to_parent`,
  `rename_variants_to_parent_in_topology`,
  `restore_variants_post_addhydrogens`) are the canonical dance for
  making OpenMM's peptide-bond inference survive a mid-chain
  LYN/HIE/HSD. `minimize/pipeline.py` calls
  `text_rename_variants_to_parent` BEFORE `PDBFile.load` so OpenMM
  can infer intra-residue bonds — its `_standardResidues` set does
  not include AMBER variants and topology has zero bonds for them
  otherwise.
- **Do not use `_classify_variant` or `_STD_PKA`** (both deleted in
  0.7.4). PROPKA-driven variant decisions belong in
  `dvbfixer.protonate.decide_protonation` — it filters by PROPKA
  group type (so N+/C- pKas don't overwrite side-chain pKas) and
  keys by `(chain, resnum, icode)`. `prep_backend.run_prep` layers
  Reduce's HID/HIE tautomer under PROPKA's HIP decision. Do NOT
  bring back the raw-pka dict path.
- **Force strip-H + `Modeller.addHydrogens(variants=[...])` when
  `amber_renames` is non-empty in minimize.** `addHydrogens` copies
  input bonds but never rebuilds missing ones from templates, and
  only ADDS missing H (never REMOVES extras). Strip-and-readd is
  the only path that produces a topology whose atoms AND bonds
  match the variant template. The keep-H fast path stays for the
  no-variant common case.
- **Cap the pH passed to `Modeller.addHydrogens` at 9.99.** OpenMM's
  `hydrogens.xml` gates HZ3 with `maxph="10.0"` and at pH > 10
  produces LYS residues that miss HZ3 (template mismatch). Variants
  from PROPKA already carry the true protonation state, so the
  auto-selection heuristic driven by pH is redundant.
- **Do not unconditionally drop SG-SG bonds in
  `_drop_spurious_inter_aa_bonds`.** The filter drops non-peptide
  inter-residue bonds between two standard AAs, but disulfides
  (SG-SG on CYS-family residues: CYS/CYX/CYM) are kept — dropping a
  *genuine* one severs the disulfide during minimize and lets it
  relax apart (fixed 0.7.4). Since 0.7.8 the function also resolves
  *spurious extra* SG-SG bonds down to one partner per atom (via
  `valid_ss_pairs` snapshot-and-restore, or nearest-distance matching
  when `positions` is given) — OpenMM's own `PDBFile.__init__` and
  `PDBFixer.addMissingAtoms()` both call `Topology.
  createDisulfideBonds()`, a pure distance-cutoff scan with no 1:1
  matching, so tightly-packed CYS clusters can otherwise give one SG
  atom two or three "partners" and crash `createSystem` with "N S
  atom(s) too many". The invariant is "every SG has at most one SG-SG
  partner, and a genuine one is never removed" — not "every SG-SG
  bond survives unconditionally."
- **Do not skip the `build_glycam_system` wrapper.** It packages
  the three-step GLYCAM ritual (`create_forcefield_with_openff` +
  `loadHydrogenDefinitions('glycam-hydrogens.xml')` +
  `add_glycam_bonds`); reproducing it inline is how prep/minimize/
  protonate got out of sync last time.
- **Do not add `# type: ignore` to files under strict mypy** without
  a concrete reason. The strict-target list in `pyproject.toml`
  (`cli.py`, `ffutils/`, `pdbutils/`, `align.py`) is enforced in CI.
- **Do not edit `docs/reference/*.md` by hand.** They're generated
  from each subcommand's `parse_args()`. Run
  `python scripts/gen_cli_reference.py` after touching any argparse.
  CI enforces this with `--check`.

## Running tests locally

```bash
pytest tests/ -q          # fast lane, ~2 s, no OpenMM/Modeller needed
ruff check src/dvbfixer   # style
mypy src/dvbfixer/cli.py src/dvbfixer/ffutils src/dvbfixer/pdbutils src/dvbfixer/align.py
```

The full test suite (Modeller-touching integration cases) needs the
scientific stack from `environment.yml`; it runs in CI's full-lane job.

## Session preferences

- **Always pass `--no-solvent`** on every dev iteration of `minimize` /
  `zbs`. Solvent-box minimize is orders of magnitude slower and eats the
  feedback loop. Add solvent back only for real evaluation runs the user
  has explicitly asked for.
- Reply in English even when the user writes in Russian.
