# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## What this project is

**dvbfixer** — a Python package providing CLI tools for preparing PDB
(Protein Data Bank) structural biology files. Installed as a single
`dvbfixer` command with 17 subcommands.

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

## Working in this codebase — hard rules

- **Default prep backend is `tleap-reduce`** (added on this branch
  `feat/tleap-reduce-backend`). `prepare` and `protonate` route
  through `dvbfixer.prep_backend.run_prep` which runs
  `tleap` (heavy atoms) + `reduce -build -nuclear` (H) + PROPKA
  (AMBER variant renames). Deterministic, L-only by construction,
  no coincident-H bug. Ships with `ambertools>=23`. The legacy
  Modeller+PDBFixer path stays behind `--backend legacy` for
  GLYCAM glycoproteins and exotic heterogens tleap rejects.
- **Do not remove `keepIds=True`** from any `PDBFile.writeFile` call.
  Losing chain IDs mid-pipeline breaks the `.dat` handoff that
  `minimize` uses for tiered restraints.
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
  regresses.
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
  LYN/HIE/HSD.
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

- Default to `--no-solvent` on dev iterations of `minimize` / `zbs` —
  solvent-box minimize is too slow for the loop. Add solvent back only
  for real evaluation runs.
- Reply in English even when the user writes in Russian.
