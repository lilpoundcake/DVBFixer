# dvbfixer 0.3.0 → 0.4.0: Phase 0–4 revision + bug fix + protonate UX + CLI reorg

Consolidates ~2 months of internal refactoring plus a hard-crash bug
fix, a new `--no-propka` escape hatch, and grouped `--help` output
across every subcommand. 25 commits, 99 files changed, +12,549/−8,642.

Version bumps: `pyproject.toml` `0.3.0 → 0.4.0`, `src/dvbfixer/__init__.py.__version__` `0.1.0 → 0.4.0` (previously drifted).

## Highlights

### Bug fix (crash on `zbs` step 5) — `bce0dda`
- `ffutils.variants.restore_variants_post_addhydrogens` used to
  crash with `ValueError: not enough values to unpack (expected 3,
  got 2)` when its fallback loop hit a 2-tuple key in `saved`. Real
  users hit this on the default `zbs` flow. Fallback now guards on
  tuple arity; 4 regression tests cover it.

### `--no-propka` in `protonate` + `zbs` — `1445c60`
- New `--propka`/`--no-propka` (BooleanOptionalAction, default ON)
  next to the existing `--protassign`. `--no-propka` skips PROPKA
  entirely and relies on MolProbity Reduce for HIS tautomers + ASN/GLN
  flips. Combining `--no-propka --no-protassign` is a hard error.
- `zbs` gets the matching passthrough. The step-5 banner now shows
  the actual engines being invoked, e.g.
  `Step 5: PROTONATE (Reduce + addHydrogens)`.

### CLI reorg — `97f6dda`
- Every subcommand's `parse_args` converted to `add_argument_group`
  so `--help` renders visually chunked. 17 tools touched, 216 flags
  regrouped. `docs/reference/*.md` regenerated; no semantic change.

### Refactor (Phases 0–4)

**Phase 0 — Guardrails (`a00dbbf`)**
- `pytest` scaffolding, `tests/` with 22 seed regression tests,
  `[tool.ruff]` / `[tool.mypy]` / `[tool.pytest.ini_options]` in
  `pyproject.toml`, `.github/workflows/ci.yml` (fast + full lanes).

**Phase 1 — Dedup shared helpers**
- `c32c5cb` `ffutils/variants.py` — consolidates AMBER/CHARMM
  variant capture/restore that `prepare` / `minimize` / `protonate` /
  `transplant` each re-implemented.
- `a77be2a` `ffutils/dat.py` — `DatRecord` dataclass codifies the
  pipeline's `.dat` schema. Six modules were hand-rolling `json.load`
  / `json.dump`; all now go through `DatRecord.load` / `.save` /
  `.merge`.
- `3ddc6c3` `ffutils.build_glycam_system` — one-call wrapper for the
  3-step GLYCAM ritual (`create_forcefield_with_openff` +
  `loadHydrogenDefinitions` + `add_glycam_bonds`).
- `5b41436` `pdbutils/` split into `inference.py` (OpenBabel-based
  CONECT perception + domain overrides) and `io.py` (line-level
  serial remap / append-before-end helpers).

**Phase 2 — God-module splits**
- `minimize` → `cli.py` / `pipeline.py` / `refine.py`
  (`7687ca1`).
- `model` → `cli.py` / `pipeline.py` / `renumber.py`
  (three-tier residue-number strategy) / `modeller_run.py`
  (`_PinnedLoopModel`, `run_modeller`, PIR helpers)
  (`b0a2655` + `a7bac6b`).
- `prepare` → `cli.py` / `pipeline.py` / `mutations.py`
  (`parse_mutations`, `apply_deletions_to_pdb_text`, SSBOND repair,
  glycan-walk BFS) / `glycan.py` (glycosylation detection, NLN/OLS/OLT
  rename, RDKit/OpenBabel heterogen-H)
  (`384e2b4` + `7511836`).
- `top` → `cli.py` / `pipeline.py` / `ff_data.py` (data-only
  constants: `PDB_TO_*`, `ION_PARAMS`, `_GLYCAN_LINKAGE_PARAMS`) /
  `writers.py` (moleculetype/PDB/posre + FF-content splicing) /
  `acpype.py` (extracted `--acpype` mode + 13 unit tests)
  (`2d37c17` + `2599a6e` + `5e07dc5` + `47ac2be`).

**Phase 3 — Types + hardening (`ff76530`)**
- Annotated `cli.py`, `ffutils/`, `pdbutils/`, `align.py` (the "API
  surface" the plan targeted for strict typing). `mypy` CI gate
  flipped from advisory to required for those modules.

**Phase 4 — Documentation revision**
- `5eaa261` `scripts/gen_cli_reference.py` generates
  `docs/reference/{cmd}.md` from each subcommand's `parse_args()`.
  CI runs `--check` on the full lane so committed reference pages
  never drift from the actual `--help`.
- `38cddb7` `CLAUDE.md` shrunk from ~85 KB to ~3.6 KB (pure index of
  the docs tree + hard rules); the body moved verbatim to
  `docs/DESIGN_NOTES.md`.
- `b8c9657` The 4 architecture-adjacent Notes sections in
  DESIGN_NOTES → `ARCHITECTURE.md`.
- `308820f` The 17 per-subcommand "algorithm" sections migrated
  under a `## How it works` heading in each
  `docs/commands/{cmd}.md`. `DESIGN_NOTES.md` reduced to a
  ~20-line stub pointer.
- `c4fdcea` + `0b943d0` `CHANGELOG.md` created and backfilled from
  git history through 0.4.0.

## Test / lint / type / docs gates

- `pytest tests/ -q` → **57 tests, all passing** (up from 0 pre-PR).
- `ruff check src/dvbfixer tests` → clean.
- `mypy src/dvbfixer/{cli.py,ffutils,pdbutils,align.py}` → clean.
- `python scripts/gen_cli_reference.py --check` → in sync.

## Non-breaking

Every subcommand keeps its name, argparse surface (all existing
flags), and I/O contract. Every `from dvbfixer.X import
main / parse_args / ...` that existed at 0.3.0 still resolves —
packages re-export through `__init__.py`.

The `docs/reference/*.md` layout DID change (new group headings)
but that file tree is machine-generated; it's not part of the API.

## What's queued but not in this PR

Follow-ups documented in each package's `__init__.py` and in
`ARCHITECTURE.md`'s "Recommended areas for future work":
- `top/rtp_build.py` — extracting `TopologyBuilder` + `build_chain`
  family (~1175 LOC still in `top/pipeline.py`, tightly coupled to
  dataclasses defined there).
- Extend `minimize --parametrize-ligands` GAFF2 path to `prepare`
  for whole-system H addition.
- `xtb v6.8+` upgrade when conda-forge ships it (fixes the `$fix`
  bug that causes ~0.1 Å drift on nominally-frozen anchors).

## Test plan

- [ ] `pytest tests/ -q` in fresh checkout on this branch.
- [ ] `ruff check src/dvbfixer tests` clean.
- [ ] `mypy src/dvbfixer/{cli.py,ffutils,pdbutils,align.py}` clean.
- [ ] `python scripts/gen_cli_reference.py --check` returns 0.
- [ ] `dvbfixer zbs test/ASN.pdb -o /tmp/out --no-solvent` succeeds
      (was the reproducer for the pre-fix crash).
- [ ] `dvbfixer protonate test/ASN.pdb --no-propka -v` skips PROPKA
      cleanly; runs Reduce.
- [ ] `dvbfixer protonate test/ASN.pdb --no-propka --no-protassign`
      hard-errors with a clear message.
- [ ] `dvbfixer <cmd> --help` for each subcommand shows grouped
      sections.
